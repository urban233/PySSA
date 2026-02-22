"""Name registry for tracking reserved names across the open project.

The ``NameRegistry`` is the single source of truth for which names are
currently unavailable for new objects (sequences, proteins, protein pairs)
within the active project.  It covers two types of reservations:

1. **Persisted objects** — names already belonging to objects stored in the
   project (proteins, protein pairs, sequences).  These are loaded via
   ``rebuild()`` whenever a project is opened.

2. **Job reservations** — names that queued or running jobs will produce on
   completion.  These are registered by the ``JobScheduler`` on ``submit()``
   and released on job completion, failure, or cancellation.

The registry is owned exclusively by ``AppState`` and exposed as a read-only
property.  Dialog controllers query it directly via
``app_state.name_registry.is_reserved(category, name)`` instead of
receiving a separate ``Watcher`` object.
"""
from __future__ import annotations

import logging
from typing import Iterable, TYPE_CHECKING

if TYPE_CHECKING:
    from src.pyssa.internal.data_structures import project

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Category constants
# ---------------------------------------------------------------------------

SEQUENCE: str = "sequence"
"""Category key for sequence names."""

PROTEIN: str = "protein"
"""Category key for protein names."""

PROTEIN_PAIR: str = "protein_pair"
"""Category key for protein pair names."""

_VALID_CATEGORIES: frozenset[str] = frozenset({SEQUENCE, PROTEIN, PROTEIN_PAIR})


class NameRegistry:
    """Tracks which names are currently reserved within the active project.

    A name becomes reserved either because it belongs to an existing project
    object (loaded via ``rebuild``) or because a queued/running job will
    produce an object with that name (managed by the ``JobScheduler``).

    All lookups are O(1) thanks to the underlying ``set`` storage.
    """

    def __init__(self) -> None:
        """Initialise an empty registry."""
        self._sets: dict[str, set[str]] = {
            SEQUENCE: set(),
            PROTEIN: set(),
            PROTEIN_PAIR: set(),
        }

    # ------------------------------------------------------------------
    # Bulk operations
    # ------------------------------------------------------------------

    def rebuild(self, a_project: "project.Project") -> None:
        """Rebuild the registry from the given open project.

        Clears all existing reservations and re-populates from:

        * ``project.sequences`` (``SeqRecord.name``)
        * ``project.proteins``  (``Protein.get_molecule_object()``)
        * ``project.protein_pairs`` (``ProteinPair.name``)

        Args:
            a_project: The currently open project domain object.

        Raises:
            ValueError: If ``a_project`` is ``None``.
        """
        if a_project is None:
            logger.error("a_project is None.")
            raise ValueError("a_project must not be None.")

        self._sets[SEQUENCE] = {seq.name for seq in a_project.sequences}
        self._sets[PROTEIN] = {
            prot.get_molecule_object() for prot in a_project.proteins
        }
        self._sets[PROTEIN_PAIR] = {pp.name for pp in a_project.protein_pairs}
        logger.info(
            "NameRegistry rebuilt — sequences=%d, proteins=%d, protein_pairs=%d.",
            len(self._sets[SEQUENCE]),
            len(self._sets[PROTEIN]),
            len(self._sets[PROTEIN_PAIR]),
        )

    def clear(self) -> None:
        """Empty all category sets.

        Call this when the active project is closed so that stale reservations
        do not leak into the next session.
        """
        for key in self._sets:
            self._sets[key].clear()
        logger.info("NameRegistry cleared.")

    # ------------------------------------------------------------------
    # Single-name operations
    # ------------------------------------------------------------------

    def reserve(self, category: str, name: str) -> None:
        """Reserve a single name in the given category.

        Args:
            category: One of ``SEQUENCE``, ``PROTEIN``, or ``PROTEIN_PAIR``.
            name: The name to reserve.

        Raises:
            ValueError: If ``category`` is unknown or ``name`` is empty.
        """
        self._validate(category, name)
        self._sets[category].add(name)
        logger.debug("Reserved '%s' in category '%s'.", name, category)

    def release(self, category: str, name: str) -> None:
        """Release a reservation so the name becomes available again.

        This operation is idempotent: releasing a name that was never reserved
        logs a warning but does not raise.

        Args:
            category: One of ``SEQUENCE``, ``PROTEIN``, or ``PROTEIN_PAIR``.
            name: The name to release.

        Raises:
            ValueError: If ``category`` is unknown or ``name`` is empty.
        """
        self._validate(category, name)
        if name in self._sets[category]:
            self._sets[category].discard(name)
            logger.debug("Released '%s' from category '%s'.", name, category)
        else:
            logger.warning(
                "Attempted to release '%s' from category '%s', but it was not reserved.",
                name,
                category,
            )

    def is_reserved(self, category: str, name: str) -> bool:
        """Check whether a name is reserved in the given category.

        Args:
            category: One of ``SEQUENCE``, ``PROTEIN``, or ``PROTEIN_PAIR``.
            name: The name to check.

        Returns:
            ``True`` if the name is currently reserved, ``False`` otherwise.

        Raises:
            ValueError: If ``category`` is unknown or ``name`` is empty.
        """
        self._validate(category, name)
        return name in self._sets[category]

    # ------------------------------------------------------------------
    # Bulk single-category operations (used by the scheduler)
    # ------------------------------------------------------------------

    def reserve_many(self, category: str, names: Iterable[str]) -> None:
        """Reserve multiple names in the given category at once.

        Args:
            category: One of ``SEQUENCE``, ``PROTEIN``, or ``PROTEIN_PAIR``.
            names: An iterable of names to reserve.

        Raises:
            ValueError: If ``category`` is unknown.
        """
        self._validate_category(category)
        for name in names:
            if name:
                self._sets[category].add(name)
                logger.debug("Reserved '%s' in category '%s'.", name, category)

    def release_many(self, category: str, names: Iterable[str]) -> None:
        """Release multiple reservations in the given category at once.

        Releases are idempotent; missing names are logged at WARNING level.

        Args:
            category: One of ``SEQUENCE``, ``PROTEIN``, or ``PROTEIN_PAIR``.
            names: An iterable of names to release.

        Raises:
            ValueError: If ``category`` is unknown.
        """
        self._validate_category(category)
        for name in names:
            if name:
                if name in self._sets[category]:
                    self._sets[category].discard(name)
                    logger.debug("Released '%s' from category '%s'.", name, category)
                else:
                    logger.warning(
                        "Attempted to release '%s' from category '%s', but it was not reserved.",
                        name,
                        category,
                    )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _validate_category(self, category: str) -> None:
        """Raise ``ValueError`` for unknown categories.

        Args:
            category: The category string to validate.

        Raises:
            ValueError: If ``category`` is not one of the valid categories.
        """
        if category not in _VALID_CATEGORIES:
            raise ValueError(
                f"Unknown category '{category}'. "
                f"Valid categories are: {sorted(_VALID_CATEGORIES)}."
            )

    def _validate(self, category: str, name: str) -> None:
        """Raise ``ValueError`` for unknown categories or empty names.

        Args:
            category: The category string to validate.
            name: The name string to validate.

        Raises:
            ValueError: If ``category`` is unknown or ``name`` is empty.
        """
        self._validate_category(category)
        if not name:
            raise ValueError("name must not be empty.")
