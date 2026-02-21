"""Slot and supporting types for handling item selection in the project QTreeView.

Architecture
------------
``on_project_tree_selection_changed`` is connected to the ``selectionChanged``
signal of the QTreeView that has ``PSAObjectsModel`` set as its model.  It
works identically for single-item and multi-item selections.

The slot's job is deliberately narrow: classify every selected index into one
of the typed buckets, assemble a ``SelectionResult``, then dispatch to a
dedicated handler for each non-empty bucket.  Handlers are independent of each
other — the sequences handler knows nothing about the residues handler — so you
can fill them in one at a time.

Connection example
------------------
::

    tree_view.setSelectionMode(
        QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection
    )
    tree_view.selectionModel().selectionChanged.connect(
        self.on_project_tree_selection_changed
    )
"""
from __future__ import annotations

import logging
import os
import pathlib
from dataclasses import dataclass, field
from typing import Optional

from src.pyssa.gui.qt import QtCore
from src.pyssa.io_pyssa import binary_data
from src.pyssa.model.protein_subtree_mixin import (
    TYPE_SECTION,
    TYPE_SEQUENCE,
    TYPE_PROTEIN,
    TYPE_PROTEIN_PAIR,
    TYPE_HEADER,
    TYPE_SCENE,
    TYPE_CHAIN,
    TYPE_RESIDUE,
    TYPE_ATOM,
    LABEL_SCENES,
    LABEL_CHAINS,
    LABEL_SEQUENCES,
    LABEL_PROTEINS,
    LABEL_PROTEIN_PAIRS,
)
from src.pyssa.util import enums, constants
from src.pyssa.logging_pyssa import log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"

_DR = QtCore.Qt.ItemDataRole.DisplayRole


# ===========================================================================
# Typed per-node selection dataclasses
# ===========================================================================

@dataclass
class SectionSelection:
    """A top-level section divider node ("Sequences", "Proteins", "Protein Pairs").

    Attributes:
        label: One of LABEL_SEQUENCES, LABEL_PROTEINS, LABEL_PROTEIN_PAIRS.
        index: Use ``index.model().rowCount(index)`` to get the direct child count.
    """

    label: str
    index: QtCore.QModelIndex


@dataclass
class SequenceSelection:
    """A sequence leaf node under the "Sequences" section.

    Attributes:
        name: Display name shown in the tree view.
        sequence_string: Full amino acid string stored in OBJECT_ROLE.
        index: QModelIndex of this node.
    """

    name: str
    sequence_string: str
    index: QtCore.QModelIndex


@dataclass
class StandaloneProteinSelection:
    """A protein node directly under the "Proteins" section.

    Attributes:
        protein_name: Molecule object name (e.g. ``"1abc"``).
        protein_object: Full ``Protein`` domain object.
        index: QModelIndex of this node.
    """

    protein_name: str
    protein_object: object
    index: QtCore.QModelIndex


@dataclass
class ProteinPairSelection:
    """A protein-pair node under the "Protein Pairs" section.

    Attributes:
        pair_name: Display name (e.g. ``"1abc_vs_2xyz"``).
        pair_object: Full ``ProteinPair`` domain object.
        index: QModelIndex of this node.
    """

    pair_name: str
    pair_object: object
    index: QtCore.QModelIndex


@dataclass
class PairMemberProteinSelection:
    """A protein node that is a child of a protein-pair node.

    Attributes:
        protein_name: Molecule object name of this member.
        protein_object: ``Protein`` instance for this member.
        pair_object: Parent ``ProteinPair``.
        index: QModelIndex of this node.
    """

    protein_name: str
    protein_object: object
    pair_object: object
    index: QtCore.QModelIndex


@dataclass
class ScenesHeaderSelection:
    """The "Scenes" header node under a protein or protein-pair.

    Attributes:
        owner_object: ``Protein`` or ``ProteinPair`` that owns these scenes.
        owner_type: ``TYPE_PROTEIN`` or ``TYPE_PROTEIN_PAIR``.
        scene_count: Number of scene nodes currently under this header.
        index: QModelIndex of this node.
    """

    owner_object: object
    owner_type: str
    scene_count: int
    index: QtCore.QModelIndex


@dataclass
class ChainsHeaderSelection:
    """The "Chains" header node under a protein node.

    Attributes:
        protein_object: The ``Protein`` that owns these chains.
        chain_count: Number of chain nodes under this header.
        parent_of_protein_type: ``TYPE_PROTEIN_PAIR`` if inside a pair,
            otherwise the enclosing section's type string.
        index: QModelIndex of this node.
    """

    protein_object: object
    chain_count: int
    parent_of_protein_type: str
    index: QtCore.QModelIndex


@dataclass
class SceneSelection:
    """A scene leaf node under a "Scenes" header.

    Attributes:
        scene_name: Name as stored in PyMOL (e.g. ``"base"``).
        owner_object: ``Protein`` or ``ProteinPair`` that owns this scene.
        owner_type: ``TYPE_PROTEIN`` or ``TYPE_PROTEIN_PAIR``.
        index: QModelIndex of this node.
    """

    scene_name: str
    owner_object: object
    owner_type: str
    index: QtCore.QModelIndex


@dataclass
class ChainSelection:
    """A chain node under a "Chains" header.

    Attributes:
        chain_letter: Single-letter chain ID (e.g. ``"A"``).
        chain_object: ``Chain`` domain object.
        chain_color: Current PyMOL colour value (from CHAIN_COLOR_ROLE).
        protein_object: The ``Protein`` that contains this chain.
        parent_of_protein_type: ``TYPE_PROTEIN_PAIR`` if inside a pair,
            otherwise the enclosing section's type string.
        index: QModelIndex of this node.
    """

    chain_letter: str
    chain_object: object
    chain_color: object
    protein_object: object
    parent_of_protein_type: str
    index: QtCore.QModelIndex


@dataclass
class ResidueSelection:
    """A residue node under a chain node.

    Attributes:
        residue_label: ``"resi - resn"`` (e.g. ``"42 - GLY"``).
        residue_object: Residue domain object.
        chain_letter: ID of the chain this residue belongs to.
        protein_object: The ``Protein`` that contains this residue.
        parent_of_protein_type: ``TYPE_PROTEIN_PAIR`` if inside a pair,
            otherwise the enclosing section's type string.
        index: QModelIndex of this node.
    """

    residue_label: str
    residue_object: object
    chain_letter: str
    protein_object: object
    parent_of_protein_type: str
    index: QtCore.QModelIndex


@dataclass
class AtomSelection:
    """An atom leaf node under a residue node.

    Attributes:
        atom_name: Atom name (e.g. ``"CA"``, ``"N"``, ``"O"``).
        residue_label: Parent residue ``"resi - resn"``.
        residue_object: Parent residue domain object.
        chain_letter: Chain this atom belongs to.
        protein_object: The ``Protein`` that contains this atom.
        parent_of_protein_type: ``TYPE_PROTEIN_PAIR`` if inside a pair,
            otherwise the enclosing section's type string.
        index: QModelIndex of this node.
    """

    atom_name: str
    residue_label: str
    residue_object: object
    chain_letter: str
    protein_object: object
    parent_of_protein_type: str
    index: QtCore.QModelIndex


# ===========================================================================
# SelectionResult
# ===========================================================================

@dataclass
class SelectionResult:
    """All currently selected nodes bucketed by their logical type.

    Produced by ``_build_selection_result`` and passed to every handler.
    All lists are independent — iterating one does not require knowledge of
    any other.
    """

    sections: list[SectionSelection] = field(default_factory=list)
    sequences: list[SequenceSelection] = field(default_factory=list)
    standalone_proteins: list[StandaloneProteinSelection] = field(default_factory=list)
    protein_pairs: list[ProteinPairSelection] = field(default_factory=list)
    pair_member_proteins: list[PairMemberProteinSelection] = field(default_factory=list)
    scenes_headers: list[ScenesHeaderSelection] = field(default_factory=list)
    chains_headers: list[ChainsHeaderSelection] = field(default_factory=list)
    scenes: list[SceneSelection] = field(default_factory=list)
    chains: list[ChainSelection] = field(default_factory=list)
    residues: list[ResidueSelection] = field(default_factory=list)
    atoms: list[AtomSelection] = field(default_factory=list)

    # --- Boolean convenience properties ---

    @property
    def is_empty(self) -> bool:
        """True when no indexes were classified."""
        return self._total_count() == 0

    @property
    def is_single(self) -> bool:
        """True when exactly one node is selected across all buckets."""
        return self._total_count() == 1

    @property
    def is_uniform(self) -> bool:
        """True when all selected nodes belong to exactly one bucket type."""
        return sum(1 for lst in self._all_lists() if lst) == 1

    @property
    def has_sections(self) -> bool:
        """True when at least one section divider node is selected."""
        return bool(self.sections)

    @property
    def has_sequences(self) -> bool:
        """True when at least one sequence node is selected."""
        return bool(self.sequences)

    @property
    def has_proteins(self) -> bool:
        """True when at least one *standalone* protein node is selected."""
        return bool(self.standalone_proteins)

    @property
    def has_protein_pairs(self) -> bool:
        """True when at least one protein-pair node is selected."""
        return bool(self.protein_pairs)

    @property
    def has_pair_member_proteins(self) -> bool:
        """True when at least one pair-member protein node is selected."""
        return bool(self.pair_member_proteins)

    @property
    def has_scenes(self) -> bool:
        """True when at least one scene node is selected."""
        return bool(self.scenes)

    @property
    def has_chains(self) -> bool:
        """True when at least one chain node is selected."""
        return bool(self.chains)

    @property
    def has_residues(self) -> bool:
        """True when at least one residue node is selected."""
        return bool(self.residues)

    @property
    def has_atoms(self) -> bool:
        """True when at least one atom node is selected."""
        return bool(self.atoms)

    # --- Cross-bucket grouping helpers ---

    def residues_grouped_by_protein(self) -> dict[str, list[ResidueSelection]]:
        """Group selected residues by their owning protein's molecule-object name.

        Returns:
            ``{protein_name: [ResidueSelection, …], …}``
        """
        groups: dict[str, list[ResidueSelection]] = {}
        for res in self.residues:
            key = res.protein_object.get_molecule_object()
            groups.setdefault(key, []).append(res)
        return groups

    def atoms_grouped_by_protein(self) -> dict[str, list[AtomSelection]]:
        """Group selected atoms by their owning protein's molecule-object name.

        Returns:
            ``{protein_name: [AtomSelection, …], …}``
        """
        groups: dict[str, list[AtomSelection]] = {}
        for atom in self.atoms:
            key = atom.protein_object.get_molecule_object()
            groups.setdefault(key, []).append(atom)
        return groups

    def chains_grouped_by_protein(self) -> dict[str, list[ChainSelection]]:
        """Group selected chains by their owning protein's molecule-object name.

        Returns:
            ``{protein_name: [ChainSelection, …], …}``
        """
        groups: dict[str, list[ChainSelection]] = {}
        for chain in self.chains:
            key = chain.protein_object.get_molecule_object()
            groups.setdefault(key, []).append(chain)
        return groups

    def sub_structure_for_protein(
        self, a_protein_name: str,
    ) -> tuple[list[ChainSelection], list[ResidueSelection], list[AtomSelection]]:
        """Return all selected chains, residues, and atoms for one protein.

        Args:
            a_protein_name: The molecule-object name to filter by.

        Returns:
            ``(chains, residues, atoms)`` — each list may be empty.
        """
        chains = [
            c for c in self.chains
            if c.protein_object.get_molecule_object() == a_protein_name
        ]
        residues = [
            r for r in self.residues
            if r.protein_object.get_molecule_object() == a_protein_name
        ]
        atoms = [
            a for a in self.atoms
            if a.protein_object.get_molecule_object() == a_protein_name
        ]
        return chains, residues, atoms

    def all_selected_protein_names(self) -> set[str]:
        """Return the set of molecule-object names touched by the selection.

        Covers chains, residues, and atoms.
        """
        names: set[str] = set()
        for chain in self.chains:
            names.add(chain.protein_object.get_molecule_object())
        for res in self.residues:
            names.add(res.protein_object.get_molecule_object())
        for atom in self.atoms:
            names.add(atom.protein_object.get_molecule_object())
        return names

    # --- Internal helpers ---

    def _all_lists(self) -> list[list]:
        """Return every bucket list for iteration."""
        return [
            self.sections, self.sequences, self.standalone_proteins,
            self.protein_pairs, self.pair_member_proteins,
            self.scenes_headers, self.chains_headers,
            self.scenes, self.chains, self.residues, self.atoms,
        ]

    def _total_count(self) -> int:
        """Return total number of classified nodes across all buckets."""
        return sum(len(lst) for lst in self._all_lists())


# ===========================================================================
# Index classifier
# ===========================================================================

def _classify_index(index: QtCore.QModelIndex) -> Optional[object]:
    """Return the typed dataclass for *index*, or ``None`` for invalid indexes.

    This function is the single place that reads raw Qt roles and converts
    them into typed dataclasses. It inspects the ``TYPE_ROLE`` and ancestor
    chain of every node to determine which bucket it belongs to.

    Args:
        index: A ``QModelIndex`` from the ``PSAObjectsModel``.

    Returns:
        A typed dataclass instance, or ``None`` if the index is invalid.
    """
    if not index.isValid():
        return None

    node_type = index.data(enums.ModelEnum.TYPE_ROLE)

    if node_type == TYPE_SECTION:
        return SectionSelection(label=index.data(_DR), index=index)

    if node_type == TYPE_SEQUENCE:
        return SequenceSelection(
            name=index.data(_DR),
            sequence_string=index.data(enums.ModelEnum.OBJECT_ROLE),
            index=index,
        )

    if node_type == TYPE_PROTEIN:
        section = _ancestor_section_label(index)
        if section == LABEL_PROTEINS:
            return StandaloneProteinSelection(
                protein_name=index.data(_DR),
                protein_object=index.data(enums.ModelEnum.OBJECT_ROLE),
                index=index,
            )
        if section == LABEL_PROTEIN_PAIRS:
            return PairMemberProteinSelection(
                protein_name=index.data(_DR),
                protein_object=index.data(enums.ModelEnum.OBJECT_ROLE),
                pair_object=index.parent().data(enums.ModelEnum.OBJECT_ROLE),
                index=index,
            )

    if node_type == TYPE_PROTEIN_PAIR:
        return ProteinPairSelection(
            pair_name=index.data(_DR),
            pair_object=index.data(enums.ModelEnum.OBJECT_ROLE),
            index=index,
        )

    if node_type == TYPE_HEADER:
        header_label = index.data(_DR)
        model = index.model()

        if header_label == LABEL_SCENES:
            return ScenesHeaderSelection(
                owner_object=index.parent().data(enums.ModelEnum.OBJECT_ROLE),
                owner_type=index.parent().data(enums.ModelEnum.TYPE_ROLE),
                scene_count=model.rowCount(index),
                index=index,
            )

        if header_label == LABEL_CHAINS:
            parent_of_protein_type = index.parent().parent().data(
                enums.ModelEnum.TYPE_ROLE,
            )
            return ChainsHeaderSelection(
                protein_object=index.parent().data(enums.ModelEnum.OBJECT_ROLE),
                chain_count=model.rowCount(index),
                parent_of_protein_type=parent_of_protein_type,
                index=index,
            )

    if node_type == TYPE_SCENE:
        owner_index = index.parent().parent()
        return SceneSelection(
            scene_name=index.data(_DR),
            owner_object=owner_index.data(enums.ModelEnum.OBJECT_ROLE),
            owner_type=owner_index.data(enums.ModelEnum.TYPE_ROLE),
            index=index,
        )

    if node_type == TYPE_CHAIN:
        protein_index = index.parent().parent()
        return ChainSelection(
            chain_letter=index.data(_DR),
            chain_object=index.data(enums.ModelEnum.OBJECT_ROLE),
            chain_color=index.data(enums.ModelEnum.CHAIN_COLOR_ROLE),
            protein_object=protein_index.data(enums.ModelEnum.OBJECT_ROLE),
            parent_of_protein_type=protein_index.parent().data(
                enums.ModelEnum.TYPE_ROLE,
            ),
            index=index,
        )

    if node_type == TYPE_RESIDUE:
        chain_index = index.parent()
        protein_index = chain_index.parent().parent()
        return ResidueSelection(
            residue_label=index.data(_DR),
            residue_object=index.data(enums.ModelEnum.OBJECT_ROLE),
            chain_letter=chain_index.data(_DR),
            protein_object=protein_index.data(enums.ModelEnum.OBJECT_ROLE),
            parent_of_protein_type=protein_index.parent().data(
                enums.ModelEnum.TYPE_ROLE,
            ),
            index=index,
        )

    if node_type == TYPE_ATOM:
        residue_index = index.parent()
        chain_index = residue_index.parent()
        protein_index = chain_index.parent().parent()
        return AtomSelection(
            atom_name=index.data(_DR),
            residue_label=residue_index.data(_DR),
            residue_object=residue_index.data(enums.ModelEnum.OBJECT_ROLE),
            chain_letter=chain_index.data(_DR),
            protein_object=protein_index.data(enums.ModelEnum.OBJECT_ROLE),
            parent_of_protein_type=protein_index.parent().data(
                enums.ModelEnum.TYPE_ROLE,
            ),
            index=index,
        )

    return None


def _build_selection_result(indexes: list[QtCore.QModelIndex]) -> SelectionResult:
    """Classify every index in *indexes* and bucket it into a ``SelectionResult``.

    Args:
        indexes: List of ``QModelIndex`` objects from the selection model.

    Returns:
        A fully populated ``SelectionResult`` with one entry per classified index.
    """
    result = SelectionResult()
    _bucket_map = {
        SectionSelection: result.sections,
        SequenceSelection: result.sequences,
        StandaloneProteinSelection: result.standalone_proteins,
        ProteinPairSelection: result.protein_pairs,
        PairMemberProteinSelection: result.pair_member_proteins,
        ScenesHeaderSelection: result.scenes_headers,
        ChainsHeaderSelection: result.chains_headers,
        SceneSelection: result.scenes,
        ChainSelection: result.chains,
        ResidueSelection: result.residues,
        AtomSelection: result.atoms,
    }
    for index in indexes:
        classified = _classify_index(index)
        if classified is None:
            continue
        bucket = _bucket_map.get(type(classified))
        if bucket is not None:
            bucket.append(classified)
    return result


# ===========================================================================
# The slot
# ===========================================================================

def on_project_tree_selection_changed(
    self,
    selected: QtCore.QItemSelection,
    deselected: QtCore.QItemSelection,
) -> None:
    """Respond to a changed item selection in the project QTreeView.

    Works identically for single-item and multi-item (Ctrl/Shift) selections.
    The full current selection is always read from the view's selection model
    (not just the ``selected`` delta) so the result is correct when the user
    adds or removes items from an existing multi-selection.

    Args:
        self: The ``MainWindowController`` instance (passed explicitly because
            this function lives outside the class).
        selected: Newly selected items in this signal emission (delta).
        deselected: Items deselected in this emission.
    """
    tree_view = self._main_window.pyssa_objects_panel.tree_view
    all_indexes = tree_view.selectionModel().selectedIndexes()

    if not all_indexes:
        _handle_empty_selection(self)
        _update_ui_state(self, SelectionResult())
        return

    result = _build_selection_result(all_indexes)

    if result.has_sections:
        _handle_sections(self, result.sections, result)
    if result.has_sequences:
        _handle_sequences(self, result.sequences, result)
    if result.has_proteins:
        _handle_standalone_proteins(self, result.standalone_proteins, result)
    if result.has_protein_pairs:
        _handle_protein_pairs(self, result.protein_pairs, result)
    if result.has_pair_member_proteins:
        _handle_pair_member_proteins(self, result.pair_member_proteins, result)
    if result.scenes_headers:
        _handle_scenes_headers(self, result.scenes_headers, result)
    if result.chains_headers:
        _handle_chains_headers(self, result.chains_headers, result)
    if result.has_scenes:
        _handle_scenes(self, result.scenes, result)
    if result.has_chains:
        _handle_chains(self, result.chains, result)
    if result.has_residues:
        _handle_residues(self, result.residues, result)
    if result.has_atoms:
        _handle_atoms(self, result.atoms, result)

    _update_ui_state(self, result)


# ===========================================================================
# Double-click handler
# ===========================================================================

def on_tree_double_clicked(self, index: QtCore.QModelIndex) -> None:
    """Handle a double-click on a node in the project QTreeView.

    Sequences: opens the sequence viewer / dialog for the full amino acid
    string.
    Proteins / protein pairs: loads the corresponding PyMOL session from
    the domain object's base64-encoded ``pymol_session`` attribute.

    Args:
        self: The ``MainWindowController`` instance.
        index: The ``QModelIndex`` that was double-clicked.
    """
    if not index.isValid():
        return

    node_type = index.data(enums.ModelEnum.TYPE_ROLE)

    if node_type == TYPE_SEQUENCE:
        _open_sequence_viewer(self, index)
        return

    if node_type == TYPE_PROTEIN:
        protein_object = index.data(enums.ModelEnum.OBJECT_ROLE)
        if protein_object is not None:
            _load_pymol_session_from_object(self, protein_object)
        return

    if node_type == TYPE_PROTEIN_PAIR:
        pair_object = index.data(enums.ModelEnum.OBJECT_ROLE)
        if pair_object is not None:
            _load_pymol_session_from_object(self, pair_object)
        return


# ===========================================================================
# Handlers
# ===========================================================================

def _handle_empty_selection(self) -> None:
    """Clear the PyMOL ``sele`` selection when nothing is selected.

    Args:
        self: The ``MainWindowController`` instance.
    """
    try:
        self._user_pymol.get_cmd_module().select("sele", "none", enable=0)
    except Exception as e:
        logger.error(f"Failed to clear PyMOL selection: {e}")


def _handle_sections(
    self,
    sections: list[SectionSelection],
    result: SelectionResult,
) -> None:
    """Handle the selection of one or more section divider nodes.

    Section nodes are purely structural — no PyMOL selection is produced.
    UI state adjustments happen in ``_update_ui_state``.

    Args:
        self: The ``MainWindowController`` instance.
        sections: Non-empty list of selected section nodes.
        result: The complete ``SelectionResult`` for cross-bucket queries.
    """
    # Sections have no content action; UI state is handled centrally.


def _handle_sequences(
    self,
    sequences: list[SequenceSelection],
    result: SelectionResult,
) -> None:
    """Handle the selection of one or more sequence nodes.

    Sequences do not map to a PyMOL 3-D object, so no ``sele`` is created.
    Context menu and double-click behaviour are wired separately.

    Args:
        self: The ``MainWindowController`` instance.
        sequences: Non-empty list of selected sequence nodes.
        result: The complete ``SelectionResult`` for cross-bucket queries.
    """
    # Sequences have no PyMOL content action; UI state is handled centrally.


def _handle_standalone_proteins(
    self,
    proteins: list[StandaloneProteinSelection],
    result: SelectionResult,
) -> None:
    """Select the entire molecule object(s) in PyMOL.

    Args:
        self: The ``MainWindowController`` instance.
        proteins: Non-empty list of selected standalone protein nodes.
        result: The complete ``SelectionResult`` for cross-bucket queries.
    """
    parts = [p.protein_name for p in proteins]
    _apply_pymol_selection(self, " or ".join(parts))


def _handle_protein_pairs(
    self,
    pairs: list[ProteinPairSelection],
    result: SelectionResult,
) -> None:
    """Select both member molecule objects for each pair in PyMOL.

    Args:
        self: The ``MainWindowController`` instance.
        pairs: Non-empty list of selected protein-pair nodes.
        result: The complete ``SelectionResult`` for cross-bucket queries.
    """
    parts: list[str] = []
    for pair in pairs:
        pair_obj = pair.pair_object
        parts.append(pair_obj.protein_1.get_molecule_object())
        parts.append(pair_obj.protein_2.get_molecule_object())
    _apply_pymol_selection(self, " or ".join(parts))


def _handle_pair_member_proteins(
    self,
    proteins: list[PairMemberProteinSelection],
    result: SelectionResult,
) -> None:
    """Select the member molecule object(s) in PyMOL.

    Args:
        self: The ``MainWindowController`` instance.
        proteins: Non-empty list of selected pair-member protein nodes.
        result: The complete ``SelectionResult`` for cross-bucket queries.
    """
    parts = [p.protein_name for p in proteins]
    _apply_pymol_selection(self, " or ".join(parts))


def _handle_scenes_headers(
    self,
    headers: list[ScenesHeaderSelection],
    result: SelectionResult,
) -> None:
    """Handle the selection of one or more "Scenes" header nodes.

    Scenes headers are purely structural — no PyMOL action is produced.

    Args:
        self: The ``MainWindowController`` instance.
        headers: Non-empty list of selected scenes-header nodes.
        result: The complete ``SelectionResult`` for cross-bucket queries.
    """
    # No content action; UI state is handled centrally.


def _handle_chains_headers(
    self,
    headers: list[ChainsHeaderSelection],
    result: SelectionResult,
) -> None:
    """Handle the selection of one or more "Chains" header nodes.

    Chains headers are purely structural — no PyMOL action is produced.

    Args:
        self: The ``MainWindowController`` instance.
        headers: Non-empty list of selected chains-header nodes.
        result: The complete ``SelectionResult`` for cross-bucket queries.
    """
    # No content action; UI state is handled centrally.


def _handle_scenes(
    self,
    scenes: list[SceneSelection],
    result: SelectionResult,
) -> None:
    """Recall the selected PyMOL scene when exactly one is selected.

    When multiple scenes are selected the recall is skipped because PyMOL
    can only display one scene at a time.

    Args:
        self: The ``MainWindowController`` instance.
        scenes: Non-empty list of selected scene nodes.
        result: The complete ``SelectionResult`` for cross-bucket queries.
    """
    if len(scenes) == 1:
        try:
            self._user_pymol.get_cmd_module().scene(
                scenes[0].scene_name, "recall",
            )
        except Exception as e:
            logger.error(f"Failed to recall scene '{scenes[0].scene_name}': {e}")


def _handle_chains(
    self,
    chains: list[ChainSelection],
    result: SelectionResult,
) -> None:
    """Select the chain(s) in PyMOL.

    Args:
        self: The ``MainWindowController`` instance.
        chains: Non-empty list of selected chain nodes.
        result: The complete ``SelectionResult`` for cross-bucket queries.
    """
    model = self._main_window.pyssa_objects_panel.tree_view.model()
    parts = [model.construct_selection_string(c.index) for c in chains]
    _apply_pymol_selection(self, " or ".join(filter(None, parts)))


def _handle_residues(
    self,
    residues: list[ResidueSelection],
    result: SelectionResult,
) -> None:
    """Select the residue(s) in PyMOL.

    Args:
        self: The ``MainWindowController`` instance.
        residues: Non-empty list of selected residue nodes.
        result: The complete ``SelectionResult`` for cross-bucket queries.
    """
    model = self._main_window.pyssa_objects_panel.tree_view.model()
    parts = [model.construct_selection_string(r.index) for r in residues]
    _apply_pymol_selection(self, " or ".join(filter(None, parts)))


def _handle_atoms(
    self,
    atoms: list[AtomSelection],
    result: SelectionResult,
) -> None:
    """Select the atom(s) in PyMOL.

    Args:
        self: The ``MainWindowController`` instance.
        atoms: Non-empty list of selected atom nodes.
        result: The complete ``SelectionResult`` for cross-bucket queries.
    """
    model = self._main_window.pyssa_objects_panel.tree_view.model()
    parts = [model.construct_selection_string(a.index) for a in atoms]
    _apply_pymol_selection(self, " or ".join(filter(None, parts)))


# ===========================================================================
# UI State
# ===========================================================================

def _update_ui_state(self, result: SelectionResult) -> None:
    """Apply enabled/visible/checked state to every context-sensitive UI element.

    Called exactly once per selection change, after all content handlers have
    run, with the complete ``SelectionResult``.

    The strategy is:
      1. Call ``self.refresh_ui()`` which resets all widgets to the canonical
         project-level state (idempotent and declarative).
      2. Apply selection-specific overrides that ``refresh_ui`` cannot know
         about because they depend on *which* tree nodes are selected.

    Args:
        self: The ``MainWindowController`` instance.
        result: The fully classified current selection.
    """
    # Step 1: reset everything to canonical project-level state.
    self.refresh_ui()

    # Step 2: selection-specific overrides.
    mw = self._main_window

    # -- Prediction sub-menu overrides --
    if result.has_sequences and result.is_uniform:
        # All selected nodes are sequences — refine monomer vs. multimer.
        has_monomer_candidate = len(result.sequences) == 1
        has_multimer_candidate = len(result.sequences) >= 2
        mw.action_predict_monomer.setEnabled(has_monomer_candidate)
        mw.action_predict_multimer.setEnabled(has_multimer_candidate)
    elif not result.has_sequences:
        # No sequences in the selection — disable prediction entirely.
        mw.action_predict_monomer.setEnabled(False)
        mw.action_predict_multimer.setEnabled(False)


# ===========================================================================
# Private helpers
# ===========================================================================

def _ancestor_section_label(index: QtCore.QModelIndex) -> str:
    """Walk up the ancestor chain and return the label of the first section node.

    Used to distinguish ``TYPE_PROTEIN`` nodes that are standalone (under
    "Proteins") from those that are pair members (under "Protein Pairs").

    Args:
        index: A ``QModelIndex`` from the ``PSAObjectsModel``.

    Returns:
        The section label string, or an empty string if no section is found.
    """
    model = index.model()
    if model is None:
        return ""
    item = model.itemFromIndex(index)
    if item is None:
        return ""
    item = item.parent()
    while item is not None:
        if item.data(enums.ModelEnum.TYPE_ROLE) == TYPE_SECTION:
            return item.text()
        item = item.parent()
    return ""


def _apply_pymol_selection(self, selection_string: str) -> None:
    """Set the PyMOL ``sele`` named selection and enable its visual indicator.

    If ``selection_string`` is empty the selection is cleared instead.

    Args:
        self: The ``MainWindowController`` instance.
        selection_string: A valid PyMOL selection expression.
    """
    if not selection_string:
        try:
            self._user_pymol.get_cmd_module().select("sele", "none", enable=0)
        except Exception as e:
            logger.error(f"Failed to clear PyMOL selection: {e}")
        return
    try:
        self._user_pymol.get_cmd_module().select(
            "sele", selection_string, enable=1,
        )
    except Exception as e:
        logger.error(f"Failed to apply PyMOL selection '{selection_string}': {e}")


def _load_pymol_session_from_object(self, domain_object: object) -> None:
    """Decode the base64 PyMOL session from a ``Protein`` or ``ProteinPair`` and load it.

    Writes the decoded session to a temporary ``.pse`` file in
    ``constants.CACHE_PYMOL_SESSION_DIR`` and then loads it into the user
    PyMOL instance.

    Args:
        self: The ``MainWindowController`` instance.
        domain_object: A ``Protein`` or ``ProteinPair`` with a ``pymol_session``
            attribute (base64-encoded string).
    """
    pymol_session = getattr(domain_object, "pymol_session", None)
    if not pymol_session:
        logger.warning("Domain object has no pymol_session data; skipping load.")
        return

    # Derive a human-readable session name for the temp file.
    session_name = getattr(domain_object, "name", None)
    if session_name is None:
        session_name = getattr(domain_object, "_pymol_molecule_object", "temp")
        if session_name is None:
            molecule_getter = getattr(domain_object, "get_molecule_object", None)
            session_name = molecule_getter() if molecule_getter else "temp"

    session_path = pathlib.Path(
        f"{constants.CACHE_PYMOL_SESSION_DIR}/session_of_{session_name}.pse",
    )

    try:
        if not os.path.exists(constants.CACHE_PYMOL_SESSION_DIR):
            os.makedirs(constants.CACHE_PYMOL_SESSION_DIR, exist_ok=True)
        binary_data.write_binary_file_from_base64_string(session_path, pymol_session)
        self._user_pymol.get_cmd_module().reinitialize()
        self._user_pymol.get_cmd_module().load(str(session_path))
        logger.info(f"Loaded PyMOL session from '{session_path}'.")
    except Exception as e:
        logger.error(f"Failed to load PyMOL session for '{session_name}': {e}")


def _open_sequence_viewer(self, index: QtCore.QModelIndex) -> None:
    """Open a read-only dialog displaying the full amino acid sequence.

    Args:
        self: The ``MainWindowController`` instance.
        index: The ``QModelIndex`` of the double-clicked sequence node.
    """
    from src.pyssa.gui.qt import QtWidgets

    name = index.data(_DR)
    sequence_string = index.data(enums.ModelEnum.OBJECT_ROLE)
    if not sequence_string:
        return

    dialog = QtWidgets.QDialog(self._main_window)
    dialog.setWindowTitle(f"Sequence — {name}")
    dialog.setMinimumSize(500, 300)

    layout = QtWidgets.QVBoxLayout(dialog)

    text_edit = QtWidgets.QPlainTextEdit(dialog)
    text_edit.setReadOnly(True)
    text_edit.setPlainText(sequence_string)
    text_edit.setLineWrapMode(QtWidgets.QPlainTextEdit.LineWrapMode.WidgetWidth)
    layout.addWidget(text_edit)

    close_button = QtWidgets.QPushButton("Close", dialog)
    close_button.clicked.connect(dialog.accept)
    layout.addWidget(close_button)

    dialog.exec()
