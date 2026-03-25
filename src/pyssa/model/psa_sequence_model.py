"""Module contains PySSASequenceModel — the sequences tree model."""
import logging

from src.pyssa.gui.qt import QtCore

from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.model import base_tree_model
from src.pyssa.model.protein_subtree_mixin import TYPE_SEQUENCE, LABEL_SEQUENCES
from src.pyssa.util import enums, exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)

__docformat__ = "google"


class PSASequenceModel(base_tree_model.BaseTreeModel):
  """Tree model that stores amino acid sequences as a flat list.

  Each sequence is one top-level node whose display text is the sequence
  name and whose ``OBJECT_ROLE`` holds the raw amino acid sequence string.
  This replaces the legacy table-based ``SequenceModel`` with a structure
  that is consistent with the rest of the project tree.

  Tree structure
  --------------
  (invisible root)
  ├── <sequence name>    [TYPE_SEQUENCE]  OBJECT_ROLE = "<amino acid string>"
  ├── <sequence name>    [TYPE_SEQUENCE]
  └── …
  """

  def __init__(self) -> None:
    """Constructor."""
    super().__init__()
    self.create_root_node()

  # ------------------------------------------------------------------
  # Public API — adding sequences
  # ------------------------------------------------------------------

  def add_sequence(self, a_sequence: str) -> None:
    """Add a sequence to the model, using the sequence string itself as the display name.

    This preserves the interface of the legacy ``SequenceModel``.  When a
    human-readable name is available prefer :meth:`add_named_sequence`.

    Args:
        a_sequence: The amino acid sequence string to add.

    Raises:
        exception.IllegalArgumentError: If ``a_sequence`` is ``None`` or empty.
    """
    if not a_sequence:
      logger.error("a_sequence is either None or an empty string.")
      raise exception.IllegalArgumentError(
        "a_sequence is either None or an empty string."
      )
    self._add_sequence_node(display_name=a_sequence, sequence=a_sequence)

  def add_named_sequence(self, a_name: str, a_sequence: str) -> None:
    """Add a sequence with a separate display name and sequence string.

    Args:
        a_name: The human-readable name shown in the tree view.
        a_sequence: The amino acid sequence string stored in ``OBJECT_ROLE``.

    Raises:
        exception.IllegalArgumentError: If ``a_name`` or ``a_sequence`` is
            ``None`` or empty.
    """
    if not a_name:
      logger.error("a_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
        "a_name is either None or an empty string."
      )
    if not a_sequence:
      logger.error("a_sequence is either None or an empty string.")
      raise exception.IllegalArgumentError(
        "a_sequence is either None or an empty string."
      )
    self._add_sequence_node(display_name=a_name, sequence=a_sequence)

  # ------------------------------------------------------------------
  # Public API — removing sequences
  # ------------------------------------------------------------------

  def remove_sequence(self, a_model_index: QtCore.QModelIndex) -> None:
    """Remove the sequence node at *a_model_index*.

    Args:
        a_model_index: Index of a sequence node.

    Raises:
        exception.IllegalArgumentError: If ``a_model_index`` is ``None``.
        ValueError: If the node at the index is not a sequence node.
    """
    if a_model_index is None:
      logger.error("a_model_index is None.")
      raise exception.IllegalArgumentError("a_model_index is None.")
    if a_model_index.data(enums.ModelEnum.TYPE_ROLE) != TYPE_SEQUENCE:
      raise ValueError("The provided index does not point to a sequence node.")

    sequence_item = self.itemFromIndex(a_model_index)
    self.removeRow(sequence_item.row())

  # ------------------------------------------------------------------
  # Public API — querying
  # ------------------------------------------------------------------

  def get_sequence_string(self, a_model_index: QtCore.QModelIndex) -> str:
    """Return the amino acid sequence string stored at *a_model_index*.

    Args:
        a_model_index: Index of a sequence node.

    Raises:
        exception.IllegalArgumentError: If ``a_model_index`` is ``None``.
        ValueError: If the node at the index is not a sequence node.
    """
    if a_model_index is None:
      logger.error("a_model_index is None.")
      raise exception.IllegalArgumentError("a_model_index is None.")
    if a_model_index.data(enums.ModelEnum.TYPE_ROLE) != TYPE_SEQUENCE:
      raise ValueError("The provided index does not point to a sequence node.")

    return a_model_index.data(enums.ModelEnum.OBJECT_ROLE)

  # ------------------------------------------------------------------
  # Private helpers
  # ------------------------------------------------------------------

  def _add_sequence_node(self, display_name: str, sequence: str) -> None:
    """Create and append a single sequence node to the root."""
    self.add_node(
      a_parent_node=self.root_node,
      an_item_name=display_name,
      an_item_type_value=TYPE_SEQUENCE,
      an_item_object_value=sequence,
    )
    logger.debug("Sequence '%s' added to the model.", display_name)
