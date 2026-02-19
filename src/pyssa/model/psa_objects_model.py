"""Module contains PySSAProjectModel — the single model for the project QTreeView.

This model combines sequences, standalone proteins, and protein pairs into one
``QStandardItemModel`` so that a single ``QTreeView`` can display all project
data.  Three permanent top-level section nodes act as visual dividers and
define the fixed display order.

Tree structure
--------------
(invisible root)
├── Sequences                          [TYPE_SECTION]
│   └── <sequence name>                [TYPE_SEQUENCE]
├── Proteins                           [TYPE_SECTION]
│   └── <protein name>                 [TYPE_PROTEIN]
│       ├── Scenes                     [TYPE_HEADER]
│       │   └── <scene name> …         [TYPE_SCENE]
│       └── Chains                     [TYPE_HEADER]
│           └── <chain id>             [TYPE_CHAIN]
│               └── <resi> - <resn>    [TYPE_RESIDUE]
│                   └── <atom name>    [TYPE_ATOM]
└── Protein Pairs                      [TYPE_SECTION]
    └── <pair name>                    [TYPE_PROTEIN_PAIR]
        ├── Scenes                     [TYPE_HEADER]
        │   └── <scene name> …         [TYPE_SCENE]
        ├── <protein 1 name>           [TYPE_PROTEIN]
        │   └── Chains                 [TYPE_HEADER]
        │       └── …
        └── <protein 2 name>           [TYPE_PROTEIN]
            └── Chains                 [TYPE_HEADER]
                └── …

Design rationale
----------------
Section nodes are created in the constructor in the required display order
(Sequences → Proteins → Protein Pairs) so the ordering is always correct
regardless of when or in what order data is added at runtime.

All per-domain logic is delegated to the three focused sub-models, keeping
this class small and easy to follow.  The ``_with_root`` helper temporarily
redirects each sub-model's root node to the matching section node so its
``add_node`` calls build inside the correct section of this model's tree
without any changes being required in the sub-model classes.
"""
import logging

import zmq
from Bio.SeqRecord import SeqRecord

from src.pyssa.gui.qt import QtGui
from src.pyssa.gui.qt import QtCore

from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.model import base_tree_model
from src.pyssa.model.protein_subtree_mixin import (
  TYPE_PROTEIN,
  TYPE_PROTEIN_PAIR,
  TYPE_SEQUENCE,
  TYPE_SECTION,
  LABEL_SEQUENCES,
  LABEL_PROTEINS,
  LABEL_PROTEIN_PAIRS,
)
from src.pyssa.model import psa_sequence_model
from src.pyssa.model import psa_protein_model
from src.pyssa.model import psa_protein_pair_model
from src.pyssa.internal.data_structures import protein, protein_pair
from src.pyssa.util import enums, exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)

__docformat__ = "google"


class PSAObjectsModel(base_tree_model.BaseTreeModel):
  """Single tree model for displaying sequences, proteins, and protein pairs in one QTreeView.

  The model owns three permanent section nodes created in the order
  Sequences → Proteins → Protein Pairs, and delegates all subtree-building
  and management operations to the three specialised sub-models.

  Typical usage
  -------------
  ::

      model = PySSAProjectModel()
      model.build_model(sequences, proteins, protein_pairs, main_socket, socket)
      tree_view.setModel(model)
  """

  def __init__(self) -> None:
    """Constructor.

    Creates the three permanent section nodes in display order and
    initialises the sub-models.
    """
    super().__init__()
    self.create_root_node()

    # Section nodes are appended in the required display order.
    self._sequences_section: QtGui.QStandardItem = self._append_section_node(LABEL_SEQUENCES)
    self._proteins_section: QtGui.QStandardItem = self._append_section_node(LABEL_PROTEINS)
    self._protein_pairs_section: QtGui.QStandardItem = self._append_section_node(LABEL_PROTEIN_PAIRS)

    # Focused sub-models that own the building and management logic.
    # Their root_node attributes are temporarily redirected to the matching
    # section node when a build/add call is made (see _with_root).
    self._sequences_model = psa_sequence_model.PSASequenceModel()
    self._proteins_model = psa_protein_model.PSAProteinModel()
    self._protein_pairs_model = psa_protein_pair_model.PSAProteinPairModel()

  # ------------------------------------------------------------------
  # Public API — building the model
  # ------------------------------------------------------------------

  def build_model(
          self,
          the_sequences: list[SeqRecord],
          the_protein_objects: list["protein.Protein"],
          the_protein_pair_objects: list["protein_pair.ProteinPair"]
  ) -> None:
    """Populate the model from all project data in one call.

    Args:
        the_sequences: Amino acid sequence SeqRecords to add under "Sequences".
        the_protein_objects: Standalone proteins to add under "Proteins".
        the_protein_pair_objects: Protein pairs to add under "Protein Pairs".

    Raises:
        exception.IllegalArgumentError: If any argument is ``None``.
    """
    if the_sequences is None:
      logger.error("the_sequences is None.")
      raise exception.IllegalArgumentError("the_sequences is None.")
    if the_protein_objects is None:
      logger.error("the_protein_objects is None.")
      raise exception.IllegalArgumentError("the_protein_objects is None.")
    if the_protein_pair_objects is None:
      logger.error("the_protein_pair_objects is None.")
      raise exception.IllegalArgumentError("the_protein_pair_objects is None.")

    for sequence in the_sequences:
      self.add_sequence(sequence.name)

    for tmp_protein in the_protein_objects:
      self.add_protein(tmp_protein)

    for tmp_pair in the_protein_pair_objects:
      self.add_protein_pair(tmp_pair)

  # ------------------------------------------------------------------
  # Public API — adding items
  # ------------------------------------------------------------------

  def add_sequence(self, a_sequence: str) -> None:
    """Add an amino acid sequence under the "Sequences" section.

    The sequence string is used as both the display name and the stored
    value.  Use :meth:`add_named_sequence` when a name is
    available.

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
    self._with_root(
      self._sequences_model,
      self._sequences_section,
      lambda: self._sequences_model.add_sequence(a_sequence),
    )

  def add_named_sequence(self, a_name: str, a_sequence: str) -> None:
    """Add a sequence with a separate display name under the "Sequences" section.

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
    self._with_root(
      self._sequences_model,
      self._sequences_section,
      lambda: self._sequences_model.add_named_sequence(a_name, a_sequence),
    )

  def add_protein(
          self,
          a_protein: "protein.Protein"
  ) -> None:
    """Add a standalone protein under the "Proteins" section.

    Args:
        a_protein: The protein to add.

    Raises:
        exception.IllegalArgumentError: If any argument is ``None``.
    """
    if a_protein is None:
      logger.error("a_protein is None.")
      raise exception.IllegalArgumentError("a_protein is None.")

    self._with_root(
      self._proteins_model,
      self._proteins_section,
      lambda: self._proteins_model.add_protein_from_protein_object(a_protein),
    )

  def add_protein_pair(
          self,
          a_protein_pair: "protein_pair.ProteinPair",
  ) -> None:
    """Add a protein pair under the "Protein Pairs" section.

    Args:
        a_protein_pair: The protein pair to add.

    Raises:
        exception.IllegalArgumentError: If any argument is ``None``.
    """
    if a_protein_pair is None:
      logger.error("a_protein_pair is None.")
      raise exception.IllegalArgumentError("a_protein_pair is None.")

    self._with_root(
      self._protein_pairs_model,
      self._protein_pairs_section,
      lambda: self._protein_pairs_model.add_protein_pair(a_protein_pair),
    )

  # ------------------------------------------------------------------
  # Public API — removing items
  # ------------------------------------------------------------------

  def remove_sequence(self, a_model_index: QtCore.QModelIndex) -> None:
    """Remove the sequence at *a_model_index* from the "Sequences" section.

    Args:
        a_model_index: Index of a sequence node.

    Raises:
        exception.IllegalArgumentError: If the argument is ``None``.
        ValueError: If the node is not a sequence node.
    """
    if a_model_index.data(enums.ModelEnum.TYPE_ROLE) != TYPE_SEQUENCE:
      raise ValueError("The provided index does not point to a sequence node.")

    sequence_item = self.itemFromIndex(a_model_index)
    sequence_item.parent().removeRow(sequence_item.row())

  def remove_protein(self, a_model_index: QtCore.QModelIndex) -> None:
    """Remove the protein at *a_model_index* from the "Proteins" section.

    Args:
        a_model_index: Index of a protein node.

    Raises:
        exception.IllegalArgumentError: If the argument is ``None``.
        ValueError: If the node is not a protein node.
    """
    if a_model_index.data(enums.ModelEnum.TYPE_ROLE) != TYPE_PROTEIN:
      raise ValueError("The provided index does not point to a protein node.")

    protein_item = self.itemFromIndex(a_model_index)
    protein_item.parent().removeRow(protein_item.row())

  def remove_protein_pair(self, a_model_index: QtCore.QModelIndex) -> None:
    """Remove the protein pair at *a_model_index* from the "Protein Pairs" section.

    Args:
        a_model_index: Index of a protein-pair node.

    Raises:
        exception.IllegalArgumentError: If the argument is ``None``.
        ValueError: If the node is not a protein-pair node.
    """
    if a_model_index.data(enums.ModelEnum.TYPE_ROLE) != TYPE_PROTEIN_PAIR:
      raise ValueError("The provided index does not point to a protein-pair node.")

    pair_item = self.itemFromIndex(a_model_index)
    pair_item.parent().removeRow(pair_item.row())

  # ------------------------------------------------------------------
  # Public API — scene management (delegates to protein/pair sub-models)
  # ------------------------------------------------------------------

  def add_scene(
          self,
          a_model_index: QtCore.QModelIndex,
          the_scene_item: QtGui.QStandardItem,
  ) -> None:
    """Add *the_scene_item* to the correct Scenes header for *a_model_index*.

    The section (protein vs protein-pair) is determined automatically.
    Sequence nodes have no scenes and will raise ``ValueError``.

    Args:
        a_model_index: Any index within a protein or protein-pair subtree.
        the_scene_item: A pre-built ``QStandardItem`` for the scene.

    Raises:
        exception.IllegalArgumentError: If any argument is ``None``.
        ValueError: If the section cannot be determined or is "Sequences".
    """
    sub_model = self._scene_capable_sub_model_for_index(a_model_index)
    sub_model.add_scene(a_model_index, the_scene_item)

  def remove_scene(self, the_model_index_of_the_scene: QtCore.QModelIndex) -> None:
    """Remove the scene at *the_model_index_of_the_scene*.

    Args:
        the_model_index_of_the_scene: Index of a scene node.

    Raises:
        exception.IllegalArgumentError: If the argument is ``None``.
        ValueError: If the node is not a scene node or is in the Sequences section.
    """
    sub_model = self._scene_capable_sub_model_for_index(the_model_index_of_the_scene)
    sub_model.remove_scene(the_model_index_of_the_scene)

  def check_if_scratch_scene_exists(self, a_model_index: QtCore.QModelIndex) -> bool:
    """Return ``True`` if a ``_scratch_`` scene exists in the relevant subtree.

    Args:
        a_model_index: Any index within a protein or protein-pair subtree.

    Raises:
        exception.IllegalArgumentError: If ``a_model_index`` is ``None``.
        ValueError: If the section cannot be determined or is "Sequences".
    """
    sub_model = self._scene_capable_sub_model_for_index(a_model_index)
    return sub_model.check_if_scratch_scene_exists(a_model_index)

  # ------------------------------------------------------------------
  # Public API — sequence data access
  # ------------------------------------------------------------------

  def get_sequence_string(self, a_model_index: QtCore.QModelIndex) -> str:
    """Return the amino acid sequence string stored at *a_model_index*.

    Args:
        a_model_index: Index of a sequence node.

    Raises:
        exception.IllegalArgumentError: If ``a_model_index`` is ``None``.
        ValueError: If the node is not a sequence node.
    """
    return self._sequences_model.get_sequence_string(a_model_index)

  # ------------------------------------------------------------------
  # Public API — section index accessors
  # ------------------------------------------------------------------

  def get_sequences_section_index(self) -> QtCore.QModelIndex:
    """Return the model index of the "Sequences" section node."""
    return self.indexFromItem(self._sequences_section)

  def get_proteins_section_index(self) -> QtCore.QModelIndex:
    """Return the model index of the "Proteins" section node."""
    return self.indexFromItem(self._proteins_section)

  def get_protein_pairs_section_index(self) -> QtCore.QModelIndex:
    """Return the model index of the "Protein Pairs" section node."""
    return self.indexFromItem(self._protein_pairs_section)

  # ------------------------------------------------------------------
  # Private helpers
  # ------------------------------------------------------------------

  def _append_section_node(self, a_label: str) -> QtGui.QStandardItem:
    """Append a permanent, non-selectable top-level section divider and return it."""
    section_item = QtGui.QStandardItem(a_label)
    section_item.setData(TYPE_SECTION, enums.ModelEnum.TYPE_ROLE)
    section_item.setFlags(
      section_item.flags() & ~QtCore.Qt.ItemFlag.ItemIsSelectable
    )
    self.root_node.appendRow(section_item)
    return section_item

  def _with_root(
          self,
          a_sub_model: base_tree_model.BaseTreeModel,
          a_section_node: QtGui.QStandardItem,
          a_callable,
  ) -> None:
    """Run *a_callable* with *a_sub_model*'s root temporarily set to *a_section_node*.

    This lets the sub-model's ``add_node`` calls build directly inside the
    correct section of this model's item tree, with no changes needed in
    the sub-model classes.  The original root is always restored, even if
    *a_callable* raises.

    Args:
        a_sub_model: The sub-model whose root to redirect.
        a_section_node: The section node to use as the temporary root.
        a_callable: Zero-argument callable to invoke with the root redirected.
    """
    original_root = a_sub_model.root_node
    a_sub_model.root_node = a_section_node
    try:
      a_callable()
    finally:
      a_sub_model.root_node = original_root

  def _scene_capable_sub_model_for_index(
          self, a_model_index: QtCore.QModelIndex
  ) -> "psa_protein_model.PSAProteinModel | psa_protein_pair_model.PSAProteinPairModel":
    """Return the protein or protein-pair sub-model that owns *a_model_index*.

    Sequences have no scenes, so this raises ``ValueError`` if the index
    belongs to the Sequences section.

    Raises:
        ValueError: If the section is "Sequences" or cannot be found in
            the ancestor chain.
    """
    item = self.itemFromIndex(a_model_index)
    while item is not None:
      if item.data(enums.ModelEnum.TYPE_ROLE) == TYPE_SECTION:
        label = item.text()
        if label == LABEL_SEQUENCES:
          raise ValueError(
            "Scene operations are not supported for the Sequences section."
          )
        if label == LABEL_PROTEINS:
          return self._proteins_model
        if label == LABEL_PROTEIN_PAIRS:
          return self._protein_pairs_model
      item = item.parent()
    raise ValueError("Could not determine the section for the given index.")
