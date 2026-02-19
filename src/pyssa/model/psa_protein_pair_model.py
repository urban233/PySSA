"""Module contains PySSAProteinPairsModel — the protein-pairs tree model."""
import logging
from typing import Optional

import zmq

from src.pyssa.gui.qt import QtGui
from src.pyssa.gui.qt import QtCore

from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.model import base_tree_model
from src.pyssa.model.protein_subtree_mixin import (
  ProteinSubtreeMixin,
  TYPE_PROTEIN,
  TYPE_PROTEIN_PAIR,
  TYPE_CHAIN,
  TYPE_HEADER,
  TYPE_SCENE,
  LABEL_SCENES,
  LABEL_CHAINS,
)
from src.pyssa.internal.data_structures import protein, protein_pair
from src.pyssa.util import enums, exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)

__docformat__ = "google"


class PSAProteinPairModel(ProteinSubtreeMixin, base_tree_model.BaseTreeModel):
  """Tree model that stores protein pairs from pair level down to atom level.

  Tree structure
  --------------
  (invisible root)
  └── <pair name>                    [TYPE_PROTEIN_PAIR]
      ├── Scenes                     [TYPE_HEADER]
      │   └── <scene name> …         [TYPE_SCENE]
      ├── <protein 1 name>           [TYPE_PROTEIN]
      │   └── Chains                 [TYPE_HEADER]
      │       └── <chain id>         [TYPE_CHAIN]
      │           └── <resi>-<resn>  [TYPE_RESIDUE]
      │               └── <atom>     [TYPE_ATOM]
      └── <protein 2 name>           [TYPE_PROTEIN]
          └── Chains                 [TYPE_HEADER]
              └── …
  """

  def __init__(self) -> None:
    """Constructor."""
    super().__init__()
    self.create_root_node()

  # ------------------------------------------------------------------
  # Public API — building the model
  # ------------------------------------------------------------------

  def build_model_from_protein_pair_objects(
          self,
          the_protein_pair_objects: list["protein_pair.ProteinPair"],
          the_main_socket: zmq.Socket,
          a_socket: zmq.Socket,
  ) -> None:
    """Build the model from a list of ProteinPair objects.

    Scenes and the full residue/atom hierarchy for both proteins in each
    pair are fetched from auxiliary PyMOL.

    Args:
        the_protein_pair_objects: Protein pairs to add to the model.
        the_main_socket: Main ZMQ socket for auxiliary PyMOL.
        a_socket: Secondary ZMQ socket for auxiliary PyMOL.

    Raises:
        exception.IllegalArgumentError: If any argument is ``None``.
    """
    if the_protein_pair_objects is None:
      logger.error("the_protein_pair_objects is None.")
      raise exception.IllegalArgumentError("the_protein_pair_objects is None.")
    if the_main_socket is None:
      logger.error("the_main_socket is None.")
      raise exception.IllegalArgumentError("the_main_socket is None.")
    if a_socket is None:
      logger.error("a_socket is None.")
      raise exception.IllegalArgumentError("a_socket is None.")

    for tmp_pair in the_protein_pair_objects:
      scenes = self._fetch_scenes(
        str(tmp_pair.pymol_session), the_main_socket, a_socket
      )
      self._add_protein_pair_node(tmp_pair, scenes, the_main_socket, a_socket)

  def add_protein_pair(
          self,
          a_protein_pair: "protein_pair.ProteinPair",
          the_main_socket: zmq.Socket,
          a_socket: zmq.Socket,
  ) -> None:
    """Add a single ProteinPair to the model with full atom-level hierarchy.

    Args:
        a_protein_pair: The protein pair to add.
        the_main_socket: Main ZMQ socket for auxiliary PyMOL.
        a_socket: Secondary ZMQ socket for auxiliary PyMOL.

    Raises:
        exception.IllegalArgumentError: If any argument is ``None``.
    """
    if a_protein_pair is None:
      logger.error("a_protein_pair is None.")
      raise exception.IllegalArgumentError("a_protein_pair is None.")
    if the_main_socket is None:
      logger.error("the_main_socket is None.")
      raise exception.IllegalArgumentError("the_main_socket is None.")
    if a_socket is None:
      logger.error("a_socket is None.")
      raise exception.IllegalArgumentError("a_socket is None.")

    scenes = self._fetch_scenes(
      str(a_protein_pair.pymol_session), the_main_socket, a_socket
    )
    self._add_protein_pair_node(a_protein_pair, scenes, the_main_socket, a_socket)

  # ------------------------------------------------------------------
  # Public API — scene management
  # ------------------------------------------------------------------

  def add_scene(
          self,
          a_model_index: QtCore.QModelIndex,
          the_scene_item: QtGui.QStandardItem,
  ) -> None:
    """Append *the_scene_item* to the Scenes header of the pair at *a_model_index*.

    *a_model_index* may point to any node within the target pair's subtree.

    Args:
        a_model_index: Any index within a protein-pair subtree.
        the_scene_item: A pre-built ``QStandardItem`` for the scene.

    Raises:
        exception.IllegalArgumentError: If any argument is ``None``.
        ValueError: If the node type cannot be resolved.
    """
    if a_model_index is None:
      logger.error("a_model_index is None.")
      raise exception.IllegalArgumentError("a_model_index is None.")
    if the_scene_item is None:
      logger.error("the_scene_item is None.")
      raise exception.IllegalArgumentError("the_scene_item is None.")

    scenes_header = self._resolve_scenes_header(a_model_index)
    scenes_header.appendRow(the_scene_item)

  def remove_scene(self, the_model_index_of_the_scene: QtCore.QModelIndex) -> None:
    """Remove the scene node at *the_model_index_of_the_scene*.

    Args:
        the_model_index_of_the_scene: Index of a scene node.

    Raises:
        exception.IllegalArgumentError: If the argument is ``None``.
        ValueError: If the node is not a scene node.
    """
    if the_model_index_of_the_scene is None:
      logger.error("the_model_index_of_the_scene is None.")
      raise exception.IllegalArgumentError("the_model_index_of_the_scene is None.")
    if self.data(the_model_index_of_the_scene, enums.ModelEnum.TYPE_ROLE) != TYPE_SCENE:
      raise ValueError("The provided index does not point to a scene node.")

    scene_item = self.itemFromIndex(the_model_index_of_the_scene)
    scene_item.parent().removeRow(scene_item.row())

  def check_if_scratch_scene_exists(self, a_model_index: QtCore.QModelIndex) -> bool:
    """Return ``True`` if a ``_scratch_`` scene exists for the pair.

    *a_model_index* may point to any node within the pair's subtree.

    Args:
        a_model_index: Any index within a protein-pair subtree.

    Raises:
        exception.IllegalArgumentError: If ``a_model_index`` is ``None``.
        ValueError: If the node type cannot be resolved.
    """
    if a_model_index is None:
      logger.error("a_model_index is None.")
      raise exception.IllegalArgumentError("a_model_index is None.")

    scenes_header = self._resolve_scenes_header(a_model_index)
    for row in range(scenes_header.rowCount()):
      if scenes_header.child(row, 0).data(QtCore.Qt.ItemDataRole.DisplayRole) == "_scratch_":
        return True
    return False

  # ------------------------------------------------------------------
  # Public API — pair management
  # ------------------------------------------------------------------

  def remove_protein_pair(
          self, the_model_index_of_the_pair: QtCore.QModelIndex
  ) -> None:
    """Remove the protein-pair subtree at *the_model_index_of_the_pair*.

    Args:
        the_model_index_of_the_pair: Index of a protein-pair node.

    Raises:
        exception.IllegalArgumentError: If the argument is ``None``.
        ValueError: If the node is not a protein-pair node.
    """
    if the_model_index_of_the_pair is None:
      logger.error("the_model_index_of_the_pair is None.")
      raise exception.IllegalArgumentError("the_model_index_of_the_pair is None.")
    if self.data(the_model_index_of_the_pair, enums.ModelEnum.TYPE_ROLE) != TYPE_PROTEIN_PAIR:
      raise ValueError("The provided index does not point to a protein-pair node.")

    pair_item = self.itemFromIndex(the_model_index_of_the_pair)
    self.removeRow(pair_item.row())

  # ------------------------------------------------------------------
  # Private helpers
  # ------------------------------------------------------------------

  def _add_protein_pair_node(
          self,
          a_protein_pair: "protein_pair.ProteinPair",
          scenes: list[str],
          the_main_socket: zmq.Socket,
          a_socket: zmq.Socket,
  ) -> QtGui.QStandardItem:
    """Create and append a full protein-pair subtree under the root node.

    Args:
        a_protein_pair: The source ``ProteinPair`` object.
        scenes: Scene names for the pair's Scenes header.
        the_main_socket: Main ZMQ socket for auxiliary PyMOL.
        a_socket: Secondary ZMQ socket for auxiliary PyMOL.

    Returns:
        The newly created protein-pair ``QStandardItem``.
    """
    pair_node = self.add_node(
      a_parent_node=self.root_node,
      an_item_name=a_protein_pair.name,
      an_item_type_value=TYPE_PROTEIN_PAIR,
      an_item_object_value=a_protein_pair,
    )

    # Scenes belong to the pair (shared session), not to individual proteins
    scenes_header = self._append_header_node(pair_node, LABEL_SCENES)
    self._append_scenes_from_list(scenes_header, scenes)

    # Each protein in the pair gets its own subtree (chains → residues → atoms)
    self._add_protein_child_node(
      pair_node, a_protein_pair.protein_1, the_main_socket, a_socket
    )
    self._add_protein_child_node(
      pair_node, a_protein_pair.protein_2, the_main_socket, a_socket
    )

    return pair_node

  def _add_protein_child_node(
          self,
          pair_node: QtGui.QStandardItem,
          a_protein: "protein.Protein",
          the_main_socket: zmq.Socket,
          a_socket: zmq.Socket,
  ) -> QtGui.QStandardItem:
    """Append a protein node (with full chain hierarchy) to *pair_node*.

    The chempy model is fetched from auxiliary PyMOL so the hierarchy is
    populated down to the atom level.  Chain colour and object data from
    *a_protein* are merged into each chain node.

    Args:
        pair_node: The parent protein-pair ``QStandardItem``.
        a_protein: The protein to add as a child.
        the_main_socket: Main ZMQ socket for auxiliary PyMOL.
        a_socket: Secondary ZMQ socket for auxiliary PyMOL.

    Returns:
        The newly created protein ``QStandardItem``.
    """
    protein_node = self.add_node(
      a_parent_node=pair_node,
      an_item_name=a_protein.get_molecule_object(),
      an_item_type_value=TYPE_PROTEIN,
      an_item_object_value=a_protein,
    )

    chempy_model = self._fetch_chempy_model_for_protein(a_protein, the_main_socket, a_socket)
    hierarchy_map = self._build_hierarchy_map(chempy_model)
    chains_for_this_protein = hierarchy_map.get(a_protein.get_molecule_object(), {})

    chains_header = self._append_header_node(protein_node, LABEL_CHAINS)
    chain_object_lookup = {chain.chain_letter: chain for chain in a_protein.chains}
    self._populate_chain_hierarchy(
      chains_header, chains_for_this_protein, chain_object_lookup=chain_object_lookup
    )

    return protein_node

  def _resolve_scenes_header(
          self, a_model_index: QtCore.QModelIndex
  ) -> QtGui.QStandardItem:
    """Return the Scenes ``QStandardItem`` for the pair that owns *a_model_index*.

    The Scenes header always belongs to the pair node, not to either
    individual protein node.

    Raises:
        ValueError: If the node type is not recognised.
    """
    node_type = a_model_index.data(enums.ModelEnum.TYPE_ROLE)

    if node_type == TYPE_PROTEIN_PAIR:
      pair_item = self.itemFromIndex(a_model_index)
    elif node_type == TYPE_SCENE:
      # scene → Scenes header → pair
      pair_item = self.itemFromIndex(a_model_index).parent().parent()
    elif node_type == TYPE_PROTEIN:
      # protein → pair
      pair_item = self.itemFromIndex(a_model_index).parent()
    elif node_type == TYPE_HEADER:
      # Either the Scenes header (parent = pair) or the Chains header (parent = protein → pair)
      header_item = self.itemFromIndex(a_model_index)
      if header_item.data(QtCore.Qt.ItemDataRole.DisplayRole) == LABEL_SCENES:
        pair_item = header_item.parent()
      else:
        # Chains header: parent = protein, grandparent = pair
        pair_item = header_item.parent().parent()
    elif node_type == TYPE_CHAIN:
      # chain → Chains header → protein → pair
      pair_item = self.itemFromIndex(a_model_index).parent().parent().parent()
    else:
      # residue or atom: walk up until we hit the pair node
      item = self.itemFromIndex(a_model_index)
      while item is not None:
        if item.data(enums.ModelEnum.TYPE_ROLE) == TYPE_PROTEIN_PAIR:
          pair_item = item
          break
        item = item.parent()
      else:
        raise ValueError(f"Could not resolve a pair node from type: '{node_type}'")

    # Scenes is always the first child of a protein-pair node
    return pair_item.child(0, 0)
