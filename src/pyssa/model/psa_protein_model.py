"""Module contains PySSAObjectsModel — the standalone-proteins tree model."""
import logging
from typing import Optional

import zmq
from chempy.models import Indexed

from src.pyssa.gui.qt import QtGui
from src.pyssa.gui.qt import QtCore

from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.model import base_tree_model
from src.pyssa.model.protein_subtree_mixin import (
  ProteinSubtreeMixin,
  TYPE_PROTEIN,
  TYPE_CHAIN,
  TYPE_RESIDUE,
  TYPE_ATOM,
  TYPE_HEADER,
  TYPE_SCENE,
  LABEL_SCENES,
  LABEL_CHAINS,
)
from src.pyssa.internal.data_structures import protein
from src.pyssa.util import enums, exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)

__docformat__ = "google"


class PSAProteinModel(ProteinSubtreeMixin, base_tree_model.BaseTreeModel):
  """Tree model that stores standalone proteins from protein level down to atom level.

  Tree structure
  --------------
  (invisible root)
  └── <protein name>                 [TYPE_PROTEIN]
      ├── Scenes                     [TYPE_HEADER]
      │   ├── base                   [TYPE_SCENE]
      │   └── <scene name> …
      └── Chains                     [TYPE_HEADER]
          └── <chain id>             [TYPE_CHAIN]
              └── <resi> - <resn>    [TYPE_RESIDUE]
                  └── <atom name>    [TYPE_ATOM]
  """

  def __init__(self) -> None:
    """Constructor."""
    super().__init__()
    self.create_root_node()

  # ------------------------------------------------------------------
  # Public API — building the model
  # ------------------------------------------------------------------

  def build_model_from_protein_objects(
          self,
          the_protein_objects: list["protein.Protein"],
          the_main_socket: zmq.Socket,
          a_socket: zmq.Socket,
  ) -> None:
    """Build the model from a list of Protein objects.

    Scenes and the full residue/atom hierarchy are both fetched from
    auxiliary PyMOL so the tree is populated down to the atom level.

    Args:
        the_protein_objects: Proteins to add to the model.
        the_main_socket: Main ZMQ socket for auxiliary PyMOL.
        a_socket: Secondary ZMQ socket for auxiliary PyMOL.

    Raises:
        exception.IllegalArgumentError: If any argument is ``None``.
    """
    if the_protein_objects is None:
      logger.error("the_protein_objects is None.")
      raise exception.IllegalArgumentError("the_protein_objects is None.")
    if the_main_socket is None:
      logger.error("the_main_socket is None.")
      raise exception.IllegalArgumentError("the_main_socket is None.")
    if a_socket is None:
      logger.error("a_socket is None.")
      raise exception.IllegalArgumentError("a_socket is None.")

    for tmp_protein in the_protein_objects:
      scenes = self._fetch_scenes_for_protein(tmp_protein, the_main_socket, a_socket)
      chempy_model = self._fetch_chempy_model_for_protein(tmp_protein, the_main_socket, a_socket)
      hierarchy_map = self._build_hierarchy_map(chempy_model)
      self._add_protein_node(tmp_protein, scenes=scenes, hierarchy_map=hierarchy_map)

  def add_protein_from_protein_object(
          self,
          a_protein: "protein.Protein"
  ) -> None:
    """Add a single Protein to the model with full atom-level hierarchy.

    A ``"base"`` scene is used as a fallback when the session contains no
    scenes yet.

    Args:
        a_protein: The protein to add.
        the_main_socket: Main ZMQ socket for auxiliary PyMOL.
        a_socket: Secondary ZMQ socket for auxiliary PyMOL.

    Raises:
        exception.IllegalArgumentError: If any argument is ``None``.
    """
    if a_protein is None:
      logger.error("a_protein is None.")
      raise exception.IllegalArgumentError("a_protein is None.")

    scenes = self._fetch_scenes_for_protein(a_protein)
    if not scenes:
      scenes = ["base"]

    chempy_model = self._fetch_chempy_model_for_protein(a_protein)
    hierarchy_map = self._build_hierarchy_map(chempy_model)
    self._add_protein_node(a_protein, scenes=scenes, hierarchy_map=hierarchy_map)

  def add_protein_from_chempy_model(self, a_chempy_protein: Indexed) -> None:
    """Add a protein from a chempy Indexed model (full hierarchy, no scenes).

    Use this when a live PyMOL connection is not available but a chempy
    model is.  Call :meth:`add_scene` afterwards to attach scenes.

    Args:
        a_chempy_protein: The PyMOL/chempy model (e.g. from ``cmd.get_model()``).

    Raises:
        exception.IllegalArgumentError: If ``a_chempy_protein`` is ``None``.
    """
    hierarchy_map = self._build_hierarchy_map(a_chempy_protein)
    existing_names = self._existing_protein_names()

    for obj_name, chains in hierarchy_map.items():
      if obj_name in existing_names:
        logger.info("Protein '%s' is already in the model – skipping.", obj_name)
        continue
      try:
        protein_node = self.add_node(
          a_parent_node=self.root_node,
          an_item_name=obj_name,
          an_item_type_value=TYPE_PROTEIN,
        )
        existing_names.add(obj_name)
        self._append_header_node(protein_node, LABEL_SCENES)  # empty; caller adds scenes
        chains_header = self._append_header_node(protein_node, LABEL_CHAINS)
        self._populate_chain_hierarchy(chains_header, chains)
        logger.info("Protein '%s' added to the model.", obj_name)
      except Exception as exc:
        logger.error("Failed to add protein '%s': %s", obj_name, exc, exc_info=True)

  # ------------------------------------------------------------------
  # Public API — scene management
  # ------------------------------------------------------------------

  def add_scene(
          self,
          a_model_index: QtCore.QModelIndex,
          the_scene_item: QtGui.QStandardItem,
  ) -> None:
    """Append *the_scene_item* to the Scenes header of the protein at *a_model_index*.

    *a_model_index* may point to any node within the target protein's subtree.

    Args:
        a_model_index: Any index within a protein subtree.
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
    """Return ``True`` if a ``_scratch_`` scene exists for the protein.

    *a_model_index* may point to any node within the target protein subtree.

    Args:
        a_model_index: Any index within a protein subtree.

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
  # Public API — protein management
  # ------------------------------------------------------------------

  def remove_protein(self, the_model_index_of_the_protein: QtCore.QModelIndex) -> None:
    """Remove the protein subtree at *the_model_index_of_the_protein*.

    Args:
        the_model_index_of_the_protein: Index of a protein node.

    Raises:
        exception.IllegalArgumentError: If the argument is ``None``.
        ValueError: If the node is not a protein node.
    """
    if the_model_index_of_the_protein is None:
      logger.error("the_model_index_of_the_protein is None.")
      raise exception.IllegalArgumentError("the_model_index_of_the_protein is None.")
    if self.data(the_model_index_of_the_protein, enums.ModelEnum.TYPE_ROLE) != TYPE_PROTEIN:
      raise ValueError("The provided index does not point to a protein node.")

    protein_item = self.itemFromIndex(the_model_index_of_the_protein)
    self.removeRow(protein_item.row())

  # ------------------------------------------------------------------
  # Public API — selection helpers
  # ------------------------------------------------------------------

  def get_selection_string(self, a_model_index: QtCore.QModelIndex) -> str:
    """Return a comma-separated selection string for *a_model_index*.

    Format: ``protein`` / ``protein,chain`` / ``protein,chain,resname,resi`` /
    ``protein,chain,resname,resi,atom``.

    Args:
        a_model_index: The index to build a selection string for.

    Raises:
        exception.IllegalArgumentError: If ``a_model_index`` is ``None``.
    """
    node_type = a_model_index.data(enums.ModelEnum.TYPE_ROLE)

    if node_type == TYPE_PROTEIN:
      return self._display(a_model_index)

    if node_type == TYPE_CHAIN:
      return f"{self._display(a_model_index.parent())},{self._display(a_model_index)}"

    if node_type == TYPE_RESIDUE:
      protein_name = self._display(a_model_index.parent().parent())
      chain_name = self._display(a_model_index.parent())
      residue_obj = a_model_index.data(enums.ModelEnum.OBJECT_ROLE)
      return f"{protein_name},{chain_name},{residue_obj.get_resname()},{residue_obj.id[1]}"

    if node_type == TYPE_ATOM:
      protein_name = self._display(a_model_index.parent().parent().parent())
      chain_name = self._display(a_model_index.parent().parent())
      residue_obj = a_model_index.parent().data(enums.ModelEnum.OBJECT_ROLE)
      atom_name = self._display(a_model_index)
      return f"{protein_name},{chain_name},{residue_obj.get_resname()},{residue_obj.id[1]},{atom_name}"

    return ""

  def construct_pymol_selection_string(self, a_model_index: QtCore.QModelIndex) -> str:
    """Return a PyMOL-format selection string for *a_model_index*.

    Examples: ``/ProteinA``, ``/ProteinA//A``, ``/ProteinA//A/42+GLY``,
    ``/ProteinA//A/42+GLY/CA``.

    Args:
        a_model_index: The index to build a PyMOL selection string for.

    Raises:
        exception.IllegalArgumentError: If ``a_model_index`` is ``None``.
    """
    hierarchy: list[str] = []
    item = self.itemFromIndex(a_model_index)
    while item is not None:
      hierarchy.insert(0, item.text())
      item = item.parent()

    if len(hierarchy) == 1:
      return f"/{hierarchy[0]}"
    if len(hierarchy) == 2:
      return f"/{hierarchy[0]}//{hierarchy[1]}"
    if len(hierarchy) == 3:
      return f"/{hierarchy[0]}//{hierarchy[1]}/{hierarchy[2].replace(' - ', '+')}"
    if len(hierarchy) == 4:
      return f"/{hierarchy[0]}//{hierarchy[1]}/{hierarchy[2].replace(' - ', '+')}/{hierarchy[3]}"
    return ""

  def create_counts_map(self) -> dict[str, dict[str, int]]:
    """Return chain, residue and atom counts for every protein in the model."""
    chain_counts: dict[str, int] = {}
    residue_counts: dict[str, int] = {}
    atom_counts: dict[str, int] = {}

    for protein_row in self.create_row_number_iterator():
      protein_index = self.get_index(protein_row)
      protein_name = self.get_display_data_of_index(protein_index)

      chains_header_index = self._find_header_index(protein_index, LABEL_CHAINS)
      if chains_header_index is None:
        chain_counts[protein_name] = 0
        continue

      chain_counts[protein_name] = self.rowCount(chains_header_index)

      for chain_row in self.create_row_number_iterator(chains_header_index):
        chain_index = self.get_index(chain_row, chains_header_index)
        chain_name = self.get_display_data_of_index(chain_index)
        residue_counts[chain_name] = self.rowCount(chain_index)

        for residue_row in self.create_row_number_iterator(chain_index):
          residue_index = self.get_index(residue_row, chain_index)
          residue_name = self.get_display_data_of_index(residue_index)
          atom_counts[residue_name] = self.rowCount(residue_index)

    return {
      "chain_counts": chain_counts,
      "residue_counts": residue_counts,
      "atom_counts": atom_counts,
    }

  # ------------------------------------------------------------------
  # Private helpers
  # ------------------------------------------------------------------

  def _add_protein_node(
          self,
          a_protein: "protein.Protein",
          scenes: list[str],
          hierarchy_map: Optional[dict] = None,
  ) -> QtGui.QStandardItem:
    """Create and append a protein subtree under the root node.

    When *hierarchy_map* is provided the tree is built down to atoms and
    chain colour/object data from *a_protein* is merged into each chain
    node.  Without it only flat chain nodes (no residues/atoms) are created.

    Args:
        a_protein: The source ``Protein`` object.
        scenes: Scene names for the Scenes header.
        hierarchy_map: Optional result of :meth:`_build_hierarchy_map`.

    Returns:
        The newly created protein ``QStandardItem``.
    """
    protein_node = self.add_node(
      a_parent_node=self.root_node,
      an_item_name=a_protein.get_molecule_object(),
      an_item_type_value=TYPE_PROTEIN,
      an_item_object_value=a_protein,
    )

    scenes_header = self._append_header_node(protein_node, LABEL_SCENES)
    self._append_scenes_from_list(scenes_header, scenes)

    chains_header = self._append_header_node(protein_node, LABEL_CHAINS)

    if hierarchy_map is not None:
      chain_object_lookup = {chain.chain_letter: chain for chain in a_protein.chains}
      chains_for_this_protein = hierarchy_map.get(a_protein.get_molecule_object(), {})
      self._populate_chain_hierarchy(
        chains_header, chains_for_this_protein, chain_object_lookup=chain_object_lookup
      )
    else:
      # Fallback: flat chain nodes only (no residue/atom depth).
      for chain in a_protein.chains:
        chain_node = self.add_node(
          a_parent_node=chains_header,
          an_item_name=chain.chain_letter,
          an_item_type_value=TYPE_CHAIN,
          an_item_object_value=chain,
        )
        chain_node.setData(
          chain.pymol_parameters[enums.PymolParameterEnum.COLOR.value],
          enums.ModelEnum.CHAIN_COLOR_ROLE,
        )

    return protein_node

  def _resolve_scenes_header(
          self, a_model_index: QtCore.QModelIndex
  ) -> QtGui.QStandardItem:
    """Return the Scenes ``QStandardItem`` for the protein that owns *a_model_index*.

    Raises:
        ValueError: If the node type is not recognised.
    """
    node_type = a_model_index.data(enums.ModelEnum.TYPE_ROLE)

    if node_type == TYPE_PROTEIN:
      protein_item = self.itemFromIndex(a_model_index)
    elif node_type == TYPE_HEADER:
      protein_item = self.itemFromIndex(a_model_index).parent()
    elif node_type == TYPE_SCENE:
      protein_item = self.itemFromIndex(a_model_index).parent().parent()
    elif node_type == TYPE_CHAIN:
      # chain → Chains header → protein
      protein_item = self.itemFromIndex(a_model_index).parent().parent()
    elif node_type == TYPE_RESIDUE:
      # residue → chain → Chains header → protein
      protein_item = self.itemFromIndex(a_model_index).parent().parent().parent()
    elif node_type == TYPE_ATOM:
      # atom → residue → chain → Chains header → protein
      protein_item = self.itemFromIndex(a_model_index).parent().parent().parent().parent()
    else:
      raise ValueError(f"Unrecognised node type: '{node_type}'")

    # Scenes is always the first child of a protein node
    return protein_item.child(0, 0)

  def _find_header_index(
          self, a_protein_index: QtCore.QModelIndex, a_label: str
  ) -> Optional[QtCore.QModelIndex]:
    """Return the index of the child header labelled *a_label*, or ``None``."""
    for row in self.create_row_number_iterator(a_protein_index):
      candidate = self.get_index(row, a_protein_index)
      if (
              candidate.data(enums.ModelEnum.TYPE_ROLE) == TYPE_HEADER
              and self.get_display_data_of_index(candidate) == a_label
      ):
        return candidate
    return None

  def _display(self, a_model_index: QtCore.QModelIndex) -> str:
    """Return the DisplayRole text for *a_model_index*."""
    return a_model_index.data(QtCore.Qt.ItemDataRole.DisplayRole)

  def _existing_protein_names(self) -> set[str]:
    """Return the set of protein names currently at the top level of the model."""
    return {self.item(row).text() for row in range(self.rowCount())}
