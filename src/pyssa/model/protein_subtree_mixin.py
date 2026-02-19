"""Mixin that provides shared protein-subtree building logic.

Both PySSAObjectsModel and PySSAProteinPairsModel need to:
  - fetch scenes from auxiliary PyMOL
  - fetch a chempy model from auxiliary PyMOL
  - turn a chempy model into a hierarchy map
  - populate a Chains header node down to the atom level
  - append typed header nodes

Placing this logic here avoids duplicating it in both model classes while
keeping each model class focused on its own domain.
"""
from typing import Optional

import zmq
from chempy.models import Indexed

from src.auxiliary_pymol import auxiliary_pymol_client
from src.pyssa.gui.qt import QtGui
from src.pyssa.internal.data_structures import protein, job
from src.pyssa.internal.pymol import pml_worker, pml_enums
from src.pyssa.util import enums

# ---------------------------------------------------------------------------
# Node-type constants (shared by all models in this package)
# ---------------------------------------------------------------------------

TYPE_PROTEIN = enums.TypesEnum.PROTEIN_TYPE
TYPE_PROTEIN_PAIR = "protein_pair"
TYPE_CHAIN = enums.TypesEnum.CHAIN_TYPE
TYPE_RESIDUE = enums.TypesEnum.RESIDUE_TYPE
TYPE_ATOM = enums.TypesEnum.ATOM_TYPE
TYPE_HEADER = "header"
TYPE_SCENE = "scene"
TYPE_SECTION = "section"    # top-level section dividers
TYPE_SEQUENCE = "sequence"  # leaf nodes under the "Sequences" section

LABEL_SCENES = "Scenes"
LABEL_CHAINS = "Chains"
LABEL_SEQUENCES = "Sequences"
LABEL_PROTEINS = "Proteins"
LABEL_PROTEIN_PAIRS = "Protein Pairs"


class ProteinSubtreeMixin:
  """Mixin that adds reusable protein-subtree helpers to a BaseTreeModel subclass.

  This mixin has no state of its own and does not call ``super().__init__``.
  Mix it in alongside ``BaseTreeModel`` (or a subclass of it):

      class MyModel(ProteinSubtreeMixin, BaseTreeModel):
          ...
  """

  # ------------------------------------------------------------------
  # Auxiliary PyMOL fetch helpers
  # ------------------------------------------------------------------

  def _fetch_scenes(
          self,
          a_pymol_session: str,
          the_main_socket: zmq.Socket,
          a_socket: zmq.Socket,
  ) -> list[str]:
    """Return the list of scene names stored in *a_pymol_session*."""
    with pml_worker.PmlWorker.session(
            pml_worker.PmlWorker.cache_session(a_pymol_session, "fetch_scenes")
    ) as tmp_pml_worker:
      return tmp_pml_worker.do(
        pml_enums.PmlCommand.GET_SCENE_LIST, sync=True
      )

    tmp_job = job.GeneralPurposeJobDescription(
      enums.JobShortDescription.GET_ALL_SCENES_OF_SESSION
    )
    tmp_job.setup_dict(
      {enums.JobDescriptionKeys.PYMOL_SESSION.value: a_pymol_session}
    )
    tmp_reply = auxiliary_pymol_client.send_request_to_auxiliary_pymol(
      the_main_socket, a_socket, tmp_job
    )
    return tmp_reply["data"]

  def _fetch_chempy_model(
          self,
          a_pymol_session: str,
          a_molecule_object: str,
          the_main_socket: zmq.Socket,
          a_socket: zmq.Socket,
  ) -> Indexed:
    """Return the chempy ``Indexed`` model for *a_molecule_object* in *a_pymol_session*."""
    with pml_worker.PmlWorker.session(
            pml_worker.PmlWorker.cache_session(a_pymol_session, "fetch_scenes")
    ) as tmp_pml_worker:
      return tmp_pml_worker.do(
        pml_enums.PmlCommand.GET_MODEL, ("all", ), sync=True
      )

    tmp_job = job.GeneralPurposeJobDescription(
      enums.JobShortDescription.GET_MODEL
    )
    tmp_job.setup_dict(
      {
        enums.JobDescriptionKeys.PYMOL_SESSION.value: a_pymol_session,
        enums.JobDescriptionKeys.MOLECULE_OBJECT.value: a_molecule_object,
      }
    )
    tmp_reply = auxiliary_pymol_client.send_request_to_auxiliary_pymol(
      the_main_socket, a_socket, tmp_job
    )
    return tmp_reply["data"]

  def _fetch_scenes_for_protein(
          self,
          a_protein: "protein.Protein",
          the_main_socket: zmq.Socket,
          a_socket: zmq.Socket,
  ) -> list[str]:
    """Convenience wrapper: fetch scenes using a ``Protein`` object."""
    return self._fetch_scenes(
      str(a_protein.pymol_session), the_main_socket, a_socket
    )

  def _fetch_chempy_model_for_protein(
          self,
          a_protein: "protein.Protein",
          the_main_socket: zmq.Socket,
          a_socket: zmq.Socket,
  ) -> Indexed:
    """Convenience wrapper: fetch chempy model using a ``Protein`` object."""
    return self._fetch_chempy_model(
      str(a_protein.pymol_session),
      str(a_protein.get_molecule_object()),
      the_main_socket,
      a_socket,
    )

  # ------------------------------------------------------------------
  # Hierarchy map builder
  # ------------------------------------------------------------------

  @staticmethod
  def _build_hierarchy_map(
          a_chempy_protein: Indexed,
  ) -> dict[str, dict[str, dict[tuple[str, str], set[str]]]]:
    """Group atoms from *a_chempy_protein* into a nested hierarchy map.

    Returns:
        A mapping of ``object_name → chain_id → (resi, resn) → {atom_names}``.
    """
    hierarchy: dict[str, dict[str, dict[tuple[str, str], set[str]]]] = {}
    for atom in getattr(a_chempy_protein, "atom", []):
      obj_name = getattr(atom, "model", None) or ""
      chain_id = getattr(atom, "chain", None) or ""
      resi = str(getattr(atom, "resi", ""))
      resn = str(getattr(atom, "resn", ""))
      atom_name = str(getattr(atom, "name", ""))
      (
        hierarchy
        .setdefault(obj_name, {})
        .setdefault(chain_id, {})
        .setdefault((resi, resn), set())
        .add(atom_name)
      )
    return hierarchy

  # ------------------------------------------------------------------
  # Node-building helpers
  # ------------------------------------------------------------------

  def _append_header_node(
          self, a_parent_node: QtGui.QStandardItem, a_label: str
  ) -> QtGui.QStandardItem:
    """Append a typed header node to *a_parent_node* and return it."""
    header_item = QtGui.QStandardItem(a_label)
    header_item.setData(TYPE_HEADER, enums.ModelEnum.TYPE_ROLE)
    a_parent_node.appendRow(header_item)
    return header_item

  def _append_scenes_from_list(
          self,
          scenes_header: QtGui.QStandardItem,
          scene_names: list[str],
  ) -> None:
    """Append one scene node per name in *scene_names* to *scenes_header*."""
    for scene_name in scene_names:
      scene_item = QtGui.QStandardItem(scene_name)
      scene_item.setData(TYPE_SCENE, enums.ModelEnum.TYPE_ROLE)
      scenes_header.appendRow(scene_item)

  def _populate_chain_hierarchy(
          self,
          chains_header: QtGui.QStandardItem,
          chains: dict[str, dict[tuple[str, str], set[str]]],
          chain_object_lookup: Optional[dict] = None,
  ) -> None:
    """Populate *chains_header* with chain → residue → atom nodes.

    When *chain_object_lookup* is provided (a ``{chain_letter: chain_object}``
    mapping from the ``Protein`` instance) each chain node is enriched with
    the chain's ``OBJECT_ROLE`` and ``CHAIN_COLOR_ROLE`` so that existing
    view delegates keep working.

    Args:
        chains_header: The Chains header ``QStandardItem`` to populate.
        chains: Per-chain data from :meth:`_build_hierarchy_map`.
        chain_object_lookup: Optional chain-letter → chain-object mapping.
    """
    for chain_id, residues in chains.items():
      chain_object = (
        chain_object_lookup.get(chain_id) if chain_object_lookup else None
      )
      # self.add_node is provided by BaseTreeModel
      chain_node = self.add_node(  # type: ignore[attr-defined]
        a_parent_node=chains_header,
        an_item_name=chain_id,
        an_item_type_value=TYPE_CHAIN,
        an_item_object_value=chain_object,
      )
      if chain_object is not None:
        chain_node.setData(
          chain_object.pymol_parameters[enums.PymolParameterEnum.COLOR.value],
          enums.ModelEnum.CHAIN_COLOR_ROLE,
        )

      sorted_residues = sorted(
        residues.items(),
        key=lambda item: self._residue_sort_key(item[0]),
      )
      for (resi, resn), atom_names in sorted_residues:
        residue_label = f"{resi} - {resn}"
        residue_node = self.add_node(  # type: ignore[attr-defined]
          a_parent_node=chain_node,
          an_item_name=residue_label,
          an_item_type_value=TYPE_RESIDUE,
        )
        for atom_name in sorted(atom_names):
          self.add_node(  # type: ignore[attr-defined]
            a_parent_node=residue_node,
            an_item_name=atom_name,
            an_item_type_value=TYPE_ATOM,
          )

  @staticmethod
  def _residue_sort_key(resi_resn_tuple: tuple[str, str]):
    """Sort residues numerically when possible, alphabetically otherwise."""
    resi_str, _ = resi_resn_tuple
    try:
      return int(resi_str)
    except ValueError:
      return resi_str
