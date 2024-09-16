import logging
import re

from PyQt6 import QtCore
from PyQt6 import QtGui
from pyssa.model.util import enums
from pyssa.model.preference import model_definitions
from pyssa.model.util import exception
from pyssa.model.qmodel import base_tree_model
from pyssa.model.central_objects import protein
from pyssa.model.data_classes import selection_filter
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)

__docformat__ = "google"


class ProteinModel(base_tree_model.BaseTreeModel):
  """Class for storing proteins."""

  def __init__(self) -> None:
    """Constructor."""
    super().__init__()

  @staticmethod
  def from_a_list_of_proteins(a_list_of_proteins: list["protein.Protein"]) -> "ProteinModel":
    """Alternative constructor.

    Notes:
      - If the list of proteins is empty, then an empty model with a root node will be created.

    Raises:
      exception.NoneValueError: If `a_list_of_proteins` is None.
    """
    # <editor-fold desc="Checks">
    if a_list_of_proteins is None:
      default_logging.append_to_log_file(logger, "a_list_of_proteins is None.", logging.ERROR)
      raise exception.NoneValueError("a_list_of_proteins is None.")
    # </editor-fold>
    tmp_protein_model = ProteinModel()
    tmp_protein_model.create_root_node()
    if len(a_list_of_proteins) == 0:
      default_logging.append_to_log_file(
        logger,
        "Initialized protein model with an empty list.",
        logging.INFO
      )
      return tmp_protein_model  # Return empty model if no proteins are given
    try:
      for tmp_protein in a_list_of_proteins:
        tmp_protein_model.add_protein(tmp_protein)
        # tmp_protein_model.add_node(
        #   a_parent_node=tmp_protein_model.root_node,
        #   an_item_name=tmp_protein.name,
        #   an_item_type_value=model_definitions.TypesEnum.PROTEIN_TYPE,
        #   an_item_object_value=tmp_protein
        # )
    except Exception as e:
      default_logging.append_to_log_file(
        logger,
        "A problem occurred while adding the protein nodes to the model. Following error was raised:\n" + str(e),
        logging.ERROR
      )
    finally:
      # Returning the model ensures at least some of the proteins are in the model
      return tmp_protein_model

  def add_protein(self, a_protein: "protein.Protein"):
    """Adds a protein to the model.

    Args:
      a_protein: The protein to add

    Raises:
      exception.NoneValueError: If `a_protein` is None.
    """
    # <editor-fold desc="Checks">
    if a_protein is None:
      default_logging.append_to_log_file(logger, "a_protein is None.", logging.ERROR)
      raise exception.NoneValueError("a_protein is None.")
    # </editor-fold>
    try:
      tmp_protein_node = self.add_node(
        a_parent_node=self.root_node,
        an_item_name=a_protein.name,
        an_item_type_value=model_definitions.TypesEnum.PROTEIN_TYPE,
        an_item_object_value=a_protein
      )
      for tmp_chain in a_protein.structure.get_chains():
        tmp_chain_node = self.add_node(
          a_parent_node=tmp_protein_node,
          an_item_name=tmp_chain.id,
          an_item_type_value=model_definitions.TypesEnum.CHAIN_TYPE,
          an_item_object_value=tmp_chain
        )
        for tmp_residue in tmp_chain.get_residues():
          tmp_residue_node = self.add_node(
            a_parent_node=tmp_chain_node,
            an_item_name=f"{tmp_residue.get_resname()}-{tmp_residue.id[1]}",
            an_item_type_value=model_definitions.TypesEnum.RESIDUE_TYPE,
            an_item_object_value=tmp_residue
          )
          for tmp_atom in tmp_residue.get_atoms():
            tmp_atom_node = self.add_node(
              a_parent_node=tmp_residue_node,
              an_item_name=tmp_atom.id,
              an_item_type_value=model_definitions.TypesEnum.ATOM_TYPE,
              an_item_object_value=tmp_atom
            )
      default_logging.append_to_log_file(logger, f"The protein {a_protein.name} was added to the protein model.",
                                         logging.INFO)
    except Exception as e:
      default_logging.append_to_log_file(
        logger,
        f"An error occurred while adding the protein {a_protein.name} to the protein model. The error message is {e.__str__()}.",
        logging.ERROR
      )

  def get_chain_index(
          self,
          a_chain_identifier: str,
          a_parent_index: QtCore.QModelIndex
  ):
    """Experimental"""
    if a_parent_index.data(enums.ModelEnum.TYPE_ROLE) != model_definitions.TypesEnum.PROTEIN_TYPE:
      raise ValueError(f"Wrong parent type: {a_parent_index.data(enums.ModelEnum.TYPE_ROLE)}")

    for tmp_chain_row in range(self.rowCount(self.index(a_parent_index.parent().row(), 0))):
      tmp_chain_index = self.index(tmp_chain_row, 0, a_parent_index.parent())
      if a_chain_identifier == tmp_chain_index.data(QtCore.Qt.ItemDataRole.DisplayRole):
        return tmp_chain_index

  def get_selection_string(self, a_model_index: QtCore.QModelIndex) -> str:
    """Gets the selection string based on the given model index.

    Args:
      a_model_index: The model index to get the selection string from.

    Returns:
      The position of the current selection as selection string.

    Raises:
      exception.NoneValueError: If `a_model_index` is None.

    """
    # <editor-fold desc="Checks">
    if a_model_index is None:
      default_logging.append_to_log_file(logger, "a_model_index is None.", logging.ERROR)
      raise exception.NoneValueError("a_model_index is None.")
    # </editor-fold>

    if a_model_index.data(enums.ModelEnum.TYPE_ROLE) == model_definitions.TypesEnum.PROTEIN_TYPE:
      return a_model_index.data(QtCore.Qt.ItemDataRole.DisplayRole)

    elif a_model_index.data(enums.ModelEnum.TYPE_ROLE) == model_definitions.TypesEnum.CHAIN_TYPE:
      tmp_protein_name = a_model_index.parent().data(QtCore.Qt.ItemDataRole.DisplayRole)
      tmp_chain_name = a_model_index.data(QtCore.Qt.ItemDataRole.DisplayRole)
      return f"{tmp_protein_name},{tmp_chain_name}"

    elif a_model_index.data(enums.ModelEnum.TYPE_ROLE) == model_definitions.TypesEnum.RESIDUE_TYPE:
      tmp_protein_name = a_model_index.parent().parent().data(QtCore.Qt.ItemDataRole.DisplayRole)
      tmp_chain_name = a_model_index.parent().data(QtCore.Qt.ItemDataRole.DisplayRole)
      tmp_residue_name = a_model_index.data(enums.ModelEnum.OBJECT_ROLE).get_resname()
      tmp_residue_number = a_model_index.data(enums.ModelEnum.OBJECT_ROLE).id[1]
      return f"{tmp_protein_name},{tmp_chain_name},{tmp_residue_name},{tmp_residue_number}"

    elif a_model_index.data(enums.ModelEnum.TYPE_ROLE) == model_definitions.TypesEnum.ATOM_TYPE:
      tmp_protein_name = a_model_index.parent().parent().parent().data(QtCore.Qt.ItemDataRole.DisplayRole)
      tmp_chain_name = a_model_index.parent().parent().data(QtCore.Qt.ItemDataRole.DisplayRole)
      tmp_residue_name = a_model_index.parent().data(enums.ModelEnum.OBJECT_ROLE).get_resname()
      tmp_residue_number = a_model_index.parent().data(enums.ModelEnum.OBJECT_ROLE).id[1]
      tmp_atom_name = a_model_index.data(QtCore.Qt.ItemDataRole.DisplayRole)
      return f"{tmp_protein_name},{tmp_chain_name},{tmp_residue_name},{tmp_residue_number},{tmp_atom_name}"

  def get_model_index_based_on_name(self, a_node_name: str):
    """WORK IN PROGRESS"""
    print(a_node_name)
    tmp_sele_filter = selection_filter.SelectionFilter(a_node_name)
    print(tmp_sele_filter)
    tmp_matches: list = []

    # Iterate over all proteins in the model
    for tmp_protein_row in self.create_row_number_iterator():
      tmp_protein_index = self.get_index(tmp_protein_row)
      if tmp_sele_filter.protein == self.get_display_data_of_index(tmp_protein_index):
        if not tmp_sele_filter.chain_part_is_empty():
          # --- Look at the chain part --- #
          # Iterate over all chains of the selection filter
          for tmp_chain_identifier in tmp_sele_filter.get_single_chain_identifiers():
            # Iterate over all chains in the model
            for tmp_chain_row in self.create_row_number_iterator(tmp_protein_index):
              tmp_chain_index = self.get_index(tmp_chain_row, tmp_protein_index)
              if tmp_chain_identifier == self.get_display_data_of_index(tmp_chain_index):
                if not tmp_sele_filter.residue_name_part_is_empty():
                  # --- Look at the residue name part --- #
                  # Iterate over all residue names of the selection filter
                  for tmp_residue_name in tmp_sele_filter.get_single_residue_names():
                    # Iterate over all residues in the model
                    for tmp_residue_row in self.create_row_number_iterator(tmp_chain_index):
                      tmp_residue_index = self.get_index(tmp_residue_row, tmp_chain_index)
                      if tmp_residue_name == self.get_object_data_of_index(tmp_residue_index).get_resname():
                        if not tmp_sele_filter.residue_number_part_is_empty():
                          # --- Look at the residue number part --- #
                          pass
                        else:
                          # --- No residue number part --- #
                          tmp_matches.append(tmp_residue_index)
                else:
                  # --- No residue name part --- #
                  tmp_matches.append(tmp_chain_index)
                  break
        else:
          # --- No chain part --- #
          tmp_matches.append(tmp_protein_index)
    return tmp_matches

  def get_info_as_table_model(self, a_index: QtCore.QModelIndex) -> QtGui.QStandardItemModel:
    # TODO: Finish this method!
    tmp_table_model = QtGui.QStandardItemModel()
    if a_index.data(model_definitions.RolesEnum.TYPE_ROLE) == model_definitions.TypesEnum.PROTEIN_TYPE:
      tmp_protein: "protein.Protein" = a_index.data(model_definitions.RolesEnum.OBJECT_ROLE)
      data = [
        ["Name", tmp_protein.name]
      ]
      for tmp_key in tmp_protein.chain_sequence_map:
        data.append([f"Chain {tmp_key}", str(tmp_protein.chain_sequence_map[tmp_key].seq)])
    else:
      data = []
    for row in data:
      items = [QtGui.QStandardItem(item) for item in row]
      tmp_table_model.appendRow(items)
    return tmp_table_model
