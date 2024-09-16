import logging
import pathlib
from typing import Optional

from PyQt6 import QtGui
from PyQt6 import QtCore

from pyssa.model.central_objects import ligand
from pyssa.model.preference import model_definitions
from pyssa.model.util import exception
from pyssa.model.qmodel import base_tree_model
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)

__docformat__ = "google"


class LigandModel(base_tree_model.BaseTreeModel):
  """Class for storing ligands of the current opened project."""

  def __init__(self) -> None:
    """Constructor."""
    super().__init__()

  @staticmethod
  def from_a_list_of_ligands(a_list_of_ligands: list["ligand.Ligand"]) -> "LigandModel":
    """Alternative constructor.

    Raises:
      exception.NoneValueError: If `a_list_of_ligands` is None.

    Notes:
      - If the list of ligand is empty, then an empty model with a root node will be created.
    """
    # <editor-fold desc="Checks">
    if a_list_of_ligands is None:
      default_logging.append_to_log_file(logger, "a_list_of_ligands is None.", logging.ERROR)
      raise exception.NoneValueError("a_list_of_ligands is None.")
    # </editor-fold>
    tmp_ligand_model = LigandModel()
    tmp_ligand_model.create_root_node()
    if len(a_list_of_ligands) == 0:
      default_logging.append_to_log_file(
        logger,
        "Initialized ligand model with an empty list.",
        logging.INFO
      )
      return tmp_ligand_model  # Return empty model if no ligands are given
    try:
      for tmp_ligand in a_list_of_ligands:
        tmp_ligand_model.add_node(
          a_parent_node=tmp_ligand_model.root_node,
          an_item_name=tmp_ligand.name,
          an_item_type_value=model_definitions.TypesEnum.LIGAND_TYPE,
          an_item_object_value=tmp_ligand
        )
    except Exception as e:
      default_logging.append_to_log_file(
        logger,
        "A problem occurred while adding the ligand nodes to the model. Following error was raised:\n" + str(e),
        logging.ERROR
      )
    finally:
      # Returning the model ensures at least some of the ligands are in the model
      return tmp_ligand_model

  def add_ligand(self, a_ligand: "ligand.Ligand") -> QtGui.QStandardItem:
    """Adds a ligand to the model.

    Args:
      a_ligand: The ligand to add

    Raises:
      exception.NoneValueError: If `a_ligand` is None.
    """
    # <editor-fold desc="Checks">
    if a_ligand is None:
      default_logging.append_to_log_file(logger, "a_ligand is None.", logging.ERROR)
      raise exception.NoneValueError("a_ligand is None.")
    # </editor-fold>
    tmp_ligand_node = self.add_node(
      a_parent_node=self.root_node,
      an_item_name=a_ligand.name,
      an_item_type_value=model_definitions.TypesEnum.LIGAND_TYPE,
      an_item_object_value=a_ligand
    )
    return tmp_ligand_node
