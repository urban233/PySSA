import copy
import logging
import pathlib
import time
from typing import Callable, Optional, TYPE_CHECKING
from Bio.PDB import PDBParser
from PyQt6 import QtWidgets
from tea.concurrent import task_result
from tea.concurrent import action

from pyssa.gui.abstract_interfaces import icomponent
from pyssa.model.central_objects import protein
from pyssa.gui.control import save_protein_controller
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)


class SaveProteinComponent(icomponent.IComponent):
  """Save protein component."""

  def __init__(self, a_callable: Callable) -> None:
    """Starts and runs the tool.

    Args:
      a_callable (Callable): The function to execute after the tool is closed.

    Raises:
      exception.NoneValueError: If `a_callable` is None.
    """
    # <editor-fold desc="Checks">
    if a_callable is None:
      default_logging.append_to_log_file(logger, "a_callable is None.", logging.ERROR)
      raise exception.NoneValueError("a_callable is None.")
    # </editor-fold>
    self._dialog = QtWidgets.QFileDialog()
    self._controller = save_protein_controller.SaveProteinController(self._dialog)
    self._controller.component_task.connect(a_callable)
    self.protein = None

  def add_protein_structure_to_tool(self, a_protein):
    """Adds a protein structure to the tool by deep copying it.

    Args:
      a_protein: A PySSA protein object
    """
    self.protein = copy.deepcopy(a_protein)

  def run_tool(
          self,
          a_callable: Optional[Callable] = None
  ) -> None:
    """Starts and runs the tool.

    Args:
      a_callable (Callable): The function to execute after the tool is closed.
    """
    self._controller.open_file_dialog()
    if self._controller.input is None or self._controller.input == "":
      self._controller.was_canceled = True
    tmp_task_result = task_result.TaskResult.from_action(
      action.Action(
        a_target=self.save_protein, args=(self.protein,)
      ),
      a_callable
    )
    if self._controller.was_canceled is False:
      self._controller.component_task.emit((True, tmp_task_result))
    else:
      self._controller.component_task.emit((False, tmp_task_result))

  def save_protein(self, a_protein: "protein.Protein") -> tuple[int]:
    """Saves an existing protein structure"""
    try:
      tmp_filepath: pathlib.Path = pathlib.Path(self._controller.input)
      a_protein.dump_to_pdb_file(tmp_filepath)
      return 0,
    except Exception as e:
      print(e)
      return -1,
