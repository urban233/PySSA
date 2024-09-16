import logging
import pathlib
import time
from typing import Callable, Optional, TYPE_CHECKING
from Bio.PDB import PDBParser
from tea.concurrent import task_result
from tea.concurrent import action

from src.python.pydd.components.internal.abstract_interfaces import icomponent
from pyssa.model.central_objects import protein
from .gui import dialog_import_protein
from .control import import_protein_controller
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)


class ImportProteinComponent(icomponent.IComponent):
  """Import protein component."""

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
    self._dialog = dialog_import_protein.DialogImportProtein()
    self._controller = import_protein_controller.ImportProteinController(self._dialog)
    self._controller.component_task.connect(a_callable)

  def run_tool(
          self,
          a_callable: Optional[Callable] = None
  ) -> None:
    """Starts and runs the tool.

    Args:
      a_callable (Callable): The function to execute after the tool is closed.
    """
    self._controller.restore_ui()
    self._dialog.exec()
    tmp_task_result = task_result.TaskResult.from_action(
      action.Action(
        a_target=self.import_protein
      ),
      a_callable
    )
    if self._controller.was_canceled is False:
      self._controller.component_task.emit((True, tmp_task_result))
    else:
      self._controller.component_task.emit((False, tmp_task_result))

  def import_protein(self) -> tuple["protein.Protein"]:
    """Imports an existing protein structure"""
    tmp_pdb_parser = PDBParser()
    tmp_filepath: pathlib.Path = pathlib.Path(self._controller.input)
    tmp_structure = tmp_pdb_parser.get_structure(tmp_filepath.name.replace(".pdb", ""), str(tmp_filepath))
    tmp_protein = protein.Protein(tmp_filepath.name.replace(".pdb", ""), a_structure=tmp_structure)
    return tmp_protein,
