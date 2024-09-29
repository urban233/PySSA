import logging
import pathlib
import time
from typing import Callable, Optional
from tea.concurrent import task_result
from tea.concurrent import action
from rdkit import Chem

from pyssa.gui.abstract_interfaces import icomponent
from pyssa.model.central_objects import ligand
from pyssa.gui.dialog import dialog_import_ligand
from pyssa.gui.control import import_ligand_controller
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)


class ImportLigandComponent(icomponent.IComponent):
  """Import protein component."""

  def __init__(self, a_callable: Callable):
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
    self._dialog = dialog_import_ligand.DialogImportLigand()
    self._controller = import_ligand_controller.ImportLigandController(self._dialog)
    self._controller.component_task.connect(a_callable)

  def run_tool(
          self,
          a_callable: Optional[Callable] = None
  ) -> None:
    """Starts and runs the tool.

    Args:
      a_callable (Callable): The function to execute after the tool is closed.
    """
    self._dialog.exec()
    tmp_task_result = task_result.TaskResult.from_action(
      action.Action(
        a_target=self.import_ligand
      ),
      a_callable
    )
    self._controller.component_task.emit((0, tmp_task_result))

  def import_ligand(self) -> tuple["ligand.Ligand"]:
    """Imports an existing ligand structure"""
    tmp_supplier = Chem.SDMolSupplier(r"C:\Users\student\Downloads\ATP_ideal.sdf")
    tmp_structure = tmp_supplier[0]
    return ligand.Ligand("ATP_ideal", tmp_structure),
