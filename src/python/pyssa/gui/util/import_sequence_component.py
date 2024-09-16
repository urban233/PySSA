import logging
from typing import Callable, Optional
from tea.concurrent import task_result
from tea.concurrent import action

from src.python.pydd.components.internal.abstract_interfaces import icomponent
from pyssa.model.central_objects import sequence
from .gui import dialog_import_sequence
from .control import import_sequence_controller
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)


class ImportSequenceComponent(icomponent.IComponent):
  """Import sequence component."""

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
    self._dialog = dialog_import_sequence.DialogImportSequence()
    self._controller = import_sequence_controller.ImportSequenceController(self._dialog)
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
        a_target=self.import_sequence
      ),
      a_callable
    )
    self._controller.component_task.emit((0, tmp_task_result))

  def import_sequence(self) -> "sequence.Sequence":
    """Imports an existing sequence structure"""
    pass
