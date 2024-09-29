import logging
from typing import Callable, Optional

from tea.concurrent import task_result
from tea.concurrent import action

from pyssa.model.central_objects import project
from pyssa.model.util import exception
from pyssa.gui.abstract_interfaces import icomponent
from pyssa.model.pyssa_logging import default_logging
from pyssa.model.qmodel import project_model
from pyssa.gui.dialog import dialog_open_project
from pyssa.gui.control import open_project_controller


logger = default_logging.setup_logger(__file__)


class OpenProjectComponent(icomponent.IComponent):
  """Create project component."""

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
    self._dialog = dialog_open_project.DialogOpenProject()
    self._controller = open_project_controller.OpenProjectController(self._dialog)
    self._controller.component_task.connect(a_callable)

  def connect_project_model_to_component(self, a_project_model: "project_model.ProjectModel") -> None:
    """Connects the project model fo the MainFrame to this component.

    Args:
      a_project_model: The project model containing all the projects of the current workspace.

    Raises:
      exception.NoneValueError: If `a_project_model` is None.

    Notes:
      The connection between any model of the MainFrame and a component is normally prohibited.
      But in this case, it is smart to use the model from the MainFrame because
      the model is not processed in any separate threads.
    """
    # <editor-fold desc="Checks">
    if a_project_model is None:
      default_logging.append_to_log_file(logger, "a_project_model is None.", logging.ERROR)
      raise exception.NoneValueError("a_project_model is None.")
    # </editor-fold>
    self._dialog.ui.projects_list_view.setModel(a_project_model)

  def run_tool(
          self,
          a_callable: Optional[Callable] = None,
  ) -> None:
    """Starts and runs the tool.

    Args:
      a_callable (Callable): The function to execute after the tool is closed.
    """
    self._controller.restore_ui()
    self._dialog.exec()
    tmp_task_result = task_result.TaskResult.from_action(
      action.Action(a_target=self.async_open_project),
      a_callable
    )
    if self._controller.was_canceled is False:
      self._controller.component_task.emit((True, tmp_task_result))
    else:
      self._controller.component_task.emit((False, tmp_task_result))

  def async_open_project(self) -> tuple["project.Project"]:
    """Opens the selected project"""
    try:
      return project.Project(self._controller.selected_project_name),
    except Exception as e:
      default_logging.append_to_log_file(logger, e.__str__(), logging.ERROR)
      raise exception.AsyncOperationFailedError(e.__str__())
    finally:
      default_logging.append_to_log_file(logger, "'open_project' method finished.", logging.DEBUG)
