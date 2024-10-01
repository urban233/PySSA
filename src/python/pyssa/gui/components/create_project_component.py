import logging
from typing import Callable, Optional

from tea.concurrent import task_result
from tea.concurrent import action

from pyssa.gui.abstract_interfaces import icomponent
from pyssa.model.central_objects import project
from pyssa.gui.dialog import dialog_create_project
from pyssa.gui.control import create_project_controller
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging
from pyssa.model.qmodel import project_model

logger = default_logging.setup_logger(__file__)

__docformat__ = 'google'


class CreateProjectComponent(icomponent.IComponent):
  """Create project component."""

  def __init__(self, a_callable: Callable) -> None:
    """Constructor.

    Args:
      a_callable: A function to run after the tools task finished

    Raises:
      exception.NoneValueError: If `a_callable` is None.
    """
    # <editor-fold desc="Checks">
    if a_callable is None:
      default_logging.append_to_log_file(logger, "a_callable is None.", logging.ERROR)
      raise exception.NoneValueError("a_callable is None.")
    # </editor-fold>
    self._dialog = dialog_create_project.DialogCreateProject()
    self._controller = create_project_controller.CreateProjectController(self._dialog)
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
    self._dialog.ui.list_create_projects_view.setModel(a_project_model)
    self._controller.project_names = self._controller.convert_model_into_set()

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
        a_target=self.async_create_project
      ),
      a_callable
    )
    if self._controller.was_canceled is False:
      self._controller.component_task.emit((True, tmp_task_result))
    else:
      self._controller.component_task.emit((False, tmp_task_result))

  def async_create_project(self) -> tuple["project.Project"]:
    """Creates the project based on the entered project name."""
    tmp_project = project.Project(self._controller.entered_project_name)
    return tmp_project,
