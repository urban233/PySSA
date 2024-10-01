import logging
from typing import TYPE_CHECKING, Optional
from PyQt6 import QtWidgets
from PyQt6 import QtCore
from pyssa.gui.base_classes import base_controller
from pyssa.gui.util import input_validator
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)

__docformat__ = "google"

if TYPE_CHECKING:
  from pyssa.gui.dialog import dialog_create_project


class CreateProjectController(base_controller.BaseController):
  """Controller for the create project dialog."""

  component_task = QtCore.pyqtSignal(tuple)

  def __init__(self, a_dialog):
    """Constructor.
    
    Args:
      a_dialog: a dialog instance to be managed by the controller

    Raises:
      exception.NoneValueError: If `a_dialog` is None.
    """
    # <editor-fold desc="Checks">
    if a_dialog is None:
      default_logging.append_to_log_file(logger, "a_dialog is None.", logging.ERROR)
      raise exception.NoneValueError("a_dialog is None.")
    # </editor-fold>
    super().__init__()
    # <editor-fold desc="Instance attributes">
    self._dialog: "dialog_create_project.DialogCreateProject" = a_dialog
    """The dialog that the controller should work with."""
    self.was_canceled: bool = False
    """Flag to indicate whether the dialog was cancelled."""
    self.entered_project_name: str = ""
    """The project name that was entered in the line edit."""
    self.project_names: Optional[set] = None
    """All project names of the current workspace."""
    # </editor-fold>
    self.connect_all_signals()

  def connect_all_signals(self):
    """Connects all signals with their appropriate slots."""
    self._dialog.dialogClosed.connect(self.set_dialog_close_as_canceled)
    self._dialog.ui.txt_new_project_name.textChanged.connect(self.__slot_check_project_name_input)
    self._dialog.ui.btn_new_create_project.clicked.connect(self.__slot_create_project)

  def get_dialog(self) -> QtWidgets.QDialog:
    """Gets the dialog of the controller."""
    return self._dialog
  
  def restore_ui(self) -> None:
    """Restores the UI to default values."""
    self._dialog.setup_ui()

  def set_dialog_close_as_canceled(self) -> None:
    """Sets the was_canceled flag to true."""
    self.was_canceled = True

  def get_project_name(self) -> str:
    """Gets the project name from the line edit of the dialog.

    Returns:
      The project name
    """
    tmp_project_name: str = ""
    try:
      tmp_project_name: str = self._dialog.ui.txt_new_project_name.text()
    except Exception as e:
      default_logging.append_to_log_file(logger, f"An error occurred while accessing the project name. {e}", logging.ERROR)
    finally:
      return tmp_project_name

  def convert_model_into_set(self) -> set:
    """Converts the model into a set of project names.

    Returns:
        A set containing the project names.
    """
    tmp_project_names = []
    for tmp_row in range(
            self._dialog.ui.list_create_projects_view.model().rowCount()
    ):
      tmp_project_names.append(
        self._dialog.ui.list_create_projects_view.model()
        .index(tmp_row, 0)
        .data(QtCore.Qt.ItemDataRole.DisplayRole),
        )
    return set(tmp_project_names)

  def __slot_check_project_name_input(self) -> None:
    """Checks if the project name is valid or not."""
    tmp_validate_flag, tmp_stylesheet_string, tmp_message = (
      input_validator.validate_input_for_project_name(
        self._dialog.ui.txt_new_project_name.text(),
        self.project_names,
      )
    )
    self._dialog.ui.txt_new_project_name.setStyleSheet(tmp_stylesheet_string)
    if tmp_validate_flag:
      self._dialog.ui.lbl_new_status_project_name.setText("")
      self._dialog.ui.btn_new_create_project.setEnabled(True)
    else:
      self._dialog.ui.lbl_new_status_project_name.setText(tmp_message)
      self._dialog.ui.btn_new_create_project.setEnabled(False)

  def __slot_create_project(self) -> None:
    """Slot method for the create project button."""
    self.entered_project_name = self._dialog.ui.txt_new_project_name.text()
    self._dialog.close()  # This triggers the close event
    self.was_canceled = False
