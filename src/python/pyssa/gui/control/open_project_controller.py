import logging
from typing import TYPE_CHECKING
from PyQt6 import QtWidgets
from PyQt6 import QtCore
from pyssa.gui.base_classes import base_controller
from pyssa.gui.util import ui_util
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging
from pyssa.model.util.gui_style import styles_utils

logger = default_logging.setup_logger(__file__)

__docformat__ = "google"

if TYPE_CHECKING:
  from pyssa.gui.dialog import dialog_open_project


class OpenProjectController(base_controller.BaseController):
  """Controller for the open project dialog."""

  component_task = QtCore.pyqtSignal(tuple)

  def __init__(self, a_dialog: "dialog_open_project.DialogOpenProject") -> None:
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
    self._dialog: "dialog_open_project.DialogOpenProject" = a_dialog
    """The dialog that the controller should work with."""
    self.was_canceled: bool = False
    """Flag to indicate whether the dialog was cancelled."""
    self.selected_project_name: str = ""
    """The name of the selected project name."""
    # </editor-fold>
    self.connect_all_signals()

  def connect_all_signals(self):
    """Connects all signals with their appropriate slots."""
    self._dialog.dialogClosed.connect(self.set_dialog_close_as_canceled)
    self._dialog.ui.txt_open_search.textChanged.connect(
      self._validate_open_search
    )
    self._dialog.ui.projects_list_view.clicked.connect(self.__slot_project_selected)
    self._dialog.ui.projects_list_view.doubleClicked.connect(self.__slot_open_project)
    self._dialog.ui.btn_open_project.clicked.connect(self.__slot_open_project)

  def get_dialog(self) -> QtWidgets.QDialog:
    """Gets the dialog of the controller."""
    return self._dialog
  
  def restore_ui(self) -> None:
    """Restores the UI to default values."""
    self._dialog.setup_ui()

  def _validate_open_search(self, the_entered_text: str) -> None:
    """Validates the input of the project name in real-time.

    Args:
      the_entered_text: The text entered by the user for the open search.

    Raises:
      exception.IllegalArgumentError: If `the_entered_text` is None.
    """
    # <editor-fold desc="Checks">
    if the_entered_text is None:
      logger.error("the_entered_text is None.")
      raise exception.IllegalArgumentError("the_entered_text is None.")
    # </editor-fold>
    ui_util.select_matching_string_in_q_list_view(
      self._dialog.ui.txt_open_search.text(),
      self._dialog.ui.projects_list_view,
      self._dialog.ui.txt_open_selected_project,
    )

  def set_dialog_close_as_canceled(self) -> None:
    """Sets the was_canceled flag to true."""
    self.was_canceled = True

  def __slot_project_selected(self) -> None:
    """Slot method for the project list clicked signal"""
    self._dialog.ui.txt_open_selected_project.setText(
      self._dialog.ui.projects_list_view.currentIndex().data(
        QtCore.Qt.ItemDataRole.DisplayRole
      )
    )
    self._dialog.ui.btn_open_project.setEnabled(True)

  def __slot_open_project(self) -> None:
    """Slot method for the open project button."""
    self.selected_project_name = self._dialog.ui.projects_list_view.currentIndex().data(
      QtCore.Qt.ItemDataRole.DisplayRole
    )
    self._dialog.close()  # This triggers the close event
    self.was_canceled = False
