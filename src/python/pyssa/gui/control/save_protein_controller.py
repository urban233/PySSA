import logging
from typing import TYPE_CHECKING
from PyQt6 import QtWidgets
from PyQt6 import QtCore
from pyssa.gui.base_classes import base_controller
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)

__docformat__ = "google"

if TYPE_CHECKING:
  from pyssa.gui.dialog import dialog_save_protein


class SaveProteinController(base_controller.BaseController):
  """Controller for the create project dialog."""

  component_task = QtCore.pyqtSignal(tuple)

  def __init__(self, a_dialog: QtWidgets.QFileDialog) -> None:
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
    self._dialog: QtWidgets.QFileDialog = a_dialog
    """The dialog that the controller should work with."""
    self.was_canceled: bool = False
    """Flag to indicate whether the dialog was cancelled."""
    self.input: str = ""
    """The filepath of the protein."""

  def get_dialog(self) -> QtWidgets.QDialog:
    """Gets the dialog of the controller."""
    return self._dialog
  
  def restore_ui(self) -> None:
    """Restores the UI to default values."""
    raise NotImplementedError("This method is not implemented because the dialog is a QFileDialog!")

  def set_dialog_close_as_canceled(self) -> None:
    """Sets the was_canceled flag to true."""
    self.was_canceled = True

  def open_file_dialog(self) -> None:
    """Open a QFileDialog for the saving process."""
    desktop_path = QtCore.QStandardPaths.standardLocations(
      QtCore.QStandardPaths.StandardLocation.DesktopLocation
    )[0]
    self._dialog.setDirectory(desktop_path)
    file_path, _ = self._dialog.getSaveFileName(
      None,
      "Save Protein Structure",
      "",
      "Protein Data Bank File (*.pdb)",
    )
    self.input = file_path
