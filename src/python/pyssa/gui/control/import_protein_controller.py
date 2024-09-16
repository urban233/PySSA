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
  from pyssa.gui.dialog import dialog_import_protein


class ImportProteinController(base_controller.BaseController):
  """Controller for the create project dialog."""

  component_task = QtCore.pyqtSignal(tuple)

  def __init__(self, a_dialog: "dialog_import_protein.DialogImportProtein") -> None:
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
    self._dialog: "dialog_import_protein.DialogImportProtein" = a_dialog
    """The dialog that the controller should work with."""
    self.was_canceled: bool = False
    """Flag to indicate whether the dialog was cancelled."""
    self.input: str = ""
    """The filepath or PDB id of the protein."""
    self.connect_all_signals()

  def connect_all_signals(self):
    """Connects all signals with their appropriate slots."""
    self._dialog.ui.btn_choose_protein.clicked.connect(self.__slot_load_protein_from_filesystem)
    self._dialog.ui.btn_add_protein.clicked.connect(self.__slot_import_protein)

  def get_dialog(self) -> QtWidgets.QDialog:
    """Gets the dialog of the controller."""
    return self._dialog
  
  def restore_ui(self) -> None:
    """Restores the UI to default values."""
    self._dialog.setup_ui()

  def set_dialog_close_as_canceled(self) -> None:
    """Sets the was_canceled flag to true."""
    self.was_canceled = True

  def __slot_load_protein_from_filesystem(self) -> None:
    """Loads a protein from the filesystem into the textbox."""
    default_logging.append_to_log_file(
      logger, "'Load protein from filesystem' button was clicked."
    )
    try:
      file_name = QtWidgets.QFileDialog.getOpenFileName(
        self._dialog,
        "Open existing protein",
        QtCore.QDir.homePath(),
        "PDB Files (*.pdb)",
      )
      if file_name == ("", ""):
        self._dialog.ui.lbl_status.setText("No file has been selected.")
      else:
        # display path in text box
        self._dialog.ui.txt_add_protein.setText(str(file_name[0]))
        self._dialog.ui.btn_add_protein.setEnabled(True)
    except FileNotFoundError:
      self._dialog.ui.lbl_status.setText("Loading the protein structure failed!")

  def __slot_import_protein(self) -> None:
    """Slot method for the open project button."""
    self.input = self._dialog.ui.txt_add_protein.text()
    self._dialog.close()  # This triggers the close event
    self.was_canceled = False
