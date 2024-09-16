import logging
from typing import TYPE_CHECKING
from PyQt6 import QtWidgets
from PyQt6 import QtCore
from pyssa.gui.base_classes import base_controller
from pyssa.model.pyssa_logging import default_logging
from pyssa.model.util import exception

from pyssa.gui.control import fasta_file_import_preview_controller

if TYPE_CHECKING:
  from pyssa.gui.dialog import dialog_import_sequence

logger = default_logging.setup_logger(__file__)

__docformat__ = "google"


class ImportSequenceController(base_controller.BaseController):
  """Controller for the create project dialog."""

  component_task = QtCore.pyqtSignal(tuple)

  def __init__(self, a_dialog: "dialog_import_sequence.DialogImportSequence") -> None:
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
    self._dialog: "dialog_import_sequence.DialogImportSequence" = a_dialog
    """The dialog that the controller should work with."""
    self.connect_all_signals()

  def connect_all_signals(self):
    """Connects all signals with their appropriate slots."""
    self._dialog.ui.btn_choose_fasta_file.clicked.connect(self.__slot_open_sequence_from_filesystem)

  def get_dialog(self) -> QtWidgets.QDialog:
    """Gets the dialog of the controller."""
    return self._dialog
  
  def restore_ui(self) -> None:
    """Restores the UI to default values."""
    self._dialog.setup_ui()

  def _open_preview(self) -> None:
    """Opens FASTA file import preview."""
    default_logging.append_to_log_file(logger, "'Preview' button was clicked.")
    self._external_controller = fasta_file_import_preview_controller.FastaFileImportPreviewController(
      self._parsed_sequences
    )
    self._external_controller.user_input.connect(self._post_open_preview)
    self._external_controller.restore_ui()
    self._external_controller.fill_sequence_table()
    self._interface_manager.get_fasta_file_import_preview_view().show()

  def _post_open_preview(self, return_value: tuple) -> None:
    """Processes the return values from the 'open_preview' method.

    Args:
        return_value: A tuple containing the return values from the 'open_preview' method.

    Raises:
        exception.IllegalArgumentError: If `return_value` is None.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      raise exception.IllegalArgumentError("return_value is None.")

    # </editor-fold>

    _, self._parsed_sequences = return_value
    self._parsed_seq_records = self._convert_seqs_to_seq_records()
    print(self._parsed_seq_records)

  def __slot_open_sequence_from_filesystem(self) -> None:
    """Loads a protein from the filesystem into the textbox."""
    default_logging.append_to_log_file(
      logger, "'Open sequence from filesystem' button was clicked."
    )
    try:
      file_name = QtWidgets.QFileDialog.getOpenFileName(
        self._dialog,
        "Open existing sequence file",
        QtCore.QDir.homePath(),
        "FASTA Files (*.fasta)",
      )
      if file_name == ("", ""):
        self._dialog.ui.lbl_status.setText("No file has been selected.")
      else:
        # display path in text box
        self._dialog.ui.txt_import_sequence.setText(str(file_name[0]))
        self._dialog.ui.btn_import_sequence.setEnabled(True)
    except FileNotFoundError:
      self._dialog.ui.lbl_status.setText("Opening the protein sequence failed!")
