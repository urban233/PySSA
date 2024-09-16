from typing import TYPE_CHECKING
from PyQt6 import QtWidgets
from PyQt6 import QtCore
from pyssa.gui.base_classes import base_controller
from pyssa.model.pyssa_logging import default_logging
from pyssa.model.util import exception


if TYPE_CHECKING:
  from pyssa.gui.dialog import dialog_fasta_file_import_preview

logger = default_logging.setup_logger(__file__)

__docformat__ = "google"


class FastaFileImportPreviewController(base_controller.BaseController):
  """Controller for the fasta file import preview dialog."""

  component_task = QtCore.pyqtSignal(tuple)

  def __init__(self, a_dialog):
    """Constructor.
    
    Args:
      a_dialog: a dialog instance to be managed by the controller
    """
    super().__init__()
    self._dialog: "dialog_fasta_file_import_preview.DialogFastaFileImportPreview" = a_dialog
    """The dialog that the controller should work with."""
    self.connect_all_signals()

  def connect_all_signals(self):
    """Connects all signals with their appropriate slots."""
    pass

  def get_dialog(self) -> QtWidgets.QDialog:
    """Gets the dialog of the controller."""
    return self._dialog
  
  def restore_ui(self) -> None:
    """Restores the UI to default values."""
    self._dialog.setup_ui()
    self._dialog.ui.lbl_description.setText("Sequences Overview")
    self._dialog.ui.lbl_description_2.setText(
      "The chain letter will be lost, if you do the import! (The chain letter is important to store the sequences in a specific order.)"
    )
    self._dialog.ui.btn_delete.setEnabled(False)
    self._dialog.ui.btn_duplicate.setEnabled(False)
    self._dialog.ui.btn_add.setEnabled(True)
    self._dialog.ui.btn_add.hide()
    self._dialog.ui.btn_insert.setEnabled(False)
    self._dialog.ui.btn_insert.hide()
    self._dialog.ui.btn_save.setEnabled(True)
    # restore table
    self._restore_table()

  def _restore_table(self) -> None:
    """Restores the table."""
    self._dialog.ui.table.clear()
    self._dialog.ui.table.setColumnCount(3)
    self._dialog.ui.table.setHorizontalHeaderLabels(["Name", "Chain", "Sequence"])

  def _restore_default_table_values(self) -> None:
    """Restores the default table values."""
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Restore' button was clicked."
    )
    self._restore_table()
    self.fill_sequence_table()

  def fill_sequence_table(self) -> None:
    """Fill the sequence table with the parsed sequences."""
    self._dialog.ui.table.setRowCount(len(self._parsed_sequences))
    i = 0
    for tmp_seq_info in self._parsed_sequences:
      tmp_name_item = QtWidgets.QTableWidgetItem(tmp_seq_info.name)
      tmp_chain_item = QtWidgets.QTableWidgetItem(tmp_seq_info.chain)
      tmp_sequence_item = QtWidgets.QTableWidgetItem(tmp_seq_info.seq)
      # tmp_sequence_item.setFlags(tmp_sequence_item.flags() & ~Qt.ItemIsEditable)
      self._dialog.ui.table.setItem(i, 0, tmp_name_item)
      self._dialog.ui.table.setItem(i, 1, tmp_chain_item)
      self._dialog.ui.table.setItem(i, 2, tmp_sequence_item)
      i += 1
    self._dialog.ui.table.resizeColumnsToContents()

  def _check_if_sequences_can_be_saved(self) -> None:
    """Checks if sequences can be saved based on the data in a table view."""
    if self._dialog.ui.table.rowCount() == 0:
      self._dialog.ui.btn_save.setEnabled(False)
      return
    else:
      self._dialog.ui.btn_save.setEnabled(True)

    try:
      for tmp_row_no in range(self._dialog.ui.table.rowCount()):
        if (
                self._dialog.ui.table.item(tmp_row_no, 1).text()
                == "Enter a chain letter!"
                or self._dialog.ui.table.item(tmp_row_no, 1).text()
                == "Invalid input!"
        ):
          self._dialog.ui.btn_save.setEnabled(False)
          return
        elif self._dialog.ui.table.item(tmp_row_no, 0).text() == "Invalid input!":
          self._dialog.ui.btn_save.setEnabled(False)
          return
    except AttributeError:
      # Is needed because the signal cellChanged triggers this function even if the item is None
      return
    else:
      self._dialog.ui.btn_save.setEnabled(True)

  # @SLOT
  def __slot_save_sequence_changes(self) -> None:
    """Save the changes made to the sequences by sending the `user_input` signal."""
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Save' button was clicked."
    )
    self._parsed_sequences = []
    for tmp_row_no in range(self._dialog.ui.table.rowCount()):
      tmp_name = self._dialog.ui.table.item(tmp_row_no, 0).text()
      tmp_chain = self._dialog.ui.table.item(tmp_row_no, 1).text()
      tmp_sequence = self._dialog.ui.table.item(tmp_row_no, 2).text()
      tmp_seq_info = basic_seq_info.BasicSeqInfo(
        tmp_name, tmp_chain, tmp_sequence
      )
      self._parsed_sequences.append(tmp_seq_info)

    self._dialog.close()
    self.user_input.emit((0, self._parsed_sequences))

  def __slot_delete_sequence_from_table(self) -> None:
    """Delete selected rows from the table."""
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Delete' button was clicked."
    )
    selected_rows = sorted(
      set(index.row() for index in self._dialog.ui.table.selectedIndexes())
    )
    # Reverse the selected rows list to ensure correct removal when removing rows
    for row in reversed(selected_rows):
      self._dialog.ui.table.removeRow(row)
    self._dialog.ui.btn_delete.setEnabled(False)
    self._check_if_sequences_can_be_saved()

  def __slot_enable_table_edit_buttons(self) -> None:
    """Enable or disable the edit buttons in the table view based on the selected items."""
    if len(self._dialog.ui.table.selectedItems()) == 0:
      self._dialog.ui.btn_delete.setEnabled(False)
      self._dialog.ui.btn_duplicate.setEnabled(False)
    else:
      self._dialog.ui.btn_delete.setEnabled(True)
      self._dialog.ui.btn_duplicate.setEnabled(True)
    self._check_if_sequences_can_be_saved()

  def __slot_duplicate_selected_sequence(self) -> None:
    """Duplicates the selected sequence in the table."""
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Duplicate' button was clicked."
    )
    tmp_name_item = QtWidgets.QTableWidgetItem(
      self._dialog.ui.table.item(self._dialog.ui.table.currentRow(), 0).text()
    )
    tmp_chain_item = QtWidgets.QTableWidgetItem("Enter a valid chain letter!")
    tmp_seq_item = QtWidgets.QTableWidgetItem(
      self._dialog.ui.table.item(self._dialog.ui.table.currentRow(), 2).text()
    )
    # tmp_seq_item.setFlags(tmp_seq_item.flags() & ~Qt.ItemIsEditable)
    self._dialog.ui.table.insertRow(self._dialog.ui.table.currentRow() + 1)
    self._dialog.ui.table.setItem(
      self._dialog.ui.table.currentRow() + 1, 0, tmp_name_item
    )
    self._dialog.ui.table.setItem(
      self._dialog.ui.table.currentRow() + 1, 1, tmp_chain_item
    )
    self._dialog.ui.table.setItem(
      self._dialog.ui.table.currentRow() + 1, 2, tmp_seq_item
    )

    self._dialog.ui.btn_duplicate.setEnabled(False)
    self._dialog.ui.btn_delete.setEnabled(False)
    self._dialog.ui.btn_save.setEnabled(False)
    self._dialog.ui.table.resizeColumnsToContents()
