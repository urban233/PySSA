#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/urban233/PySSA>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""Module for the import sequence view controller."""
import logging
from collections import defaultdict
from Bio.SeqRecord import SeqRecord
from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.qt import QtWidgets
from src.pyssa.gui.ui.views import import_sequence_view
from src.pyssa.controller import fasta_file_import_preview_view_controller
from src.pyssa.internal.data_structures.data_classes import basic_seq_info
from src.pyssa.logging_pyssa import log_levels, log_handlers
from src.pyssa.util import exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class ImportSequenceViewController(QtCore.QObject):
  """Class for the ImportSequenceViewController class."""



  def __init__(
      self, the_app_state: "app_state.AppState", a_parent=None
  ):
    """Constructor.

    Args:
        the_app_state (app_state.AppState): The AppState object.
        a_parent: Parent widget to pass to the view.

    Raises:
        exception.IllegalArgumentError: If `the_app_state` is None.
    """
    # <editor-fold desc="Checks">
    if the_app_state is None:
      logger.error("the_app_state is None.")
      raise exception.IllegalArgumentError("the_app_state is None.")

    # </editor-fold>

    super().__init__()
    self._app_state = the_app_state
    self._view = import_sequence_view.ImportSequenceView(a_parent)
    self._external_controller = None
    self._parsed_sequences: list[basic_seq_info.BasicSeqInfo] = []
    self._parsed_seq_records = []
    self._connect_all_ui_elements_to_slot_functions()

  def get_view(self):
    return self._view

  def _connect_all_ui_elements_to_slot_functions(self) -> None:
    """Connects all UI elements to their corresponding slot functions in the class."""
    self._view.ui.btn_choose_fasta_file.clicked.connect(self.load_model)
    self._view.ui.btn_preview.clicked.connect(self._open_preview)
    self._view.ui.btn_import_sequence.clicked.connect(self.import_sequence)
    self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)

  def restore_ui(self) -> None:
    """Restores the UI."""
    self._view.ui.txt_import_sequence.setText("")
    self._view.ui.lbl_status.setText("")
    self._view.ui.btn_preview.setEnabled(False)
    self._view.ui.btn_import_sequence.setEnabled(False)

  def _open_help_for_dialog(self) -> None:
    """Opens the help dialog for the corresponding dialog."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Help' button was clicked."
    )
    # self._interface_manager.help_manager.open_sequence_import_page()

  # @SLOT
  def _open_preview(self) -> None:
    """Opens FASTA file import preview."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Preview' button was clicked."
    )
    self._external_controller = fasta_file_import_preview_view_controller.FastaFileImportPreviewViewController(
        the_app_state=self._app_state,
        the_parsed_sequences=self._parsed_sequences,
        on_save_callback=self._post_open_preview,
        a_parent=self._view.window()
    )
    self._external_controller.restore_ui()
    self._external_controller.fill_sequence_table()
    self._external_controller.get_view().show()

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

  def load_model(self) -> None:
    """Loads a protein from the filesystem into the textbox."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Load fasta file' button was clicked.",
    )
    try:
      # open file dialog
      file_name = QtWidgets.QFileDialog.getOpenFileName(
          self._view,
          "Import existing sequence",
          QtCore.QDir.homePath(),
          "FASTA Files (*.fasta)",
      )
      if file_name == ("", ""):
        self._view.ui.lbl_status.setText("No file has been selected.")
      else:
        # display path in text box
        self._view.ui.txt_import_sequence.setText(str(file_name[0]))
        self._parsed_sequences: list[basic_seq_info.BasicSeqInfo] = (
            self.parse_fasta(self._view.ui.txt_import_sequence.text())
        )
        self._parsed_seq_records = self._convert_seqs_to_seq_records()
        self._view.ui.btn_preview.setEnabled(True)
        self._view.ui.btn_import_sequence.setEnabled(True)
    except FileNotFoundError:
      self._view.ui.lbl_status.setText("Loading the protein sequence failed!")
      self._view.ui.btn_preview.setEnabled(False)
      self._view.ui.btn_import_sequence.setEnabled(False)

  def parse_fasta(self, file_path: str) -> list["basic_seq_info.BasicSeqInfo"]:
    """Parses a FASTA file and returns a list of BasicSeqInfo objects.

    Args:
        file_path (str): The path to the FASTA file.

    Returns:
        A list of BasicSeqInfo objects, each representing a sequence entry in the FASTA file.

    Raises:
        exception.IllegalArgumentError: If `file_path` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if file_path is None or file_path == "":
      logger.error("file_path is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "file_path is either None or an empty string."
      )
    # </editor-fold>

    sequences = []
    with open(file_path, "r") as file:
      tmp_chains = []
      tmp_sequence = ""
      tmp_name = ""
      for line in file:
        line = line.strip()
        if line.startswith(">"):
          for tmp_chain in tmp_chains:
            tmp_seq_info = basic_seq_info.BasicSeqInfo(
                tmp_name, tmp_chain, tmp_sequence
            )
            sequences.append(tmp_seq_info)
          tmp_chains = []
          tmp_sequence = ""
          tmp_name = line.split("|")[0].replace(">", "")

          if line.find("Chains") != -1:
            chain_info = line.split("|")[1].replace("Chains ", "")
            tmp_chains = chain_info.split(", ")  # Extracting chain letters
          elif line.find("Chain") != -1:
            chain_info = line.split("|")[1].replace("Chain ", "")
            tmp_chains = chain_info.split(", ")  # Extracting chain letters
          else:
            tmp_chains = ["A"]
        else:
          tmp_sequence += line
      # for the last fasta file entry
      for tmp_chain in tmp_chains:
        tmp_seq_info = basic_seq_info.BasicSeqInfo(
            tmp_name, tmp_chain, tmp_sequence
        )
        sequences.append(tmp_seq_info)
    return sequences

  def merge_sequences(self, seq_infos: list) -> list:
    """Merges sequences from multiple sequence info objects.

    Args:
        seq_infos (list): A list of sequence info objects.

    Returns:
        A list of merged sequence info objects.

    Raises:
        exception.IllegalArgumentError: If `seq_infos` is None.
    """
    # <editor-fold desc="Checks">
    if seq_infos is None:
      logger.error("seq_infos is None.")
      raise exception.IllegalArgumentError("seq_infos is None.")

    # </editor-fold>

    merged_seqs = defaultdict(lambda: {"name": "", "chain": "", "seq": ""})

    for seq_info in seq_infos:
      name_base = seq_info.name

      if merged_seqs[name_base]["name"] == "":
        merged_seqs[name_base]["name"] = name_base

      if merged_seqs[name_base]["chain"] == "":
        merged_seqs[name_base]["chain"] = seq_info.chain
      else:
        merged_seqs[name_base]["chain"] += seq_info.chain

      if merged_seqs[name_base]["seq"] == "":
        merged_seqs[name_base]["seq"] = seq_info.seq
      else:
        merged_seqs[name_base]["seq"] += "," + seq_info.seq

    merged_seq_infos = [
        basic_seq_info.BasicSeqInfo(
            name=val["name"], chain=val["chain"], seq=val["seq"]
        )
        for val in merged_seqs.values()
    ]
    return merged_seq_infos

  def _convert_seqs_to_seq_records(self) -> list:
    """Converts parsed sequences into sequence records.

    Returns:
        list: A list of SeqRecord objects representing the parsed sequences.
    """
    merged_seq_infos = self.merge_sequences(self._parsed_sequences)

    tmp_seq_records = []
    for seq_info in merged_seq_infos:
      tmp_seq_record = SeqRecord(seq_info.seq, name=seq_info.name)
      tmp_seq_records.append(tmp_seq_record)
    return tmp_seq_records

  def import_sequence(self) -> None:
    """Imports a sequence into the project."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Import' button was clicked."
    )

    from src.pyssa.internal.thread.thread_api import thread_runtime
    
    if not self._parsed_seq_records:
        logger.error("No sequences parsed to import.")
        QtWidgets.QMessageBox.critical(self._view, "Import Error", "No sequences to import.")
        return

    self._view.ui.btn_import_sequence.setEnabled(False)

    def import_task(progress_callback, is_cancelled):
        for tmp_seq_record in self._parsed_seq_records:
            logger.info(
              f"Adding new sequence {tmp_seq_record.name} with {tmp_seq_record.seq} to the current project."
            )
            self._app_state.hot_db.insert_sequence(
              str(tmp_seq_record.id),
              str(tmp_seq_record.seq),
              tmp_seq_record.name,
              self._app_state.project.get_id()
            )

    def on_success(result):
        for tmp_seq_record in self._parsed_seq_records:
            self._app_state.project.sequences.append(tmp_seq_record)
            
        self._app_state.pyssa_objects_model.build_model(self._app_state.project)
        self._app_state.status_bar_manager.show_permanent_message("", False)
        self._app_state.status_bar_manager.show_temporary_message("Sequence(s) imported.")
        
    def on_error(exc):
        logger.exception("Failed to insert sequence into database.", exc_info=exc)
        QtWidgets.QMessageBox.critical(self._view, "Import Error", f"An error occurred: {exc}")
        self._view.ui.btn_import_sequence.setEnabled(True)
        self._app_state.status_bar_manager.show_error_message(
          "Failed to import sequence(s)", True
        )

    thread_runtime.get_singleton_thread_runtime().run(import_task).on_success(on_success).on_error(on_error)
    self._app_state.status_bar_manager.show_permanent_message(
      "Importing sequence(s) ...", True
    )
    self._view.close()
