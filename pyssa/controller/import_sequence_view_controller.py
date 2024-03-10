#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2022
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
"""Module for the FastaFileImportPreviewView Dialog."""
import os
import pymol
from pymol import cmd
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from pyssa.controller import interface_manager, fasta_file_import_preview_view_controller
from pyssa.internal.portal import pymol_io
from pyssa.internal.thread import tasks
from pyssa.internal.thread.async_pyssa import validate_async
from pyssa.io_pyssa import bio_data
from pyssa.util import constants


class ImportSequenceViewController(QtCore.QObject):
    """Class for the ImportSequenceViewController class."""
    user_input = QtCore.pyqtSignal(tuple)

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager"):
        super().__init__()
        self._interface_manager = the_interface_manager
        self._view = the_interface_manager.get_import_sequence_view()
        self._external_controller = None
        self._parsed_sequences: dict[str, tuple[str, str]] = {}
        self._parsed_seq_records = []
        self._connect_all_ui_elements_to_slot_functions()

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.btn_choose_fasta_file.clicked.connect(self.load_model)
        self._view.ui.btn_preview.clicked.connect(self._open_preview)
        self._view.ui.btn_import_sequence.clicked.connect(self.import_sequence)

    def restore_ui(self):
        self._view.ui.txt_import_sequence.setText("")
        self._view.ui.lbl_status.setText("")
        self._view.ui.btn_preview.setEnabled(False)
        self._view.ui.btn_import_sequence.setEnabled(False)

    # @SLOT
    def _open_preview(self):
        self._external_controller = fasta_file_import_preview_view_controller.FastaFileImportPreviewViewController(
            self._interface_manager, self._parsed_sequences)
        self._external_controller.user_input.connect(self._post_open_preview)
        self._external_controller.restore_ui()
        self._external_controller.fill_sequence_table()
        self._interface_manager.get_fasta_file_import_preview_view().show()

    def _post_open_preview(self):
        pass

    def load_model(self) -> None:
        """Loads a protein from the filesystem into the textbox."""
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
                self._parsed_sequences: dict = self.parse_fasta(self._view.ui.txt_import_sequence.text())
                print(self._parsed_sequences)
                self._convert_seqs_to_seq_records()
                self._view.ui.btn_preview.setEnabled(True)
                self._view.ui.btn_import_sequence.setEnabled(True)
        except FileNotFoundError:
            self._view.ui.lbl_status.setText("Loading the protein sequence failed!")
            self._view.ui.btn_preview.setEnabled(False)
            self._view.ui.btn_import_sequence.setEnabled(False)

    def parse_fasta(self, file_path):
        sequences = {}
        with open(file_path, 'r') as file:
            tmp_chains = []
            tmp_sequence = ""
            tmp_name = ""
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    for tmp_chain in tmp_chains:
                        sequences[tmp_chain] = (tmp_name, tmp_sequence)
                    tmp_chains = []
                    tmp_sequence = ""
                    tmp_name = line.split('|')[0].replace(">", "")
                    chain_info = line.split('|')[1].replace("Chains ", "")
                    tmp_chains = chain_info.split(", ")  # Extracting chain letters
                else:
                    tmp_sequence += line
            # for the last fasta file entry
            for tmp_chain in tmp_chains:
                sequences[tmp_chain] = (tmp_name, tmp_sequence)
        return sequences

    def _convert_seqs_to_seq_records(self):
        # Create SeqRecord objects
        from Bio.SeqRecord import SeqRecord
        seq_records = []
        for chain_letter, (name, sequence) in sorted(self._parsed_sequences.items()):
            name = name.split('_')[0]  # Extracting the name without _1 or _2
            seq_record = SeqRecord(sequence, id=name, name=name, description="")
            seq_records.append(seq_record)

        # Join SeqRecord objects
        joined_seq = ','.join([str(record.seq) for record in seq_records])

        print(joined_seq)

        # for tmp_key in self._parsed_sequences:
        #     tmp_name = self._parsed_sequences[tmp_key][0]
        #     tmp_chain_letter = tmp_key
        #     tmp_seq = self._parsed_sequences[tmp_key][1]

    def import_sequence(self) -> None:
        """Adds a protein to the global variable and closes the dialog."""
        self._view.close()
        self.user_input.emit((self._view.ui.txt_import_sequence.text(), len(self._view.ui.txt_import_sequence.text())))
