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
from PyQt5.QtWidgets import QHBoxLayout
from PyQt5 import QtCore

from pyssa.model import sequence_model
from pyssa.util import constants
from pyssa.gui.ui.forms.auto_generated.auto_dialog_sequence_viewer import Ui_Dialog
from PyQt5.QtWidgets import *
import PyQt5
from PyQt5 import QtGui


class SequenceViewer(QtWidgets.QDialog):
    """This class opens a settings customization dialog."""

    def __init__(self, sequence, filename, parent=None):
        """Constructor.

        Args:
            args
            kwargs
        """
        QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)

        self.sequence_as_list = list(sequence[0].sequence)
        tmp_model = sequence_model.SequenceModel()
        tmp_model.add_sequence(sequence[0].sequence)
        self.protein_filenames = [filename.replace(".pdb", "")]
        self.ui.table_seq_viewer = QTableView(self)
        self.ui.table_seq_viewer.setModel(tmp_model)
        # # fill table
        # self.ui.table_seq_viewer.setRowCount(len(self.protein_filenames))
        # self.ui.table_seq_viewer.setColumnCount(len(self.sequence_as_list))
        # i = 0
        # for amino_acid in self.sequence_as_list:
        #     self.ui.table_seq_viewer.setItem(0, i, QTableWidgetItem(amino_acid))
        #     i += 1
        # i = 0
        # for protein_filename in self.protein_filenames:
        #     self.ui.table_seq_viewer.setVerticalHeaderItem(i, QTableWidgetItem(protein_filename))

        self.ui.table_seq_viewer.setSizeAdjustPolicy(PyQt5.QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.ui.table_seq_viewer.resizeColumnsToContents()
        # setting context menu policy on my table, "self.ui.tableWidgetGraph"
        self.ui.table_seq_viewer.setContextMenuPolicy(PyQt5.QtCore.Qt.CustomContextMenu)
        # setting context menu request  by calling a function,"self.on_context_menu"
        self.ui.table_seq_viewer.customContextMenuRequested.connect(self.on_context_menu)
        self.ui.btn_duplicate_seq.clicked.connect(self.duplicate_sequence)
        self.ui.btn_remove_seq.clicked.connect(self.remove_sequence)

        stylesheet = """
        QDialog {background-color: #F6F4F8;}
        QTableWidget {background-color: white;}
        """
        self.setStyleSheet(stylesheet)
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        self.setWindowTitle("Sequence Viewer")
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
        self.show()

    def on_context_menu(self, pos):
        menu = PyQt5.QtWidgets.QMenu()
        color_action = menu.addAction("Color residue")
        color_action.triggered.connect(self.color_residue)
        menu.exec_(self.ui.table_seq_viewer.viewport().mapToGlobal(pos))

    def color_residue(self):
        item = self.ui.table_seq_viewer.currentItem()
        self.ui.table_seq_viewer.currentItem().setBackground(PyQt5.QtGui.QColor(193, 20, 26))

    def duplicate_sequence(self):
        items = self.ui.table_seq_viewer.selectedItems()
        row_count = self.ui.table_seq_viewer.rowCount()
        self.ui.table_seq_viewer.setRowCount(row_count + 1)
        i = 0
        for amino_acid in items:
            new_item = QTableWidgetItem(amino_acid.text())
            self.ui.table_seq_viewer.setItem(row_count, i, new_item)
            i += 1

    def remove_sequence(self):
        self.ui.table_seq_viewer.removeRow(self.ui.table_seq_viewer.currentRow())
