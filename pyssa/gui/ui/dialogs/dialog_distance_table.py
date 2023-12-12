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
from pymol import Qt
from PyQt5.QtWidgets import (QHBoxLayout)
from PyQt5 import QtCore
from pyssa.util import constants
from pyssa.gui.ui.forms.auto_generated.auto_dialog_sequence_viewer import Ui_Dialog
from PyQt5.QtWidgets import *
import PyQt5


class DialogDistanceTable(Qt.QtWidgets.QDialog):
    """This class opens a settings customization dialog."""

    def __init__(self, path, parent=None):
        """Constructor.

        Args:
            args
            kwargs
        """
        Qt.QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
        self.labels = [
            "Residue pair no.", "Protein 1 Chain", "Protein 1 Position", "Protein 1 Residue",
            "Protein 2 Chain", "Protein 2 Position", "Protein 2 Residue", "Distance in Å",
        ]
        self.standard_item_list = []
        with open(path, 'r', encoding="utf-8") as csv_file:
            i = 0
            for line in csv_file:
                tmp_list = line.split(",")
                tmp_list.pop(0)
                if i == 0:
                    tmp_list_1 = []
                    for tmp in tmp_list:
                        tmp_list_1.append(tmp)
                    self.standard_item_list.append(tmp_list_1)
                else:
                    pair_no: int = int(tmp_list[0])
                    ref_chain: str = tmp_list[1]
                    ref_pos: int = int(tmp_list[2])
                    ref_resi: str = tmp_list[3]
                    model_chain: str = tmp_list[4]
                    model_pos: int = int(tmp_list[5])
                    model_resi: str = tmp_list[6]
                    distance: float = float(tmp_list[7])
                    self.standard_item_list.append([pair_no, ref_chain, ref_pos, ref_resi, model_chain, model_pos, model_resi, distance])
                    # self.standard_item_list.append(pair_no)
                    # self.standard_item_list.append(ref_chain)
                    # self.standard_item_list.append(ref_pos)
                    # self.standard_item_list.append(ref_resi)
                    # self.standard_item_list.append(model_chain)
                    # self.standard_item_list.append(model_pos)
                    # self.standard_item_list.append(model_resi)
                    # self.standard_item_list.append(distance)
                i += 1
            csv_file.close()
        # fill table
        self.ui.table_seq_viewer.setRowCount(len(self.standard_item_list))
        self.ui.table_seq_viewer.setColumnCount(i)

        i = 0
        for distance_table_line in self.standard_item_list:
            j = 0
            for line_content in distance_table_line:
                tmp_item = QTableWidgetItem()
                tmp_item.setData(QtCore.Qt.DisplayRole, line_content)
                #self.ui.table_seq_viewer.setItem(i, j, QTableWidgetItem(line_content))
                j += 1
            i += 1

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
        self.show()

    def on_context_menu(self, pos):
        menu = PyQt5.QtWidgets.QMenu()
        color_action = menu.addAction('Color residue')
        color_action.triggered.connect(self.color_residue)
        menu.exec_(self.ui.table_seq_viewer.viewport().mapToGlobal(pos))

    def color_residue(self):
        item = self.ui.table_seq_viewer.currentItem()
        self.ui.table_seq_viewer.currentItem().setBackground(PyQt5.QtGui.QColor(193, 20, 26))

    def duplicate_sequence(self):
        items = self.ui.table_seq_viewer.selectedItems()
        row_count = self.ui.table_seq_viewer.rowCount()
        self.ui.table_seq_viewer.setRowCount(row_count+1)
        i = 0
        for amino_acid in items:
            new_item = QTableWidgetItem(amino_acid.text())
            self.ui.table_seq_viewer.setItem(row_count, i, new_item)
            i += 1

    def remove_sequence(self):
        self.ui.table_seq_viewer.removeRow(self.ui.table_seq_viewer.currentRow())
