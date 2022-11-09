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

from pyssa.gui.ui.forms.auto_generated.auto_dialog_add_models import Ui_Dialog

global_var_pdb_files = []


class DialogAddModels(Qt.QtWidgets.QDialog):

    def __init__(self, parent=None):
        """Constructor

        Args:
            args
            kwargs
        """
        Qt.QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        self.ui.btn_create_projects.setEnabled(False)
        self.ui.btn_remove.setEnabled(False)
        #self.ui.list_models.setSelectionMode(PyQt5.QtWidgets.QAbstractItemView.ExtendedSelection)

        self.ui.btn_cancel.clicked.connect(self.close_dialog)
        self.ui.btn_add.clicked.connect(self.add_model_to_list)
        self.ui.btn_remove.clicked.connect(self.remove_pdb_from_list)
        self.ui.btn_create_projects.clicked.connect(self.create_projects)
        self.ui.list_models.currentItemChanged.connect(self.activate_remove_button)

        self.setWindowTitle("Add multiple models")

    # @SLOT
    def close_dialog(self):
        self.close()

    def add_model_to_list(self):
        try:
            # open file dialog
            file_names = Qt.QtWidgets.QFileDialog.getOpenFileNames(self, "Open Model",
                                                                  Qt.QtCore.QDir.homePath(),
                                                                  "PDB Files (*.pdb)")
            if file_names == ([], ""):
                print("No file has been selected.")
                return
            # display path in text box
            for file in file_names[0]:
                self.ui.list_models.addItem(str(file))
            self.ui.btn_create_projects.setEnabled(True)
        except FileNotFoundError:
            self.status_bar.showMessage("Loading the model failed!")

    def activate_remove_button(self):
        self.ui.btn_remove.setEnabled(True)

    def remove_pdb_from_list(self):
        self.ui.list_models.takeItem(self.ui.list_models.currentRow())

    def create_projects(self):
        for row_index in range(len(self.ui.list_models)):
            global global_var_pdb_files
            global_var_pdb_files.append(self.ui.list_models.item(row_index).text())
        self.close()
