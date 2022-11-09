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
from pyssa.gui.ui.forms.auto_generated.auto_dialog_add_model import Ui_Dialog


global_var_add_model = ("", False)


class DialogAddModel(Qt.QtWidgets.QDialog):

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
        self.ui.btn_add_model.setEnabled(False)
        self.ui.btn_choose_model.clicked.connect(self.load_model)
        self.ui.btn_cancel.clicked.connect(self.close_dialog)
        self.ui.btn_add_model.clicked.connect(self.add_model)
        self.setWindowTitle("Add a model to the current project")

    # @SLOT
    def close_dialog(self):
        self.close()

    def load_model(self):
        try:
            # open file dialog
            file_name = Qt.QtWidgets.QFileDialog.getOpenFileName(self, "Open Model",
                                                                 Qt.QtCore.QDir.homePath(),
                                                                 "PDB Files (*.pdb)")
            if file_name == ("", ""):
                print("No file has been selected.")
            # display path in text box
            self.ui.txt_model.setText(str(file_name[0]))
            self.ui.btn_add_model.setEnabled(True)
        except FileNotFoundError:
            self.status_bar.showMessage("Loading the model failed!")

    def add_model(self):
        global global_var_add_model
        global_var_add_model = (self.ui.txt_model.text(), True)
        self.close()
