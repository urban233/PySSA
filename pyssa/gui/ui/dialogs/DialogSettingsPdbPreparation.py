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

import pyssa.gui.data_structures.settings
from pyssa.gui.ui.forms.auto_generated.auto_DialogSettingsPdbPreparation import Ui_Dialog


class DialogSettingsPdbPreparation(Qt.QtWidgets.QDialog):
    SETTINGS_FULL_FILENAME = project_constants.SETTINGS_DIR
    PDB_STORAGE_PATH_NODE = project_constants.PDB_STORAGE_PATH_TAG
    ZIP_STORAGE_PATH_NODE = project_constants.ZIP_STORAGE_PATH_TAG
    ATTRIBUTE = project_constants.ATTRIBUTE

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

        self.ui.txtPdbStorageDir.setEnabled(False)
        self.ui.txtZipStorageDir.setEnabled(False)
        # sets default values

        self.xmlObj = pyssa.gui.data_structures.settings.SettingsXml(self.SETTINGS_FULL_FILENAME)
        self.xmlFile = self.xmlObj.load_xml_in_memory()
        self.ui.txtPdbStorageDir.setText(self.xmlObj.get_path(self.xmlFile,
                                                              self.PDB_STORAGE_PATH_NODE,
                                                              self.ATTRIBUTE))
        self.ui.txtZipStorageDir.setText(self.xmlObj.get_path(self.xmlFile,
                                                              self.ZIP_STORAGE_PATH_NODE,
                                                              self.ATTRIBUTE))

        self.ui.btnPdbStorageDir.clicked.connect(self.choosePdbStorageDir)
        self.ui.btnZipStorageDir.clicked.connect(self.chooseZipStorageDir)
        self.ui.btnCancel.clicked.connect(self.cancelDialog)
        self.ui.btnOk.clicked.connect(self.okDialog)

        self.setWindowTitle("PdbPreparation Settings")


    def getDirectoryPath(self):
        tmpDialog = Qt.QtWidgets.QFileDialog()
        tmpDialog.setFileMode(Qt.QtWidgets.QFileDialog.Directory)
        tmpDialog.setOption(Qt.QtWidgets.QFileDialog.ShowDirsOnly)
        tmpDialog.exec_()
        return tmpDialog.directory().path()

    # @SLOT()
    def choosePdbStorageDir(self):
        self.ui.txtPdbStorageDir.setText(self.getDirectoryPath())

    def chooseZipStorageDir(self):
        self.ui.txtZipStorageDir.setText(self.getDirectoryPath())

    def cancelDialog(self):
        self.close()

    def okDialog(self):
        self.xmlObj.set_value(self.xmlFile, self.PDB_STORAGE_PATH_NODE,
                              self.ATTRIBUTE, self.ui.txtPdbStorageDir.text())
        self.xmlObj.set_value(self.xmlFile, self.ZIP_STORAGE_PATH_NODE,
                              self.ATTRIBUTE, self.ui.txtZipStorageDir.text())
        self.xmlObj.save_xml_file(self.xmlFile)
        self.close()
