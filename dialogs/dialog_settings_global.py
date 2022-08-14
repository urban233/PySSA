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

from utils import constants, tools
from uiForms.auto.auto_dialog_settings_global import Ui_Dialog


class DialogSettingsGlobal(Qt.QtWidgets.QDialog):
    # define const var for further gui
    SETTINGS_DIR = constants.SETTINGS_DIR
    WORKSPACE_PATH_NODE = constants.WORKSPACE_PATH_NODE
    PDB_STORAGE_PATH_NODE = constants.PDB_STORAGE_PATH_NODE
    ZIP_STORAGE_PATH_NODE = constants.ZIP_STORAGE_PATH_NODE
    CYCLES_VALUE_NODE = constants.CYCLES_VALUE_NODE
    CUTOFF_VALUE_NODE = constants.CUTOFF_VALUE_NODE
    ATTRIBUTE = constants.ATTRIBUTE

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

        self.ui.txt_workspace_dir.setEnabled(False)
        self.ui.txt_pdb_storage_dir.setEnabled(False)
        self.ui.txt_zip_storage_dir.setEnabled(False)

        # sets default values
        self.xmlObj = tools.SettingsXml(self.SETTINGS_DIR)
        self.xmlFile = self.xmlObj.load_xml_in_memory()
        self.ui.txt_workspace_dir.setText(self.xmlObj.get_path(self.xmlFile,
                                                               self.WORKSPACE_PATH_NODE,
                                                               self.ATTRIBUTE))
        self.ui.txt_pdb_storage_dir.setText(self.xmlObj.get_path(self.xmlFile,
                                                                 self.PDB_STORAGE_PATH_NODE,
                                                                 self.ATTRIBUTE))
        self.ui.txt_zip_storage_dir.setText(self.xmlObj.get_path(self.xmlFile,
                                                                 self.ZIP_STORAGE_PATH_NODE,
                                                                 self.ATTRIBUTE))
        self.ui.spb_cycles.setValue(int(self.xmlObj.get_path(self.xmlFile,
                                                             self.CYCLES_VALUE_NODE,
                                                             self.ATTRIBUTE)))
        self.ui.dspb_cutoff.setValue(float(self.xmlObj.get_path(self.xmlFile,
                                                                self.CUTOFF_VALUE_NODE,
                                                                self.ATTRIBUTE)))

        # customize spin boxes
        self.ui.spb_cycles.setMinimum(0)
        # self.ui.spb_cycles.setMaximum(20) # is a maximum needed?
        self.ui.spb_cycles.setSingleStep(1)
        self.ui.dspb_cutoff.setMinimum(0.00)
        self.ui.dspb_cutoff.setMaximum(20.00)
        self.ui.dspb_cutoff.setSingleStep(0.1)

        # connect elements with function
        self.ui.btn_workspace_dir.clicked.connect(self.chooseWorkspaceDir)
        self.ui.btn_pdb_storage_dir.clicked.connect(self.choosePdbStorageDir)
        self.ui.btn_zip_storage_dir.clicked.connect(self.chooseZipStorageDir)
        self.ui.btn_cancel.clicked.connect(self.cancelDialog)
        self.ui.btn_ok.clicked.connect(self.okDialog)

        self.setWindowTitle("Global Settings")

    # @SLOT()
    def chooseWorkspaceDir(self):
        currentFilePath = self.ui.txt_workspace_dir.text()
        newFilePath = Qt.QtWidgets.QFileDialog.getExistingDirectory(
            self, "Open Directory",currentFilePath,
            options=Qt.QtWidgets.QFileDialog.ShowDirsOnly)
        if newFilePath == "":
            self.ui.txt_workspace_dir.setText(currentFilePath)
        else:
            self.ui.txt_workspace_dir.setText(newFilePath)

    def choosePdbStorageDir(self):
        currentFilePath = self.ui.txt_pdb_storage_dir.text()
        newFilePath = Qt.QtWidgets.QFileDialog.getExistingDirectory(
            self, "Open Directory", currentFilePath,
            options=Qt.QtWidgets.QFileDialog.ShowDirsOnly)
        if newFilePath == "":
            self.ui.txt_pdb_storage_dir.setText(currentFilePath)
        else:
            self.ui.txt_pdb_storage_dir.setText(newFilePath)

    def chooseZipStorageDir(self):
        currentFilePath = self.ui.txt_zip_storage_dir.text()
        newFilePath = Qt.QtWidgets.QFileDialog.getExistingDirectory(
            self, "Open Directory", currentFilePath,
            options=Qt.QtWidgets.QFileDialog.ShowDirsOnly)
        if newFilePath == "":
            self.ui.txt_zip_storage_dir.setText(currentFilePath)
        else:
            self.ui.txt_zip_storage_dir.setText(newFilePath)

    def cancelDialog(self):
        self.close()

    def okDialog(self):
        self.xmlObj.set_value(self.xmlFile, self.WORKSPACE_PATH_NODE,
                              self.ATTRIBUTE,
                              self.ui.txt_workspace_dir.text())
        self.xmlObj.set_value(self.xmlFile, self.PDB_STORAGE_PATH_NODE,
                              self.ATTRIBUTE,
                              self.ui.txt_pdb_storage_dir.text())
        self.xmlObj.set_value(self.xmlFile, self.ZIP_STORAGE_PATH_NODE,
                              self.ATTRIBUTE,
                              self.ui.txt_zip_storage_dir.text())
        self.xmlObj.set_value(self.xmlFile, self.CYCLES_VALUE_NODE,
                              self.ATTRIBUTE,
                              str(self.ui.spb_cycles.value()))
        self.xmlObj.set_value(self.xmlFile, self.CUTOFF_VALUE_NODE,
                              self.ATTRIBUTE,
                              str(self.ui.dspb_cutoff.value()))
        self.xmlObj.save_xml_file(self.xmlFile)
        self.close()
