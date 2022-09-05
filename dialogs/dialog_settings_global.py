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
import logging
from pymol import Qt
from utils import constants, tools, gui_utils
from uiForms.auto.auto_dialog_settings_global import Ui_Dialog

# setup logger
logging.basicConfig(level=logging.DEBUG)


class DialogSettingsGlobal(Qt.QtWidgets.QDialog):
    """This class opens a settings customization dialog.

    """
    """This variable is for controlling wether the dialog opens or not"""
    ERROR = False

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

        self.xmlObj = tools.SettingsXml(constants.SETTINGS)
        try:
            self.xmlFile = self.xmlObj.load_xml_in_memory()
        except IsADirectoryError:
            logging.error("There is only a directory and not a file.")
            gui_utils.error_dialog("There is only a directory.", "More information in the log.")
            # is used to stop opening this dialog
            self.ERROR = True
        except FileNotFoundError:
            logging.error("The settings.xml file was not found and could therefore not be opened.")
            gui_utils.error_dialog("The settings.xml file was not found and could therefore not be opened.", "")
            # is used to stop opening this dialog
            self.ERROR = True

        logging.info("Settings dialog was opened.")
        # loading information from the settings.xml
        self.ui.txt_workspace_dir.setText(self.xmlObj.get_path(self.xmlFile,
                                                               constants.WORKSPACE_PATH_TAG,
                                                               constants.ATTRIBUTE))
        self.ui.txt_pdb_storage_dir.setText(self.xmlObj.get_path(self.xmlFile,
                                                                 constants.PDB_STORAGE_PATH_TAG,
                                                                 constants.ATTRIBUTE))
        self.ui.txt_zip_storage_dir.setText(self.xmlObj.get_path(self.xmlFile,
                                                                 constants.ZIP_STORAGE_PATH_TAG,
                                                                 constants.ATTRIBUTE))
        self.ui.spb_cycles.setValue(int(self.xmlObj.get_path(self.xmlFile,
                                                             constants.CYCLES_VALUE_TAG,
                                                             constants.ATTRIBUTE)))
        self.ui.dspb_cutoff.setValue(float(self.xmlObj.get_path(self.xmlFile,
                                                                constants.CUTOFF_VALUE_TAG,
                                                                constants.ATTRIBUTE)))
        logging.info("Loading values from settings.xml was successful.")
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
        gui_utils.choose_directory(self, self.ui.txt_workspace_dir)

    def choosePdbStorageDir(self):
        gui_utils.choose_directory(self, self.ui.txt_pdb_storage_dir)

    def chooseZipStorageDir(self):
        gui_utils.choose_directory(self, self.ui.txt_zip_storage_dir)

    def cancelDialog(self):
        self.close()

    def okDialog(self):
        self.xmlObj.set_value(self.xmlFile, constants.WORKSPACE_PATH_TAG,
                              constants.ATTRIBUTE,
                              self.ui.txt_workspace_dir.text())
        self.xmlObj.set_value(self.xmlFile, constants.PDB_STORAGE_PATH_TAG,
                              constants.ATTRIBUTE,
                              self.ui.txt_pdb_storage_dir.text())
        self.xmlObj.set_value(self.xmlFile, constants.ZIP_STORAGE_PATH_TAG,
                              constants.ATTRIBUTE,
                              self.ui.txt_zip_storage_dir.text())
        self.xmlObj.set_value(self.xmlFile, constants.CYCLES_VALUE_TAG,
                              constants.ATTRIBUTE,
                              str(self.ui.spb_cycles.value()))
        self.xmlObj.set_value(self.xmlFile, constants.CUTOFF_VALUE_TAG,
                              constants.ATTRIBUTE,
                              str(self.ui.dspb_cutoff.value()))
        self.xmlObj.save_xml_file(self.xmlFile)
        self.close()
