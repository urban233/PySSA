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
import pathlib
import subprocess
import os
import logging
from pyssa.gui.utilities import gui_utils
from pyssa.gui.ui.forms.auto_generated.auto_dialog_settings_global import Ui_Dialog
from pyssa.gui.ui.dialogs import dialog_message_wsl
from pyssa.gui.ui.dialogs import dialog_message_local_colabfold
from pyssa.gui.data_structures import settings
from pyssa.gui.utilities import constants
from pymol import Qt
from PyQt5.QtWidgets import QMessageBox
from pyssa.gui.utilities import tools

# setup logger
logging.basicConfig(level=logging.DEBUG)


class DialogSettingsGlobal(Qt.QtWidgets.QDialog):
    """This class opens a settings customization dialog.

    """
    """This variable is for controlling whether the dialog opens or not"""
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
        self.tmp_settings = settings.Settings(constants.SETTINGS_DIR, constants.SETTINGS_FILENAME)
        self.settings = self.tmp_settings.deserialize_settings()

        self.ui.txt_workspace_dir.setEnabled(False)
        self.ui.txt_zip_storage_dir.setEnabled(False)

        # try:
        #     self.xmlFile = self.xmlObj.load_xml_in_memory()
        # except IsADirectoryError:
        #     logging.error("There is only a directory and not a file.")
        #     gui_utils.error_dialog("There is only a directory.", "More information in the log.")
        #     # is used to stop opening this dialog
        #     self.ERROR = True
        # except FileNotFoundError:
        #     logging.error("The settings.xml file was not found and could therefore not be opened.")
        #     gui_utils.error_dialog("The settings.xml file was not found and could therefore not be opened.", "")
        #     # is used to stop opening this dialog
        #     self.ERROR = True

        logging.info("Settings dialog was opened.")
        # loading information from the settings.xml
        self.ui.txt_workspace_dir.setText(str(self.settings.get_workspace_path()))
        self.ui.txt_zip_storage_dir.setText(str(self.settings.get_prediction_path()))
        self.ui.spb_cycles.setValue(int(self.settings.get_cycles()))
        self.ui.dspb_cutoff.setValue(float(self.settings.get_cutoff()))
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
        self.ui.btn_zip_storage_dir.clicked.connect(self.chooseZipStorageDir)
        self.ui.btn_cancel.clicked.connect(self.cancelDialog)
        self.ui.btn_ok.clicked.connect(self.okDialog)
        self.ui.btn_install_local_prediction.clicked.connect(self.install_local_colabfold)
        self.ui.btn_install_wsl2.clicked.connect(self.install_wsl)
        
        self.setWindowTitle("Global Settings")

    # @SLOT()
    def chooseWorkspaceDir(self):
        gui_utils.choose_directory(self, self.ui.txt_workspace_dir)

    def chooseZipStorageDir(self):
        gui_utils.choose_directory(self, self.ui.txt_zip_storage_dir)

    def cancelDialog(self):
        self.close()

    def okDialog(self):
        self.settings.set_workspace_path(self.ui.txt_workspace_dir.text())
        self.settings.set_prediction_path(self.ui.txt_zip_storage_dir.text())
        self.settings.set_cycles(str(self.ui.spb_cycles.value()))
        self.settings.set_cutoff(str(self.ui.dspb_cutoff.value()))
        self.settings.serialize_settings()
        self.close()

    def install_local_colabfold(self):
        home_path_wsl = r"\\wsl$\Ubuntu\home"
        colabfold_username = os.listdir(r"\\wsl$\Ubuntu\home")
        colabbatch_path = r"\.pyssa\colabfold_batch\bin\colabfold_batch"
        path_colabfold = home_path_wsl  + "\\" + colabfold_username[0] + colabbatch_path
        if os.path.exists(path_colabfold):
            subprocess.run(["wsl", "rm", "-r", "/home/$USER/.pyssa"])
            self.ui.btn_install_local_prediction.setText("Install")
        else:
            dialog_message_local_colabfold.installation_local_colabfold_accept(self)
            dialog_message_local_colabfold.installation_local_colabfold_progress(self)
            dialog_message_local_colabfold.installation_local_colabfold_end(self)
            self.ui.btn_install_local_prediction.setText("Uninstall")

    def install_wsl(self):
        dialog = dialog_message_wsl.DialogMessageWsl()
        dialog.exec_()
