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
        try:
            self.settings = self.tmp_settings.deserialize_settings()
        except ValueError:
            logging.error("Settings dialog cannot be opened due to an error.")
            return
        logging.info("Loading values from settings.xml was successful.")
        self.ui.txt_workspace_dir.setEnabled(False)
        self.ui.txt_zip_storage_dir.setEnabled(False)

        logging.info("Settings dialog was opened.")
        # loading information from the settings.xml
        self.ui.txt_workspace_dir.setText(str(self.settings.get_workspace_path()))
        self.ui.txt_zip_storage_dir.setText(str(self.settings.get_prediction_path()))
        self.ui.spb_cycles.setValue(int(self.settings.get_cycles()))
        self.ui.dspb_cutoff.setValue(float(self.settings.get_cutoff()))
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
        logging.info("Settings were successfully saved.")
        self.close()

    def install_local_colabfold(self):
        home_path_wsl = pathlib.Path("//wsl$/Ubuntu/home//")
        colabfold_username = os.getlogin()
        # colabfold_username = os.listdir(r"\\wsl$\Ubuntu\home")
        colabbatch_path = str(pathlib.Path("/.pyssa/colabfold_batch/bin/colabfold_batch"))
        path_colabfold = pathlib.Path("home_path_wsl + colabfold_username[0] + colabbatch_path")
        if os.path.exists(path_colabfold):
            msg = QMessageBox()
            msg.setIcon(QMessageBox.Question)
            msg.setText("Are you sure that you want to remove Local Colabfold from your system?")
            msg.setWindowTitle("Remove Local Colabfold")
            msg.setStandardButtons(QMessageBox.Yes)
            msg.setStandardButtons(QMessageBox.No)
            msg.exec_()
            # --- Answer ---
            # Here you will need a message box, which asks the user, if he would like to remove local colabfold from
            # his system.
            # The logical message of this if-statement is: the local colabfold is installed

            msg = QMessageBox()
            msg.setIcon(QMessageBox.Question)
            msg.setText("Do you want to remove Local Colabfold from your system and and remove the folder?")
            msg.setWindowTitle("Remove Local Colabfold")
            msg.setStandardButtons(QMessageBox.Yes)
            msg.setStandardButtons(QMessageBox.No)
            msg.exec_()
            # here you will need to implement a message box which asks the user if he wants to remove local colabfold
            # and evaluate the return value as it's done in line 139

            # ???
            # Is it right? Why we ask the same question in every message box?

            # subprocess.run(["wsl", "rm", "-r", "/home/$USER/.pyssa"])
            self.ui.btn_install_local_prediction.setText("Install")
        else:
            if dialog_message_local_colabfold.installation_local_colabfold_accept() is True:
                # logical message: the user wants to install local colabfold
                # substitute "pass" with the action which needs to be done for installing local colabfold
                pass
            else:
                # logical message: the user does NOT want to install local colabfold
                # substitute "pass" with the action which needs to be done for aborting the local colabfold installation
                pass

            if dialog_message_local_colabfold.installation_local_colabfold_progress() is True:
            # user_name = os.getlogin()
            # subprocess.run(["C:\\Windows\\System32\\WindowsPowerShell\\v1.0\\powershell.exe", pathlib.Path(f"{os.path.expanduser('~')}/github_repos/tmpPySSA/pyssa/scripts/convert_dos_to_unix.ps1")])
            # subprocess.run(["wsl", "mkdir", "/home/$USER/.pyssa"])
            # subprocess.run(
            #     ["wsl", f"/mnt/c/Users/{user_name}/github_repos/tmpPySSA/pyssa/scripts/installation_colabfold.sh"])
            # subprocess.run(
            #     ["wsl", "cd", "/home/$USER/.pyssa", "&&", "wget", "-q", "-P", ".",
            #      "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"])
            # subprocess.run(
            #     ["wsl", "cd", "/home/$USER/.pyssa", "&&", "./install_colabbatch_linux.sh"])
            # subprocess.run(["wsl", "cd", "/home/$USER/.pyssa", "&&", "./post_colabfold_installation.sh"])
            # subprocess.run(["wsl", "cd", "/home/$USER/.pyssa", "&&", "./update.sh"])
                pass
            else:
                pass

            if dialog_message_local_colabfold.installation_local_colabfold_end() is True:
                pass
            else:
                pass

            self.ui.btn_install_local_prediction.setText("Uninstall")

    def install_wsl(self):
        dialog = dialog_message_wsl.DialogMessageWsl()
        dialog.exec_()
