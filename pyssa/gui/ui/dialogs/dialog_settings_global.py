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
import time

from pyssa.gui.utilities import gui_utils
from pyssa.gui.utilities import styles
from pyssa.gui.ui.forms.auto_generated.auto_dialog_settings_global import Ui_Dialog
from pyssa.gui.data_structures import settings
from pyssa.gui.utilities import constants
from pymol import Qt
from PyQt5.QtWidgets import QMessageBox
import PyQt5
from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.gui.ui.messageboxes import settings_boxes

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
        logging.info("Loading values from settings.json was successful.")
        self.ui.txt_workspace_dir.setEnabled(False)
        self.ui.txt_zip_storage_dir.setEnabled(False)

        logging.info("Settings dialog was opened.")
        # loading information from the settings.json
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

        if self.settings.wsl_install == 0:
            # wsl is not installed
            self.ui.btn_install_wsl2.setText("Install")
        elif self.settings.wsl_install == 1:
            # wsl is installed
            self.ui.btn_install_wsl2.setText("Uninstall")
        else:
            gui_utils.error_dialog_settings("The settings are corrupted, please restore the settings!", "", self.settings)
        if self.settings.local_colabfold == 0:
            # local colabfold is not installed
            self.ui.btn_install_local_prediction.setText("Install")
        elif self.settings.local_colabfold == 1:
            # local colabfold is installed
            self.ui.btn_install_local_prediction.setText("Uninstall")
        else:
            gui_utils.error_dialog_settings("The settings are corrupted, please restore the settings!", "", self.settings)

        # connect elements with function
        self.ui.btn_workspace_dir.clicked.connect(self.chooseWorkspaceDir)
        self.ui.btn_zip_storage_dir.clicked.connect(self.chooseZipStorageDir)
        self.ui.btn_cancel.clicked.connect(self.cancelDialog)
        self.ui.btn_ok.clicked.connect(self.okDialog)
        self.ui.btn_install_local_prediction.clicked.connect(self.install_local_colabfold)
        self.ui.btn_install_wsl2.clicked.connect(self.install_wsl)
        styles.set_stylesheet(self)
        self.setWindowIcon(PyQt5.QtGui.QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
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
        if self.settings.local_colabfold == 1:
            # colabfold installed on system, user wants to uninstall local colabfold
            if basic_boxes.yes_or_no("Remove Local Colabfold", "Are you sure that you want to remove Local Colabfold from your system?", QMessageBox.Question):
                try:
                    subprocess.run(["C:\\Windows\\System32\\WindowsPowerShell\\v1.0\\powershell.exe", pathlib.Path(
                        f"{os.path.expanduser('~')}/github_repos/tmpPySSA/pyssa/scripts/remove_wsl_env.ps1")])
                except:
                    basic_boxes.ok("Local Colabfold removal", "The uninstallation failed. Please re-run the process or consult the documentation.",
                                   QMessageBox.Critical)
                    return
                self.ui.btn_install_local_prediction.setText("Install")
                self.settings.local_colabfold = 0
            else:
                return
        else:
            if basic_boxes.yes_or_no("Local Colabfold installation", "Are you sure that you want to install Local Colabfold?", QMessageBox.Question) is True:
                # logical message: the user wants to install local colabfold
                user_name = os.getlogin()
                try:
                    if not os.path.exists(pathlib.Path(f"C:/Users/{os.getlogin()}/.pyssa/wsl/")):
                        os.mkdir(pathlib.Path(f"C:/Users/{os.getlogin()}/.pyssa/wsl/"))
                    if not os.path.exists(constants.WSL_STORAGE_PATH):
                        os.mkdir(constants.WSL_STORAGE_PATH)
                    subprocess.run(["C:\\Windows\\System32\\WindowsPowerShell\\v1.0\\powershell.exe", pathlib.Path(
                        f"{os.path.expanduser('~')}/github_repos/tmpPySSA/pyssa/scripts/convert_dos_to_unix.ps1")])
                    subprocess.run(["wsl", "--import", constants.WSL_DISTRO_NAME, str(constants.WSL_STORAGE_PATH), str(constants.WSL_DISTRO_IMPORT_PATH)])
                    subprocess.run(["wsl", "-s", constants.WSL_DISTRO_NAME])
                    subprocess.run(
                        ["wsl",
                         "cp", f"/mnt/c/Users/{user_name}/github_repos/tmpPySSA/pyssa/wsl_extras/wsl.conf", "/etc"])
                    subprocess.run(["wsl", "--shutdown"])
                    print("Shutting down all WSL distros ...")
                    time.sleep(10)
                    subprocess.run(["wsl", "mkdir", "/home/$USER/.pyssa"])
                    subprocess.run(
                        ["wsl", f"/mnt/c/Users/{user_name}/github_repos/tmpPySSA/pyssa/scripts/installation_colabfold.sh"])
                    subprocess.run(
                        ["wsl", "cd", "/home/$USER/.pyssa", "&&", "wget", "-q", "-P", ".",
                         "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"])
                    subprocess.run(
                        ["wsl", "cd", "/home/$USER/.pyssa", "&&", "./install_colabbatch_linux.sh"])
                    subprocess.run(["wsl", "cd", "/home/$USER/.pyssa", "&&", "./post_colabfold_installation.sh"])
                    subprocess.run(["wsl", "cd", "/home/$USER/.pyssa", "&&", "./update.sh"])
                except:
                    # this subprocess cleans the installation directory in case of an error
                    subprocess.run(["wsl", "rm", "-r", "/home/$USER/.pyssa"])
                    basic_boxes.ok("Local Colabfold installation", "Installation failed. Please re-run the process.", QMessageBox.Critical)
                    return
                self.ui.btn_install_local_prediction.setText("Uninstall")
                self.settings.local_colabfold = 1
                basic_boxes.ok("Local Colabfold installation", "Installation is finished!", QMessageBox.Information)
            else:
                # logical message: the user does NOT want to install local colabfold
                basic_boxes.ok("Local Colabfold installation", "Installation process aborted.", QMessageBox.Information)
                return

    def install_wsl(self):
        if self.settings.wsl_install == 1:
            # WSL is installed on system, user wants to uninstall WSL2
            if basic_boxes.yes_or_no("Remove WSL2", "Are you sure that you want to remove WSL2 from your system?", QMessageBox.Question):
                # TODO: Open the documentation page, because the WSL2 does require more user input than the installation!
                self.ui.btn_install_wsl2.setText("Install")
                self.settings.wsl_install = 0
            else:
                return
        elif self.settings.wsl_install == 0:
            if basic_boxes.yes_or_no("WSL2 installation", "Are you sure that you want to install WSL2?", QMessageBox.Question) is True:
                # the user wants to install WSL2
                try:
                    subprocess.run("wsl --install")
                except:
                    basic_boxes.ok("WSL2 installation", "Installation failed. Please re-run the process or look in the documentation.", QMessageBox.Critical)
                    return
                self.ui.btn_install_wsl2.setText("Uninstall")
                self.settings.wsl_install = 1
                if settings_boxes.restart_now_later():
                    os.system("shutdown /r")
            else:
                # the user does NOT want to install WSL
                basic_boxes.ok("WSL2 installation", "Installation process aborted.", QMessageBox.Information)
                return
        else:
            # unexpected case
            pass
