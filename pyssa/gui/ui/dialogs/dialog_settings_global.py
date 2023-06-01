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
import os
import logging
import subprocess
from pyssa.gui.ui.forms.auto_generated.auto_dialog_settings_global import Ui_Dialog
from pyssa.internal.data_structures import settings
from pyssa.util import constants, gui_utils
from pyssa.gui.ui.styles import styles
from pyssa.internal.thread import workers
from PyQt5.QtWidgets import QMessageBox
import PyQt5
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.gui.ui.messageboxes import settings_boxes

# setup logger
logging.basicConfig(level=logging.DEBUG)


def is_wsl2_installed():
    output = subprocess.run("wsl --list --verbose")
    if output.returncode == 0:
        return True
    else:
        return False


def is_local_colabfold_installed():
    return os.path.exists(constants.WSL_DISK_PATH)


class DialogSettingsGlobal(QtWidgets.QDialog):
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
        QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)

        # <editor-fold desc="Class attributes">
        self.tmp_settings = settings.Settings(constants.SETTINGS_DIR, constants.SETTINGS_FILENAME)
        try:
            self.settings = self.tmp_settings.deserialize_settings()
        except ValueError:
            logging.error("Settings dialog cannot be opened due to an error.")
            return
        self.threadpool = QtCore.QThreadPool()
        self.colabfold_installer_worker = workers.ColabfoldInstallerWorkerPool(True)
        self.wsl_installer_worker = workers.WslInstallerWorkerPool(True)
        self.block_box_colabfold = basic_boxes.no_buttons("Process running", "A process is currently running.", QMessageBox.Information)
        self.block_box_wsl = basic_boxes.no_buttons("Process running", "A process is currently running.", QMessageBox.Information)

        # </editor-fold>

        # <editor-fold desc="Check WSL status">
        if is_wsl2_installed():
            self.settings.wsl_install = 1
            self.ui.btn_install_wsl2.setText("Uninstall")
            self.ui.btn_install_local_prediction.setEnabled(True)
        else:
            self.settings.wsl_install = 0
            self.ui.btn_install_wsl2.setText("Install")
            self.ui.btn_install_local_prediction.setEnabled(False)

        # </editor-fold>

        # <editor-fold desc="Check local colabfold status">
        if is_local_colabfold_installed():
            self.settings.local_colabfold = 1
            self.ui.btn_install_local_prediction.setText("Uninstall")
        else:
            self.settings.local_colabfold = 0
            self.ui.btn_install_local_prediction.setText("Install")

        # </editor-fold>

        self._connect_all_gui_elements()

        # <editor-fold desc="Set up defaults">
        self.ui.txt_workspace_dir.setEnabled(False)
        self.ui.txt_workspace_dir.setText(str(self.settings.get_workspace_path()))
        self.ui.spb_cycles.setValue(int(self.settings.get_cycles()))
        self.ui.dspb_cutoff.setValue(float(self.settings.get_cutoff()))
        # customize spin boxes
        self.ui.spb_cycles.setMinimum(0)
        # self.ui.spb_cycles.setMaximum(20) # fixme: is a maximum needed?
        self.ui.spb_cycles.setSingleStep(1)
        self.ui.dspb_cutoff.setMinimum(0.00)
        self.ui.dspb_cutoff.setMaximum(20.00)
        self.ui.dspb_cutoff.setSingleStep(0.1)

        # </editor-fold>

        styles.set_stylesheet(self)
        self.setWindowIcon(PyQt5.QtGui.QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
        self.setWindowTitle("Global Settings")

    # @SLOT()
    def choose_workspace_dir(self):
        gui_utils.choose_directory(self, self.ui.txt_workspace_dir)

    def _connect_all_gui_elements(self):
        self.ui.btn_workspace_dir.clicked.connect(self.choose_workspace_dir)
        self.ui.btn_cancel.clicked.connect(self.cancel_dialog)
        self.ui.btn_ok.clicked.connect(self.ok_dialog)
        self.ui.btn_install_local_prediction.clicked.connect(self.install_local_colabfold)
        self.ui.btn_install_wsl2.clicked.connect(self.install_wsl)
        # thread-related
        self.colabfold_installer_worker.signals.finished.connect(self.post_colabfold_install)
        self.wsl_installer_worker.signals.finished.connect(self.post_wsl_install)

    def cancel_dialog(self):
        self.close()

    def ok_dialog(self):
        self.settings.set_workspace_path(self.ui.txt_workspace_dir.text())
        #self.settings.set_prediction_path(self.ui.txt_zip_storage_dir.text())
        self.settings.set_cycles(str(self.ui.spb_cycles.value()))
        self.settings.set_cutoff(str(self.ui.dspb_cutoff.value()))
        self.settings.serialize_settings()
        logging.info("Settings were successfully saved.")
        self.close()

    def post_colabfold_install(self):
        self.block_box_colabfold.destroy(True)
        if self.colabfold_installer_worker.install is False:
            self.ui.btn_install_local_prediction.setText("Install")
            self.settings.local_colabfold = 0
            self.settings.serialize_settings()
            basic_boxes.ok("Local Colabfold", "Removing local colabfold is finished!", QMessageBox.Information)
        else:
            self.ui.btn_install_local_prediction.setText("Uninstall")
            self.settings.local_colabfold = 1
            self.settings.serialize_settings()
            basic_boxes.ok("Local Colabfold installation", "Installation is finished!", QMessageBox.Information)
            if settings_boxes.restart_now_later():
                os.system("shutdown /r")

    def install_local_colabfold(self):
        if self.settings.local_colabfold == 1:
            # colabfold installed on system, user wants to uninstall local colabfold
            if basic_boxes.yes_or_no("Remove Local Colabfold", "Are you sure that you want to remove Local Colabfold from your system?", QMessageBox.Question):
                self.colabfold_installer_worker.install = False
                self.threadpool.start(self.colabfold_installer_worker)
                self.block_box_colabfold.exec_()
            else:
                return
        else:
            if basic_boxes.yes_or_no("Local Colabfold installation", "Are you sure that you want to install Local Colabfold?", QMessageBox.Question) is True:
                self.colabfold_installer_worker.install = True
                self.threadpool.start(self.colabfold_installer_worker)
                self.block_box_colabfold.exec_()
            else:
                # logical message: the user does NOT want to install local colabfold
                basic_boxes.ok("Local Colabfold installation", "Installation process aborted.", QMessageBox.Information)
                return

    def post_wsl_install(self):
        self.block_box_wsl.destroy(True)
        if self.wsl_installer_worker.install is False:
            #self.ui.btn_install_wsl2.setText("Install")
            self.settings.wsl_install = 0
        else:
            basic_boxes.ok("WSL2", "Installation of WSL2 is finished!", QMessageBox.Information)
            self.ui.btn_install_wsl2.setText("Uninstall")
            self.settings.wsl_install = 1
            self.settings.serialize_settings()
            if settings_boxes.restart_now_later():
                os.system("shutdown /r")

    def install_wsl(self):
        if self.settings.wsl_install == 1:
            # WSL is installed on system, user wants to uninstall WSL2
            if basic_boxes.yes_or_no("Remove WSL2", "Are you sure that you want to remove WSL2 from your system?", QMessageBox.Question):
                basic_boxes.ok("Uninstall WSL2", "Please refer to the documentation under Help, to uninstall the WSL2.",
                               QMessageBox.Information)
                #self.wsl_installer_worker.install = False
                #self.threadpool.start(self.wsl_installer_worker)
                #self.block_box_wsl.exec_()
            else:
                return
        elif self.settings.wsl_install == 0:
            if basic_boxes.yes_or_no("WSL2 installation", "Are you sure that you want to install WSL2?", QMessageBox.Question) is True:
                # the user wants to install WSL2
                self.wsl_installer_worker.install = True
                self.threadpool.start(self.wsl_installer_worker)
                self.block_box_wsl.exec_()
            else:
                # the user does NOT want to install WSL
                basic_boxes.ok("WSL2 installation", "Installation process aborted.", QMessageBox.Information)
                return
