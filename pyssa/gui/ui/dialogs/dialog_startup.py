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
import pathlib
import sys
import zipfile
import PyQt5.QtWidgets
from urllib import request
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication
from pyssa.util import constants
from pyssa.gui.ui.forms.auto_generated.auto_dialog_startup import Ui_Dialog

global_var_startup_workspace = ""
global_var_terminate_app = 0


class DialogStartup(QtWidgets.QDialog):

    def __init__(self, parent=None):
        """Constructor.

        Args:
            args
            kwargs
        """
        QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        self.workspace_dir = ""

        self.ui.btn_cancel.clicked.connect(self.close_dialog)
        self.ui.btn_choose_workspace.clicked.connect(self.choose_workspace)
        self.ui.btn_launch.clicked.connect(self.launch_app)

        path_list = [
            f"{os.path.expanduser('~')}/Documents",
            f"{os.path.expanduser('~')}/Documents",
            f"{os.path.expanduser('~')}\\Documents",
        ]
        # appends the os specific python path
        if sys.platform.startswith("darwin"):
            # macOS path
            self.ui.txt_workspace.setText(path_list[1])
        elif sys.platform.startswith("linux"):
            # Linux path
            self.ui.txt_workspace.setText(path_list[0])
        elif sys.platform.startswith("win32"):
            # Windows path
            self.ui.txt_workspace.setText(str(constants.DEFAULT_WORKSPACE_PATH))
        self.setWindowFlag(QtCore.Qt.WindowCloseButtonHint, True)
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
        self.setWindowTitle("PySSA Launch")

    # @SLOT
    def close_dialog(self):
        global global_var_terminate_app
        global_var_terminate_app = 1
        self.close()

    def choose_workspace(self):
        """This function opens a file dialog to choose a workspace directory."""
        self.workspace_dir = PyQt5.QtWidgets.QFileDialog.getExistingDirectory(
            self, "Open Workspace Directory", QtCore.QDir.homePath(),
            PyQt5.QtWidgets.QFileDialog.ShowDirsOnly | PyQt5.QtWidgets.QFileDialog.DontResolveSymlinks)
        if self.workspace_dir != "":
            self.ui.txt_workspace.setText(self.workspace_dir)

    def launch_app(self):
        """This function launches the plugin."""
        QApplication.setOverrideCursor(Qt.WaitCursor)
        global global_var_startup_workspace
        global_var_startup_workspace = self.ui.txt_workspace.text()

        url = constants.UNIX_SCRIPTS_SCIEBO_URL  # Replace with the actual file URL
        destination = pathlib.Path(f"{constants.SETTINGS_DIR}/unix.zip")  # Replace with the desired file path and name
        if not os.path.exists(pathlib.Path(constants.SETTINGS_DIR)):
            os.mkdir(pathlib.Path(constants.SETTINGS_DIR))
        if not os.path.exists(pathlib.Path(f"{constants.SETTINGS_DIR}/scripts")):
            os.mkdir(pathlib.Path(f"{constants.SETTINGS_DIR}/scripts"))
        if not os.path.exists(pathlib.Path(f"{constants.SETTINGS_DIR}/scripts/unix")):
            os.mkdir(pathlib.Path(f"{constants.SETTINGS_DIR}/scripts/unix"))

        request.urlretrieve(url, str(destination))

        zip_file_path = str(destination)  # Replace with the actual path to your downloaded zip file
        destination_folder = str(pathlib.Path(f"{constants.SETTINGS_DIR}/scripts"))  # Replace with the desired folder to extract the contents
        with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
            zip_ref.extractall(destination_folder)

        if not os.path.exists(pathlib.Path(f"{constants.SCRATCH_DIR}")):
            os.mkdir(pathlib.Path(f"{constants.SCRATCH_DIR}"))
        if not os.path.exists(pathlib.Path(f"{constants.CACHE_DIR}")):
            os.mkdir(pathlib.Path(f"{constants.CACHE_DIR}"))
        request.urlretrieve(constants.DEMO_PROJECT_SCIEBO_URL, str(pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects.zip")))

        self.close()
