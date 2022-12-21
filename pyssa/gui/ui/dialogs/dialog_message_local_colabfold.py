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
import PyQt5.QtWidgets
from pymol import Qt
from datetime import datetime
from pyssa.gui.ui.forms.auto_generated.auto_dialog_message_local_colabfold import Ui_Dialog
from pyssa.gui.utilities import gui_utils
import subprocess
import os

class DialogMessageLocalColabfold(Qt.QtWidgets.QDialog):
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

        self.setWindowTitle("Local Colabfold installation")
        self.ui.lbl_message_localcolabfold.setText("Are you sure that you want the Local Colabfold installation?")
        # btn
        self.ui.btn_message_localcolabfold_ok.show()
        self.ui.btn_message_localcolabfold_cancel.show()
        self.ui.btn_message_localcolabfold_restart.hide()
        self.ui.btn_message_localcolabfold_restart_later.hide()
        self.ui.btn_message_localcolabfold_ok_2.hide()

        # btn connections
        self.ui.btn_message_localcolabfold_ok.clicked.connect(self.installation_in_progress)
        self.ui.btn_message_localcolabfold_cancel.clicked.connect(self.cancel_installation)
        self.ui.btn_message_localcolabfold_ok.clicked.connect(self.installation_is_finished)
        self.ui.btn_message_localcolabfold_ok_2.clicked.connect(self.close_dlg_installation_interface)
        # self.ui.btn_message_localcolabfold_restart_later.clicked.connect(self.restart_later)
        # self.ui.btn_message_localcolabfold_restart.clicked.connect(self.restart_system)

    def installation_in_progress(self):
        self.ui.lbl_message_localcolabfold.setText("Don't close the window and wait!")
        self.ui.btn_message_localcolabfold_ok.hide()
        self.ui.btn_message_localcolabfold_ok_2.hide()
        self.ui.btn_message_localcolabfold_cancel.hide()
        self.ui.btn_message_localcolabfold_restart.hide()
        self.ui.btn_message_localcolabfold_restart_later.hide()
        user_name = os.getlogin()
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

    def cancel_installation(self):
        self.close()

    def installation_is_finished(self):
        self.ui.lbl_message_localcolabfold.setText("Installation is finished!")
        self.ui.btn_message_localcolabfold_ok.hide()
        self.ui.btn_message_localcolabfold_ok_2.show()
        self.ui.btn_message_localcolabfold_cancel.hide()
        self.ui.btn_message_localcolabfold_restart.hide()
        self.ui.btn_message_localcolabfold_restart_later.hide()

    def close_dlg_installation_interface(self):
        self.close()

    def restart_later(self):
        self.close()

    # def restart_system(self):
    #     os.system("shutdown /r")