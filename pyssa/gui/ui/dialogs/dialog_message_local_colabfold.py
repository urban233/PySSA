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

import PyQt5.QtWidgets
from pymol import Qt
from datetime import datetime
from pyssa.gui.ui.forms.auto_generated.auto_dialog_message_local_colabfold import Ui_Dialog
from pyssa.gui.utilities import gui_utils
import subprocess
import os
import sys
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QMessageBox
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import pyqtSlot


def installation_local_colabfold_accept() -> bool:
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Question)
    msg.setText("Are you sure that you want to install Local Colabfold?")
    msg.setWindowTitle("Local Colabfold installation")
    btn_installation_local_colabfold_accept_yes = msg.addButton("Yes", QMessageBox.ActionRole)
    btn_installation_local_colabfold_accept_no = msg.addButton("No", QMessageBox.ActionRole)
    msg.exec_()

     # def def get_a_random_number() -> int:
     #   ilcp = installation_local_colabfold_progress()
     #  return ???
    # with this line, you will save the return value from "installation_local_colabfold_progress" into ilcp but
    # the function does not have any return value
    # Extra tip: annotate your function heads like this: def get_a_random_number() -> int:
    # this helps to identify if a function has a return value or not

    # ??? Question ???
    # I think, that I don't understand to 100% how this function works. So I understand, that we need a return value, but when the return value must be an int,
    # it don't can be True or False. But when I took e.g. a 1 as a return value, a lot of mistakes appears. I play with the function, but I don't find a solution.
    # So I need your help again. Thanks a lot for your patience!!!
    # --- Answer ---
    # The function "get_a_random_number()" was just an example to show you how you can annotate your function heads.
    # You do NOT need this function here or elsewhere.
    # For better visualization of function annotation, I added the annotation
    # to the function "installation_local_colabfold_accept()"

    if msg.clickedButton() == btn_installation_local_colabfold_accept_yes:
        return True
        # you can better return a True, which gets evaluated in the dialog_settings_global.py through an
        # if/else statement
    else:
        msg.close()
        return False
        # closing the message box is the right way, but you are missing a return value like False


def installation_local_colabfold_progress():
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Warning)
    msg.setText("Don't close the window and wait!")
    msg.setWindowTitle("Local Colabfold installation")
    # fixme: hide msg.setStandardButtons -> read source from Todoist
    msg.exec_()


def installation_local_colabfold_end():
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Information)
    msg.setText("Installation is finished!")
    msg.setWindowTitle("Local Colabfold installation")
    msg.setStandardButtons(QMessageBox.Ok)


# garbage
# class DialogMessageLocalColabfold(QWidget):
#
#     def __init__(self):
#         super().__init__()
#         self.title = 'PyQt5 messagebox - pythonspot.com'
#         self.left = 10
#         self.top = 10
#         self.width = 320
#         self.height = 200
#         self.initUI()
#
#     def initUI(self):
#         self.setWindowTitle(self.title)
#         self.setGeometry(self.left, self.top, self.width, self.height)
#
#         askInstallation = QMessageBox.question(self, 'Local Colabfold installation', "Are you sure that you want the Local Colabfold installation?",
#                                            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
#
#         installationProgress = QMessageBox.warning(self, 'Local Colabfold installation',  "Don't close the window and wait!")
#
#         installationEnd = QMessageBox.information(self, 'Local Colabfold installation', "Installation is finished!",
#                                            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
#
#         if askInstallation == QMessageBox.Yes:
#             # installation begins
#             print(installationProgress)
#             # user_name = os.getlogin()
#             # subprocess.run(["C:\\Windows\\System32\\WindowsPowerShell\\v1.0\\powershell.exe",pathlib.Path(f"{os.path.expanduser('~')}/github_repos/tmpPySSA/pyssa/scripts/convert_dos_to_unix.ps1")])
#             # subprocess.run(["wsl", "mkdir", "/home/$USER/.pyssa"])
#             # subprocess.run(
#             #     ["wsl", f"/mnt/c/Users/{user_name}/github_repos/tmpPySSA/pyssa/scripts/installation_colabfold.sh"])
#             # subprocess.run(
#             #     ["wsl", "cd", "/home/$USER/.pyssa", "&&", "wget", "-q", "-P", ".",
#             #      "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"])
#             # subprocess.run(
#             #     ["wsl", "cd", "/home/$USER/.pyssa", "&&", "./install_colabbatch_linux.sh"])
#             # subprocess.run(["wsl", "cd", "/home/$USER/.pyssa", "&&", "./post_colabfold_installation.sh"])
#             # subprocess.run(["wsl", "cd", "/home/$USER/.pyssa", "&&", "./update.sh"])
#
#             # installation ends
#             print(installationEnd)
#
#
#         else:
#             sys.exit()


# if __name__ == '__main__':
#     app = QApplication(sys.argv)
#     ex = DialogMessageLocalColabfold()
#     sys.exit(app.exec_())

# more btns
# buttonReply = QMessageBox.question(self, 'PyQt5 message', "Do you want to save?", QMessageBox.Yes | QMessageBox.No | QMessageBox.Cancel, QMessageBox.Cancel)
# print(int(buttonReply))
# if buttonReply == QMessageBox.Yes:
#     print('Yes clicked.')
# if buttonReply == QMessageBox.No:
#     print('No clicked.')
# if buttonReply == QMessageBox.Cancel:
#     print('Cancel')
# class DialogMessageLocalColabfold(Qt.QtWidgets.QDialog):
#     def __init__(self, parent=None):
#         """Constructor
#
#         Args:
#             args
#             kwargs
#         """
#         Qt.QtWidgets.QDialog.__init__(self, parent)
#         # build ui object
#         self.ui = Ui_Dialog()
#         self.ui.setupUi(self)
#
#         self.setWindowTitle("Local Colabfold installation")
#         self.ui.lbl_message_localcolabfold.setText("Are you sure that you want the Local Colabfold installation?")
#         # btn
#         self.ui.btn_message_localcolabfold_ok.show()
#         self.ui.btn_message_localcolabfold_cancel.show()
#         self.ui.btn_message_localcolabfold_restart.hide()
#         self.ui.btn_message_localcolabfold_restart_later.hide()
#         self.ui.btn_message_localcolabfold_ok_2.hide()
#
#         # btn connections
#         self.ui.btn_message_localcolabfold_ok.clicked.connect(self.installation_localcolabfold)
#         self.ui.btn_message_localcolabfold_cancel.clicked.connect(self.cancel_installation)
#         self.ui.btn_message_localcolabfold_ok_2.clicked.connect(self.close_dlg_installation_interface)
#
#     def installation_localcolabfold(self):
#         # installation is started
#         self.ui.lbl_message_localcolabfold.setText("Don't close the window and wait!")
#         self.ui.btn_message_localcolabfold_ok.hide()
#         self.ui.btn_message_localcolabfold_ok_2.hide()
#         self.ui.btn_message_localcolabfold_cancel.hide()
#         self.ui.btn_message_localcolabfold_restart.hide()
#         self.ui.btn_message_localcolabfold_restart_later.hide()
#         # user_name = os.getlogin()
#         # subprocess.run(["C:\\Windows\\System32\\WindowsPowerShell\\v1.0\\powershell.exe",pathlib.Path(f"{os.path.expanduser('~')}/github_repos/tmpPySSA/pyssa/scripts/convert_dos_to_unix.ps1")])
#         # subprocess.run(["wsl", "mkdir", "/home/$USER/.pyssa"])
#         # subprocess.run(
#         #     ["wsl", f"/mnt/c/Users/{user_name}/github_repos/tmpPySSA/pyssa/scripts/installation_colabfold.sh"])
#         # subprocess.run(
#         #     ["wsl", "cd", "/home/$USER/.pyssa", "&&", "wget", "-q", "-P", ".",
#         #      "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"])
#         # subprocess.run(
#         #     ["wsl", "cd", "/home/$USER/.pyssa", "&&", "./install_colabbatch_linux.sh"])
#         # subprocess.run(["wsl", "cd", "/home/$USER/.pyssa", "&&", "./post_colabfold_installation.sh"])
#         # subprocess.run(["wsl", "cd", "/home/$USER/.pyssa", "&&", "./update.sh"])
#
#         # installation is finished
#         # self.ui.lbl_message_localcolabfold.setText("Installation is finished!")
#         # self.ui.btn_message_localcolabfold_ok.hide()
#         # self.ui.btn_message_localcolabfold_ok_2.show()
#         # self.ui.btn_message_localcolabfold_cancel.hide()
#         # self.ui.btn_message_localcolabfold_restart.hide()
#         # self.ui.btn_message_localcolabfold_restart_later.hide()
#
#     def cancel_installation(self):
#         self.close()
#
#     # def installation_is_finished(self):
#
#
#     def close_dlg_installation_interface(self):
#         self.close()
#
#     def restart_later(self):
#         self.close()
#
#     # def restart_system(self):
#     #     os.system("shutdown /r")
#
#     def install_localcolabfold_command(self):
#         time.sleep(5)