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
from pymol import Qt
from pyssa.gui.ui.forms.auto_generated.auto_dialog_message_wsl import Ui_Dialog
from pyssa.util import gui_utils


class DialogMessageWsl(Qt.QtWidgets.QDialog):

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
        # self.status = ""

        self.setWindowTitle("WSL2 installation")
        self.ui.lbl_message_wsl.setText("Are you sure that you want the WSL2 installation?")
        # btn
        self.ui.btn_message_wsl_ok.show()
        self.ui.btn_message_wsl_cancel.show()
        self.ui.btn_message_wsl_restart.hide()
        self.ui.btn_message_wsl_restart_later.hide()

        # btn connections
        self.ui.btn_message_wsl_ok.clicked.connect(self.installation_wsl)
        self.ui.btn_message_wsl_cancel.clicked.connect(self.cancel_installation)
        self.ui.btn_message_wsl_restart_later.clicked.connect(self.restart_later)
        self.ui.btn_message_wsl_restart.clicked.connect(self.restart_system)

    def installation_wsl (self):
        # installation is started
        self.ui.lbl_message_wsl.setText("Don't close the window and wait!")
        gui_elements_to_hide = [
            self.ui.btn_message_wsl_ok,
            self.ui.btn_message_wsl_cancel,
            self.ui.btn_message_wsl_restart,
            self.ui.btn_message_wsl_restart_later,
        ]
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        print("go_finished")
        # subprocess.run("wsl --install")

        # installation is finished
        # self.ui.lbl_message_wsl.setText("Installation is finished! Restart is necessary.")
        # gui_elements_to_hide = [
        #     self.ui.btn_message_wsl_ok,
        #     self.ui.btn_message_wsl_cancel,
        # ]
        # gui_utils.hide_gui_elements(gui_elements_to_hide)
        # self.ui.btn_message_wsl_restart.show()
        # self.ui.btn_message_wsl_restart_later.show()

    def cancel_installation(self):
        self.close()

    # def installation_is_finished(self):
    #     self.ui.lbl_message_wsl.setText("Installation is finished! Restart is necessary.")
    #     self.ui.btn_message_wsl_ok.hide()
    #     self.ui.btn_message_wsl_cancel.hide()
    #     self.ui.btn_message_wsl_restart.show()
    #     self.ui.btn_message_wsl_restart_later.show()

    def close_dlg_installation_interface(self):
        self.close()

    def restart_later(self):
        self.close()

    def restart_system(self):
        os.system("shutdown /r")
