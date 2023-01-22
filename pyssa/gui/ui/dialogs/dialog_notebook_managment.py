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
from datetime import datetime
from pyssa.gui.ui.forms.auto_generated.auto_dialog_notebook_managment import Ui_Dialog
from util import gui_utils

global_var_startup_workspace = ""
global_var_terminate_app = 0


class DialogNotebookManagment(Qt.QtWidgets.QDialog):

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
        # vars
        self.web_page = None
        self.check_status = None
        self.interface = None
        # changes of some gui elements
        self.ui.lbl_title.setText("A prediction is currently running, you can check the status or abort.")
        self.ui.lbl_current_status.setText("")
        self.ui.btn_disable_gpu.hide()
        self.ui.btn_reconnect.hide()
        # connections
        self.ui.btn_disable_gpu.clicked.connect(self.disable_gpu)
        self.ui.btn_check_status.clicked.connect(self.check_current_status)
        self.ui.btn_abort.clicked.connect(self.close_interface)
        self.ui.btn_reconnect.clicked.connect(self.reconnect_session)
        self.setWindowTitle("Warning")


    # @SLOT
    def disable_gpu(self):
        self.web_page.runJavaScript("document.getElementsByTagName('paper-button').item(2).click()")
        self.ui.lbl_current_status.setText("Prediction runs without GPU acceleration.")
        self.ui.btn_disable_gpu.hide()
        self.ui.btn_check_status.setEnabled(True)

    def reconnect_session(self):
        self.web_page.runJavaScript("document.getElementsByTagName('paper-button').item(2).click()")
        self.ui.lbl_current_status.setText("Google Colab session will be reconnected.")
        self.ui.btn_reconnect.hide()
        self.ui.btn_check_status.setEnabled(True)

    def check_current_status(self):
        current_time = datetime.now().strftime("%H:%M:%S")
        if self.check_status() == 2:
            self.ui.lbl_current_status.setText(f"Predictions runs normal. Last checked: {current_time}")
        elif self.check_status() == 3:
            self.ui.lbl_current_status.setText(f"The GPU cannot be used, due to an account timeout. "
                                               f"You can still run the prediction without a GPU."
                                               f"Last checked: {current_time}")
            self.ui.btn_disable_gpu.show()
            self.ui.btn_check_status.setEnabled(False)
        elif self.check_status() == 4:
            self.ui.lbl_current_status.setText("Please reconnect your session.")
            self.ui.btn_reconnect.show()
            self.ui.btn_check_status.setEnabled(False)
        else:
            self.close_interface()

    def close_interface(self):
        self.interface.close()
        self.close()
        gui_utils.error_prediction_progress_lost()
