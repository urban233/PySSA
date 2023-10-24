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
from PyQt5 import QtGui
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
    """This class opens a settings customization dialog."""
    """This variable is for controlling whether the dialog opens or not"""
    ERROR = False

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

        # <editor-fold desc="Info button changes">
        pixmapi = QtWidgets.QStyle.SP_MessageBoxInformation
        icon = self.style().standardIcon(pixmapi)
        self.ui.btn_info.setIcon(icon)
        self.ui.btn_info.setText("")
        self.ui.btn_info.setFixedWidth(50)

        # </editor-fold>

        self.ui.label.hide()
        self.ui.lbl_color_vision_mode.hide()
        self.ui.cb_color_vision_mode.hide()

        # <editor-fold desc="Class attributes">
        self.tmp_settings = settings.Settings(constants.SETTINGS_DIR, constants.SETTINGS_FILENAME)
        try:
            self.settings = self.tmp_settings.deserialize_settings()
        except ValueError:
            logging.error("Settings dialog cannot be opened due to an error.")
            return

        # </editor-fold>

        # <editor-fold desc="Check WSL status">
        if is_wsl2_installed():
            self.settings.wsl_install = 1
        else:
            self.settings.wsl_install = 0

        # </editor-fold>

        # <editor-fold desc="Check local colabfold status">
        if is_local_colabfold_installed():
            self.settings.local_colabfold = 1
        else:
            self.settings.local_colabfold = 0

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
        item_list = [
            "normal",
            "Red-green (green weak, deuteranopia)",
            "Red-green (red weak, protanopia)",
            "Blue-yellow (tritanopia)",
        ]
        gui_utils.fill_combo_box(self.ui.cb_color_vision_mode, item_list)

        # </editor-fold>
        if self.settings.color_vision_mode == constants.CVM_NORMAL:
            self.ui.cb_color_vision_mode.setCurrentIndex(0)
        elif self.settings.color_vision_mode == constants.CVM_DEUTERANOPIA:
            self.ui.cb_color_vision_mode.setCurrentIndex(1)
        elif self.settings.color_vision_mode == constants.CVM_PROTANOPIA:
            self.ui.cb_color_vision_mode.setCurrentIndex(2)
        elif self.settings.color_vision_mode == constants.CVM_TRITANOPIA:
            self.ui.cb_color_vision_mode.setCurrentIndex(3)

        styles.set_stylesheet(self)
        self.setWindowIcon(PyQt5.QtGui.QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
        self.setWindowTitle("Global Settings")
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)

    # @SLOT()
    def choose_workspace_dir(self):
        gui_utils.choose_directory(self, self.ui.txt_workspace_dir)

    def set_color_vision_mode(self):
        pass

    def _connect_all_gui_elements(self):
        self.ui.cb_color_vision_mode.currentIndexChanged.connect(self.set_color_vision_mode)
        self.ui.btn_workspace_dir.clicked.connect(self.choose_workspace_dir)
        self.ui.btn_cancel.clicked.connect(self.cancel_dialog)
        self.ui.btn_ok.clicked.connect(self.ok_dialog)
        self.ui.btn_info.clicked.connect(self.open_page_information)

    def cancel_dialog(self):
        self.close()

    def ok_dialog(self):
        self.settings.set_workspace_path(self.ui.txt_workspace_dir.text())
        self.settings.set_cycles(str(self.ui.spb_cycles.value()))
        self.settings.set_cutoff(str(self.ui.dspb_cutoff.value()))
        self.settings.color_vision_mode = self.ui.cb_color_vision_mode.currentText()
        self.settings.serialize_settings()
        logging.info("Settings were successfully saved.")
        self.close()

    def open_page_information(self) -> None:
        """Opens the message box, to display extra information based on the page."""
        msg = QtWidgets.QMessageBox()
        msg.setIcon(QtWidgets.QMessageBox.Information)
        msg.setWindowIcon(QtGui.QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
        styles.set_stylesheet(msg)
        msg.setWindowTitle("Information")
        msg.setStyleSheet("QLabel{font-size: 11pt;}")

        msg.setText(
            "Global Settings\n\nCycles: Maximum number of outlier rejection cycles\n\n"
            "Cutoff: Outlier rejection cutoff for sequence alignment\n\n"
            "Note: Cutoff value is neglected unless numbers of cycles are more than 0.",
        )
        msg.exec_()
        return
