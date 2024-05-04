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
"""Module for the global settings dialog."""
import logging
import subprocess
from pyssa.gui.ui.forms.auto_generated.auto_dialog_settings_global import Ui_Dialog
from pyssa.internal.data_structures import settings
from pyssa.internal.thread import tasks
from pyssa.internal.thread.async_pyssa import util_async
from pyssa.util import constants, gui_utils
from pyssa.gui.ui.styles import styles
from pyssa.util import globals
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from PyQt5 import QtGui

# setup logger
logging.basicConfig(level=logging.DEBUG)


def is_wsl2_installed() -> bool:
    """Checks if the WSL2 is installed."""
    output = subprocess.run("wsl --list --verbose")
    if output.returncode == 0:
        return True
    return False


def is_local_colabfold_installed() -> bool:
    """Checks if the local colabfold is installed."""
    powershell_results = subprocess.run(["wsl", "-d", "almaColabfold9", "ls"])
    if powershell_results.returncode == 0:
        subprocess.run(["wsl", "--shutdown"])
        return True
    subprocess.run(["wsl", "--shutdown"])
    return False


class DialogSettingsGlobal(QtWidgets.QDialog):
    """This class opens a settings customization dialog."""

    """This variable is for controlling whether the dialog opens or not"""
    ERROR = False

    def __init__(self, the_interface_manager, parent=None) -> None:  # noqa: ANN001
        """Constructor.

        Args:
            parent: the parent
        """
        QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        self._interface_manager = the_interface_manager

        # <editor-fold desc="Info button changes">
        self.ui.btn_help.setIcon(QtGui.QIcon(":/icons/help_w200.png"))
        self.ui.btn_help.setIconSize(self.ui.btn_help.icon().actualSize(QtCore.QSize(30, 30)))
        self.ui.btn_help.setText("")

        # </editor-fold>

        # <editor-fold desc="Class attributes">
        self.tmp_settings = settings.Settings(constants.SETTINGS_DIR, constants.SETTINGS_FILENAME)
        try:
            self.settings = self.tmp_settings.deserialize_settings()
        except ValueError:
            logging.error("Settings dialog cannot be opened due to an error.")
            return

        # </editor-fold>

        # # <editor-fold desc="Check WSL status">
        # if globals.g_os == "win32":
        #     if is_wsl2_installed():
        #         self.settings.wsl_install = 1
        #     else:
        #         self.settings.wsl_install = 0
        # else:
        #     self.settings.wsl_install = 1
        # # </editor-fold>
        #
        # # <editor-fold desc="Check local colabfold status">
        # if is_local_colabfold_installed():
        #     self.settings.local_colabfold = 1
        # else:
        #     self.settings.local_colabfold = 0
        #
        # # </editor-fold>

        self._connect_all_gui_elements()

        self._initialize_ui()
        self.setMinimumWidth(450)
        self.resize(450, 450)
        styles.set_stylesheet(self)
        styles.color_bottom_frame_button(self.ui.btn_ok)
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        self.setWindowTitle("Global Settings")
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)

    def _initialize_ui(self):
        # combo box BgColor
        item_list_bg_color = [
            "black",
            "white",
        ]
        gui_utils.fill_combo_box(self.ui.box_bg_color, item_list_bg_color)
        # combo box Renderer
        item_list_renderer = [
            "default renderer",
            "PyMOL internal renderer",
        ]
        gui_utils.fill_combo_box(self.ui.box_renderer, item_list_renderer)
        # combo box RayTraceMode
        item_list_ray_trace_mode = [
            "normal color",
            "normal color + black outline",
            "black outline only",
            "quantized color + black outline",
        ]
        gui_utils.fill_combo_box(self.ui.box_ray_trace_mode, item_list_ray_trace_mode)
        # combo box Ray Texture
        item_list_ray_texture = [
            "None",
            "Matte 1",
            "Matte 2",
            "Swirl 1",
            "Fiber",
        ]
        gui_utils.fill_combo_box(self.ui.box_ray_texture, item_list_ray_texture)
        item_list = [
            "normal",
            "Red-green (green weak, deuteranopia)",
            "Red-green (red weak, protanopia)",
            "Blue-yellow (tritanopia)",
        ]
        gui_utils.fill_combo_box(self.ui.cb_color_vision_mode, item_list)
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
        self.ui.lbl_color_vision_mode.hide()
        self.ui.cb_color_vision_mode.hide()

        self.ui.cb_color_vision_mode.setCurrentIndex(
            self.ui.cb_color_vision_mode.findText(self.settings.color_vision_mode)
        )
        self.ui.box_bg_color.setCurrentIndex(
            self.ui.box_bg_color.findText(self.settings.image_background_color)
        )
        if self.settings.start_help_at_startup == 1:
            self.ui.check_box_start_help.setChecked(True)
        else:
            self.ui.check_box_start_help.setChecked(False)
        if self.settings.image_renderer == "0":
            self.ui.box_renderer.setCurrentIndex(0)
        else:
            self.ui.box_renderer.setCurrentIndex(1)

        self.ui.box_ray_trace_mode.setCurrentIndex(self.settings.image_ray_trace_mode)
        self.ui.box_ray_texture.setCurrentIndex(self.settings.image_ray_texture)

    # @SLOT()
    def choose_workspace_dir(self) -> None:
        """Opens a QFileDialog to choose a workspace directory."""
        gui_utils.choose_directory(self, self.ui.txt_workspace_dir)

    def set_color_vision_mode(self) -> None:
        """Sets the color vision mode."""
        raise NotImplementedError("This function is not yet implemented!")

    def _connect_all_gui_elements(self) -> None:
        """Connects all dialog gui elements."""
        #self.ui.cb_color_vision_mode.currentIndexChanged.connect(self.set_color_vision_mode)
        self.ui.btn_workspace_dir.clicked.connect(self.choose_workspace_dir)
        self.ui.btn_ok.clicked.connect(self.ok_dialog)
        self.ui.btn_help.clicked.connect(self.open_page_information)

    def cancel_dialog(self) -> None:
        """Closes the dialog."""
        self.close()

    def ok_dialog(self) -> None:
        """Sets all settings from the gui elements into the settings object and closes the dialog window."""
        self.settings.set_workspace_path(self.ui.txt_workspace_dir.text())
        self.settings.set_cycles(self.ui.spb_cycles.value())
        self.settings.set_cutoff(self.ui.dspb_cutoff.value())
        self.settings.color_vision_mode = self.ui.cb_color_vision_mode.currentText()
        self.settings.image_background_color = self.ui.box_bg_color.currentText()
        if self.ui.check_box_start_help.isChecked():
            self.settings.start_help_at_startup = 1
        else:
            self.settings.start_help_at_startup = 0
        if self.ui.box_renderer.currentText() == "default renderer":
            self.settings.image_renderer = "0"
        else:
            self.settings.image_renderer = "-1"
        self.settings.image_ray_trace_mode = self.ui.box_ray_trace_mode.currentIndex()
        self.settings.image_ray_texture = self.ui.box_ray_texture.currentIndex()

        self.settings.serialize_settings()
        logging.info("Settings were successfully saved.")
        self.close()

    def open_help(self, a_page_name: str):
        """Opens the pyssa documentation window if it's not already open.

        Args:
            a_page_name (str): a name of a documentation page to display
        """
        self._active_task = tasks.Task(
            target=util_async.open_documentation_on_certain_page,
            args=(a_page_name, self._interface_manager.documentation_window),
            post_func=self.__await_open_help,
        )
        self._active_task.start()

    def __await_open_help(self):
        subprocess.run([constants.HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH])

    def open_page_information(self) -> None:
        """Opens the message box, to display extra information based on the page."""
        self.open_help("help/settings/pyssa_settings/")
