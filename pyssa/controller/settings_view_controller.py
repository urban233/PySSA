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
"""Module for the Settings Dialog."""
import os
import subprocess
import logging

from PyQt5 import QtCore
from pyssa.controller import interface_manager
from pyssa.gui.ui.custom_dialogs import custom_message_box
from pyssa.internal.thread import tasks
from pyssa.internal.thread.async_pyssa import util_async
from pyssa.util import gui_utils
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class SettingsViewController(QtCore.QObject):
    """Class for the Create Project View Controller."""
    user_input = QtCore.pyqtSignal(tuple)

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager") -> None:
        super().__init__()
        self._interface_manager = the_interface_manager
        self._settings_manager = the_interface_manager.get_settings_manager()
        self._view = the_interface_manager.get_settings_view()
        self._initialize_ui()
        self.restore_ui()
        self._connect_all_ui_elements_to_slot_functions()

    # <editor-fold desc="Util methods">
    def open_help(self, a_page_name: str):
        """Opens the pyssa documentation window if it's not already open.

        Args:
            a_page_name (str): a name of a documentation page to display
        """
        self._interface_manager.update_status_bar("Opening help center ...")
        self._active_task = tasks.Task(
            target=util_async.open_documentation_on_certain_page,
            args=(a_page_name, self._interface_manager.documentation_window),
            post_func=self.__await_open_help,
        )
        self._active_task.start()

    def __await_open_help(self, return_value):
        self._interface_manager.documentation_window = return_value[2]
        if not os.path.exists(constants.HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH):
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "The script for bringing the documentation window in front could not be found!", "Documentation",
                custom_message_box.CustomMessageBoxIcons.ERROR.value
            )
            tmp_dialog.exec_()
        else:
            self._interface_manager.documentation_window.restore()
            subprocess.run([constants.HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH])
            self._interface_manager.update_status_bar("Opening help center finished.")

    def _open_help_for_dialog(self):
        self.open_help("help/settings/pyssa_settings/")

    def restore_ui(self):
        """Restores the ui."""
        self._view.ui.tabWidget.setCurrentIndex(0)

    # </editor-fold>
    
    def _initialize_ui(self):
        # combo box BgColor
        item_list_bg_color = [
            "black",
            "white",
        ]
        gui_utils.fill_combo_box(self._view.ui.box_bg_color, item_list_bg_color)
        # combo box Renderer
        item_list_renderer = [
            "default renderer",
            "PyMOL internal renderer",
        ]
        gui_utils.fill_combo_box(self._view.ui.box_renderer, item_list_renderer)
        # combo box RayTraceMode
        item_list_ray_trace_mode = [
            "normal color",
            "normal color + black outline",
            "black outline only",
            "quantized color + black outline",
        ]
        gui_utils.fill_combo_box(self._view.ui.box_ray_trace_mode, item_list_ray_trace_mode)
        # combo box Ray Texture
        item_list_ray_texture = [
            "None",
            "Matte 1",
            "Matte 2",
            "Swirl 1",
            "Fiber",
        ]
        gui_utils.fill_combo_box(self._view.ui.box_ray_texture, item_list_ray_texture)
        item_list_proteins_tab_representation_gui_elements = [
            "Toggle Buttons",
            "Checkboxes",
        ]
        gui_utils.fill_combo_box(self._view.ui.box_proteins_tab_toggle_checkbox,
                                 item_list_proteins_tab_representation_gui_elements)
        gui_utils.fill_combo_box(self._view.ui.box_protein_pairs_tab_toggle_checkbox,
                                 item_list_proteins_tab_representation_gui_elements)
        item_list_proteins_tab_color_gui_elements = [
            "Color Grid",
            "Combobox",
        ]
        gui_utils.fill_combo_box(self._view.ui.box_proteins_tab_box_color,
                                 item_list_proteins_tab_color_gui_elements)
        gui_utils.fill_combo_box(self._view.ui.box_protein_pairs_tab_box_color,
                                 item_list_proteins_tab_color_gui_elements)
        # item_list = [
        #     "normal",
        #     "Red-green (green weak, deuteranopia)",
        #     "Red-green (red weak, protanopia)",
        #     "Blue-yellow (tritanopia)",
        # ]
        # gui_utils.fill_combo_box(self._view.ui.cb_color_vision_mode, item_list)
        self._view.ui.txt_workspace_dir.setEnabled(False)
        self._view.ui.txt_workspace_dir.setText(str(self._settings_manager.settings.get_workspace_path()))
        self._view.ui.spb_cycles.setValue(int(self._settings_manager.settings.get_cycles()))
        self._view.ui.dspb_cutoff.setValue(float(self._settings_manager.settings.get_cutoff()))
        # customize spin boxes
        self._view.ui.spb_cycles.setMinimum(0)
        # self._view.ui.spb_cycles.setMaximum(20) # fixme: is a maximum needed?
        self._view.ui.spb_cycles.setSingleStep(1)
        self._view.ui.dspb_cutoff.setMinimum(0.00)
        self._view.ui.dspb_cutoff.setMaximum(20.00)
        self._view.ui.dspb_cutoff.setSingleStep(0.1)
        # self._view.ui.lbl_color_vision_mode.hide()
        # self._view.ui.cb_color_vision_mode.hide()
        #
        # self._view.ui.cb_color_vision_mode.setCurrentIndex(
        #     self._view.ui.cb_color_vision_mode.findText(self._settings_manager.settings.color_vision_mode)
        # )
        self._view.ui.box_bg_color.setCurrentIndex(
            self._view.ui.box_bg_color.findText(self._settings_manager.settings.image_background_color)
        )
        if self._settings_manager.settings.start_help_at_startup == 1:
            self._view.ui.check_box_start_help.setChecked(True)
        else:
            self._view.ui.check_box_start_help.setChecked(False)
        if self._settings_manager.settings.image_renderer == "0":
            self._view.ui.box_renderer.setCurrentIndex(0)
        else:
            self._view.ui.box_renderer.setCurrentIndex(1)
        self._view.ui.box_ray_trace_mode.setCurrentIndex(self._settings_manager.settings.image_ray_trace_mode)
        self._view.ui.box_ray_texture.setCurrentIndex(self._settings_manager.settings.image_ray_texture)

        if self._settings_manager.settings.proteins_tab_use_toggle == 1:
            self._view.ui.box_proteins_tab_toggle_checkbox.setCurrentIndex(0)
        else:
            self._view.ui.box_proteins_tab_toggle_checkbox.setCurrentIndex(1)
        if self._settings_manager.settings.proteins_tab_use_combobox_for_colors == 1:
            self._view.ui.box_proteins_tab_box_color.setCurrentIndex(1)
        else:
            self._view.ui.box_proteins_tab_box_color.setCurrentIndex(0)
        if self._settings_manager.settings.protein_pairs_tab_use_toggle == 1:
            self._view.ui.box_protein_pairs_tab_toggle_checkbox.setCurrentIndex(0)
        else:
            self._view.ui.box_protein_pairs_tab_toggle_checkbox.setCurrentIndex(1)
        if self._settings_manager.settings.protein_pairs_tab_use_combobox_for_colors == 1:
            self._view.ui.box_protein_pairs_tab_box_color.setCurrentIndex(1)
        else:
            self._view.ui.box_protein_pairs_tab_box_color.setCurrentIndex(0)
    
    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.btn_workspace_dir.clicked.connect(self.choose_workspace_dir)
        self._view.ui.btn_ok.clicked.connect(self.ok_dialog)
        self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)
        
    def choose_workspace_dir(self) -> None:
        """Opens a QFileDialog to choose a workspace directory."""
        gui_utils.choose_directory(self, self._view.ui.txt_workspace_dir)
    
    def ok_dialog(self) -> None:
        """Sets all settings from the gui elements into the settings object and closes the dialog window."""
        self._settings_manager.settings.set_workspace_path(self._view.ui.txt_workspace_dir.text())
        self._settings_manager.settings.set_cycles(self._view.ui.spb_cycles.value())
        self._settings_manager.settings.set_cutoff(self._view.ui.dspb_cutoff.value())
        # self._settings_manager.settings.color_vision_mode = self._view.ui.cb_color_vision_mode.currentText()
        self._settings_manager.settings.image_background_color = self._view.ui.box_bg_color.currentText()
        if self._view.ui.check_box_start_help.isChecked():
            self._settings_manager.settings.start_help_at_startup = 1
        else:
            self._settings_manager.settings.start_help_at_startup = 0
        if self._view.ui.box_renderer.currentText() == "default renderer":
            self._settings_manager.settings.image_renderer = "0"
        else:
            self._settings_manager.settings.image_renderer = "-1"
        self._settings_manager.settings.image_ray_trace_mode = self._view.ui.box_ray_trace_mode.currentIndex()
        self._settings_manager.settings.image_ray_texture = self._view.ui.box_ray_texture.currentIndex()

        if self._view.ui.box_proteins_tab_box_color.currentIndex() == 0:
            self._settings_manager.settings.proteins_tab_use_combobox_for_colors = 0
        else:
            self._settings_manager.settings.proteins_tab_use_combobox_for_colors = 1
        if self._view.ui.box_proteins_tab_toggle_checkbox.currentIndex() == 0:
            self._settings_manager.settings.proteins_tab_use_toggle = 1
        else:
            self._settings_manager.settings.proteins_tab_use_toggle = 0
        if self._view.ui.box_protein_pairs_tab_box_color.currentIndex() == 0:
            self._settings_manager.settings.protein_pairs_tab_use_combobox_for_colors = 0
        else:
            self._settings_manager.settings.protein_pairs_tab_use_combobox_for_colors = 1
        if self._view.ui.box_protein_pairs_tab_toggle_checkbox.currentIndex() == 0:
            self._settings_manager.settings.protein_pairs_tab_use_toggle = 1
        else:
            self._settings_manager.settings.protein_pairs_tab_use_toggle = 0

        self._settings_manager.settings.serialize_settings()
        logging.info("Settings were successfully saved.")
        self._view.close()
        self.user_input.emit((0, True))
