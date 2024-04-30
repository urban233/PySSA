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
"""Module for the Open Dialog."""
import logging
import os
import subprocess

import pygetwindow
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from pyssa.controller import interface_manager
from pyssa.gui.ui.custom_dialogs import custom_message_box
from pyssa.internal.data_structures.data_classes import prediction_configuration
from pyssa.internal.thread import tasks
from pyssa.internal.thread.async_pyssa import validate_async, util_async
from pyssa.io_pyssa import bio_data
from pyssa.util import constants, gui_utils
from pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class AdvancedPredictionConfigurationsViewController(QtCore.QObject):
    """Class for the RenameProteinViewController class"""
    user_input = QtCore.pyqtSignal(tuple)

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager", a_prediction_configuration):
        super().__init__()
        self._interface_manager = the_interface_manager
        self._view = the_interface_manager.get_advanced_prediction_configurations_view()
        self.prediction_config: prediction_configuration.PredictionConfiguration = a_prediction_configuration
        self._connect_all_ui_elements_to_slot_functions()
        self._view.ui.cb_amber.setChecked(self.prediction_config.amber_force_field)
        item_list_templates = [
            "none",
            "pdb70",
            # "custom", TODO: implement a way to add a custom MSA
        ]
        gui_utils.fill_combo_box(self._view.ui.combo_box_template, item_list_templates)
        self._view.ui.combo_box_template.setCurrentIndex(
            self._view.ui.combo_box_template.findText(self.prediction_config.templates),
        )

    def open_help(self, a_page_name: str):
        """Opens the pyssa documentation window if it's not already open.

        Args:
            a_page_name (str): a name of a documentation page to display
        """
        self._interface_manager.status_bar_manager.show_temporary_message("Opening help center ...")
        if len(pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)) != 1:
            self._interface_manager.documentation_window = None
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
            self._interface_manager.status_bar_manager.show_temporary_message("Opening help center finished.")

    def _open_help_for_dialog(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Help' button was clicked.")
        self.open_help("help/protein_structure_prediction/advanced_prediction_configurations/")

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)
        self._view.ui.btn_ok.clicked.connect(self.save_config)

    def restore_ui(self):
        self._view.setMinimumWidth(500)

    def save_config(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'OK' button was clicked.")
        self.prediction_config.amber_force_field = self._view.ui.cb_amber.isChecked()
        self.prediction_config.templates = self._view.ui.combo_box_template.currentText()
        self._view.close()
        self.user_input.emit((0, self.prediction_config))
