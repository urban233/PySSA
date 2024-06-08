#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/zielesny/PySSA>
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
"""Module for the advanced prediction configuration view controller."""
import logging
import os
import subprocess

import pygetwindow
from PyQt5 import QtCore
from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.internal.data_structures.data_classes import prediction_configuration
from src.pyssa.internal.thread import tasks
from src.pyssa.internal.thread.async_pyssa import util_async
from src.pyssa.util import constants, exception
from src.pyssa.util import gui_utils
from src.pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class AdvancedPredictionConfigurationsViewController(QtCore.QObject):
  """Class for the AdvancedPredictionConfigurationsViewController."""

  user_input = QtCore.pyqtSignal(tuple)
  """Singal used to transfer data back to the previous window."""

  def __init__(
      self,
      the_interface_manager: "interface_manager.InterfaceManager",
      a_prediction_configuration,
  ) -> None:
    """Constructor.

    Args:
        the_interface_manager (interface_manager.InterfaceManager): The InterfaceManager object.
        a_prediction_configuration: The configuration for the prediction.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_interface_manager is None:
      logger.error("the_interface_manager is None.")
      raise exception.IllegalArgumentError("the_interface_manager is None.")
    if a_prediction_configuration is None:
      logger.error("a_prediction_configuration is None.")
      raise exception.IllegalArgumentError(
          "a_prediction_configuration is None."
      )

    # </editor-fold>

    super().__init__()
    self._interface_manager = the_interface_manager
    self._view = (
        the_interface_manager.get_advanced_prediction_configurations_view()
    )
    self.prediction_config: prediction_configuration.PredictionConfiguration = (
        a_prediction_configuration
    )
    self._connect_all_ui_elements_to_slot_functions()
    self._view.ui.cb_amber.setChecked(self.prediction_config.amber_force_field)
    item_list_templates = [
        "none",
        "pdb70",
        # "custom", TODO: implement a way to add a custom MSA
    ]
    gui_utils.fill_combo_box(
        self._view.ui.combo_box_template, item_list_templates
    )
    self._view.ui.combo_box_template.setCurrentIndex(
        self._view.ui.combo_box_template.findText(
            self.prediction_config.templates
        ),
    )

  def _open_help_for_dialog(self) -> None:
    """Opens the help dialog for the corresponding dialog."""
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Help' button was clicked."
    )
    self._interface_manager.help_manager.open_advanced_prediction_configuration_page()

  def _connect_all_ui_elements_to_slot_functions(self) -> None:
    """Connects all UI elements to their corresponding slot functions in the class."""
    self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)
    self._view.ui.btn_ok.clicked.connect(self.save_config)

  def restore_ui(self) -> None:
    """Restores the UI."""
    self._view.setMinimumWidth(500)

  def save_config(self) -> None:
    """Saves the configuration settings from the view by sending the `user_input` signal and closing the dialog window."""
    logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'OK' button was clicked.")
    self.prediction_config.amber_force_field = (
        self._view.ui.cb_amber.isChecked()
    )
    self.prediction_config.templates = (
        self._view.ui.combo_box_template.currentText()
    )
    self._view.close()
    self.user_input.emit((0, self.prediction_config))
