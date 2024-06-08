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
"""Module for the distance analysis view controller."""
import logging
import os
import subprocess

import pygetwindow
from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal

from src.pyssa.controller import add_protein_pair_view_controller
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


class DistanceAnalysisViewController(QtCore.QObject):
  """Class for the DistanceAnalysisViewController."""

  job_input = pyqtSignal(tuple)
  """Singal used to transfer data back to the previous window."""

  def __init__(
      self,
      the_interface_manager: "interface_manager.InterfaceManager",
      the_watcher: "watcher.Watcher",
  ) -> None:
    """Constructor.

    Args:
        the_interface_manager (interface_manager.InterfaceManager): The InterfaceManager object.
        the_watcher (watcher.Watcher): The Watcher object.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_interface_manager is None:
      logger.error("the_interface_manager is None.")
      raise exception.IllegalArgumentError("the_interface_manager is None.")
    if the_watcher is None:
      logger.error("the_watcher is None.")
      raise exception.IllegalArgumentError("the_watcher is None.")

    # </editor-fold>

    super().__init__()
    self._interface_manager = the_interface_manager
    self._watcher = the_watcher
    self._view: "distance_analysis_view.DistanceAnalysisView" = (
        the_interface_manager.get_distance_analysis_view()
    )
    self.prediction_configuration = (
        prediction_configuration.PredictionConfiguration(True, "pdb70")
    )

    self._connect_all_ui_elements_to_slot_functions()
    self.display_distance_analysis()

  def _open_help_for_dialog(self) -> None:
    """Opens the help dialog window."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Help' button was clicked."
    )
    self._interface_manager.help_manager.open_distance_analysis_page()

  def _connect_all_ui_elements_to_slot_functions(self) -> None:
    """Connects all UI elements to their corresponding slot functions in the class."""
    self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)
    self._view.ui.btn_distance_analysis_add.clicked.connect(
        self._add_protein_pair
    )
    self._view.ui.btn_distance_analysis_remove.clicked.connect(
        self.remove_analysis_run
    )
    self._view.ui.btn_distance_analysis_start.clicked.connect(
        self.start_process_batch
    )
    self._view.ui.list_distance_analysis_overview.clicked.connect(
        self.structure_analysis_overview_clicked
    )

  def display_distance_analysis(self) -> None:
    """Displays the job analysis work area."""
    gui_elements_to_show = [
        self._view.ui.btn_distance_analysis_add,
        self._view.ui.lbl_distance_analysis_overview,
        self._view.ui.list_distance_analysis_overview,
    ]
    gui_elements_to_hide = [
        self._view.ui.btn_distance_analysis_remove,
        self._view.ui.btn_distance_analysis_start,
    ]
    gui_utils.show_gui_elements(gui_elements_to_show)
    gui_utils.hide_gui_elements(gui_elements_to_hide)

    self._view.ui.list_distance_analysis_overview.clear()

  def start_process_batch(self) -> None:
    """Starts the process batch based on selected items in the distance analysis overview."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Start' button was clicked."
    )
    tmp_raw_analysis_run_names: list = []
    for row_no in range(self._view.ui.list_distance_analysis_overview.count()):
      tmp_raw_analysis_run_names.append(
          self._view.ui.list_distance_analysis_overview.item(row_no).text()
      )
    self._view.close()
    self.job_input.emit(("job_input", tmp_raw_analysis_run_names, False))

  def _get_all_current_analysis_runs(self) -> list[str]:
    """Retrieves a list of all current analysis runs.

    Returns:
        A list of strings representing the current analysis runs.
    """
    tmp_analysis_runs = []
    for tmp_row in range(self._view.ui.list_distance_analysis_overview.count()):
      tmp_analysis_runs.append(
          self._view.ui.list_distance_analysis_overview.item(tmp_row).text()
      )
    return tmp_analysis_runs

  def _get_all_current_protein_pair_names(self) -> list[str]:
    """Retrieves the names of all current protein pairs.

    Returns:
        A list of strings representing the names of all current protein pairs.
    """
    tmp_protein_pair_names = []
    for (
        tmp_protein_pair
    ) in self._interface_manager.get_current_project().protein_pairs:
      tmp_protein_pair_names.append(tmp_protein_pair.name)
    return tmp_protein_pair_names

  def _add_protein_pair(self) -> None:
    """Instantiates the AddProteinPairViewController class and shows the 'AddProteinPair' view."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Add' button was clicked."
    )
    self._external_controller = (
        add_protein_pair_view_controller.AddProteinPairViewController(
            self._interface_manager,
            self._watcher,
            self._get_all_current_analysis_runs(),
            self._get_all_current_protein_pair_names(),
        )
    )
    self._external_controller.user_input.connect(self._post_add_protein_pair)
    self._interface_manager.get_add_protein_pair_view().show()

  def _post_add_protein_pair(self, return_value: tuple) -> None:
    """Adds the protein pair to the list widget and updates the UI.

    Args:
        return_value (tuple): The return value from the method.

    Raises:
        exception.IllegalArgumentError: If `return_value` is None.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      raise exception.IllegalArgumentError("return_value is None.")

    # </editor-fold>

    tmp_item, _ = return_value
    self._view.ui.list_distance_analysis_overview.addItem(tmp_item)
    self._view.ui.btn_distance_analysis_remove.show()
    self._view.ui.btn_distance_analysis_remove.setEnabled(False)
    self._view.ui.btn_distance_analysis_start.show()
    self._view.ui.btn_distance_analysis_start.setEnabled(True)

  def structure_analysis_overview_clicked(self) -> None:
    """Enables the remove button."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "A distance analysis run from the list was clicked.",
    )
    self._view.ui.btn_distance_analysis_remove.setEnabled(True)

  def remove_analysis_run(self) -> None:
    """Removes the selected protein pair from the list of protein pairs to analyze."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Remove' button was clicked."
    )
    self._view.ui.list_distance_analysis_overview.takeItem(
        self._view.ui.list_distance_analysis_overview.currentRow()
    )
    if self._view.ui.list_distance_analysis_overview.count() == 0:
      gui_elements_to_show = [
          self._view.ui.lbl_distance_analysis_overview,
          self._view.ui.list_distance_analysis_overview,
          self._view.ui.btn_distance_analysis_add,
      ]

      gui_elements_to_hide = [
          self._view.ui.btn_distance_analysis_remove,
          self._view.ui.btn_distance_analysis_start,
      ]

      gui_utils.show_gui_elements(gui_elements_to_show)
      gui_utils.hide_gui_elements(gui_elements_to_hide)
      self._view.ui.btn_distance_analysis_remove.hide()
    else:
      if self._view.ui.list_distance_analysis_overview.count() > 0:
        try:
          self._view.ui.list_distance_analysis_overview.currentItem().setSelected(
              False
          )
        except AttributeError:
          constants.PYSSA_LOGGER.debug(
              "No selection in struction analysis overview."
          )
    self._view.ui.btn_distance_analysis_remove.setEnabled(False)
