#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
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
"""Module for the results view controller."""
import logging
import os

from typing import TYPE_CHECKING

from src.pyssa.gui.qt import QtWidgets
from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.qt import Qt
from src.pyssa.controller import database_manager
from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.gui.ui.views import results_view, help_view
from src.pyssa.gui.ui.views import plot_view
from src.pyssa.internal.thread import tasks
from src.pyssa.internal.thread.async_pyssa import protein_pair_async
from src.pyssa.util import constants, enums, exception
from src.pyssa.util import gui_utils
from src.pyssa.logging_pyssa import log_levels, log_handlers

if TYPE_CHECKING:
  from src.pyssa.gui import app_state, user_pymol
  from src.pyssa.internal.data_structures import protein_pair

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class ResultsViewController(QtCore.QObject):
  """Class for the ResultsViewController."""

  def __init__(
      self,
      the_protein_pair: "protein_pair.ProteinPair",
      the_app_state: "app_state.AppState",
      the_user_pymol: "user_pymol.UserPyMOL",
      a_parent=None
  ) -> None:
    """Constructor.

    Args:
        the_protein_pair (protein_pair.ProteinPair): An instance of the protein pair class.
        the_pymol_session_manager (pymol_session_manager.PymolSessionManager): An instance of the Pymol session manager class.
        the_app_state (app_state.AppState): The AppState object.
        a_parent: Parent widget to pass to the view.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_app_state is None:
      logger.error("the_app_state is None.")
      raise exception.IllegalArgumentError("the_app_state is None.")
    if the_protein_pair is None:
      logger.error("the_protein_pair is None.")
      raise exception.IllegalArgumentError("the_protein_pair is None.")

    # </editor-fold>

    super().__init__()
    self._app_state = the_app_state
    self._user_pymol = the_user_pymol
    self._protein_pair = the_protein_pair
    self._color_configuration_protein_pair: tuple[list, list] = ([], [])
    self._view: "results_view.ResultsView" = results_view.ResultsView(a_parent)
    self._distance_data_visualizer = None
    self.cb_protein_pair_color = QtWidgets.QComboBox()
    self._connect_all_ui_elements_to_slot_functions()
    self._build_table_widget()
    self._fill_table_widget()
    self._check_results()
    self._view.setWindowTitle("Results Summary")

  def get_view(self):
    return self._view

  def restore_default_view(self):
    if self._user_pymol.get_currently_loaded_object() == self._protein_pair:
      self._view.ui.btn_color_by_rmsd.setEnabled(True)
    else:
      self._view.ui.btn_color_by_rmsd.setEnabled(False)

  def _open_help_for_dialog(self) -> None:
    """Opens the help dialog."""
    tmp_dialog = help_view.HelpView(
      constants.HELP_TEXT_MAP["ResultsSummaryDialog"]
    )
    tmp_dialog.exec()

  def _connect_all_ui_elements_to_slot_functions(self) -> None:
    """Connects all UI elements to their corresponding slot functions in the class."""
    self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)
    self._view.ui.btn_view_plots.clicked.connect(self._open_plot_view)
    self._view.ui.btn_color_by_rmsd.clicked.connect(
        self.__slot_color_protein_pair
    )
    self._view.ui.btn_export_data.clicked.connect(self.__slot_export_data)

  # <editor-fold desc="Util methods">
  def _build_table_widget(self) -> None:
    """Builds the table widget for displaying results."""
    self._view.ui.table_widget_results.clear()
    self._view.ui.table_widget_results.setRowCount(4)
    self._view.ui.table_widget_results.setColumnCount(2)
    self._view.ui.table_widget_results.verticalHeader().setVisible(False)
    self._view.ui.table_widget_results.setHorizontalHeaderLabels(
        ["Name", "Value"]
    )
    self._view.ui.table_widget_results.horizontalHeader().setSectionResizeMode(
        QtWidgets.QHeaderView.ResizeMode.ResizeToContents
    )

  def _fill_table_widget(self) -> None:
    """Fill the table widget with data and configure it."""
    # RMSD value item for table widget
    tmp_rmsd_label_item = QtWidgets.QTableWidgetItem("RMSD (Å)")
    tmp_result_rmsd_item = QtWidgets.QTableWidgetItem(
        str(self._protein_pair.distance_analysis.analysis_results.rmsd)
    )
    self._view.ui.table_widget_results.setItem(0, 0, tmp_rmsd_label_item)
    self._view.ui.table_widget_results.setItem(0, 1, tmp_result_rmsd_item)
    # Aligned residues value item for table widget
    tmp_aligned_aa_label_item = QtWidgets.QTableWidgetItem("Aligned Residues")
    tmp_result_aligned_aa_item = QtWidgets.QTableWidgetItem(
        str(self._protein_pair.distance_analysis.analysis_results.aligned_aa)
    )
    self._view.ui.table_widget_results.setItem(1, 0, tmp_aligned_aa_label_item)
    self._view.ui.table_widget_results.setItem(1, 1, tmp_result_aligned_aa_item)
    # Color information about protein 1 & 2
    tmp_color_description_protein_1_item = QtWidgets.QTableWidgetItem(
        f"Color {self._protein_pair.protein_1.get_molecule_object()}"
    )
    self._view.ui.table_widget_results.setItem(
        2, 0, tmp_color_description_protein_1_item
    )
    tmp_color_protein_1_item = QtWidgets.QTableWidgetItem("green (default)")
    self._view.ui.table_widget_results.setItem(2, 1, tmp_color_protein_1_item)
    tmp_color_description_protein_2_item = QtWidgets.QTableWidgetItem(
        f"Color {self._protein_pair.protein_2.get_molecule_object()}"
    )
    self._view.ui.table_widget_results.setItem(
        3, 0, tmp_color_description_protein_2_item
    )
    tmp_color_protein_2_item = QtWidgets.QTableWidgetItem("blue (default)")
    self._view.ui.table_widget_results.setItem(3, 1, tmp_color_protein_2_item)
    # Set editing flags
    tmp_color_description_protein_1_item.setFlags(
        tmp_color_description_protein_1_item.flags() & ~Qt.ItemIsEditable
    )
    tmp_color_protein_1_item.setFlags(
        tmp_color_protein_1_item.flags() & ~Qt.ItemIsEditable
    )
    tmp_color_description_protein_2_item.setFlags(
        tmp_color_description_protein_2_item.flags() & ~Qt.ItemIsEditable
    )
    tmp_color_protein_2_item.setFlags(
        tmp_color_protein_2_item.flags() & ~Qt.ItemIsEditable
    )
    tmp_result_rmsd_item.setFlags(
        tmp_result_rmsd_item.flags() & ~Qt.ItemIsEditable
    )
    tmp_rmsd_label_item.setFlags(
        tmp_rmsd_label_item.flags() & ~Qt.ItemIsEditable
    )
    tmp_aligned_aa_label_item.setFlags(
        tmp_aligned_aa_label_item.flags() & ~Qt.ItemIsEditable
    )
    tmp_result_aligned_aa_item.setFlags(
        tmp_result_aligned_aa_item.flags() & ~Qt.ItemIsEditable
    )
    # Resize
    self._view.ui.table_widget_results.resizeColumnToContents(0)
    self._view.ui.table_widget_results.resizeColumnToContents(1)

  def _check_results(self) -> None:
    """Checks the analysis results and enable/disable the "View Plots" button accordingly."""
    if self._protein_pair.distance_analysis.analysis_results.rmsd == 0:
      self._view.ui.btn_view_plots.setEnabled(False)
      self._view.ui.btn_view_plots.setToolTip(
          "The RMSD value is exact 0, therefore no plots can be displayed."
      )
    else:
      self._view.ui.btn_view_plots.setEnabled(True)
      self._view.ui.btn_view_plots.setToolTip("")

  # </editor-fold>

  # <editor-fold desc="Methods for coloring protein pair">
  def __slot_color_protein_pair(self) -> None:
    """Colors the protein pair based on the selection made in the 'By RMSD' checkbox."""
    self._user_pymol.color_protein_pair_by_rmsd(self._protein_pair)
  # </editor-fold>

  def _open_plot_view(self) -> None:
    """Opens the plot view."""
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'View plots' button was clicked."
    )
    self._distance_data_visualizer = plot_view.PlotView(
      self._protein_pair,
      self._app_state.project,
      self._protein_pair,
      None,
      parent=self._view,
    )
    self._distance_data_visualizer.show()

  def __slot_export_data(self) -> None:
    """Handles the export data functionality."""
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
      "'Export data' button was clicked.",
    )
    file_dialog = QtWidgets.QFileDialog()
    desktop_path = QtCore.QStandardPaths.standardLocations(
      QtCore.QStandardPaths.StandardLocation.DesktopLocation
    )[0]
    file_dialog.setDirectory(desktop_path)
    file_path, _ = file_dialog.getSaveFileName(
      self._view,
      "Export Distance Data",
      "",
      "Comma-Separated Values (*.csv)",
    )
    if file_path:
      self._protein_pair.distance_analysis.analysis_results.export_distance_data_as_csv(
        file_path
      )
      if os.path.exists(file_path):
        tmp_dialog = custom_message_box.CustomMessageBoxOk(
          "Export data as .csv file finished.",
          "Export Data",
          custom_message_box.CustomMessageBoxIcons.INFORMATION.value,
        )
      else:
        tmp_dialog = custom_message_box.CustomMessageBoxOk(
          "Export data as .csv file failed!",
          "Export Data",
          custom_message_box.CustomMessageBoxIcons.ERROR.value,
        )
      tmp_dialog.exec()
