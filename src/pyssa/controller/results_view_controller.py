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
"""Module for the results view controller."""
import logging
import os
import subprocess

import pygetwindow
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from src.pyssa.controller import database_manager
from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.gui.ui.views import plot_view
from src.pyssa.internal.thread import tasks
from src.pyssa.internal.thread.async_pyssa import util_async, protein_pair_async
from src.pyssa.util import constants, enums, exception
from src.pyssa.util import gui_utils
from src.pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class ResultsViewController(QtCore.QObject):
  """Class for the ResultsViewController."""

  def __init__(
      self,
      the_interface_manager: "interface_manager.InterfaceManager",
      the_protein_pair: "protein_pair.ProteinPair",
      the_pymol_session_manager: "pymol_session_manager.PymolSessionManager",
  ) -> None:
    """Constructor.

    Args:
        the_interface_manager (interface_manager.InterfaceManager): An instance of the interface manager class.
        the_protein_pair (protein_pair.ProteinPair): An instance of the protein pair class.
        the_pymol_session_manager (pymol_session_manager.PymolSessionManager): An instance of the Pymol session manager class.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_interface_manager is None:
      logger.error("the_interface_manager is None.")
      raise exception.IllegalArgumentError("the_interface_manager is None.")
    if the_protein_pair is None:
      logger.error("the_protein_pair is None.")
      raise exception.IllegalArgumentError("the_protein_pair is None.")
    if the_pymol_session_manager is None:
      logger.error("the_pymol_session_manager is None.")
      raise exception.IllegalArgumentError("the_pymol_session_manager is None.")

    # </editor-fold>

    super().__init__()
    self._interface_manager = the_interface_manager
    self._pymol_session_manager = the_pymol_session_manager
    self._protein_pair = the_protein_pair
    self._color_configuration_protein_pair: tuple[list, list] = ([], [])
    self._view: "results_view.ResultsView" = (
        the_interface_manager.get_results_view()
    )
    self._distance_data_visualizer = None
    self.cb_protein_pair_color = QtWidgets.QComboBox()

    self._connect_all_ui_elements_to_slot_functions()
    self._build_table_widget()
    self._fill_table_widget()
    
    if not self._pymol_session_manager.is_the_current_protein_pair_in_session(
        self._interface_manager.get_current_active_protein_pair_object().name
    ):
      self.cb_protein_pair_color.setEnabled(False)
    else:
      self._check_coloring()
      
    self._check_results()
    
    self._view.setWindowTitle(f"Results Summary Of {self._protein_pair.name}")

  def _open_help_for_dialog(self) -> None:
    """Opens the help dialog."""
    self._interface_manager.help_manager.open_results_summary_page()

  def _connect_all_ui_elements_to_slot_functions(self) -> None:
    """Connects all UI elements to their corresponding slot functions in the class."""
    self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)
    self._view.ui.btn_view_plots.clicked.connect(self._open_plot_view)
    self.cb_protein_pair_color.currentIndexChanged.connect(
        self.__slot_color_protein_pair
    )
    self._view.ui.btn_export_data.clicked.connect(self.__slot_export_data)
    self._view.ui.table_widget_results.clicked.connect(
        self.__slot_highlight_protein
    )

  # <editor-fold desc="Util methods">
  def _build_table_widget(self) -> None:
    """Builds the table widget for displaying results."""
    self.cb_protein_pair_color.currentIndexChanged.disconnect(
        self.__slot_color_protein_pair
    )
    self._view.ui.table_widget_results.clear()
    self._view.ui.table_widget_results.setRowCount(5)
    self._view.ui.table_widget_results.setColumnCount(2)
    self._view.ui.table_widget_results.verticalHeader().setVisible(False)
    self._view.ui.table_widget_results.setHorizontalHeaderLabels(
        ["Name", "Value"]
    )
    gui_utils.fill_combo_box(self.cb_protein_pair_color, ["Default", "By RMSD"])
    self.cb_protein_pair_color.adjustSize()
    self.cb_protein_pair_color.currentIndexChanged.connect(
        self.__slot_color_protein_pair
    )

  def _fill_table_widget(self) -> None:
    """Fill the table widget with data and configure it."""
    self._get_color_configuration_of_protein_pair()
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
    # Color protein pair value item for table widget
    tmp_color_protein_pair_label_item = QtWidgets.QTableWidgetItem(
        "Color Protein Pair"
    )
    self._view.ui.table_widget_results.setItem(
        4, 0, tmp_color_protein_pair_label_item
    )
    self._view.ui.table_widget_results.setCellWidget(
        4, 1, self.cb_protein_pair_color
    )
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
    tmp_color_protein_pair_label_item.setFlags(
        tmp_color_protein_pair_label_item.flags() & ~Qt.ItemIsEditable
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
  
  def _check_coloring(self):
    """Checks if the protein pair is already colored by RMSD and sets combo box index accordingly."""
    tmp_chain_letter = None
    for tmp_chain in self._protein_pair.protein_1.chains:
      if tmp_chain.chain_type == "protein_chain":
        tmp_chain_letter = tmp_chain.chain_letter
    if tmp_chain_letter is None:
      return
    self._protein_pair.protein_1.pymol_selection.set_selection_for_a_single_chain(
      tmp_chain_letter
    )
    tmp_chain_color = self._pymol_session_manager.get_chain_color(
      self._protein_pair.protein_1.pymol_selection.selection_string, tmp_chain_letter
    )
    # if constants.PYMOL_COLORS_WITH_INDICES["rmsd_color_1"] in tmp_chain_color[1]:
    #   self.cb_protein_pair_color.setCurrentIndex(1)
    # elif constants.PYMOL_COLORS_WITH_INDICES["rmsd_color_2"] in tmp_chain_color[1]:
    #   self.cb_protein_pair_color.setCurrentIndex(1)
    # elif constants.PYMOL_COLORS_WITH_INDICES["rmsd_color_3"] in tmp_chain_color[1]:
    #   self.cb_protein_pair_color.setCurrentIndex(1)
    # elif constants.PYMOL_COLORS_WITH_INDICES["rmsd_color_4"] in tmp_chain_color[1]:
    #   self.cb_protein_pair_color.setCurrentIndex(1)
    # elif constants.PYMOL_COLORS_WITH_INDICES["rmsd_color_5"] in tmp_chain_color[1]:
    #   self.cb_protein_pair_color.setCurrentIndex(1)
    if "By RMSD" in tmp_chain_color[1]:
      self.cb_protein_pair_color.setCurrentIndex(1)
    else:
      self.cb_protein_pair_color.setCurrentIndex(0)
    print(tmp_chain_color)
    
  # </editor-fold>

  def __slot_highlight_protein(self) -> None:
    # fixme: This could be important if its needed
    # if self._view.ui.table_widget_results.currentItem().text().find("green (default)") != -1:
    #     tmp_pymol_selection = selection.Selection(self._protein_pair.protein_1.get_molecule_object())
    #     tmp_pymol_selection.selection_string = f"/{tmp_pymol_selection.molecule_object}"
    #     cmd.col(tmp_pymol_selection.selection_string)
    # elif self._view.ui.table_widget_results.currentItem().text().find("blue (default)") != -1:
    #     tmp_pymol_selection = selection.Selection(self._protein_pair.protein_2.get_molecule_object())
    #     tmp_pymol_selection.selection_string = f"/{tmp_pymol_selection.molecule_object}"
    #     cmd.select(tmp_pymol_selection.selection_string)
    # else:
    #     cmd.select(name="", selection="none")
    return
    raise NotImplementedError()

  def _open_plot_view(self) -> None:
    """Opens the plot view."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'View plots' button was clicked."
    )
    self._distance_data_visualizer = plot_view.PlotView(
        self._protein_pair,
        self._interface_manager.get_current_project(),
        self._protein_pair,
        self._interface_manager.help_manager
    )
    self._distance_data_visualizer.show()

  # <editor-fold desc="Methods for coloring protein pair">
  def __slot_color_protein_pair(self) -> None:
    """Colors the protein pair based on the selection made in the 'By RMSD' checkbox."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'By RMSD' checkbox was clicked."
    )
    if self.cb_protein_pair_color.currentText().find("By RMSD") != -1:
      self._color_protein_pair_by_rmsd()
    else:
      self._color_protein_back_to_normal()

  def _get_color_configuration_of_protein_pair(self) -> None:
    """Retrieves the color configuration for the protein pair from the database."""
    tmp_database_filepath: str = str(
        self._interface_manager.get_current_project().get_database_filepath()
    )
    with database_manager.DatabaseManager(tmp_database_filepath) as db_manager:
      for tmp_chain in self._protein_pair.protein_1.chains:
        self._color_configuration_protein_pair[0].append(
            db_manager.get_pymol_parameter_for_certain_protein_chain_in_protein_pair(
                self._protein_pair.get_id(),
                self._protein_pair.protein_1.get_id(),
                tmp_chain.chain_letter,
                enums.PymolParameterEnum.COLOR.value,
            )
        )
      for tmp_chain in self._protein_pair.protein_2.chains:
        self._color_configuration_protein_pair[1].append(
            db_manager.get_pymol_parameter_for_certain_protein_chain_in_protein_pair(
                self._protein_pair.get_id(),
                self._protein_pair.protein_2.get_id(),
                tmp_chain.chain_letter,
                enums.PymolParameterEnum.COLOR.value,
            )
        )

  def _color_protein_pair_by_rmsd(self) -> None:
    """Colors the residues in 5 colors depending on their distance to the reference."""
    self._active_task = tasks.LegacyTask(
        target=protein_pair_async.color_protein_pair_by_rmsd_value,
        args=(
            self._protein_pair,
            self._pymol_session_manager,
        ),
        post_func=self.__await_color_protein_pair_by_rmsd,
    )
    self._active_task.start()
    # hide unnecessary representations
    # fixme: it might be a problem to hide any representation at this point
    # cmd.hide("cartoon", tmp_protein_pair.protein_1.get_molecule_object())
    # cmd.hide("cartoon", f"{tmp_protein_pair.protein_2.get_molecule_object()}")
    # cmd.hide("cartoon", f"{tmp_protein_pair.protein_2.get_molecule_object()}")

  def __await_color_protein_pair_by_rmsd(self, result: tuple) -> None:
    """Checks if the coloring process failed.

    Args:
        result: A tuple containing the result of the protein pair coloring process.

    Raises:
        exception.IllegalArgumentError: If `result` is None.
    """
    # <editor-fold desc="Checks">
    if result is None:
      logger.error("result is None.")
      raise exception.IllegalArgumentError("result is None.")

    # </editor-fold>

    if result[0] == "":
      self._interface_manager.status_bar_manager.show_error_message(
          "Coloring the protein pair failed!"
      )

  def _color_protein_back_to_normal(self) -> None:
    """Reverts the color of the protein back to its original state."""
    self._get_color_configuration_of_protein_pair()
    i = 0
    tmp_protein_1_colors: list = self._color_configuration_protein_pair[0]
    for tmp_chain in self._protein_pair.protein_1.chains:
      self._protein_pair.protein_1.pymol_selection.set_selection_for_a_single_chain(
          tmp_chain.chain_letter
      )
      self._interface_manager.pymol_session_manager.color_protein(
          "green",  # tmp_protein_1_colors[i],
          self._protein_pair.protein_1.pymol_selection.selection_string,
      )
      # self._protein_pair.protein_1.pymol_selection.color_selection(tmp_protein_1_colors[i])
      i += 1

    i = 0
    tmp_protein_2_colors: list = self._color_configuration_protein_pair[1]
    for tmp_chain in self._protein_pair.protein_2.chains:
      self._protein_pair.protein_2.pymol_selection.set_selection_for_a_single_chain(
          tmp_chain.chain_letter
      )
      self._interface_manager.pymol_session_manager.color_protein(
          "blue",  # tmp_protein_2_colors[i],
          self._protein_pair.protein_2.pymol_selection.selection_string,
      )
      # self._protein_pair.protein_2.pymol_selection.color_selection(tmp_protein_2_colors[i])
      i += 1

  # </editor-fold>

  def __slot_export_data(self) -> None:
    """Handles the export data functionality."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Export data' button was clicked.",
    )
    file_dialog = QtWidgets.QFileDialog()
    desktop_path = QtCore.QStandardPaths.standardLocations(
        QtCore.QStandardPaths.DesktopLocation
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
      tmp_dialog.exec_()
