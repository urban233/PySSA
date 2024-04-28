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
"""Module for the Results Dialog."""
import logging
import os
import subprocess

import pygetwindow
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from PyQt5 import QtGui
from pymol import cmd
from pyssa.controller import interface_manager, database_manager, pymol_session_manager
from pyssa.gui.ui.custom_dialogs import custom_message_box
from pyssa.gui.ui.dialogs import dialog_advanced_prediction_configurations
from pyssa.gui.ui.styles import styles
from pyssa.gui.ui.views import results_view, plot_view
from pyssa.internal.data_structures import protein_pair, selection
from pyssa.internal.data_structures.data_classes import prediction_protein_info, prediction_configuration
from pyssa.internal.thread import tasks
from pyssa.internal.thread.async_pyssa import util_async, protein_pair_async
from pyssa.io_pyssa import safeguard
from pyssa.presenter import main_presenter_async
from pyssa.util import gui_utils, tools, constants, exit_codes, prediction_util, enums
from pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class ResultsViewController(QtCore.QObject):

    def __init__(self,
                 the_interface_manager: "interface_manager.InterfaceManager",
                 the_protein_pair: "protein_pair.ProteinPair",
                 the_pymol_session_manager: "pymol_session_manager.PymolSessionManager"):
        super().__init__()
        self._interface_manager = the_interface_manager
        self._pymol_session_manager = the_pymol_session_manager
        self._protein_pair = the_protein_pair
        self._color_configuration_protein_pair: tuple[list, list] = ([], [])
        self._view: "results_view.ResultsView" = the_interface_manager.get_results_view()
        self._distance_data_visualizer = None
        self.cb_protein_pair_color = QtWidgets.QComboBox()

        if not self._pymol_session_manager.is_the_current_protein_pair_in_session(self._interface_manager.get_current_active_protein_pair_object().name):
            self.cb_protein_pair_color.setEnabled(False)

        self._connect_all_ui_elements_to_slot_functions()
        self._build_table_widget()
        self._fill_table_widget()
        self._check_results()
        self._view.setWindowTitle(f"Results Summary Of {self._protein_pair.name}")

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
        self.open_help("help/results/summary/")

    def _connect_all_ui_elements_to_slot_functions(self):
        self._view.ui.btn_view_plots.clicked.connect(self._open_plot_view)
        self.cb_protein_pair_color.currentIndexChanged.connect(self.__slot_color_protein_pair)
        self._view.ui.btn_export_data.clicked.connect(self.__slot_export_data)
        self._view.ui.table_widget_results.clicked.connect(self.__slot_highlight_protein)

    # <editor-fold desc="Util methods">
    def _build_table_widget(self):
        self.cb_protein_pair_color.currentIndexChanged.disconnect(self.__slot_color_protein_pair)
        self._view.ui.table_widget_results.clear()
        self._view.ui.table_widget_results.setRowCount(5)
        self._view.ui.table_widget_results.setColumnCount(2)
        self._view.ui.table_widget_results.verticalHeader().setVisible(False)
        self._view.ui.table_widget_results.setHorizontalHeaderLabels(["Name", "Value"])
        gui_utils.fill_combo_box(self.cb_protein_pair_color, ["Normal", "By RMSD"])
        self.cb_protein_pair_color.adjustSize()
        self.cb_protein_pair_color.currentIndexChanged.connect(self.__slot_color_protein_pair)

    def _fill_table_widget(self):
        self._get_color_configuration_of_protein_pair()
        # RMSD value item for table widget
        tmp_rmsd_label_item = QtWidgets.QTableWidgetItem("RMSD (Å)")
        tmp_result_rmsd_item = QtWidgets.QTableWidgetItem(
            str(self._protein_pair.distance_analysis.analysis_results.rmsd))
        self._view.ui.table_widget_results.setItem(0, 0, tmp_rmsd_label_item)
        self._view.ui.table_widget_results.setItem(0, 1, tmp_result_rmsd_item)
        # Aligned residues value item for table widget
        tmp_aligned_aa_label_item = QtWidgets.QTableWidgetItem("Aligned Residues")
        tmp_result_aligned_aa_item = QtWidgets.QTableWidgetItem(
            str(self._protein_pair.distance_analysis.analysis_results.aligned_aa))
        self._view.ui.table_widget_results.setItem(1, 0, tmp_aligned_aa_label_item)
        self._view.ui.table_widget_results.setItem(1, 1, tmp_result_aligned_aa_item)
        # Color information about protein 1 & 2
        tmp_color_description_protein_1_item = QtWidgets.QTableWidgetItem(
            f"Color {self._protein_pair.protein_1.get_molecule_object()}")
        self._view.ui.table_widget_results.setItem(2, 0, tmp_color_description_protein_1_item)
        tmp_color_protein_1_item = QtWidgets.QTableWidgetItem("green (default)")
        self._view.ui.table_widget_results.setItem(2, 1, tmp_color_protein_1_item)
        tmp_color_description_protein_2_item = QtWidgets.QTableWidgetItem(
            f"Color {self._protein_pair.protein_2.get_molecule_object()}")
        self._view.ui.table_widget_results.setItem(3, 0, tmp_color_description_protein_2_item)
        tmp_color_protein_2_item = QtWidgets.QTableWidgetItem("blue (default)")
        self._view.ui.table_widget_results.setItem(3, 1, tmp_color_protein_2_item)
        # Color protein pair value item for table widget
        tmp_color_protein_pair_label_item = QtWidgets.QTableWidgetItem("Color Protein Pair")
        self._view.ui.table_widget_results.setItem(4, 0, tmp_color_protein_pair_label_item)
        self._view.ui.table_widget_results.setCellWidget(4, 1, self.cb_protein_pair_color)
        # Set editing flags
        tmp_color_description_protein_1_item.setFlags(tmp_color_description_protein_1_item.flags() & ~Qt.ItemIsEditable)
        tmp_color_protein_1_item.setFlags(tmp_color_protein_1_item.flags() & ~Qt.ItemIsEditable)
        tmp_color_description_protein_2_item.setFlags(tmp_color_description_protein_2_item.flags() & ~Qt.ItemIsEditable)
        tmp_color_protein_2_item.setFlags(tmp_color_protein_2_item.flags() & ~Qt.ItemIsEditable)
        tmp_result_rmsd_item.setFlags(tmp_result_rmsd_item.flags() & ~Qt.ItemIsEditable)
        tmp_rmsd_label_item.setFlags(tmp_rmsd_label_item.flags() & ~Qt.ItemIsEditable)
        tmp_aligned_aa_label_item.setFlags(tmp_aligned_aa_label_item.flags() & ~Qt.ItemIsEditable)
        tmp_result_aligned_aa_item.setFlags(tmp_result_aligned_aa_item.flags() & ~Qt.ItemIsEditable)
        tmp_color_protein_pair_label_item.setFlags(tmp_color_protein_pair_label_item.flags() & ~Qt.ItemIsEditable)
        # Resize
        self._view.ui.table_widget_results.resizeColumnToContents(0)
        self._view.ui.table_widget_results.resizeColumnToContents(1)

    def _check_results(self):
        if self._protein_pair.distance_analysis.analysis_results.rmsd == 0:
            self._view.ui.btn_view_plots.setEnabled(False)
            self._view.ui.btn_view_plots.setToolTip("The RMSD value is exact 0, therefore no plots can be displayed.")
        else:
            self._view.ui.btn_view_plots.setEnabled(True)
            self._view.ui.btn_view_plots.setToolTip("")

    # </editor-fold>

    def __slot_highlight_protein(self):
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
        pass

    def _open_plot_view(self) -> None:
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'View plots' button was clicked.")
        self._distance_data_visualizer = plot_view.PlotView(self._protein_pair,
                                                            self._interface_manager.get_current_project(),
                                                            self._protein_pair)
        self._distance_data_visualizer.show()

    # <editor-fold desc="Methods for coloring protein pair">
    def __slot_color_protein_pair(self) -> None:
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'By RMSD' checkbox was clicked.")
        if self.cb_protein_pair_color.currentText().find("By RMSD") != -1:
            self._color_protein_pair_by_rmsd()
        else:
            self._color_protein_back_to_normal()

    def _get_color_configuration_of_protein_pair(self):
        tmp_database_filepath: str = str(self._interface_manager.get_current_project().get_database_filepath())
        with database_manager.DatabaseManager(tmp_database_filepath) as db_manager:
            for tmp_chain in self._protein_pair.protein_1.chains:
                self._color_configuration_protein_pair[0].append(
                    db_manager.get_pymol_parameter_for_certain_protein_chain_in_protein_pair(
                        self._protein_pair.get_id(),
                        self._protein_pair.protein_1.get_id(),
                        tmp_chain.chain_letter,
                        enums.PymolParameterEnum.COLOR.value
                    )[0]  # the [0] is needed because the db return the color as tuple e.g. ('green',)
                )
            for tmp_chain in self._protein_pair.protein_2.chains:
                self._color_configuration_protein_pair[1].append(
                    db_manager.get_pymol_parameter_for_certain_protein_chain_in_protein_pair(
                        self._protein_pair.get_id(),
                        self._protein_pair.protein_2.get_id(),
                        tmp_chain.chain_letter,
                        enums.PymolParameterEnum.COLOR.value
                    )[0]  # the [0] is needed because the db return the color as tuple e.g. ('green',)
                )

    def _color_protein_pair_by_rmsd(self) -> None:
        """Colors the residues in 5 colors depending on their distance to the reference."""
        self._active_task = tasks.Task(
            target=protein_pair_async.color_protein_pair_by_rmsd_value,
            args=(
                self._protein_pair,
                0
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
        pass

    def _color_protein_back_to_normal(self):
        self._get_color_configuration_of_protein_pair()
        i = 0
        tmp_protein_1_colors: list = self._color_configuration_protein_pair[0]
        for tmp_chain in self._protein_pair.protein_1.chains:
            self._protein_pair.protein_1.pymol_selection.set_selection_for_a_single_chain(tmp_chain.chain_letter)
            self._protein_pair.protein_1.pymol_selection.color_selection(tmp_protein_1_colors[i])
            i += 1

        i = 0
        tmp_protein_2_colors: list = self._color_configuration_protein_pair[1]
        for tmp_chain in self._protein_pair.protein_2.chains:
            self._protein_pair.protein_2.pymol_selection.set_selection_for_a_single_chain(tmp_chain.chain_letter)
            self._protein_pair.protein_2.pymol_selection.color_selection(tmp_protein_2_colors[i])
            i += 1

    # </editor-fold>

    def __slot_export_data(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Export data' button was clicked.")
        file_dialog = QtWidgets.QFileDialog()
        desktop_path = QtCore.QStandardPaths.standardLocations(QtCore.QStandardPaths.DesktopLocation)[0]
        file_dialog.setDirectory(desktop_path)
        file_path, _ = file_dialog.getSaveFileName(
            self._view,
            "Export Distance Data",
            "",
            "Comma-Separated Values (*.csv)",
        )
        if file_path:
            self._protein_pair.distance_analysis.analysis_results.export_distance_data_as_csv(file_path)
            if os.path.exists(file_path):
                tmp_dialog = custom_message_box.CustomMessageBoxOk(
                    "Export data as .csv file finished.", "Export Data",
                    custom_message_box.CustomMessageBoxIcons.INFORMATION.value
                )
            else:
                tmp_dialog = custom_message_box.CustomMessageBoxOk(
                    "Export data as .csv file failed!", "Export Data",
                    custom_message_box.CustomMessageBoxIcons.ERROR.value
                )
            tmp_dialog.exec_()
