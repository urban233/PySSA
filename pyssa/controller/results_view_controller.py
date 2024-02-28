import os
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from PyQt5 import QtGui
from pymol import cmd
from pyssa.controller import interface_manager, database_manager, pymol_session_manager
from pyssa.gui.ui.custom_dialogs import custom_message_box
from pyssa.gui.ui.dialogs import dialog_advanced_prediction_configurations
from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.gui.ui.styles import styles
from pyssa.gui.ui.views import results_view, plot_view
from pyssa.internal.data_structures import protein_pair, selection
from pyssa.internal.data_structures.data_classes import prediction_protein_info, prediction_configuration
from pyssa.internal.thread import tasks
from pyssa.io_pyssa import safeguard
from pyssa.presenter import main_presenter_async
from pyssa.util import gui_utils, tools, constants, exit_codes, prediction_util, enums


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

        if not self._pymol_session_manager.is_the_current_protein_pair_in_session():
            self.cb_protein_pair_color.setEnabled(False)

        self._connect_all_ui_elements_to_slot_functions()
        self._build_table_widget()
        self._fill_table_widget()
        self._view.setWindowTitle(f"Results Summary Of {self._protein_pair.name}")

    def _connect_all_ui_elements_to_slot_functions(self):
        self._view.ui.btn_view_plots.clicked.connect(self._open_plot_view)
        self.cb_protein_pair_color.currentIndexChanged.connect(self.__slot_color_protein_pair)
        self._view.ui.btn_save_data.clicked.connect(self.__slot_export_data)
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
        self._distance_data_visualizer = plot_view.PlotView(self._protein_pair,
                                        self._interface_manager.get_current_project(),
                                        self._protein_pair)
        self._distance_data_visualizer.show()

    # <editor-fold desc="Methods for coloring protein pair">
    def __slot_color_protein_pair(self) -> None:
        if self.cb_protein_pair_color.currentText().find("By RMSD") != -1:
            self._color_protein_pair_by_rmsd()
        else:
            self._color_protein_back_to_normal()

    def _get_color_configuration_of_protein_pair(self):
        tmp_database_filepath: str = str(self._interface_manager.get_current_project().get_database_filepath())
        with database_manager.DatabaseManager(tmp_database_filepath) as db_manager:
            db_manager.open_project_database()
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
            db_manager.close_project_database()

    def _color_protein_pair_by_rmsd(self) -> None:
        """Colors the residues in 5 colors depending on their distance to the reference."""
        #self._view.wait_spinner.start()
        self._active_task = tasks.Task(
            target=main_presenter_async.color_protein_pair_by_rmsd_value,
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
        #self._view.wait_spinner.stop()

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
        file_dialog = QtWidgets.QFileDialog()
        desktop_path = QtCore.QStandardPaths.standardLocations(QtCore.QStandardPaths.DesktopLocation)[0]
        file_dialog.setDirectory(desktop_path)
        file_path, _ = file_dialog.getSaveFileName(
            self._view,
            "Save distance data",
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
