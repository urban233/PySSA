import os
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from PyQt5 import QtGui
from pyssa.controller import interface_manager
from pyssa.gui.ui.dialogs import dialog_advanced_prediction_configurations
from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.gui.ui.styles import styles
from pyssa.gui.ui.views import results_view, plot_view
from pyssa.internal.data_structures import protein_pair
from pyssa.internal.data_structures.data_classes import prediction_protein_info, prediction_configuration
from pyssa.internal.thread import tasks
from pyssa.io_pyssa import safeguard
from pyssa.presenter import main_presenter_async
from pyssa.util import gui_utils, tools, constants, exit_codes, prediction_util


class ResultsViewController(QtCore.QObject):

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager", the_protein_pair: protein_pair.ProteinPair):
        super().__init__()
        self._interface_manager = the_interface_manager
        self._protein_pair = the_protein_pair
        self._view: "results_view.ResultsView" = the_interface_manager.get_results_view()
        self.cb_protein_pair_color = QtWidgets.QComboBox()
        self._connect_all_ui_elements_to_slot_functions()
        self._build_table_widget()
        self._fill_table_widget()
        self._view.setWindowTitle(f"Results Summary Of {self._protein_pair.name}")

    def _connect_all_ui_elements_to_slot_functions(self):
        self._view.ui.btn_view_plots.clicked.connect(self._open_plot_view)
        #self._view.ui.btn_save_data.clicked.connect()

    def _build_table_widget(self):
        self._view.ui.table_widget_results.clear()
        self._view.ui.table_widget_results.setRowCount(3)
        self._view.ui.table_widget_results.setColumnCount(2)
        self._view.ui.table_widget_results.verticalHeader().setVisible(False)
        self._view.ui.table_widget_results.setHorizontalHeaderLabels(["Name", "Value"])
        gui_utils.fill_combo_box(self.cb_protein_pair_color, ["Normal", "By RMSD"])
        self.cb_protein_pair_color.adjustSize()

    def _fill_table_widget(self):
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
        # Color protein pair value item for table widget
        tmp_color_protein_pair_label_item = QtWidgets.QTableWidgetItem("Color Protein Pair")
        self._view.ui.table_widget_results.setItem(2, 0, tmp_color_protein_pair_label_item)
        self._view.ui.table_widget_results.setCellWidget(2, 1, self.cb_protein_pair_color)
        # Set editing flags
        tmp_result_rmsd_item.setFlags(tmp_result_rmsd_item.flags() & ~Qt.ItemIsEditable)
        tmp_rmsd_label_item.setFlags(tmp_rmsd_label_item.flags() & ~Qt.ItemIsEditable)
        tmp_aligned_aa_label_item.setFlags(tmp_aligned_aa_label_item.flags() & ~Qt.ItemIsEditable)
        tmp_result_aligned_aa_item.setFlags(tmp_result_aligned_aa_item.flags() & ~Qt.ItemIsEditable)
        tmp_color_protein_pair_label_item.setFlags(tmp_color_protein_pair_label_item.flags() & ~Qt.ItemIsEditable)
        # Resize
        self._view.ui.table_widget_results.resizeColumnToContents(0)
        self._view.ui.table_widget_results.resizeColumnToContents(1)

    def _open_plot_view(self) -> None:
        tmp_dialog = plot_view.PlotView(self._protein_pair,
                                        self._interface_manager.get_current_project(),
                                        self._protein_pair)
        tmp_dialog.exec_()
