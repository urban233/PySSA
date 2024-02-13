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
"""Module for the distance plot dialog."""
import copy
import csv
import os
import pathlib

import numpy as np
from PyQt5.QtCore import Qt
from PyQt5 import QtGui
from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt import NavigationToolbar2QT
from pymol import cmd
import pyssa.gui.ui.styles.styles as custom_pyssa_styles
from PyQt5 import QtCore

from pyssa.gui.ui.styles import styles
from pyssa.gui.ui.views import histogram_properties_view
from pyssa.internal.data_structures import protein_pair, selection
from pyssa.internal.thread import tasks
from pyssa.internal.thread.async_pyssa import util_async
from pyssa.util import pyssa_keys, enums
from pyssa.util import constants
from PyQt5.QtWidgets import QVBoxLayout, QWidget, QScrollArea
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backend_bases import MouseButton
from matplotlib.figure import Figure
from matplotlib import ticker


class PlotWidget(QWidget):
    """Class for a custom QWidget for plotting."""

    def __init__(self, parent=None) -> None:  # noqa: ANN001
        """Constructor.

        Args:
            parent: the parent
        """
        super(PlotWidget, self).__init__(parent)
        self.figure = Figure(figsize=(9, 5.25))
        self.canvas = FigureCanvas(self.figure)
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)

    def set_figure_size(self, width: float, height: float):
        self.figure.set_size_inches(width, height)


class PlotView(QtWidgets.QDialog):
    def __init__(self, protein_pair_from_project: "protein_pair.ProteinPair", a_project, the_protein_pair,
                 parent=None) -> None:  # noqa: ANN001
        """Constructor.

        Args:
            protein_pair_from_project: the protein pair which should be used for the distance histogram.
            parent: the parent.
        """
        QtWidgets.QDialog.__init__(self, parent)
        self._protein_pair = the_protein_pair
        self._current_project = a_project
        self.clicked_point_scatter = None
        self.highlighted_bin_index = None
        self.active_row_information = None
        self._sync_with_pymol_flag = False
        self.bars = None
        self._histogram_properties = {
            enums.HistogramPropertiesEnum.X_AXIS_UNITS: 10,
            enums.HistogramPropertiesEnum.DISTANCE_INTERVAL: 1.0,
        }

        # self.resizeEvent = self.handle_resize
        # Create a timer for delayed updates
        self.resize_timer = QtCore.QTimer(self)
        self.resize_timer.setInterval(250)  # Set the interval in milliseconds
        self.resize_timer.setSingleShot(True)

        self.protein_pair_for_analysis: protein_pair.ProteinPair = protein_pair_from_project
        #custom_pyssa_styles.set_stylesheet(self)

        self._initialize_ui()

        self.setModal(False)
        # --BEGIN

        # <editor-fold desc="Distance table logic">
        self.csv_model = QtGui.QStandardItemModel()
        self.csv_model.setColumnCount(7)
        labels = [
            "Residue pair no.",
            "Protein 1 Chain",
            "Protein 1 Position",
            "Protein 1 Residue",
            "Protein 2 Chain",
            "Protein 2 Position",
            "Protein 2 Residue",
            "Distance in Å",
        ]
        self.csv_model.setHorizontalHeaderLabels(labels)
        self.table_view.setModel(self.csv_model)

        csv_filepath = pathlib.Path(f"{constants.CACHE_CSV_DIR}/{self._protein_pair.name}.csv")
        if not os.path.exists(constants.CACHE_CSV_DIR):
            os.mkdir(constants.CACHE_CSV_DIR)

        distance_data = self._protein_pair.distance_analysis.analysis_results.distance_data
        distance_data_array = np.array(
            [
                distance_data[pyssa_keys.ARRAY_DISTANCE_INDEX],
                distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN],
                distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION],
                distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI],
                distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN],
                distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION],
                distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI],
                distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES],
            ],
        )
        distance_data_array_transpose = distance_data_array.transpose()
        with open(csv_filepath, mode="w", newline="") as file:
            writer = csv.writer(file, delimiter=",")
            writer.writerows(distance_data_array_transpose)

        file.close()

        with open(csv_filepath, "r", encoding="utf-8") as csv_file:
            i = 0
            for line in csv_file:
                tmp_list = line.split(",")
                # tmp_list.pop(0)
                standard_item_list = []
                pair_no_item = QtGui.QStandardItem()
                pair_no_item.setData(int(tmp_list[0]), role=QtCore.Qt.DisplayRole)
                ref_chain_item = QtGui.QStandardItem()
                ref_chain_item.setData(str(tmp_list[1]), role=QtCore.Qt.DisplayRole)
                ref_pos_item = QtGui.QStandardItem()
                ref_pos_item.setData(int(tmp_list[2]), role=QtCore.Qt.DisplayRole)
                ref_resi_item = QtGui.QStandardItem()
                ref_resi_item.setData(str(tmp_list[3]), role=QtCore.Qt.DisplayRole)
                model_chain_item = QtGui.QStandardItem()
                model_chain_item.setData(str(tmp_list[4]), role=QtCore.Qt.DisplayRole)
                model_pos_item = QtGui.QStandardItem()
                model_pos_item.setData(int(tmp_list[5]), role=QtCore.Qt.DisplayRole)
                model_resi_item = QtGui.QStandardItem()
                model_resi_item.setData(str(tmp_list[6]), role=QtCore.Qt.DisplayRole)
                distance_item = QtGui.QStandardItem()
                distance_item.setData(float(tmp_list[7]), role=QtCore.Qt.DisplayRole)
                standard_item_list.append(pair_no_item)
                standard_item_list.append(ref_chain_item)
                standard_item_list.append(ref_pos_item)
                standard_item_list.append(ref_resi_item)
                standard_item_list.append(model_chain_item)
                standard_item_list.append(model_pos_item)
                standard_item_list.append(model_resi_item)
                standard_item_list.append(distance_item)
                self.csv_model.insertRow(i, standard_item_list)
            i += 1
        csv_file.close()
        self.csv_model.removeRow(0)
        self.table_view.setAlternatingRowColors(True)
        self.table_view.resizeColumnsToContents()
        self.table_view.verticalHeader().setVisible(False)
        self.table_view.setSortingEnabled(True)
        self.table_view.sortByColumn(0, QtCore.Qt.AscendingOrder)
        self.table_view.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)  # disables editing of cells
        # </editor-fold>

        # --END

        self.create_all_graphics()
        self.delayed_resize()
        self.setup_context_menu()
        self._connect_all_signals()

    def _initialize_ui(self):
        # <editor-fold desc="Define basic ui elements">
        self.scroll_area = QScrollArea()
        self.plot_widget_dplot = PlotWidget()
        self.plot_widget_dhistogram = PlotWidget()
        self.menubar = QtWidgets.QMenuBar(self)
        self.lbl_status = QtWidgets.QLabel()
        self.lbl_status.setText("Status: ")

        # </editor-fold>

        # <editor-fold desc="Create a menu and add it to the menu bar">
        plots_menu = QtWidgets.QMenu('Views', self)
        table_menu = QtWidgets.QMenu('Table', self)
        hide_column_menu = QtWidgets.QMenu('Hide Column', self)
        show_column_menu = QtWidgets.QMenu('Show Column', self)
        help_menu = QtWidgets.QMenu('Help', self)
        self.menubar.addMenu(plots_menu)
        self.menubar.addMenu(table_menu)
        table_menu.addMenu(hide_column_menu)
        table_menu.addMenu(show_column_menu)
        self.menubar.addMenu(help_menu)

        # </editor-fold>

        # <editor-fold desc="Create actions and add them to the menu">
        self.action_plot = QtWidgets.QAction('Plot', self)
        self.action_plot.setCheckable(True)
        self.action_plot.setChecked(True)
        self.action_histogram = QtWidgets.QAction('Histogram', self)
        self.action_histogram.setCheckable(True)
        self.action_histogram.setChecked(True)

        self.action_table = QtWidgets.QAction('Table', self)
        self.action_table.setCheckable(True)
        self.action_table.setChecked(True)

        self.action_sync_with_pymol = QtWidgets.QAction('Sync With PyMOL', self)
        self.action_sync_with_pymol.setCheckable(True)
        self.action_sync_with_pymol.setChecked(False)
        
        plots_menu.addAction(self.action_plot)
        plots_menu.addAction(self.action_histogram)
        plots_menu.addAction(self.action_table)
        plots_menu.addAction(self.action_sync_with_pymol)

        self.action_hide_residue_pair_no = QtWidgets.QAction("Residue Pair No.")
        self.action_hide_protein_1_chain = QtWidgets.QAction("Protein 1 Chain")
        self.action_hide_protein_1_position = QtWidgets.QAction("Protein 1 Position")
        self.action_hide_protein_1_residue = QtWidgets.QAction("Protein 1 Residue")
        self.action_hide_protein_2_chain = QtWidgets.QAction("Protein 2 Chain")
        self.action_hide_protein_2_position = QtWidgets.QAction("Protein 2 Position")
        self.action_hide_protein_2_residue = QtWidgets.QAction("Protein 2 Residue")
        self.action_hide_distance = QtWidgets.QAction("Distance")

        hide_column_menu.addAction(self.action_hide_residue_pair_no)
        hide_column_menu.addAction(self.action_hide_protein_1_chain)
        hide_column_menu.addAction(self.action_hide_protein_1_position)
        hide_column_menu.addAction(self.action_hide_protein_1_residue)
        hide_column_menu.addAction(self.action_hide_protein_2_chain)
        hide_column_menu.addAction(self.action_hide_protein_2_position)
        hide_column_menu.addAction(self.action_hide_protein_2_residue)
        hide_column_menu.addAction(self.action_hide_distance)

        self.action_show_all_columns = QtWidgets.QAction("All Columns")
        self.action_show_residue_pair_no = QtWidgets.QAction("Residue Pair No.")
        self.action_show_protein_1_chain = QtWidgets.QAction("Protein 1 Chain")
        self.action_show_protein_1_position = QtWidgets.QAction("Protein 1 Position")
        self.action_show_protein_1_residue = QtWidgets.QAction("Protein 1 Residue")
        self.action_show_protein_2_chain = QtWidgets.QAction("Protein 2 Chain")
        self.action_show_protein_2_position = QtWidgets.QAction("Protein 2 Position")
        self.action_show_protein_2_residue = QtWidgets.QAction("Protein 2 Residue")
        self.action_show_distance = QtWidgets.QAction("Distance")

        show_column_menu.addAction(self.action_show_all_columns)
        show_column_menu.addSeparator()
        show_column_menu.addAction(self.action_show_residue_pair_no)
        show_column_menu.addAction(self.action_show_protein_1_chain)
        show_column_menu.addAction(self.action_show_protein_1_position)
        show_column_menu.addAction(self.action_show_protein_1_residue)
        show_column_menu.addAction(self.action_show_protein_2_chain)
        show_column_menu.addAction(self.action_show_protein_2_position)
        show_column_menu.addAction(self.action_show_protein_2_residue)
        show_column_menu.addAction(self.action_show_distance)

        self.action_docs = QtWidgets.QAction('PySSA Documentation', self)
        self.action_docs.triggered.connect(self._open_help_center)
        self.action_help = QtWidgets.QAction('Get Help', self)
        self.action_help.triggered.connect(self._open_distance_data_visualizer_help)

        help_menu.addAction(self.action_docs)
        help_menu.addAction(self.action_help)

        # </editor-fold>

        # self.toolbar = NavigationToolbar2QT(self.plot_widget.canvas, self)
        # Find the action you want to remove
        # items_to_remove = ["Home", "Back", "Forward", "Pan", "Zoom", "Subplots"]
        # for tmp_item in items_to_remove:
        #     for action in self.toolbar.actions():
        #         print(action.text())
        #         if action.text() == tmp_item:
        #             self.toolbar.removeAction(action)

        # <editor-fold desc="Set layouts">
        self.scroll_area.setWidget(self.plot_widget_dhistogram)

        # Create labels
        self.lbl_status1 = QtWidgets.QLabel(f"Protein 1: {self.protein_pair_for_analysis.protein_1.get_molecule_object()}")
        self.lbl_status2 = QtWidgets.QLabel(f"Protein 2: {self.protein_pair_for_analysis.protein_2.get_molecule_object()}")

        # Create a QTableView
        self.table_view = QtWidgets.QTableView()
        self.table_view.setMinimumWidth(450)

        # Create a QWidget to hold labels and QTableView
        self.container_widget = QtWidgets.QWidget()
        container_layout = QtWidgets.QVBoxLayout(self.container_widget)
        container_layout.addWidget(self.lbl_status1)
        container_layout.addWidget(self.lbl_status2)
        container_layout.addWidget(self.table_view)

        stylesheet = """
                    QLabel {
                        background-color: white;
                        font-size: 12px;
                        padding: 5px;
                        border-style: solid;
                        border-width: 2px;
                        border-radius: 6px;
                        border-color: #DCDBE3;
                    }
                    """
        self.container_widget.setStyleSheet(stylesheet)
        # # Create a QHBoxLayout for the scroll area
        # self.scroll_area_layout = QtWidgets.QHBoxLayout()
        #
        # # Assuming self.plot_widget and self.scroll_area are already defined
        # self.scroll_area_layout.addWidget(self.plot_widget)
        # self.scroll_area_layout.addWidget(self.container_widget)
        #
        # # Set up scroll area
        # self.scroll_area = QtWidgets.QScrollArea()
        # self.scroll_area.setLayout(self.scroll_area_layout)
        #
        # # Create main layout
        # self.main_Layout = QtWidgets.QVBoxLayout()
        # self.main_Layout.addWidget(self.scroll_area)
        # self.main_Layout.addWidget(self.lbl_status)
        # self.main_Layout.setMenuBar(self.menubar)
        # self.setLayout(self.main_Layout)

        # Create the first splitter to divide the dialog into two sections
        self.vertical_splitter = QtWidgets.QSplitter()
        self.main_layout = QtWidgets.QVBoxLayout(self)
        self.main_layout.addWidget(self.vertical_splitter)

        # Left side (plot area)
        plot_area = QtWidgets.QFrame()
        self.vertical_splitter.addWidget(plot_area)
        # Right side (table)
        self.vertical_splitter.addWidget(self.container_widget)
        # Second splitter within the plot area to split it horizontally
        self.horizontal_splitter = QtWidgets.QSplitter()
        plot_area_layout = QVBoxLayout()
        plot_area.setLayout(plot_area_layout)
        plot_area_layout.addWidget(self.horizontal_splitter)
        # Left part of the plot area
        self.horizontal_splitter.addWidget(self.plot_widget_dplot)
        # Right part of the plot area
        self.horizontal_splitter.addWidget(self.scroll_area)
        self.horizontal_splitter.setOrientation(0)  # Set orientation to horizontal

        self.vertical_splitter.setCollapsible(0, False)
        self.main_layout.setMenuBar(self.menubar)
        self.main_layout.addWidget(self.lbl_status)
        self.setLayout(self.main_layout)

        # </editor-fold>

        # <editor-fold desc="Setup subplots">
        self._ax_plot = self.plot_widget_dplot.figure.add_subplot()
        self._ax_hist = self.plot_widget_dhistogram.figure.add_subplot()

        # </editor-fold>

        # <editor-fold desc="Set window styles">
        styles.set_stylesheet(self)
        stylesheet = """
                QDialog {background-color: #F6F4F8;}
                QTableWidget {background-color: white;}
                """
        self.setStyleSheet(stylesheet)
        self.resize(1400, 800)
        self.setWindowFlag(QtCore.Qt.WindowMaximizeButtonHint, True)
        self.setWindowFlag(QtCore.Qt.WindowCloseButtonHint, True)
        #self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        self.setWindowTitle("Distance Data Visualizer")

        # </editor-fold>

    # <editor-fold desc="Help related methods">
    def open_help(self, a_page_name: str):
        """Opens the pyssa documentation window if it's not already open.

        Args:
            a_page_name (str): a name of a documentation page to display
        """
        self._active_task = tasks.Task(
            target=util_async.open_documentation_on_certain_page,
            args=(a_page_name, 0),
            post_func=self.__await_open_help,
        )
        self._active_task.start()

    @staticmethod
    def __await_open_help():
        constants.PYSSA_LOGGER.info("Opening help center finished.")

    def _open_help_center(self):
        self.open_help("help/")

    def _open_distance_data_visualizer_help(self):
        self.open_help("help/results/distance_data_visualizer")
    # </editor-fold>

    def _connect_all_signals(self):
        self.action_table.triggered.connect(self.hide_distance_table)

        self.action_hide_residue_pair_no.triggered.connect(self.__action_hide_residue_pair_no)
        self.action_hide_protein_1_chain.triggered.connect(self.__action_hide_protein_1_chain)
        self.action_hide_protein_1_position.triggered.connect(self.__action_hide_protein_1_position)
        self.action_hide_protein_1_residue.triggered.connect(self.__action_hide_protein_1_residue)
        self.action_hide_protein_2_chain.triggered.connect(self.__action_hide_protein_2_chain)
        self.action_hide_protein_2_position.triggered.connect(self.__action_hide_protein_2_position)
        self.action_hide_protein_2_residue.triggered.connect(self.__action_hide_protein_2_residue)
        self.action_hide_distance.triggered.connect(self.__action_hide_distances)

        self.action_show_residue_pair_no.triggered.connect(self.__action_show_residue_pair_no)
        self.action_show_protein_1_chain.triggered.connect(self.__action_show_protein_1_chain)
        self.action_show_protein_1_position.triggered.connect(self.__action_show_protein_1_position)
        self.action_show_protein_1_residue.triggered.connect(self.__action_show_protein_1_residue)
        self.action_show_protein_2_chain.triggered.connect(self.__action_show_protein_2_chain)
        self.action_show_protein_2_position.triggered.connect(self.__action_show_protein_2_position)
        self.action_show_protein_2_residue.triggered.connect(self.__action_show_protein_2_residue)
        self.action_show_distance.triggered.connect(self.__action_show_distances)
        self.action_show_all_columns.triggered.connect(self.__action_show_all_columns)

        self.table_view.clicked.connect(self.highlight_table_selection_in_plot)
        self.action_plot.triggered.connect(self.move_horizontal_splitter)
        self.action_histogram.triggered.connect(self.move_horizontal_splitter)
        self.action_sync_with_pymol.triggered.connect(self.__slot_sync_with_pymol)

        self.resize_timer.timeout.connect(self.actual_resize)
        self.vertical_splitter.splitterMoved.connect(self.delayed_resize)
        self.horizontal_splitter.splitterMoved.connect(self.delayed_resize)

        self.plot_widget_dplot.canvas.mpl_connect('button_press_event', self.on_canvas_click)
        #self.plot_widget_dhistogram.canvas.mpl_connect('button_press_event', self.__slot_open_context_menu_for_histogram)
        #self.plot_widget_dplot.canvas.mpl_connect('motion_notify_event', self._on_move)

    # <editor-fold desc="QSplitter related methods">
    def hide_distance_table(self):
        if self.action_table.isChecked():
            tmp_index = self.vertical_splitter.indexOf(self.container_widget)
            tmp_lowest, tmp_highest = self.vertical_splitter.getRange(tmp_index)
            self.vertical_splitter.moveSplitter(tmp_highest - 450, tmp_index)
        else:
            tmp_index = self.vertical_splitter.indexOf(self.container_widget)
            tmp_lowest, tmp_highest = self.vertical_splitter.getRange(tmp_index)
            self.vertical_splitter.moveSplitter(tmp_highest, tmp_index)

    def move_horizontal_splitter(self):
        if self.action_plot.isChecked() and self.action_histogram.isChecked():
            # Both should be displayed
            tmp_index = self.horizontal_splitter.indexOf(self.plot_widget_dhistogram)
            tmp_lowest, tmp_highest = self.horizontal_splitter.getRange(tmp_index)
            self.horizontal_splitter.moveSplitter(tmp_lowest + 400, tmp_index)
        elif self.action_plot.isChecked() and not self.action_histogram.isChecked():
            # Only the plot should be displayed
            tmp_index = self.horizontal_splitter.indexOf(self.plot_widget_dhistogram)
            tmp_lowest, tmp_highest = self.horizontal_splitter.getRange(tmp_index)
            self.horizontal_splitter.moveSplitter(tmp_highest, tmp_index)
        elif not self.action_plot.isChecked() and self.action_histogram.isChecked():
            # Only the histogram should be displayed
            tmp_index = self.horizontal_splitter.indexOf(self.plot_widget_dhistogram)
            tmp_lowest, tmp_highest = self.horizontal_splitter.getRange(tmp_index)
            self.horizontal_splitter.moveSplitter(tmp_lowest, tmp_index)
        elif not self.action_plot.isChecked() and not self.action_histogram.isChecked() and self.action_table.isChecked():
            # No plots should be displayed
            tmp_index = self.vertical_splitter.indexOf(self.container_widget)
            tmp_lowest, tmp_highest = self.vertical_splitter.getRange(tmp_index)
            self.vertical_splitter.moveSplitter(tmp_lowest, tmp_index)
        elif not self.action_plot.isChecked() and not self.action_histogram.isChecked() and not self.action_table.isChecked():
            tmp_index = self.vertical_splitter.indexOf(self.container_widget)
            tmp_lowest, tmp_highest = self.vertical_splitter.getRange(tmp_index)
            self.vertical_splitter.moveSplitter(tmp_highest, tmp_index)

    def delayed_resize(self):
        # Start or restart the timer when the splitter is moved
        self.plot_widget_dplot.figure.clear()
        self.plot_widget_dhistogram.figure.clear()
        self.resize_timer.start()

    def actual_resize(self):
        """Perform the actual resizing operation. This method will be called after the timer interval has elapsed"""
        print("Resizing starts now ...")
        if self.container_widget.size().width() == 0:
            self.action_table.setChecked(False)
        else:
            self.action_table.setChecked(True)

        if self.plot_widget_dplot.size().height() == 0:
            self.action_plot.setChecked(False)
        else:
            self.action_plot.setChecked(True)

        if self.plot_widget_dhistogram.size().height() == 0:
            self.action_histogram.setChecked(False)
        elif self.plot_widget_dhistogram.size().width() == 0:
            self.action_plot.setChecked(False)
            self.action_histogram.setChecked(False)
        else:
            self.action_histogram.setChecked(True)

        self.toggle_graphics_visibility()

    def toggle_graphics_visibility(self):
        self.plot_widget_dplot.figure.clear()
        self.plot_widget_dhistogram.figure.clear()
        if self.action_plot.isChecked() and self.action_histogram.isChecked():
            print(self.scroll_area.size())
            print(self.scroll_area.width() / 100)
            tmp_histogram_width = self.scroll_area.width() - 20
            tmp_histogram_height = ((5/6) * len(self.bars)) * 100
            self.plot_widget_dhistogram.resize(tmp_histogram_width, tmp_histogram_height)
            self.plot_widget_dhistogram.set_figure_size(tmp_histogram_width / 100, tmp_histogram_height/ 100)

            self.plot_widget_dplot.show()
            self.plot_widget_dhistogram.show()

            self._ax_plot = self.plot_widget_dplot.figure.subplots()
            self._ax_hist = self.plot_widget_dhistogram.figure.subplots()
            self.create_distance_plot()
            self.setup_plot_defaults()
            self.create_distance_histogram()
            self.setup_histogram_defaults()

            self.plot_widget_dplot.figure.tight_layout()
            self.plot_widget_dplot.canvas.draw()
            try:  # TODO: this is not an ideal way, but I didn't find anything better
                self.plot_widget_dhistogram.figure.tight_layout()
                self.plot_widget_dhistogram.canvas.draw()
            except np.linalg.LinAlgError:
                print("A layout cannot be applied to the histogram.")

        elif self.action_plot.isChecked() and not self.action_histogram.isChecked():
            self.plot_widget_dplot.show()

            self._ax_plot = self.plot_widget_dplot.figure.subplots()
            self.create_distance_plot()
            self.setup_plot_defaults()

            self.plot_widget_dplot.figure.tight_layout()
            self.plot_widget_dplot.canvas.draw()

        elif not self.action_plot.isChecked() and self.action_histogram.isChecked():
            self.plot_widget_dhistogram.show()
            self._ax_hist = self.plot_widget_dhistogram.figure.subplots()
            self.create_distance_histogram()
            self.setup_histogram_defaults()

            try:  # TODO: this is not an ideal way, but I didn't find anything better
                self.plot_widget_dhistogram.figure.tight_layout()
                self.plot_widget_dhistogram.canvas.draw()
            except np.linalg.LinAlgError:
                print("A layout cannot be applied to the histogram.")
        else:
            tmp_index = self.vertical_splitter.indexOf(self.container_widget)
            tmp_lowest, tmp_highest = self.vertical_splitter.getRange(tmp_index)
            self.vertical_splitter.moveSplitter(tmp_lowest, tmp_index)
            self.plot_widget_dplot.hide()
            self.plot_widget_dhistogram.hide()

    # </editor-fold>

    def __slot_sync_with_pymol(self):
        if self.action_sync_with_pymol.isChecked():
            self._sync_with_pymol_flag = True
        else:
            self._sync_with_pymol_flag = False

    # <editor-fold desc="Plotting related methods">
    def create_all_graphics(self):
        self.plot_widget_dplot.show()
        self.plot_widget_dhistogram.show()
        #num_subplots = 2  # Adjust this number based on your requirements
        #self.plot_widget_dplot.figure.clear()  # Clear existing figure
        #self.plot_widget_dhistogram.figure.clear()

        # # Adjust figure layout based on the number of subplots
        # self.plot_widget.figure.subplots(num_subplots, 1)
        #
        # # Create and set up the subplots
        # self._ax_plot = self.plot_widget.figure.axes[0]  # First subplot
        # self.create_distance_plot()
        # self.setup_plot_defaults()
        #
        # self._ax_hist = self.plot_widget.figure.axes[1]  # Second subplot
        # self.create_distance_histogram()
        # self.setup_histogram_defaults()
        #
        # # Adjust layout
        # self.plot_widget.figure.tight_layout()
        # self.plot_widget.canvas.draw()

        self.create_distance_plot()
        self.setup_plot_defaults()
        self.create_distance_histogram()
        self.setup_histogram_defaults()
        self.plot_widget_dplot.figure.tight_layout()
        self.plot_widget_dhistogram.figure.tight_layout()
        self.plot_widget_dplot.canvas.draw()
        self.plot_widget_dhistogram.canvas.draw()

    def setup_plot_defaults(self):
        # Set axis labels
        self._ax_plot.set_xlabel("Residue Pair No.")
        self._ax_plot.set_ylabel("Distance In Å")
        self._ax_plot.set_title("Distance Plot")
        # Set the x-axis limits with minimum value set to 0
        self._ax_plot.set_xlim(0)
        # Set the number of ticks for the x and y axes
        self._ax_plot.xaxis.set_major_locator(ticker.MultipleLocator(10))

    def setup_histogram_defaults(self):
        # Move the entire x-axis to the top
        self._ax_hist.xaxis.tick_top()
        self._ax_hist.xaxis.set_label_position("top")
        self._ax_hist.set_title("Distance Histogram")
        # Set axis labels
        self._ax_hist.set_xlabel("Count")
        self._ax_hist.set_ylabel("Bins")
        # Invert the y-axis
        self._ax_hist.invert_yaxis()
        # Remove the spines where no ax is are present
        self._ax_hist.spines["right"].set_visible(False)
        self._ax_hist.spines["bottom"].set_visible(False)
        # Remove the spines where no axis are present
        self._ax_hist.spines["right"].set_visible(False)
        self._ax_hist.spines["bottom"].set_visible(False)
        # Set axis labels
        self._ax_hist.set_xlabel("Frequency Of C-α Distances")
        self._ax_hist.set_ylabel("Distance Interval In Å")
        self._ax_hist.xaxis.set_major_locator(ticker.MultipleLocator(self._histogram_properties[enums.HistogramPropertiesEnum.X_AXIS_UNITS]))

    def create_distance_plot(self):
        # data for actual distance plot line
        distance_data = self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data
        distance_list = distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES]
        # data for cutoff line
        # cutoff_line = []
        # for i in range(len(distance_list)):
        #     cutoff_line.append(self.protein_pair_for_analysis.distance_analysis.cutoff)
        self._ax_plot.plot(distance_list, color="#367AF6")

    def create_distance_histogram(self):
        distance_data: dict[
            str,
            np.ndarray,
        ] = self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data
        distance_list = copy.deepcopy(distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES])
        distance_list.sort()
        length = len(distance_list)
        max_distance = distance_list[length - 1]
        self.bins = np.arange(0, max_distance + 1, self._histogram_properties[enums.HistogramPropertiesEnum.DISTANCE_INTERVAL])
        hist, bin_edges = np.histogram(distance_list, bins=self.bins)
        
        frequencies = hist.tolist()
        self.freqs_without_zeros = []
        self.bins_without_zeros_label = []
        self.bins_without_zeros = []
        for i, freq in enumerate(frequencies):
            if freq != 0:
                self.freqs_without_zeros.append(freq)
                tmp_str_bin = str(float(bin_edges[i]))
                tmp_str_bin_2 = str(float(bin_edges[i]) + self._histogram_properties[enums.HistogramPropertiesEnum.DISTANCE_INTERVAL])
                self.bins_without_zeros.append(bin_edges[i])
                self.bins_without_zeros_label.append(f"[{tmp_str_bin},{tmp_str_bin_2}]")

        self.bars = self._ax_hist.barh(self.bins_without_zeros_label,
                                       self.freqs_without_zeros,
                                       color="#367AF6",
                                       height=0.6)
        self._ax_hist.bar_label(self.bars, padding=4)

        # bar_width = 1.0 / (num_bars + 1)  # Adjusting bar width based on number of bars
        #
        # fig_width = max(10, 2 * num_bars)  # Minimum figure width to ensure bars are visible
        # fig_height = 6  # Adjust this as needed
        #
        # fig, ax = plt.subplots(figsize=(fig_width, fig_height))  # Adjust figure size
        #
        # # Calculate x positions for bars
        # x_positions = [i * (1 + bar_width) for i in range(num_bars)]
        #
        # # Plot bars
        # ax.bar(x_positions, data, width=bar_width)
        #
        # ax.set_xticks([i + 0.5 * bar_width for i in x_positions])  # Adjusting x ticks position
        #
        # plt.show()

    def create_distance_histogram_old(self):
        distance_data: dict[
            str,
            np.ndarray,
        ] = self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data
        distance_list = copy.deepcopy(distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES])
        distance_list.sort()
        length = len(distance_list)
        max_distance = distance_list[length - 1]

        n, self.bins, self.patches = self._ax_hist.hist(
            distance_list,
            bins=np.arange(0, max_distance + 0.25, 0.25),
            orientation="horizontal",
            rwidth=0.7,
            color="#4B91F7",
        )

        # Add labels to the non-zero frequency histogram bars
        for bin_value, patch in zip(n, self.patches):
            x = patch.get_x() + patch.get_width() + 0.1
            y = patch.get_y() + patch.get_height() / 2
            self._ax_hist.annotate(
                f"{bin_value}",
                xy=(x, y),
                xycoords="data",
                xytext=(3, 0),
                textcoords="offset points",
                ha="left",
                va="center",
            )

            # Calculate the midpoints between bin edges
            bin_midpoints = (self.bins[:-1] + self.bins[1:]) / 2
            # Set y-ticks at the bin midpoints
            self._ax_hist.set_yticks(bin_midpoints)
            # Set custom tick labels
            custom_labels = [f"{bin_start} - {bin_end}" for bin_start, bin_end in zip(self.bins[:-1], self.bins[1:])]
            self._ax_hist.set_yticklabels(custom_labels)

    # </editor-fold>
    # def _on_move(self, event):
    #     if event.inaxes:
    #         print(f'data coords {event.xdata} {event.ydata},',
    #               f'pixel coords {event.x} {event.y}')

    def on_canvas_click(self, event):
        if event.button is MouseButton.LEFT:
            x_clicked, y_clicked = event.xdata, event.ydata
            # Clear the previous clicked point
            if self.clicked_point_scatter is not None:
                try:
                    self.clicked_point_scatter.remove()
                except Exception as e:
                    constants.PYSSA_LOGGER.warning(f"Something went wrong: {e}")

            # Find the nearest point on the line
            distance_data = self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES]
            x_line = range(len(distance_data))
            y_line = distance_data
            distances = np.sqrt((x_line - x_clicked) ** 2 + (y_line - y_clicked) ** 2)
            index_nearest = np.argmin(distances)
            x_nearest, y_nearest = x_line[index_nearest], y_line[index_nearest]
            self.clicked_point_scatter = self._ax_plot.scatter(x_nearest, y_nearest, color='red', marker='o', label='Clicked Point')
            self.plot_widget_dplot.canvas.draw()
            self.active_row_information = self.search_point_by_id(x_nearest)
            if self.active_row_information is not None:
                tmp_prot_1_name = self.protein_pair_for_analysis.protein_1.get_molecule_object()
                tmp_prot_2_name = self.protein_pair_for_analysis.protein_2.get_molecule_object()
                msg: str = (f"Status: Residue pair no. = {x_nearest}, Distance (Å) = {y_nearest} Between {tmp_prot_1_name} "
                            f"Chain {self.active_row_information[1]} Position {self.active_row_information[2]} Residue {self.active_row_information[3]} and "
                            f"{tmp_prot_2_name} Chain {self.active_row_information[4]} Position {self.active_row_information[5]} Residue {self.active_row_information[6]}")
                self.lbl_status.setText(msg)
                if self._sync_with_pymol_flag:
                    self.highlight_residue_in_pymol(x_nearest, tmp_prot_1_name)
            else:
                self.lbl_status.setText(f"Status: Residue pair no. = {x_nearest}, Distance (Å) = {y_nearest}")
            # Change the selection color
            #self.change_selection_color(QtGui.QColor(75, 145, 247, 200))
            self.highlight_histogram_bar(y_nearest)
        elif event.button is MouseButton.RIGHT:
            print("Open Context Menu ...")

    def change_selection_color(self, color):
        palette = self.table_view.palette()
        palette.setColor(QtGui.QPalette.Highlight, color)
        self.table_view.setPalette(palette)

    def search_point_by_id(self, target_id):
        for row in range(self.csv_model.rowCount()):
            id_item = self.csv_model.item(row, 0)  # Assuming running id is in the first column
            if id_item and id_item.text() == str(target_id):
                # Select the entire row in the QTableView
                self.table_view.selectRow(row)
                self.active_row_information = (
                    self.table_view.model().index(row, 0).data(Qt.DisplayRole),
                    self.table_view.model().index(row, 1).data(Qt.DisplayRole),
                    self.table_view.model().index(row, 2).data(Qt.DisplayRole),
                    self.table_view.model().index(row, 3).data(Qt.DisplayRole),
                    self.table_view.model().index(row, 4).data(Qt.DisplayRole),
                    self.table_view.model().index(row, 5).data(Qt.DisplayRole),
                    self.table_view.model().index(row, 6).data(Qt.DisplayRole),
                )
                return self.active_row_information  # Return True if the id is found
        return None

    def highlight_table_selection_in_plot(self):
        tmp_x_value = self.table_view.model().index(
            self.table_view.currentIndex().row(), 0
        ).data(Qt.DisplayRole)
        tmp_y_value = self.table_view.model().index(
            self.table_view.currentIndex().row(), 7
        ).data(Qt.DisplayRole)

        # Clear the previous clicked point
        if self.clicked_point_scatter is not None:
            try:
                self.clicked_point_scatter.remove()
            except Exception as e:
                constants.PYSSA_LOGGER.warning(f"Something went wrong: {e}")
        self.clicked_point_scatter = self._ax_plot.scatter(tmp_x_value, tmp_y_value, color='red', marker='o',
                                                           label='Clicked Point')
        self.plot_widget_dplot.canvas.draw()

        tmp_row = self.table_view.currentIndex().row()
        self.table_view.selectRow(tmp_row)
        self.active_row_information = (
            self.table_view.model().index(tmp_row, 0).data(Qt.DisplayRole),
            self.table_view.model().index(tmp_row, 1).data(Qt.DisplayRole),
            self.table_view.model().index(tmp_row, 2).data(Qt.DisplayRole),
            self.table_view.model().index(tmp_row, 3).data(Qt.DisplayRole),
            self.table_view.model().index(tmp_row, 4).data(Qt.DisplayRole),
            self.table_view.model().index(tmp_row, 5).data(Qt.DisplayRole),
            self.table_view.model().index(tmp_row, 6).data(Qt.DisplayRole),
        )
        if self.active_row_information is not None:
            tmp_prot_1_name = self.protein_pair_for_analysis.protein_1.get_molecule_object()
            tmp_prot_2_name = self.protein_pair_for_analysis.protein_2.get_molecule_object()
            msg: str = (f"Status: Residue pair no. = {tmp_x_value}, Distance (Å) = {tmp_y_value} Between {tmp_prot_1_name} "
                        f"Chain {self.active_row_information[1]} Position {self.active_row_information[2]} Residue {self.active_row_information[3]} and "
                        f"{tmp_prot_2_name} Chain {self.active_row_information[4]} Position {self.active_row_information[5]} Residue {self.active_row_information[6]}")
            self.lbl_status.setText(msg)
        else:
            self.lbl_status.setText(f"Status: Residue pair no. = {tmp_x_value}, Distance (Å) = {tmp_y_value}")
        self.highlight_histogram_bar(tmp_y_value)

    def highlight_histogram_bar(self, point_to_highlight):
        # Highlight a specific data point (example: point_to_highlight)
        if self.highlighted_bin_index is not None:
            self.bars[self.highlighted_bin_index].set_facecolor('#367AF6')

        for tmp_bin in self.bins:
            tmp_lower_bin = tmp_bin - 1
            tmp_upper_bin = tmp_bin
            if point_to_highlight >= tmp_lower_bin and point_to_highlight <= tmp_upper_bin:
                self.highlighted_bin_index = int(tmp_lower_bin)
                break
        self.highlighted_bin_index = self.bins_without_zeros.index(self.highlighted_bin_index)  # TODO: Could be cleaner
        if 0 <= self.highlighted_bin_index < len(self.bars):
            bar = self.bars[self.highlighted_bin_index]
            bar.set_facecolor("#2D5794")

        # Update the plot
        self.plot_widget_dhistogram.canvas.draw()

    def highlight_residue_in_pymol(self, the_current_id, tmp_prot_1_name):
        self.table_view.currentIndex()
        tmp_id, tmp_chain_1, tmp_pos_1, tmp_residue_1, _, _, _ = self.search_point_by_id(the_current_id)
        tmp_pymol_selection = selection.Selection(tmp_prot_1_name)
        tmp_pymol_selection.set_single_selection("", tmp_chain_1, tmp_residue_1, "")
        cmd.zoom(tmp_pymol_selection.selection_string)
        cmd.show("sticks", tmp_pymol_selection.selection_string)

    # <editor-fold desc="Distance table related methods for show and hide colums">
    def __action_hide_residue_pair_no(self):
        self.table_view.hideColumn(0)

    def __action_hide_protein_1_chain(self):
        self.table_view.hideColumn(1)

    def __action_hide_protein_1_position(self):
        self.table_view.hideColumn(2)

    def __action_hide_protein_1_residue(self):
        self.table_view.hideColumn(3)

    def __action_hide_protein_2_chain(self):
        self.table_view.hideColumn(4)

    def __action_hide_protein_2_position(self):
        self.table_view.hideColumn(5)

    def __action_hide_protein_2_residue(self):
        self.table_view.hideColumn(6)

    def __action_hide_distances(self):
        self.table_view.hideColumn(7)

    def __action_show_residue_pair_no(self):
        self.table_view.showColumn(0)

    def __action_show_protein_1_chain(self):
        self.table_view.showColumn(1)

    def __action_show_protein_1_position(self):
        self.table_view.showColumn(2)

    def __action_show_protein_1_residue(self):
        self.table_view.showColumn(3)

    def __action_show_protein_2_chain(self):
        self.table_view.showColumn(4)

    def __action_show_protein_2_position(self):
        self.table_view.showColumn(5)

    def __action_show_protein_2_residue(self):
        self.table_view.showColumn(6)

    def __action_show_distances(self):
        self.table_view.showColumn(7)

    def __action_show_all_columns(self):
        self.table_view.showColumn(0)
        self.table_view.showColumn(1)
        self.table_view.showColumn(2)
        self.table_view.showColumn(3)
        self.table_view.showColumn(4)
        self.table_view.showColumn(5)
        self.table_view.showColumn(6)
        self.table_view.showColumn(7)
    # </editor-fold>

    def __slot_open_context_menu_for_histogram(self):
        pass

    def setup_context_menu(self) -> None:
        self.context_menu = QtWidgets.QMenu()
        self.hide_selected_column = QtWidgets.QAction(self.tr("Hide Column"))
        self.context_menu.addAction(self.hide_selected_column)
        self.hide_selected_column.triggered.connect(self._hide_selected_column)

        # Set the context menu for the buttons
        self.table_view.setContextMenuPolicy(3)
        self.table_view.customContextMenuRequested.connect(self._show_context_menu_for_seq_list)

        # <editor-fold desc="Context menu setup for histogram">
        # for matplotlib histogram
        self.context_menu_hist = QtWidgets.QMenu()
        self.properties_hist = QtWidgets.QAction(self.tr("Properties"))
        self.context_menu_hist.addAction(self.properties_hist)
        self.properties_hist.triggered.connect(self._open_properties_view)

        self.properties_hist_restore = QtWidgets.QAction(self.tr("Restore Defaults"))
        self.context_menu_hist.addAction(self.properties_hist_restore)
        self.properties_hist_restore.triggered.connect(self._restore_default_histogram_properties)

        # Set the context menu for the buttons
        self.plot_widget_dhistogram.setContextMenuPolicy(3)
        self.plot_widget_dhistogram.customContextMenuRequested.connect(self._show_context_menu_for_histogram)
        # </editor-fold>

    def _show_context_menu_for_seq_list(self, a_point):
        self.hide_selected_column.triggered.disconnect()
        self.hide_selected_column.triggered.connect(self._hide_selected_column)
        self.context_menu.exec_(self.table_view.mapToGlobal(a_point))

    def _show_context_menu_for_histogram(self, a_point):
        self.properties_hist.triggered.disconnect()
        self.properties_hist.triggered.connect(self._open_properties_view)
        self.properties_hist_restore.triggered.disconnect()
        self.properties_hist_restore.triggered.connect(self._restore_default_histogram_properties)
        # add here more action connections
        self.context_menu_hist.exec_(self.plot_widget_dhistogram.mapToGlobal(a_point))

    def _open_properties_view(self):
        # self.properties_view =
        self._histogram_properties_view = histogram_properties_view.HistogramPropertiesView(self._histogram_properties)
        self._histogram_properties_view.new_properties.connect(self.post_open_properties_view)
        self._histogram_properties_view.show()

    def post_open_properties_view(self, new_properties: tuple):
        tmp_current_x_axis_units = self._histogram_properties[enums.HistogramPropertiesEnum.X_AXIS_UNITS]
        tmp_current_distance_interval = self._histogram_properties[enums.HistogramPropertiesEnum.DISTANCE_INTERVAL]
        tmp_x_axis_units = int(new_properties[0])
        tmp_distance_interval = float(new_properties[1])
        if tmp_current_x_axis_units != tmp_x_axis_units or tmp_current_distance_interval != tmp_distance_interval:
            # at least one property changed
            self._histogram_properties[enums.HistogramPropertiesEnum.X_AXIS_UNITS] = tmp_x_axis_units
            self._histogram_properties[enums.HistogramPropertiesEnum.DISTANCE_INTERVAL] = tmp_distance_interval

            self.actual_resize()

    def _restore_default_histogram_properties(self):
        self._histogram_properties[enums.HistogramPropertiesEnum.X_AXIS_UNITS] = constants.DEFAULT_HISTOGRAM_PROPERTIES[enums.HistogramPropertiesEnum.X_AXIS_UNITS]
        self._histogram_properties[enums.HistogramPropertiesEnum.DISTANCE_INTERVAL] = constants.DEFAULT_HISTOGRAM_PROPERTIES[enums.HistogramPropertiesEnum.DISTANCE_INTERVAL]
        self.actual_resize()

    def _hide_selected_column(self):
        self.table_view.hideColumn(self.table_view.currentIndex().column())
