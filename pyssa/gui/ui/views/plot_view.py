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
from pyssa.internal.data_structures import protein_pair
from pyssa.internal.thread import tasks
from pyssa.internal.thread.async_pyssa import util_async
from pyssa.util import pyssa_keys
from pyssa.util import constants
from PyQt5.QtWidgets import QVBoxLayout, QWidget, QScrollArea
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
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
        self.figure = Figure(figsize=(18, 7.25))
        self.canvas = FigureCanvas(self.figure)
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)


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

        # self.resizeEvent = self.handle_resize
        # Create a timer for delayed updates
        self.resize_timer = QtCore.QTimer(self)
        #self.resize_timer.timeout.connect(self.handle_resize_timeout)
        self.protein_pair_for_analysis: protein_pair.ProteinPair = protein_pair_from_project
        #custom_pyssa_styles.set_stylesheet(self)

        self._initialize_ui()

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
        self._connect_all_signals()

    def _initialize_ui(self):
        # <editor-fold desc="Define basic ui elements">
        self.scroll_area = QScrollArea()
        self.plot_widget = PlotWidget()
        self.menubar = QtWidgets.QMenuBar(self)
        self.lbl_status = QtWidgets.QLabel()
        self.lbl_status.setText("Status: ")

        # </editor-fold>

        # <editor-fold desc="Create a menu and add it to the menu bar">
        plots_menu = QtWidgets.QMenu('Plots', self)
        help_menu = QtWidgets.QMenu('Help', self)
        self.menubar.addMenu(plots_menu)
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
        self.action_show_all = QtWidgets.QAction('Show all', self)

        plots_menu.addAction(self.action_plot)
        plots_menu.addAction(self.action_histogram)
        plots_menu.addAction(self.action_table)
        plots_menu.addAction(self.action_show_all)

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
        # self.main_Layout = QtWidgets.QVBoxLayout()
        # self.main_Layout.addWidget(self.scroll_area)
        # self.main_Layout.addWidget(self.lbl_status)
        # self.scroll_area_layout = QtWidgets.QHBoxLayout()
        # self.scroll_area_layout.addWidget(self.plot_widget)
        # self.table_view = QtWidgets.QTableView()
        # self.table_view.setMinimumWidth(450)
        # self.scroll_area_layout.addWidget(self.table_view)
        # self.scroll_area.setLayout(self.scroll_area_layout)
        # self.main_Layout.setMenuBar(self.menubar)
        # self.setLayout(self.main_Layout)

        # Create labels
        self.lbl_status1 = QtWidgets.QLabel(f"Protein 1: {self.protein_pair_for_analysis.protein_1.get_molecule_object()}")
        self.lbl_status2 = QtWidgets.QLabel(f"Protein 2: {self.protein_pair_for_analysis.protein_2.get_molecule_object()}")

        # Create a QTableView
        self.table_view = QtWidgets.QTableView()
        self.table_view.setMinimumWidth(450)

        # Create a QWidget to hold labels and QTableView
        container_widget = QtWidgets.QWidget()
        container_layout = QtWidgets.QVBoxLayout(container_widget)
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
        container_widget.setStyleSheet(stylesheet)
        # Create a QHBoxLayout for the scroll area
        self.scroll_area_layout = QtWidgets.QHBoxLayout()

        # Assuming self.plot_widget and self.scroll_area are already defined
        self.scroll_area_layout.addWidget(self.plot_widget)
        self.scroll_area_layout.addWidget(container_widget)

        # Set up scroll area
        self.scroll_area = QtWidgets.QScrollArea()
        self.scroll_area.setLayout(self.scroll_area_layout)

        # Create main layout
        self.main_Layout = QtWidgets.QVBoxLayout()
        self.main_Layout.addWidget(self.scroll_area)
        self.main_Layout.addWidget(self.lbl_status)
        self.main_Layout.setMenuBar(self.menubar)
        self.setLayout(self.main_Layout)

        # </editor-fold>

        # <editor-fold desc="Setup subplots">
        self._ax_plot = self.plot_widget.figure.add_subplot(211)
        self._ax_hist = self.plot_widget.figure.add_subplot(212)

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

    def _connect_all_signals(self):
        self.action_table.triggered.connect(self.hide_distance_table)
        self.table_view.clicked.connect(self.highlight_table_selection_in_plot)
        self.action_plot.triggered.connect(self.toggle_graphics_visablity)
        self.action_histogram.triggered.connect(self.toggle_graphics_visablity)
        self.action_show_all.triggered.connect(self.show_all)

        self.plot_widget.canvas.mpl_connect('button_press_event', self.on_canvas_click)

    def hide_distance_table(self):
        if self.action_table.isChecked():
            self.table_view.show()
        else:
            self.table_view.hide()

    def toggle_graphics_visablity(self):
        self.plot_widget.figure.clear()
        self.plot_widget.show()
        if self.action_plot.isChecked() and self.action_histogram.isChecked():
            self._ax_plot, self._ax_hist = self.plot_widget.figure.subplots(2)
            self.create_distance_plot()
            self.setup_plot_defaults()
            self.create_distance_histogram()
            self.setup_histogram_defaults()
        elif self.action_plot.isChecked() and not self.action_histogram.isChecked():
            self._ax_plot = self.plot_widget.figure.subplots(1)
            self.create_distance_plot()
            self.setup_plot_defaults()
        elif not self.action_plot.isChecked() and self.action_histogram.isChecked():
            self._ax_hist = self.plot_widget.figure.subplots(1)
            self.create_distance_histogram()
            self.setup_histogram_defaults()
        else:
            self.plot_widget.hide()
        self.plot_widget.figure.tight_layout()
        self.plot_widget.canvas.draw()

    def show_all(self):
        self.plot_widget.figure.clear()
        self.plot_widget.show()
        self.table_view.show()
        self.create_all_graphics()

    def create_all_graphics(self):
        self.create_distance_plot()
        self.setup_plot_defaults()
        self.create_distance_histogram()
        self.setup_histogram_defaults()
        self.plot_widget.figure.tight_layout()
        self.plot_widget.canvas.draw()

    def setup_plot_defaults(self):
        # Set axis labels
        self._ax_plot.set_xlabel("Residue pair no.")
        self._ax_plot.set_ylabel("Distance in Å")
        self._ax_plot.set_title("Distance plot")
        # Set the x-axis limits with minimum value set to 0
        self._ax_plot.set_xlim(0)
        # Set the number of ticks for the x and y axes
        self._ax_plot.xaxis.set_major_locator(ticker.MultipleLocator(5))

    def setup_histogram_defaults(self):
        # Move the entire x-axis to the top
        self._ax_hist.xaxis.tick_top()
        self._ax_hist.xaxis.set_label_position("top")
        self._ax_hist.set_title("Distance histogram")
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
        self._ax_hist.set_xlabel("Frequency of C-α distances")
        self._ax_hist.set_ylabel("Distance in Å")

    def create_distance_plot(self):
        # data for actual distance plot line
        distance_data = self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data
        distance_list = distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES]
        # data for cutoff line
        # cutoff_line = []
        # for i in range(len(distance_list)):
        #     cutoff_line.append(self.protein_pair_for_analysis.distance_analysis.cutoff)
        self._ax_plot.plot(distance_list)

    def create_distance_histogram(self):
        distance_data: dict[
            str,
            np.ndarray,
        ] = self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data
        distance_list = copy.deepcopy(distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES])
        distance_list.sort()
        length = len(distance_list)
        max_distance = distance_list[length - 1]

        n, bins, patches = self._ax_hist.hist(
            distance_list,
            bins=np.arange(0, max_distance + 0.25, 0.25),
            orientation="horizontal",
            rwidth=0.7,
            color="#4B91F7",
        )

        # Add labels to the non-zero frequency histogram bars
        for bin_value, patch in zip(n, patches):
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
            bin_midpoints = (bins[:-1] + bins[1:]) / 2
            # Set y-ticks at the bin midpoints
            self._ax_hist.set_yticks(bin_midpoints)
            # Set custom tick labels
            custom_labels = [f"{bin_start} - {bin_end}" for bin_start, bin_end in zip(bins[:-1], bins[1:])]
            self._ax_hist.set_yticklabels(custom_labels)

    def on_canvas_click(self, event):
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
        self.plot_widget.canvas.draw()
        tmp_row_information = self.search_point_by_id(x_nearest)
        if tmp_row_information is not None:
            tmp_prot_1_name = self.protein_pair_for_analysis.protein_1.get_molecule_object()
            tmp_prot_2_name = self.protein_pair_for_analysis.protein_2.get_molecule_object()
            msg: str = (f"Status: Residue pair no. = {x_nearest}, Distance (Å) = {y_nearest} Between {tmp_prot_1_name} "
                        f"Chain {tmp_row_information[1]} Position {tmp_row_information[2]} Residue {tmp_row_information[3]} and "
                        f"{tmp_prot_2_name} Chain {tmp_row_information[4]} Position {tmp_row_information[5]} Residue {tmp_row_information[6]}")
            self.lbl_status.setText(msg)
        else:
            self.lbl_status.setText(f"Status: Residue pair no. = {x_nearest}, Distance (Å) = {y_nearest}")
        # Change the selection color
        self.change_selection_color(QtGui.QColor(75, 145, 247, 200))

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
                tmp_row_information = (
                    self.table_view.model().index(row, 0).data(Qt.DisplayRole),
                    self.table_view.model().index(row, 1).data(Qt.DisplayRole),
                    self.table_view.model().index(row, 2).data(Qt.DisplayRole),
                    self.table_view.model().index(row, 3).data(Qt.DisplayRole),
                    self.table_view.model().index(row, 4).data(Qt.DisplayRole),
                    self.table_view.model().index(row, 5).data(Qt.DisplayRole),
                    self.table_view.model().index(row, 6).data(Qt.DisplayRole),
                )
                return tmp_row_information  # Return True if the id is found
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
        self.plot_widget.canvas.draw()

        tmp_row = self.table_view.currentIndex().row()
        self.table_view.selectRow(tmp_row)
        tmp_row_information = (
            self.table_view.model().index(tmp_row, 0).data(Qt.DisplayRole),
            self.table_view.model().index(tmp_row, 1).data(Qt.DisplayRole),
            self.table_view.model().index(tmp_row, 2).data(Qt.DisplayRole),
            self.table_view.model().index(tmp_row, 3).data(Qt.DisplayRole),
            self.table_view.model().index(tmp_row, 4).data(Qt.DisplayRole),
            self.table_view.model().index(tmp_row, 5).data(Qt.DisplayRole),
            self.table_view.model().index(tmp_row, 6).data(Qt.DisplayRole),
        )
        if tmp_row_information is not None:
            tmp_prot_1_name = self.protein_pair_for_analysis.protein_1.get_molecule_object()
            tmp_prot_2_name = self.protein_pair_for_analysis.protein_2.get_molecule_object()
            msg: str = (f"Status: Residue pair no. = {tmp_x_value}, Distance (Å) = {tmp_y_value} Between {tmp_prot_1_name} "
                        f"Chain {tmp_row_information[1]} Position {tmp_row_information[2]} Residue {tmp_row_information[3]} and "
                        f"{tmp_prot_2_name} Chain {tmp_row_information[4]} Position {tmp_row_information[5]} Residue {tmp_row_information[6]}")
            self.lbl_status.setText(msg)
        else:
            self.lbl_status.setText(f"Status: Residue pair no. = {tmp_x_value}, Distance (Å) = {tmp_y_value}")

    # <editor-fold desc="Experimental!">
    # def resizeEvent(self, event) -> None:  # noqa: N802, ANN001
    #     """The actual resize event overwritten from PyQt5."""
    #     # Let the base class handle the event
    #     super().resizeEvent(event)
    #
    #     # Start or restart the timer when resizing
    #     self.resize_timer.start(200)  # Adjust the timeout as needed

    def handle_resize_timeout(self) -> None:
        """Handles the resize process."""
        # Handle the resize event after the timeout
        size = self.scroll_area.size()
        self.plot_widget.setFixedSize(size.width() - 10, size.height() - 10)
        self.plot_widget.figure.tight_layout()
        self.plot_widget.canvas.draw()
    # </editor-fold>
