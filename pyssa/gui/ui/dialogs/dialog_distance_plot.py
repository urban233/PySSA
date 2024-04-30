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
import pyssa.gui.ui.styles.styles as custom_pyssa_styles
from PyQt5 import QtCore
from pyssa.internal.data_structures import protein_pair
from pyssa.util import pyssa_keys
from pyssa.gui.ui.forms.auto_generated.auto_dialog_distance_plot import Ui_Dialog
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


class DialogDistancePlot(QtWidgets.QDialog):
    """Class for a distance plot dialog."""

    def __init__(self, protein_pair_from_project: "protein_pair.ProteinPair", a_project, a_results_name, parent=None) -> None:  # noqa: ANN001
        """Constructor.

        Args:
            protein_pair_from_project: the protein pair which should be used for the distance histogram.
            parent: the parent.
        """
        QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)

        # self.resizeEvent = self.handle_resize
        # Create a timer for delayed updates
        self.resize_timer = QtCore.QTimer(self)
        self.resize_timer.timeout.connect(self.handle_resize_timeout)

        # self.graph_widget = pg.PlotWidget()
        self.protein_pair_for_analysis: protein_pair.ProteinPair = protein_pair_from_project
        custom_pyssa_styles.set_stylesheet(self)
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)

        self.scroll_area = QScrollArea()
        self.plot_widget = PlotWidget()
        self.toolbar = NavigationToolbar2QT(self.plot_widget.canvas, self)
        # Find the action you want to remove
        items_to_remove = ["Home", "Back", "Forward", "Pan", "Zoom", "Subplots"]
        for tmp_item in items_to_remove:
            for action in self.toolbar.actions():
                print(action.text())
                if action.text() == tmp_item:
                    self.toolbar.removeAction(action)

        self.lbl_status = QtWidgets.QLabel()
        self.lbl_status.setText("Status: ")
        self.ui.main_Layout.addWidget(self.toolbar)
        self.ui.main_Layout.addWidget(self.scroll_area)
        self.ui.main_Layout.addWidget(self.lbl_status)
        #self.ax = self.plot_widget.figure.add_subplot(111)
        self.ax = self.plot_widget.figure.subplots(2)
        self._current_project = a_project
        self.results_name = a_results_name
        self.scroll_area_layout = QtWidgets.QHBoxLayout()
        self.scroll_area_layout.addWidget(self.plot_widget)
        self.table_view = QtWidgets.QTableView()
        self.table_view.setMinimumWidth(300)
        self.scroll_area_layout.addWidget(self.table_view)
        self.scroll_area.setLayout(self.scroll_area_layout)

        # --BEGIN

        # <editor-fold desc="Distance table logic">
        csv_model = QtGui.QStandardItemModel()
        csv_model.setColumnCount(7)
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
        csv_model.setHorizontalHeaderLabels(labels)
        self.table_view.setModel(csv_model)

        tmp_protein_pair = self._current_project.search_protein_pair(
            self.results_name,
        )
        csv_filepath = pathlib.Path(f"{constants.CACHE_CSV_DIR}/{tmp_protein_pair.name}.csv")
        if not os.path.exists(constants.CACHE_CSV_DIR):
            os.mkdir(constants.CACHE_CSV_DIR)
        tmp_protein_pair = self._current_project.search_protein_pair(self.results_name)

        distance_data = tmp_protein_pair.distance_analysis.analysis_results.distance_data
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
                csv_model.insertRow(i, standard_item_list)
            i += 1
        csv_file.close()
        csv_model.removeRow(0)
        self.table_view.setAlternatingRowColors(True)
        self.table_view.resizeColumnsToContents()
        self.table_view.verticalHeader().setVisible(False)
        self.table_view.setSortingEnabled(True)
        self.table_view.sortByColumn(0, QtCore.Qt.AscendingOrder)
        self.table_view.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)  # disables editing of cells
        # </editor-fold>

        # --END

        #self.scroll_area.setWidget(self.plot_widget)
        self.plot_distance_data()
        self.toolbar.show()
        # styles
        stylesheet = """
        QDialog {background-color: #F6F4F8;}
        QTableWidget {background-color: white;}
        """
        self.setStyleSheet(stylesheet)

        self.ui.btn_distance_plot_save.hide()
        self.ui.btn_distance_plot_update.clicked.connect(self.update_plot)
        # self.ui.btn_distance_plot_save.clicked.connect(self.save_plot_to_file)
        self.ui.btn_distance_plot_reset.clicked.connect(self.reset_distance_plot)
        self.ui.cb_turn_on_grid.stateChanged.connect(self.turn_gird_on_off)

        # Connect the mouse click event
        self.plot_widget.canvas.mpl_connect('button_press_event', self.on_canvas_click)

        self.setWindowFlag(QtCore.Qt.WindowMaximizeButtonHint, True)
        self.setWindowFlag(QtCore.Qt.WindowCloseButtonHint, True)
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        self.setWindowTitle("Distance Plot")

    def keyPressEvent(self, event) -> None:  # noqa: ANN001, N802
        """Event handler for a key press event."""
        if event.key() == Qt.Key_Return or event.key() == Qt.Key_Enter:
            self.update_plot()

    def on_canvas_click(self, event):
        x_clicked, y_clicked = event.xdata, event.ydata
        # Find the nearest point on the line
        distance_data = self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES]
        x_line = range(len(distance_data))
        y_line = distance_data
        distances = np.sqrt((x_line - x_clicked) ** 2 + (y_line - y_clicked) ** 2)
        index_nearest = np.argmin(distances)
        x_nearest, y_nearest = x_line[index_nearest], y_line[index_nearest]
        self.lbl_status.setText(f"Status: Residue pair no. = {x_nearest}, Distance (Å) = {y_nearest}")
        # Highlight the nearest point with a red dot using scatter
        # Clear any existing plot
        self.plot_widget.figure.clear()
        # Create an axis and plot some data
        ax = self.plot_widget.figure.subplots(2)
        # data for actual distance plot line
        distance_data = self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data
        distance_list = distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES]
        # data for cutoff line
        cutoff_line = []
        for i in range(len(distance_list)):
            cutoff_line.append(self.protein_pair_for_analysis.distance_analysis.cutoff)
        #ax1.plot(distance_list)
        #ax.plot(cutoff_line)
        #ax1.plot(x_nearest, y_nearest, 'ro', markersize=6)
        # Set axis labels
        ax.set_xlabel("Residue pair no.")
        ax.set_ylabel("Distance in Å")
        # Adjust subplot parameters to reduce white space
        self.plot_widget.figure.tight_layout()
        # Set the x-axis limits with minimum value set to 0
        ax.set_xlim(0)
        # Set the number of ticks for the x and y axes
        ax.xaxis.set_major_locator(ticker.MultipleLocator(5))

        distance_data: dict[
            str,
            np.ndarray,
        ] = self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data
        distance_list = copy.deepcopy(distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES])
        distance_list.sort()
        length = len(distance_list)
        max_distance = distance_list[length - 1]

        # Create an axis and plot a histogram
        ax = self.plot_widget.figure.add_subplot(111)
        n, bins, patches = ax.hist(
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
            ax.annotate(
                f"{bin_value}",
                xy=(x, y),
                xycoords="data",
                xytext=(3, 0),
                textcoords="offset points",
                ha="left",
                va="center",
            )

        # <editor-fold desc="Histogram logic">
        # Move the entire x-axis to the top
        ax[1].xaxis.tick_top()
        ax.xaxis.set_label_position("top")

        # Calculate the midpoints between bin edges
        bin_midpoints = (bins[:-1] + bins[1:]) / 2
        # Set y-ticks at the bin midpoints
        ax[1].set_yticks(bin_midpoints)
        # Set custom tick labels
        custom_labels = [f"{bin_start} - {bin_end}" for bin_start, bin_end in zip(bins[:-1], bins[1:])]
        ax[1].set_yticklabels(custom_labels)
        # Set axis labels
        ax[1].set_xlabel("Count")
        ax[1].set_ylabel("Bins")
        # Invert the y-axis
        ax[1].invert_yaxis()
        # Remove the spines where no ax[1]is are present
        ax[1].spines["right"].set_visible(False)
        ax[1].spines["bottom"].set_visible(False)

        # Remove the spines where no axis are present
        ax[1].spines["right"].set_visible(False)
        ax[1].spines["bottom"].set_visible(False)

        # Set axis labels
        ax[1].set_xlabel("Frequency of C-α distances")
        ax[1].set_ylabel("Distance in Å")
        # </editor-fold>

        # Refresh the canvas
        self.plot_widget.canvas.draw()

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

    def plot_distance_data(self) -> None:
        """Plots the distance data."""
        # Clear any existing plot
        self.plot_widget.figure.clear()
        # Create an axis and plot some data
        ax1, ax2 = self.plot_widget.figure.subplots(2)
        # data for actual distance plot line
        distance_data = self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data
        distance_list = distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES]
        # data for cutoff line
        cutoff_line = []
        for i in range(len(distance_list)):
            cutoff_line.append(self.protein_pair_for_analysis.distance_analysis.cutoff)
        ax1.plot(distance_list)
        # Set axis labels
        ax1.set_xlabel("Residue pair no.")
        ax1.set_ylabel("Distance in Å")
        ax1.set_title("Distance plot")
        # Adjust subplot parameters to reduce white space
        self.plot_widget.figure.tight_layout()
        # Set the x-axis limits with minimum value set to 0
        ax1.set_xlim(0)
        # Set the number of ticks for the x and y axes
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(5))

        distance_data: dict[
            str,
            np.ndarray,
        ] = self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data
        distance_list = copy.deepcopy(distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES])
        distance_list.sort()
        length = len(distance_list)
        max_distance = distance_list[length - 1]

        # Create an axis and plot a histogram
        #ax = self.plot_widget.figure.add_subplot(111)
        n, bins, patches = ax2.hist(
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
            ax2.annotate(
                f"{bin_value}",
                xy=(x, y),
                xycoords="data",
                xytext=(3, 0),
                textcoords="offset points",
                ha="left",
                va="center",
            )

        # <editor-fold desc="Histogram logic">
        # Move the entire x-axis to the top
        ax2.xaxis.tick_top()
        ax2.xaxis.set_label_position("top")
        ax2.set_title("Distance histogram")
        # Calculate the midpoints between bin edges
        bin_midpoints = (bins[:-1] + bins[1:]) / 2
        # Set y-ticks at the bin midpoints
        ax2.set_yticks(bin_midpoints)
        # Set custom tick labels
        custom_labels = [f"{bin_start} - {bin_end}" for bin_start, bin_end in zip(bins[:-1], bins[1:])]
        ax2.set_yticklabels(custom_labels)
        # Set axis labels
        ax2.set_xlabel("Count")
        ax2.set_ylabel("Bins")
        # Invert the y-axis
        ax2.invert_yaxis()
        # Remove the spines where no ax is are present
        ax2.spines["right"].set_visible(False)
        ax2.spines["bottom"].set_visible(False)

        # Remove the spines where no axis are present
        ax2.spines["right"].set_visible(False)
        ax2.spines["bottom"].set_visible(False)

        # Set axis labels
        ax2.set_xlabel("Frequency of C-α distances")
        ax2.set_ylabel("Distance in Å")
        # </editor-fold>

        self.plot_widget.figure.tight_layout()
        # Refresh the canvas
        self.plot_widget.canvas.draw()

        # # Clear any existing plot
        # self.plot_widget.figure.clear()
        # # Create an axis and plot some data
        # ax = self.plot_widget.figure.add_subplot(111)
        # # data for actual distance plot line
        # distance_data = self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data
        # distance_list = distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES]
        # # data for cutoff line
        # cutoff_line = []
        # for i in range(len(distance_list)):
        #     cutoff_line.append(self.protein_pair_for_analysis.distance_analysis.cutoff)
        # ax.plot(distance_list)
        # ax.plot(cutoff_line)
        # # Set axis labels
        # ax.set_xlabel("Residue pair no.")
        # ax.set_ylabel("Distance in Å")
        # # Adjust subplot parameters to reduce white space
        # self.plot_widget.figure.tight_layout()
        # # Set the x-axis limits with minimum value set to 0
        # ax.set_xlim(0)
        # # Set the number of ticks for the x and y axes
        # ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
        # # Refresh the canvas
        # self.plot_widget.canvas.draw()

    def update_plot(self) -> None:
        """Updates the distance plot."""
        self.update_distance_plot(self.ui.cb_turn_on_grid.checkState())

    def update_distance_plot(self, grid: bool) -> None:
        """This function updates the distance plot."""
        from_aa = int(self.ui.sp_distance_plot_from.text())
        to_aa = int(self.ui.sp_distance_plot_to.text())
        from_range = float(self.ui.dsp_distance_plot_from_range.text().replace(",", "."))
        to_range = float(self.ui.dsp_distance_plot_to_range.text().replace(",", "."))

        # Clear any existing plot
        self.plot_widget.figure.clear()
        # Create an axis and plot some data
        ax = self.plot_widget.figure.add_subplot(111)

        # data for actual distance plot line
        self.distance_data = self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data
        self.distance_list = self.distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES]

        # data for cutoff line
        self.cutoff_line = []
        for i in range(len(self.distance_list)):
            self.cutoff_line.append(self.protein_pair_for_analysis.distance_analysis.cutoff)

        ax.plot(self.distance_list)
        ax.plot(self.cutoff_line)

        # Set axis labels
        ax.set_xlabel("Residue pair no.")
        ax.set_ylabel("Distance in Å")
        # Adjust subplot parameters to reduce white space
        self.plot_widget.figure.tight_layout()

        # Set the x-axis limits with minimum value set to 0
        ax.set_xlim(0)
        # Set the number of ticks for the x and y axes
        x_locator = ticker.MultipleLocator(1)
        ax.xaxis.set_major_locator(x_locator)
        ax.xaxis.set_tick_params(which="major", pad=15)
        # Set the number of ticks for the y-axis and adjust spacing
        y_locator = ticker.MultipleLocator(1)
        ax.yaxis.set_major_locator(y_locator)
        ax.yaxis.set_tick_params(which="major", pad=15)
        # Set axis labels
        ax.set_xlabel("Residue pair no.")
        ax.set_ylabel("Distance in Å")
        # Zoom into the plot by setting the limits for x and y axes
        ax.set_xlim(from_aa, to_aa)  # Set the desired x-axis limits
        ax.set_ylim(from_range, to_range)  # Set the desired y-axis limits
        # Add grid to the plot
        ax.grid(grid)
        # Refresh the canvas
        self.plot_widget.canvas.draw()

        if self.ui.cb_sync_with_pymol.isChecked():
            self.highlight_selection_in_pymol(from_aa, to_aa)

    def highlight_selection_in_pymol(self, from_aa: int, to_aa: int) -> None:
        """Highlights the selection in PyMOL.

        Args:
            from_aa: the index of the first residue of the selection.
            to_aa: the index of the last residue of the selection.
        """
        pass
        # zoom_selection = f"/{self.protein_pair_for_analysis.protein_1.get_molecule_object()}///{from_aa}-{to_aa}/CA"
        # cmd.select("zoom_sele", zoom_selection)
        # cmd.show("spheres", "zoom_sele")
        # cmd.alter("zoom_sele", "vdw=0.7")
        # cmd.rebuild()
        # cmd.zoom("zoom_sele")

    def hide_highlight_selection_in_pymol(self, from_aa: int, to_aa: int) -> None:
        """Hides the selection in PyMOL.

        Args:
            from_aa: the index of the first residue of the selection.
            to_aa: the index of the last residue of the selection.
        """
        pass
        # zoom_selection = f"/{self.protein_pair_for_analysis.protein_1.get_molecule_object()}///{from_aa}-{to_aa}/CA"
        # cmd.select("zoom_sele", zoom_selection)
        # cmd.hide("spheres", "zoom_sele")

    def reset_distance_plot(self) -> None:
        """Resets the view of the distance plot."""
        from_aa = int(self.ui.sp_distance_plot_from.text())
        to_aa = int(self.ui.sp_distance_plot_to.text())
        self.ui.cb_turn_on_grid.setChecked(False)
        self.plot_distance_data()
        self.hide_highlight_selection_in_pymol(from_aa, to_aa)

    def turn_gird_on_off(self) -> None:
        """Turns gird on/off based on a checkbox state."""
        if self.ui.cb_turn_on_grid.isChecked():
            self.update_distance_plot(True)
        else:
            self.update_distance_plot(False)
