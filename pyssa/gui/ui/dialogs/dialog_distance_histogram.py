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
"""Module for the distance histogram dialog."""
import copy
import numpy as np
from PyQt5 import QtGui
from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt import NavigationToolbar2QT

from pyssa.internal.data_structures import protein_pair
import pyssa.gui.ui.styles.styles as custom_pyssa_styles
from pyssa.gui.ui.forms.auto_generated.auto_dialog_distance_histogram import Ui_Dialog
from pyssa.util import pyssa_keys
from pyssa.util import constants
from pyssa.util import gui_utils
from PyQt5 import QtCore
from PyQt5.QtWidgets import QVBoxLayout, QWidget, QScrollArea
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class PlotWidget(QWidget):
    """Class for a custom QWidget for plotting."""

    def __init__(self, figure_size: tuple[float, float], parent=None) -> None:  # noqa: ANN001
        """Constructor.

        Args:
            figure_size: the size of the figure.
            parent: the parent
        """
        super(PlotWidget, self).__init__(parent)
        self.figure = Figure(figsize=(figure_size[0], figure_size[1]))
        self.canvas = FigureCanvas(self.figure)
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)


class DialogDistanceHistogram(QtWidgets.QDialog):
    """Class for a distance histogram dialog."""

    def __init__(self, protein_pair_from_project: "protein_pair.ProteinPair", parent=None) -> None:  # noqa: ANN001
        """Constructor.

        Args:
            protein_pair_from_project: the protein pair which should be used for the distance histogram.
            parent: the parent.
        """
        QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        self.protein_pair_for_analysis: protein_pair.ProteinPair = protein_pair_from_project
        self.histogram_pos = (0, 2.5)
        self.scroll_pos = ""
        custom_pyssa_styles.set_stylesheet(self)
        self.ui.btn_scroll_up.hide()
        self.ui.btn_scroll_down.hide()
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)

        # <editor-fold desc="Matplotlib">
        self.scroll_area = QScrollArea()
        self.plot_widget = PlotWidget((19, 20))
        self.toolbar = NavigationToolbar2QT(self.plot_widget.canvas, self)
        # Find the action you want to remove
        items_to_remove = ["Home", "Back", "Forward", "Pan", "Zoom", "Subplots", "Customize"]
        for tmp_item in items_to_remove:
            for action in self.toolbar.actions():
                print(action.text())
                if action.text() == tmp_item:
                    self.toolbar.removeAction(action)
        self.ui.main_Layout.addWidget(self.toolbar)
        self.ui.main_Layout.addWidget(self.scroll_area)
        self.scroll_area.setWidget(self.plot_widget)
        self.plot_histogram()

        # </editor-fold>

        # styles
        stylesheet = """
        QDialog {background-color: #F6F4F8;}
        QTableWidget {background-color: white;}
        """
        self.setStyleSheet(stylesheet)

        items = [
            "Very large",
            "Large",
            "Medium",
            "Small",
        ]
        gui_utils.fill_combo_box(self.ui.cb_bar_size, items)
        self.ui.cb_bar_size.setCurrentIndex(2)
        self.ui.btn_scroll_up.clicked.connect(self.scroll_up)
        self.ui.btn_scroll_down.clicked.connect(self.scroll_down)
        self.ui.cb_bar_size.currentTextChanged.connect(self.update_bar_width)
        self.setWindowFlag(QtCore.Qt.WindowMaximizeButtonHint, True)
        self.setWindowFlag(QtCore.Qt.WindowCloseButtonHint, True)
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        self.setWindowTitle("Distance Histogram")

    def update_bar_width(self) -> None:
        """Updates the width of the bars."""
        if self.ui.cb_bar_size.currentText() == "Very large":
            self.setup_histogram_view((19, 50))
        elif self.ui.cb_bar_size.currentText() == "Large":
            self.setup_histogram_view((19, 30))
        elif self.ui.cb_bar_size.currentText() == "Medium":
            self.setup_histogram_view((19, 20))
        elif self.ui.cb_bar_size.currentText() == "Small":
            self.setup_histogram_view((19, 10))

        self.plot_histogram()

    def setup_histogram_view(self, figure_size: tuple[int, int]) -> None:
        """Sets up the histogram view."""
        self.ui.main_Layout.removeWidget(self.toolbar)
        self.ui.main_Layout.removeWidget(self.scroll_area)
        self.scroll_area = QScrollArea()
        self.plot_widget = PlotWidget(figure_size)
        self.toolbar = NavigationToolbar2QT(self.plot_widget.canvas, self)
        # Find the action you want to remove
        items_to_remove = ["Home", "Back", "Forward", "Pan", "Zoom", "Subplots", "Customize"]
        for tmp_item in items_to_remove:
            for action in self.toolbar.actions():
                print(action.text())
                if action.text() == tmp_item:
                    self.toolbar.removeAction(action)
        self.ui.main_Layout.addWidget(self.toolbar)
        self.ui.main_Layout.addWidget(self.scroll_area)
        self.scroll_area.setWidget(self.plot_widget)
        self.plot_histogram()

    def plot_histogram(self) -> None:
        """Plots the histogram."""
        # Clear any existing plot
        self.plot_widget.figure.clear()

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

        # Move the entire x-axis to the top
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position("top")

        # Calculate the midpoints between bin edges
        bin_midpoints = (bins[:-1] + bins[1:]) / 2
        # Set y-ticks at the bin midpoints
        ax.set_yticks(bin_midpoints)
        # Set custom tick labels
        custom_labels = [f"{bin_start} - {bin_end}" for bin_start, bin_end in zip(bins[:-1], bins[1:])]
        ax.set_yticklabels(custom_labels)
        # Set axis labels
        ax.set_xlabel("Count")
        ax.set_ylabel("Bins")
        # Invert the y-axis
        ax.invert_yaxis()
        # Remove the spines where no axis are present
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)

        # Remove the spines where no axis are present
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_visible(False)

        # Set axis labels
        ax.set_xlabel("Frequency of C-α distances")
        ax.set_ylabel("Distance in Å")

        # Adjust subplot parameters to reduce white space
        self.plot_widget.figure.tight_layout()

        # Refresh the canvas
        self.plot_widget.canvas.draw()

    def scroll_down(self) -> None:
        """Scroll down the histogram."""
        if self.scroll_pos == "up":
            self.histogram_pos = (self.histogram_pos[0] + 0.3, self.histogram_pos[1] + 0.3)
        self.view_box.setRange(yRange=[self.histogram_pos[0], self.histogram_pos[1]])
        self.histogram_pos = (self.histogram_pos[0] + 0.3, self.histogram_pos[1] + 0.3)
        self.scroll_pos = "down"

    def scroll_up(self) -> None:
        """Scroll up the histogram."""
        if self.scroll_pos == "down":
            self.histogram_pos = (self.histogram_pos[0] - 0.6, self.histogram_pos[1] - 0.6)
        else:
            self.histogram_pos = (self.histogram_pos[0] - 0.3, self.histogram_pos[1] - 0.3)
        self.view_box.setRange(yRange=[self.histogram_pos[0], self.histogram_pos[1]])
        self.scroll_pos = "up"
