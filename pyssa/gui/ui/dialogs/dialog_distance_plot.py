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
from PyQt5.QtCore import Qt
from PyQt5 import QtGui
from PyQt5 import QtWidgets
from matplotlib.backends.backend_qt import NavigationToolbar2QT
from pymol import cmd
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
    def __init__(self, parent=None) -> None:  # noqa: ANN001
        super(PlotWidget, self).__init__(parent)
        self.figure = Figure(figsize=(18, 7.25))
        self.canvas = FigureCanvas(self.figure)
        layout = QVBoxLayout()
        layout.addWidget(self.canvas)
        self.setLayout(layout)


class DialogDistancePlot(QtWidgets.QDialog):
    def __init__(self, protein_pair_from_project, parent=None):
        """Constructor.

        Args:
            args
            kwargs
        """
        QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        # self.graph_widget = pg.PlotWidget()
        self.protein_pair_for_analysis: protein_pair.ProteinPair = protein_pair_from_project
        custom_pyssa_styles.set_stylesheet(self)
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
        # <editor-fold desc="PyQtGraph">
        # creates actual distance plot line
        # distance_data = self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data
        # distance_list = distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES]
        # self.graph_widget.plotItem.plot(distance_list, pen=pg.mkPen(color="#4B91F7", width=6),
        #                            symbol="o", symbolSize=10, symbolBrush=('b'))
        # self.view_box = self.graph_widget.plotItem.getViewBox()
        # # creates cutoff line
        # cutoff_line = []
        # for i in range(len(distance_list)):
        #     cutoff_line.append(self.protein_pair_for_analysis.distance_analysis.cutoff)
        # self.graph_widget.plotItem.plot(cutoff_line, pen=pg.mkPen(color="#f83021", width=6))
        # # styling the plot
        # self.graph_widget.setBackground('w')
        # self.graph_widget.setTitle(f"Distance Plot of {self.protein_pair_for_analysis.name}", size="23pt")
        # styles = {'font-size': '14px'}
        # ax_label_y = "Distance in Å"
        # self.graph_widget.setLabel('left', ax_label_y, **styles)
        # self.graph_widget.setLabel('bottom', "Residue pair no.", **styles)
        # self.graph_widget.plotItem.showGrid(x=True, y=True)
        # self.ui.main_Layout.addWidget(self.graph_widget)
        # </editor-fold>

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

        self.ui.main_Layout.addWidget(self.toolbar)
        self.ui.main_Layout.addWidget(self.scroll_area)
        self.ax = self.plot_widget.figure.add_subplot(111)
        self.scroll_area.setWidget(self.plot_widget)
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
        self.setWindowFlag(QtCore.Qt.WindowMaximizeButtonHint, True)
        self.setWindowFlag(QtCore.Qt.WindowCloseButtonHint, True)
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        self.setWindowTitle("Distance Plot")

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_Return or event.key() == Qt.Key_Enter:
            self.update_plot()

    def plot_distance_data(self):
        # Clear any existing plot
        self.plot_widget.figure.clear()
        # Create an axis and plot some data
        ax = self.plot_widget.figure.add_subplot(111)
        # data for actual distance plot line
        distance_data = self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data
        distance_list = distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES]
        # data for cutoff line
        cutoff_line = []
        for i in range(len(distance_list)):
            cutoff_line.append(self.protein_pair_for_analysis.distance_analysis.cutoff)
        ax.plot(distance_list)
        ax.plot(cutoff_line)
        # Set axis labels
        ax.set_xlabel("Residue pair no.")
        ax.set_ylabel("Distance in Å")
        # Adjust subplot parameters to reduce white space
        self.plot_widget.figure.tight_layout()
        # Set the x-axis limits with minimum value set to 0
        ax.set_xlim(0)
        # Set the number of ticks for the x and y axes
        ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
        # Refresh the canvas
        self.plot_widget.canvas.draw()

    def update_plot(self):
        self.update_distance_plot(self.ui.cb_turn_on_grid.checkState())

    def update_distance_plot(self, grid: bool):
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
        distance_data = self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data
        distance_list = distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES]

        # data for cutoff line
        cutoff_line = []
        for i in range(len(distance_list)):
            cutoff_line.append(self.protein_pair_for_analysis.distance_analysis.cutoff)

        ax.plot(distance_list)
        ax.plot(cutoff_line)

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

    def highlight_selection_in_pymol(self, from_aa: int, to_aa: int):
        print("Sync is active.")
        zoom_selection = f"/{self.protein_pair_for_analysis.protein_1.get_molecule_object()}///{from_aa}-{to_aa}/CA"
        cmd.select("zoom_sele", zoom_selection)
        cmd.show("spheres", "zoom_sele")
        cmd.alter("zoom_sele", "vdw=0.7")
        cmd.rebuild()
        cmd.zoom("zoom_sele")

    def hide_highlight_selection_in_pymol(self, from_aa: int, to_aa: int):
        zoom_selection = f"/{self.protein_pair_for_analysis.protein_1.get_molecule_object()}///{from_aa}-{to_aa}/CA"
        cmd.select("zoom_sele", zoom_selection)
        cmd.hide("spheres", "zoom_sele")

    def reset_distance_plot(self):
        from_aa = int(self.ui.sp_distance_plot_from.text())
        to_aa = int(self.ui.sp_distance_plot_to.text())
        self.ui.cb_turn_on_grid.setChecked(False)
        self.plot_distance_data()
        self.hide_highlight_selection_in_pymol(from_aa, to_aa)

    def turn_gird_on_off(self):
        if self.ui.cb_turn_on_grid.isChecked():
            self.update_distance_plot(True)
        else:
            self.update_distance_plot(False)
