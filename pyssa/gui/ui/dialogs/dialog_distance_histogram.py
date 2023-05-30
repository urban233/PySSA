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
import copy
import numpy as np
from PyQt5.QtGui import QIcon
from PyQt5 import QtWidgets
from PyQt5 import QtGui
import pyqtgraph as pg
from pyssa.internal.data_structures import protein_pair
from pyssa.gui.ui.forms.auto_generated.auto_dialog_distance_histogram import Ui_Dialog
from pyssa.util import pyssa_keys
from pyssa.util import constants
from PyQt5 import QtCore


class DialogDistanceHistogram(QtWidgets.QDialog):

    def __init__(self, protein_pair_from_project, parent=None):
        """Constructor

        Args:
            args
            kwargs
        """
        QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        self.protein_pair_for_analysis: protein_pair.ProteinPair = protein_pair_from_project
        self.histogram_pos = (0, 2.5)
        self.scroll_pos = ""

        self.graph_widget = pg.PlotWidget()
        self.view_box = self.graph_widget.plotItem.getViewBox()
        # read csv file
        # TODO: Needs correction for xml model
        # path = pathlib.Path(f"{self.protein_pair_for_analysis.results_dir}/distance_csv/distances.csv")
        # distance_list = []
        # with open(path, 'r', encoding="utf-8") as csv_file:
        #     i = 0
        #     for line in csv_file:
        #         cleaned_line = line.replace("\n", "")
        #         if cleaned_line.split(",")[8] != 'distance':
        #             distance_list.append(float(cleaned_line.split(",")[8]))
        distance_data: dict[str, np.ndarray] = self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data
        distance_list = copy.deepcopy(distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES])
        distance_list.sort()
        length = len(distance_list)
        max_distance = distance_list[length - 1]
        x, y = np.histogram(distance_list, bins=np.arange(0, max_distance, 0.25))
        if x.size != y.size:
            x = np.resize(x, (1, y.size))
        # this conversion is needed for the pyqtgraph library!
        x = x.tolist()
        try:
            x = x[0]
        except IndexError:
            # Error got raised where the distances where all 0
            # tools.quick_log_and_display("error", "The histogram could not be created.",
            #                             self.status_bar, "The histogram could not be created. "
            #                                              " Check the distance table!")
            return

        y = y.tolist()
        color = QtGui.QColor.fromRgb(255, 128, 128)
        # creates bar chart item
        graph_bar_item = pg.BarGraphItem(x0=0, y=y, height=0.2, width=x,
                                         pen=pg.mkPen(color="#4B91F7"), brush=pg.mkBrush(color="#4B91F7"))
        # creates y-labels for bar chart
        y_labels = []
        for i in range(len(y)):
            try:
                label = f"{y[i]} - {y[i + 1]}"
            except IndexError:
                # detects if a last label is necessary
                label = f"{y[i]} - {y[i] + 0.25}"
            y_labels.append(label)
        y_values = y
        ticks = []
        for i, item in enumerate(y_labels):
            ticks.append((y_values[i], item))
        ticks = [ticks]

        # styling the plot
        self.graph_widget.setBackground('w')
        self.graph_widget.setTitle(f"Distance Histogram of {self.protein_pair_for_analysis.name}", size="23pt")
        styles = {'font-size': '14px'}
        ax_label_x = "Distance in Å"
        self.graph_widget.setLabel('left', ax_label_x, **styles)
        self.graph_widget.setLabel('bottom', "Frequency of C-α distances", **styles)
        self.graph_widget.addItem(graph_bar_item)
        bar_ax = self.graph_widget.getAxis('left')
        bar_ax.setTicks(ticks)
        self.view_box.invertY(True)
        # text = pg.TextItem(text='', border='w', fill=(0, 0, 255, 100), anchor=(0, 0))
        # self.graph_widget.addItem(text)
        self.ui.main_Layout.addWidget(self.graph_widget)

        # styles
        stylesheet = """
        QDialog {background-color: #F6F4F8;}
        QTableWidget {background-color: white;}
        """
        self.setStyleSheet(stylesheet)

        self.ui.btn_scroll_up.clicked.connect(self.scroll_up)
        self.ui.btn_scroll_down.clicked.connect(self.scroll_down)
        self.setWindowFlag(QtCore.Qt.WindowMaximizeButtonHint, True)
        self.setWindowFlag(QtCore.Qt.WindowCloseButtonHint, True)
        self.setWindowIcon(QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
        self.setWindowTitle("Distance Histogram")

    def scroll_down(self):
        if self.scroll_pos == "up":
            self.histogram_pos = (self.histogram_pos[0] + 0.3, self.histogram_pos[1] + 0.3)
        self.view_box.setRange(yRange=[self.histogram_pos[0], self.histogram_pos[1]])
        self.histogram_pos = (self.histogram_pos[0] + 0.3, self.histogram_pos[1] + 0.3)
        self.scroll_pos = "down"

    def scroll_up(self):
        if self.scroll_pos == "down":
            self.histogram_pos = (self.histogram_pos[0] - 0.6, self.histogram_pos[1] - 0.6)
        else:
            self.histogram_pos = (self.histogram_pos[0] - 0.3, self.histogram_pos[1] - 0.3)
        self.view_box.setRange(yRange=[self.histogram_pos[0], self.histogram_pos[1]])
        self.scroll_pos = "up"
