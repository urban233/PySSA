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

import os.path
import pathlib

from PyQt5.QtGui import QIcon
from pymol import Qt
from pymol import cmd
import pyqtgraph as pg
import pyqtgraph.exporters
from PyQt5 import QtCore
from internal.data_structures import protein_pair
from pyssa.util import pyssa_keys
from pyssa.gui.ui.forms.auto_generated.auto_dialog_distance_plot import Ui_Dialog
from util import constants


class DialogDistancePlot(Qt.QtWidgets.QDialog):

    def __init__(self, protein_pair_from_project, parent=None):
        """Constructor

        Args:
            args
            kwargs
        """
        Qt.QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        self.graph_widget = pg.PlotWidget()
        self.protein_pair_for_analysis: protein_pair.ProteinPair = protein_pair_from_project
        # read csv file
        # path = pathlib.Path(f"{self.protein_pair_for_analysis.results_dir}/distance_csv/distances.csv")
        # distance_list = []
        # cutoff_line = []
        # with open(path, 'r', encoding="utf-8") as csv_file:
        #     for line in csv_file:
        #         cleaned_line = line.replace("\n", "")
        #         if cleaned_line.split(",")[8] != 'distance':
        #             distance_list.append(float(cleaned_line.split(",")[8]))
        #             cutoff_line.append(float(self.protein_pair_for_analysis.cutoff))
        # creates actual distance plot line
        distance_data = self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data
        distance_list = distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES]
        self.graph_widget.plotItem.plot(distance_list, pen=pg.mkPen(color="#4B91F7", width=6),
                                   symbol="o", symbolSize=10, symbolBrush=('b'))
        self.view_box = self.graph_widget.plotItem.getViewBox()
        # creates cutoff line
        cutoff_line = []
        for i in range(len(distance_list)):
            cutoff_line.append(self.protein_pair_for_analysis.distance_analysis.cutoff)
        self.graph_widget.plotItem.plot(cutoff_line, pen=pg.mkPen(color="#f83021", width=6))
        # styling the plot
        self.graph_widget.setBackground('w')
        self.graph_widget.setTitle(f"Distance Plot of {self.protein_pair_for_analysis.name}", size="23pt")
        styles = {'font-size': '14px'}
        ax_label_y = "Distance in Å"
        self.graph_widget.setLabel('left', ax_label_y, **styles)
        self.graph_widget.setLabel('bottom', "Residue pair no.", **styles)
        self.graph_widget.plotItem.showGrid(x=True, y=True)
        self.ui.main_Layout.addWidget(self.graph_widget)

        # styles
        stylesheet = """
        QDialog {background-color: #F6F4F8;}
        QTableWidget {background-color: white;}
        """
        self.setStyleSheet(stylesheet)

        self.ui.btn_distance_plot_update.clicked.connect(self.update_distance_plot)
        self.ui.btn_distance_plot_save.clicked.connect(self.save_plot_to_file)
        self.setWindowFlag(QtCore.Qt.WindowMaximizeButtonHint, True)
        self.setWindowFlag(QtCore.Qt.WindowCloseButtonHint, True)
        self.setWindowIcon(QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
        self.setWindowTitle("Distance Plot")

    def update_distance_plot(self):
        """This function updates the distance plot

        """
        from_aa = int(self.ui.sp_distance_plot_from.text())
        to_aa = int(self.ui.sp_distance_plot_to.text())
        from_range = float(self.ui.dsp_distance_plot_from_range.text().replace(",", "."))
        to_range = float(self.ui.dsp_distance_plot_to_range.text().replace(",", "."))
        self.view_box.setRange(xRange=[from_aa, to_aa], yRange=[from_range, to_range])
        if self.ui.cb_sync_with_pymol.isChecked():
            print("Sync is active.")
            zoom_selection = f"/{self.protein_pair_for_analysis.protein_1.get_molecule_object()}///{from_aa}-{to_aa}/CA"
            cmd.select("zoom_sele", zoom_selection)
            cmd.zoom("zoom_sele")

    def save_plot_to_file(self):
        # create an exporter instance, as an argument give it
        # the item you wish to export
        exporter = pg.exporters.ImageExporter(self.graph_widget.plotItem)

        # set export parameters if needed
        exporter.parameters()['width'] = 1000  # (note this also affects height parameter)
        exporter.parameters()['height'] = 1440
        # save to file
        # TODO: save file into xml or png with file dialog
        results_image_path = pathlib.Path(
            f"{self.protein_pair_for_analysis.results_dir}/user_images")
        if not os.path.exists(results_image_path):
            os.mkdir(results_image_path)
        file_path = Qt.QtWidgets.QFileDialog.getSaveFileName(self, "Save Plot as Image", str(results_image_path), "Portable Network Graphic (.png)")
        if file_path == ("", ""):
            return
        exporter.export(f'{file_path[0]}.png')
