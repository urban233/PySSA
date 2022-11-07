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
from PyQt5.QtWidgets import QHBoxLayout
from pymol import Qt
import pyqtgraph as pg
from uiForms.auto.auto_dialog_distance_plot import Ui_Dialog
from utils import global_utils


class DialogDistancePlot(Qt.QtWidgets.QDialog):

    def __init__(self, parent=None):
        """Constructor

        Args:
            args
            kwargs
        """
        Qt.QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        graph_widget = pg.PlotWidget()

        # read csv file
        file_path = global_utils.global_var_tmp_project_info[0]
        model_name = global_utils.global_var_tmp_project_info[1]

        path = f"{file_path}/distance_csv/distances.csv"
        distance_list = []
        cutoff_line = []
        with open(path, 'r', encoding="utf-8") as csv_file:
            for line in csv_file:
                cleaned_line = line.replace("\n", "")
                if cleaned_line.split(",")[8] != 'distance':
                    distance_list.append(float(cleaned_line.split(",")[8]))
                    cutoff_line.append(global_utils.global_var_settings_obj.get_cutoff())
        # creates actual distance plot line
        graph_widget.plotItem.plot(distance_list, pen=pg.mkPen(color="#4B91F7", width=6),
                                   symbol="o", symbolSize=10, symbolBrush=('b'))
        self.view_box = graph_widget.plotItem.getViewBox()
        # creates cutoff line
        graph_widget.plotItem.plot(cutoff_line, pen=pg.mkPen(color="#f83021", width=6))
        # styling the plot
        graph_widget.setBackground('w')
        graph_widget.setTitle(f"Distance Plot of {model_name}", size="23pt")
        styles = {'font-size': '14px'}
        ax_label_y = "Distance in Å"
        graph_widget.setLabel('left', ax_label_y, **styles)
        graph_widget.setLabel('bottom', "Residue pair no.", **styles)
        graph_widget.plotItem.showGrid(x=True, y=True)
        self.ui.main_Layout.addWidget(graph_widget)

        self.ui.btn_distance_plot_update.clicked.connect(self.update_distance_plot)

        self.setWindowTitle("Distance Plot")

    def update_distance_plot(self):
        """This function updates the distance plot

        """
        from_aa = int(self.ui.sp_distance_plot_from.text())
        to_aa = int(self.ui.sp_distance_plot_to.text())
        from_range = float(self.ui.dsp_distance_plot_from_range.text())
        to_range = float(self.ui.dsp_distance_plot_to_range.text())
        self.view_box.setRange(xRange=[from_aa, to_aa], yRange=[from_range, to_range])
        # if self.ui.cb_sync_with_pymol.isChecked():
        #     print("Sync is active.")
        #     # TODO: handle objects better
        #     pdb_list = os.listdir(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb")
        #     pdb_name = pdb_list[0].replace(".pdb", "")
        #     zoom_selection = f"/{pdb_name}///{from_aa}-{to_aa}/CA"
        #     cmd.select("zoom_sele", zoom_selection)
        #     cmd.zoom("zoom_sele")
