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
"""Module for the main view of the pyssa plugin."""
import pathlib

from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from pyqtspinner import spinner

from pyssa.gui.ui.forms.auto_generated import auto_main_view
from pyssa.gui.ui.styles import styles
from pyssa.internal.data_structures import project
from pyssa.util import constants, gui_utils
from pyssa.model import application_model
from pyssa.util import enums


class MainView(QtWidgets.QMainWindow):
    """Class representing the main view of PySSA."""

    """
    The spinner that shows up if something runs as Task
    """
    wait_spinner: spinner.WaitingSpinner

    def __init__(self) -> None:
        """Constructor."""
        super().__init__()
        # build ui object
        self.ui = auto_main_view.Ui_MainWindow()
        self.ui.setupUi(self)
        self.status_bar = QtWidgets.QStatusBar()
        self.wait_spinner = spinner.WaitingSpinner(
            parent=self,
            center_on_parent=True,
            disable_parent_when_spinning=True,
            modality=Qt.ApplicationModal,
            roundness=100.0,
            fade=45.0,
            radius=14,
            lines=8,
            line_length=17,
            line_width=10,
            speed=1.25,
            color=QtGui.QColor(75, 145, 247),
        )
        self.initialize_ui()

    def initialize_ui(self) -> None:
        """Initialize the UI elements."""
        self.ui.lbl_project_name.hide()
        self.ui.project_tab_widget.hide()

        self.ui.action_close_project.setEnabled(False)
        self.ui.action_predict_monomer.setEnabled(False)
        self.ui.action_predict_multimer.setEnabled(False)
        self.ui.action_distance_analysis.setEnabled(False)
        self.ui.action_results_summary.setEnabled(False)
        self.ui.actionPreview.setEnabled(False)
        self.ui.actionCreate_ray_ta.setEnabled(False)
        self.ui.actionCreate_drawn.setEnabled(False)
        self.ui.actionPlaceholder.setEnabled(False)

        self.setMinimumWidth(700)
        self.setMinimumHeight(700)
        pixmap = QtGui.QPixmap(str(constants.PLUGIN_LOGO_WITH_FONT_FILEPATH))
        # Resize the pixmap
        scaled_pixmap = pixmap.scaled(700, 700, aspectRatioMode=Qt.KeepAspectRatio, transformMode=Qt.SmoothTransformation)
        # Set the scaled pixmap to the QLabel
        self.ui.lbl_logo.setPixmap(scaled_pixmap)
        self.ui.lbl_logo.setAlignment(Qt.AlignCenter)
        styles.set_stylesheet(self)
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        self.setWindowTitle("PySSA")
        # constants.PYSSA_LOGGER.info(f"PySSA started with version {constants.VERSION_NUMBER}.")
        # constants.PYSSA_LOGGER.info("Successful initialization of basic UI.")

    def refresh_ui_based_on_app_model(self, an_app_model: "application_model.ApplicationModel"):
        """Modifies the UI of the main view based on an app model."""
        tmp_project: "project.Project" = an_app_model.application_state[enums.ApplicationModelEnum.PROJECT]
        self.ui.lbl_logo.hide()
        if tmp_project.get_project_name() != "":
            self.ui.lbl_project_name.show()
            self.ui.lbl_project_name.setText(f"Project name: {tmp_project.get_project_name()}")
            self.ui.project_tab_widget.show()
            self.ui.action_new_project.setEnabled(False)
            self.ui.action_open_project.setEnabled(False)
            self.ui.action_delete_project.setEnabled(False)
            self.ui.actionImport.setEnabled(False)
            self.ui.actionExport.setEnabled(False)
            self.ui.action_close_project.setEnabled(True)
            self.ui.actionPreview.setEnabled(True)
            self.ui.actionCreate_ray_ta.setEnabled(True)
            self.ui.actionCreate_drawn.setEnabled(True)
            self.ui.actionPlaceholder.setEnabled(True)
            self.ui.project_tab_widget.setCurrentIndex(0)
        else:
            self.ui.lbl_project_name.hide()
            self.ui.project_tab_widget.hide()
            self.ui.lbl_logo.show()
            self.ui.action_new_project.setEnabled(True)
            self.ui.action_open_project.setEnabled(True)
            self.ui.action_delete_project.setEnabled(True)
            self.ui.actionImport.setEnabled(True)
            self.ui.actionExport.setEnabled(True)
            self.ui.action_close_project.setEnabled(False)
            self.ui.actionPreview.setEnabled(False)
            self.ui.actionCreate_ray_ta.setEnabled(False)
            self.ui.actionCreate_drawn.setEnabled(False)
            self.ui.actionPlaceholder.setEnabled(False)

        if len(tmp_project.proteins) > 0:
            list_model = QtGui.QStandardItemModel()
            for protein in tmp_project.proteins:
                item = QtGui.QStandardItem(protein.get_molecule_object())
                list_model.appendRow(item)
            self.ui.proteins_list_view.setModel(list_model)

            # Set row and column count
            self.ui.proteins_table_widget.setRowCount(len(tmp_project.proteins))
            self.ui.proteins_table_widget.setColumnCount(2)  # Assuming two columns: attribute name and attribute value

            # Set header labels
            self.ui.proteins_table_widget.setHorizontalHeaderLabels(["Parameter Name", "Parameter Value"])

            item_names = ["color", "name"]
            for row, protein in enumerate(tmp_project.proteins):
                # Set attribute name in the first column
                item_name = QtWidgets.QTableWidgetItem(item_names[row])
                self.ui.proteins_table_widget.setItem(row, 0, item_name)

                # Set attribute value in the second column
                item_value = QtWidgets.QTableWidgetItem(str(getattr(protein, list(protein.__annotations__.keys())[row])))
                if row == 0:
                    self.ui.proteins_table_widget.setItem(row, 1, QtWidgets.QTableWidgetItem("yellow"))
                else:
                    self.ui.proteins_table_widget.setItem(row, 1, item_value)

            combo = QtWidgets.QComboBox()
            gui_utils.fill_combo_box(combo, ["red", "green", "yellow"])
            self.ui.proteins_table_widget.setCellWidget(0, 1, combo)
            combo.adjustSize()

            for row in range(self.ui.proteins_table_widget.rowCount()):
                item = self.ui.proteins_table_widget.item(row, 0)
                if item:
                    item.setFlags(item.flags() & ~Qt.ItemIsEditable)

            self.ui.proteins_table_widget.verticalHeader().setVisible(False)
            self.ui.proteins_table_widget.resizeColumnsToContents()

    #
    # def _create_all_tooltips(self) -> None:
    #     """Creates all tooltips for the gui elements."""
    #     self.status_bar.setToolTip("Status information: Current process")
    #     # new project page
    #     self.ui.btn_new_choose_reference.setToolTip("Click to add a .pdb file")
    #     # sidebar
    #     self.ui.lbl_current_project_name.setToolTip("Name of the current project")
    #     # edit page
    #     self.ui.btn_edit_project_save.setToolTip("Save as a .pdb file")
    #     # view page
    #     # fixme: is this important? self.ui.list_view_project_proteins.setToolTip("Proteins of the current project")
    #     self.ui.txtedit_view_sequence.setToolTip("Protein sequence of the selected protein")
    #     # use page
    #     self.ui.txt_use_search.setToolTip("Enter a protein name to search in your current workspace")
    #     # prediction Monomer
    #     self.ui.table_pred_mono_prot_to_predict.setToolTip("Protein monomers which get predicted")
    #     self.ui.btn_pred_mono_seq_to_predict.setToolTip("Set up a protein which can be used for a prediction")
    #     self.ui.table_pred_analysis_mono_prot_to_predict.setToolTip("Protein monomers which get predicted")
    #     self.ui.list_pred_analysis_mono_overview.setToolTip("Protein pairs which get analyzed")
    #     self.ui.btn_pred_analysis_mono_seq_to_predict.setToolTip("Set up a protein which can be used for a prediction")
    #     # prediction Multimer
    #     self.ui.table_pred_multi_prot_to_predict.setToolTip("Protein multimers which get predicted")
    #     self.ui.btn_pred_multi_prot_to_predict_add.setToolTip("Set up a protein which can be used for a prediction")
    #     self.ui.table_pred_analysis_multi_prot_to_predict.setToolTip("Protein multimers which get predicted")
    #     self.ui.list_pred_analysis_multi_overview.setToolTip("Protein pairs which get analyzed")
    #     self.ui.btn_pred_analysis_multi_prot_to_predict_add.setToolTip(
    #         "Set up a protein which can be used for a prediction",
    #     )
    #     # image page
    #     self.ui.btn_save_scene.setToolTip("Create new PyMOL scene")
    #     self.ui.btn_update_scene.setToolTip("Overwrite current scene")
    #     self.ui.btn_preview_image.setToolTip("Preview current viewpoint")
    #     self.ui.btn_save_image.setToolTip("Save current viewpoint as png file")
    #     self.ui.cb_ray_tracing.setToolTip("Enable ray-tracing")
    #     self.ui.cb_transparent_bg.setToolTip("Enable transparent background")
    #     self.ui.box_representation.setToolTip("Choose a representation")
    #     self.ui.box_bg_color.setToolTip("Choose a background color")
    #     self.ui.box_renderer.setToolTip("Choose a ray-tracing renderer")
    #     self.ui.box_ray_trace_mode.setToolTip("Choose a ray-trace mode")
    #
    # def quit_app(self) -> None:
    #     """Closes the entire plugin."""
    #     self.close()
