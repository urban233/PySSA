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
from pyssa.gui.ui.custom_widgets import custom_tree_view, custom_line_edit
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
        self.ui.action_use_project.setEnabled(False)
        self.ui.action_export_project.setEnabled(False)
        self.ui.action_close_project.setEnabled(False)
        self.ui.action_predict_monomer.setEnabled(False)
        self.ui.action_predict_multimer.setEnabled(False)
        self.ui.action_distance_analysis.setEnabled(False)
        self.ui.action_results_summary.setEnabled(False)
        self.ui.action_preview_image.setEnabled(False)
        self.ui.action_ray_tracing_image.setEnabled(False)
        self.ui.action_simple_image.setEnabled(False)
        self.ui.action_protein_regions.setEnabled(False)

        # Sequences tab
        self.ui.btn_save_sequence.setEnabled(False)
        self.ui.btn_delete_sequence.setEnabled(False)

        # Proteins tab
        self.ui.btn_save_protein.setEnabled(False)
        self.ui.btn_delete_protein.setEnabled(False)
        self.ui.btn_open_protein_session.setEnabled(False)
        self.ui.btn_create_protein_scene.setEnabled(False)
        self.ui.btn_update_protein_scene.setEnabled(False)

        # Protein Pairs tab
        self.ui.btn_delete_protein_pair.setEnabled(False)
        self.ui.btn_open_protein_pair_session.setEnabled(False)
        self.ui.btn_create_protein_pair_scene.setEnabled(False)
        self.ui.btn_update_protein_pair_scene.setEnabled(False)

        # Extra UI elements
        self.cb_chain_color = QtWidgets.QComboBox()
        self.cb_chain_representation = QtWidgets.QComboBox()
        self.cb_chain_color_protein_pair = QtWidgets.QComboBox()
        self.cb_chain_representation_protein_pair = QtWidgets.QComboBox()
        self.line_edit_seq_name = custom_line_edit.CustomLineEdit()
        self.build_sequence_table()
        self.build_proteins_table()
        self.build_protein_pairs_table()

        self.ui.proteins_tree_view.setContextMenuPolicy(Qt.CustomContextMenu)
        self.setMinimumWidth(700)
        self.setMinimumHeight(700)

        pixmapi = QtWidgets.QStyle.SP_MessageBoxQuestion
        icon = self.style().standardIcon(pixmapi)
        self.ui.btn_help.setIcon(icon)
        self.ui.btn_help.setText("")
        self.ui.btn_help_2.setIcon(icon)
        self.ui.btn_help_2.setText("")
        self.ui.btn_help_3.setIcon(icon)
        self.ui.btn_help_3.setText("")

        self._create_all_tooltips()
        pixmap = QtGui.QPixmap(str(constants.PLUGIN_LOGO_WITH_FONT_FILEPATH))
        # Resize the pixmap
        scaled_pixmap = pixmap.scaled(700, 700, aspectRatioMode=Qt.KeepAspectRatio, transformMode=Qt.SmoothTransformation)
        # Set the scaled pixmap to the QLabel
        self.ui.lbl_logo.setPixmap(scaled_pixmap)
        self.ui.lbl_logo.setAlignment(Qt.AlignCenter)
        styles.set_stylesheet(self)
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        self.setWindowTitle("PySSA")
        constants.PYSSA_LOGGER.info(f"PySSA started with version {constants.VERSION_NUMBER}.")
        constants.PYSSA_LOGGER.info("Successful initialization of basic UI.")

    def build_sequence_table(self):
        self.ui.seqs_table_widget.verticalHeader().setVisible(False)
        self.ui.seqs_table_widget.setColumnCount(2)
        self.ui.seqs_table_widget.setHorizontalHeaderLabels(["Name", "Value"])

    def build_proteins_table(self):
        self.ui.proteins_table_widget.verticalHeader().setVisible(False)
        self.ui.proteins_table_widget.setColumnCount(2)
        self.ui.proteins_table_widget.setHorizontalHeaderLabels(["Name", "Value"])
        gui_utils.fill_combo_box(self.cb_chain_color, constants.PYMOL_COLORS)
        self.cb_chain_color.adjustSize()
        gui_utils.fill_combo_box(self.cb_chain_representation, constants.PYMOL_REPRESENTATIONS)
        self.cb_chain_representation.adjustSize()

    def build_protein_pairs_table(self):
        self.ui.protein_pairs_table_widget.setColumnCount(2)
        self.ui.protein_pairs_table_widget.verticalHeader().setVisible(False)
        self.ui.protein_pairs_table_widget.setHorizontalHeaderLabels(["Name", "Value"])
        gui_utils.fill_combo_box(self.cb_chain_color_protein_pair, constants.PYMOL_COLORS)
        self.cb_chain_color_protein_pair.adjustSize()
        gui_utils.fill_combo_box(self.cb_chain_representation_protein_pair, constants.PYMOL_REPRESENTATIONS)
        self.cb_chain_representation_protein_pair.adjustSize()

    def setup_sequences_table(self, row_count):
        self.line_edit_seq_name.setStyleSheet("QLineEdit { background-color: white; border-radius: 0; }")
        self.ui.seqs_table_widget.setRowCount(row_count)
        self.ui.seqs_table_widget.setCellWidget(0, 1, self.line_edit_seq_name)

    def setup_proteins_table(self, row_count):
        self.ui.proteins_table_widget.setRowCount(row_count)
        self.ui.proteins_table_widget.setCellWidget(0, 1, self.cb_chain_color)
        self.ui.proteins_table_widget.setCellWidget(1, 1, self.cb_chain_representation)

    def setup_protein_pairs_table(self, row_count):
        self.ui.protein_pairs_table_widget.setRowCount(row_count)
        self.ui.protein_pairs_table_widget.setCellWidget(0, 1, self.cb_chain_color_protein_pair)
        self.ui.protein_pairs_table_widget.setCellWidget(1, 1, self.cb_chain_representation_protein_pair)

    def _create_all_tooltips(self) -> None:
        """Creates all tooltips for the gui elements."""
        self.ui.seqs_list_view.setToolTip("A list of all sequences in the project")
        self.ui.seqs_table_widget.setToolTip("A table with additional information about the selected sequence")
        self.ui.btn_import_seq.setToolTip("Click to import an existing .fasta file")
        self.ui.btn_add_sequence.setToolTip("Click to add a sequence by pasting the sequence string")
        self.ui.btn_save_sequence.setToolTip("Click to save the selected sequence as .fasta file")
        self.ui.btn_delete_sequence.setToolTip("Click to delete the selected sequence from the project")

        self.ui.proteins_tree_view.setToolTip("A tree of all proteins in the project")
        self.ui.proteins_table_widget.setToolTip(
            "A table with changeable PyMOL parameters for the currently active session"
        )
        self.ui.btn_import_protein.setToolTip("Click to import an existing .pdb file")
        self.ui.btn_save_protein.setToolTip("Click to save the selected protein as .pdb file")
        self.ui.btn_delete_protein.setToolTip("Click to delete the selected protein from the project")
        self.ui.btn_open_protein_session.setToolTip("Click to open protein PyMOL session")
        self.ui.btn_create_protein_scene.setToolTip("Click to create a new PyMOL scene")
        self.ui.btn_update_protein_scene.setToolTip("Click to update the current scene in PyMOL")

        self.ui.protein_pairs_tree_view.setToolTip("A tree of all protein pairs in the project")
        self.ui.protein_pairs_table_widget.setToolTip(
            "A table with changeable PyMOL parameters for the currently active session"
        )
        self.ui.btn_delete_protein_pair.setToolTip("Click to delete the selected protein pair from the project")
        self.ui.btn_open_protein_pair_session.setToolTip("Click to open protein pair PyMOL session")
        self.ui.btn_create_protein_pair_scene.setToolTip("Click to create a new PyMOL scene")
        self.ui.btn_update_protein_pair_scene.setToolTip("Click to update the current scene in PyMOL")

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
