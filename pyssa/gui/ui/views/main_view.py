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
from pyssa.gui.ui.custom_widgets import custom_line_edit
from pyssa.gui.ui.forms.auto_generated import auto_main_view
from pyssa.gui.ui.styles import styles
from pyssa.gui.ui import icon_resources  # this import is used for the icons! DO NOT DELETE THIS
from pyssa.util import constants, gui_utils


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
        self.ui.lbl_session_name.hide()
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

        # fixme: the tutorial action of the menubar IS HIDDEN!!!
        self.ui.action_tutorials.setVisible(False)

        # Sequences tab
        self.ui.btn_save_sequence.setEnabled(False)
        self.ui.btn_delete_sequence.setEnabled(False)
        self.ui.seqs_list_view.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)

        # Proteins tab
        self.ui.btn_save_protein.setEnabled(False)
        self.ui.btn_delete_protein.setEnabled(False)
        self.ui.btn_open_protein_session.setEnabled(False)
        self.ui.btn_create_protein_scene.setEnabled(False)
        self.ui.btn_update_protein_scene.setEnabled(False)
        self.ui.proteins_tree_view.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)

        # Protein Pairs tab
        self.ui.btn_delete_protein_pair.setEnabled(False)
        self.ui.btn_open_protein_pair_session.setEnabled(False)
        self.ui.btn_create_protein_pair_scene.setEnabled(False)
        self.ui.btn_update_protein_pair_scene.setEnabled(False)
        self.ui.protein_pairs_tree_view.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)

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

        # <editor-fold desc="Set icons">
        pixmapi = QtWidgets.QStyle.SP_MessageBoxQuestion
        # <editor-fold desc="Help">
        icon = self.style().standardIcon(pixmapi)
        self.ui.btn_help.setIcon(QtGui.QIcon(":/icons/help_w200.svg"))
        self.ui.btn_help.setIconSize(self.ui.btn_help.icon().actualSize(QtCore.QSize(30, 30)))
        self.ui.btn_help.setText("")
        self.ui.btn_help_2.setIcon(QtGui.QIcon(":/icons/help_w200.svg"))
        self.ui.btn_help_2.setIconSize(self.ui.btn_help_2.icon().actualSize(QtCore.QSize(30, 30)))
        self.ui.btn_help_2.setText("")
        self.ui.btn_help_3.setIcon(QtGui.QIcon(":/icons/help_w200.svg"))
        self.ui.btn_help_3.setIconSize(self.ui.btn_help_3.icon().actualSize(QtCore.QSize(30, 30)))
        self.ui.btn_help_3.setText("")
        # </editor-fold>
        # Fixme: The icon size is not adjustable!!!
        # <editor-fold desc="Sequence">
        # add
        add_sequence_icon = QtGui.QIcon(QtGui.QPixmap(":icons/note_add_w200.svg"))
        add_sequence_icon.addPixmap(QtGui.QPixmap(":icons/note_add_disabled_w200.svg"),
                                    mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_add_sequence.setIcon(add_sequence_icon)
        self.ui.btn_add_sequence.setText("")
        self.ui.btn_add_sequence.setIconSize(add_sequence_icon.actualSize(QtCore.QSize(30, 30)))
        # self.ui.btn_add_sequence.setIcon(
        #     QtGui.QIcon(":/icons/note_add_w200.svg")
        # )
        # self.ui.btn_add_sequence.setText("")
        # self.ui.btn_add_sequence.setIconSize(self.ui.btn_add_sequence.icon().actualSize(QtCore.QSize(30, 30)))

        # import
        import_sequence_icon = QtGui.QIcon(QtGui.QPixmap(":icons/upload_file_w200.svg"))
        import_sequence_icon.addPixmap(QtGui.QPixmap(":icons/upload_file_disabled_w200.svg"),
                                       mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_import_seq.setIcon(import_sequence_icon)
        self.ui.btn_import_seq.setText("")
        self.ui.btn_import_seq.setIconSize(import_sequence_icon.actualSize(QtCore.QSize(30, 30)))
        # self.ui.btn_import_seq.setIcon(
        #     QtGui.QIcon(":/icons/upload_file_w200.svg")
        # )
        # self.ui.btn_import_seq.setText("")
        # self.ui.btn_import_seq.setIconSize(self.ui.btn_import_seq.icon().actualSize(QtCore.QSize(30, 30)))

        # save
        save_seq_icon = QtGui.QIcon(QtGui.QPixmap(":icons/file_save_w200.svg"))
        save_seq_icon.addPixmap(QtGui.QPixmap(":icons/file_save_disabled_w200.svg"), mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_save_sequence.setIcon(save_seq_icon)
        self.ui.btn_save_sequence.setText("")
        self.ui.btn_save_sequence.setIconSize(save_seq_icon.actualSize(QtCore.QSize(30, 30)))
        # save_icon = QtGui.QIcon(":/icons/file_save_w200.svg")
        # save_icon.addPixmap(QtGui.QPixmap(":icons/file_save_disabled_w200.svg"), mode=QtGui.QIcon.Mode.Disabled)
        # self.ui.btn_save_sequence.setIcon(save_icon)
        # self.ui.btn_save_sequence.setText("")
        # self.ui.btn_save_sequence.setIconSize(self.ui.btn_save_sequence.icon().actualSize(QtCore.QSize(30, 30)))
        # self.ui.btn_save_sequence.setIcon(QtGui.QIcon(":/icons/file_save_w200.svg"))

        # delete
        delete_seq_icon = QtGui.QIcon(QtGui.QPixmap(":icons/scan_delete_w200.svg"))
        delete_seq_icon.addPixmap(QtGui.QPixmap(":icons/scan_delete_disabled_w200.svg"), mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_delete_sequence.setIcon(delete_seq_icon)
        self.ui.btn_delete_sequence.setText("")
        self.ui.btn_delete_sequence.setIconSize(delete_seq_icon.actualSize(QtCore.QSize(30, 30)))
        # self.ui.btn_delete_sequence.setIcon(
        #     QtGui.QIcon(":/icons/scan_deletew200.svg"))
        # self.ui.btn_delete_sequence.setText("")
        # self.ui.btn_delete_sequence.setIconSize(self.ui.btn_delete_sequence.icon().actualSize(QtCore.QSize(30, 30)))
        # </editor-fold>

        # <editor-fold desc="Protein">
        # import
        import_protein_icon = QtGui.QIcon(QtGui.QPixmap(":icons/upload_file_w200.svg"))
        import_protein_icon.addPixmap(QtGui.QPixmap(":icons/upload_file_disabled_w200.svg"),
                                      mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_import_protein.setIcon(import_protein_icon)
        self.ui.btn_import_protein.setText("")
        self.ui.btn_import_protein.setIconSize(import_protein_icon.actualSize(QtCore.QSize(30, 30)))
        # self.ui.btn_import_protein.setIcon(
        #     QtGui.QIcon(":/icons/upload_file_w200.svg")
        # )
        # self.ui.btn_import_protein.setText("")
        # self.ui.btn_import_protein.setIconSize(self.ui.btn_import_protein.icon().actualSize(QtCore.QSize(30, 30)))

        # save
        save_protein_icon = QtGui.QIcon(QtGui.QPixmap(":icons/file_save_w200.svg"))
        save_protein_icon.addPixmap(QtGui.QPixmap(":icons/file_save_disabled_w200.svg"), mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_save_protein.setIcon(save_protein_icon)
        self.ui.btn_save_protein.setText("")
        self.ui.btn_save_protein.setIconSize(save_protein_icon.actualSize(QtCore.QSize(30, 30)))
        # self.ui.btn_save_protein.setIcon(
        #     QtGui.QIcon(":/icons/file_save_w200.svg")
        # )
        # self.ui.btn_save_protein.setText("")
        # self.ui.btn_save_protein.setIconSize(self.ui.btn_save_protein.icon().actualSize(QtCore.QSize(30, 30)))

        # delete
        delete_protein_icon = QtGui.QIcon(QtGui.QPixmap(":icons/scan_delete_w200.svg"))
        delete_protein_icon.addPixmap(QtGui.QPixmap(":icons/scan_delete_disabled_w200.svg"), mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_delete_protein.setIcon(delete_protein_icon)
        self.ui.btn_delete_protein.setText("")
        self.ui.btn_delete_protein.setIconSize(delete_protein_icon.actualSize(QtCore.QSize(30, 30)))
        # self.ui.btn_delete_protein.setIcon(
        #     QtGui.QIcon(":/icons/scan_delete_w200.svg"))
        # self.ui.btn_delete_protein.setText("")
        # self.ui.btn_delete_protein.setIconSize(
        #     self.ui.btn_delete_protein.icon().actualSize(QtCore.QSize(30, 30)))
        # </editor-fold>

        # <editor-fold desc="Protein Session">
        # open
        open_protein_session_icon = QtGui.QIcon(QtGui.QPixmap(":/icons/open_in_new_w200.svg"))
        open_protein_session_icon.addPixmap(QtGui.QPixmap(":/icons/open_in_new_disabled_w200.svg"),
                                            mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_open_protein_session.setIcon(open_protein_session_icon)
        self.ui.btn_open_protein_session.setText("")
        self.ui.btn_open_protein_session.setIconSize(open_protein_session_icon.actualSize(QtCore.QSize(30, 30)))
        # self.ui.btn_open_protein_session.setIcon(
        #     QtGui.QIcon(":/icons/open_in_new_w200.svg"))
        # self.ui.btn_open_protein_session.setText("")
        # self.ui.btn_open_protein_session.setIconSize(
        #     self.ui.btn_open_protein_session.icon().actualSize(QtCore.QSize(30, 30)))

        # create
        create_protein_session_icon = QtGui.QIcon(QtGui.QPixmap(":/icons/add_circle_w200.svg"))
        create_protein_session_icon.addPixmap(QtGui.QPixmap(":/icons/add_circle_disabled_w200.svg"),
                                            mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_create_protein_scene.setIcon(create_protein_session_icon)
        self.ui.btn_create_protein_scene.setText("")
        self.ui.btn_create_protein_scene.setIconSize(create_protein_session_icon.actualSize(QtCore.QSize(30, 30)))
        # self.ui.btn_create_protein_scene.setIcon(QtGui.QIcon(":/icons/add_circle_w200.svg"))
        # self.ui.btn_create_protein_scene.setText("")
        # self.ui.btn_create_protein_scene.setIconSize(
        #     self.ui.btn_create_protein_scene.icon().actualSize(QtCore.QSize(30, 30)))

        # update
        update_protein_session_icon = QtGui.QIcon(QtGui.QPixmap(":/icons/change_circle_w200.svg"))
        update_protein_session_icon.addPixmap(QtGui.QPixmap(":/icons/change_circle_disabled_w200.svg"),
                                              mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_update_protein_scene.setIcon(update_protein_session_icon)
        self.ui.btn_update_protein_scene.setText("")
        self.ui.btn_update_protein_scene.setIconSize(update_protein_session_icon.actualSize(QtCore.QSize(30, 30)))
        # self.ui.btn_update_protein_scene.setIcon(
        #     QtGui.QIcon(":/icons/change_circle_w200.svg"))
        # self.ui.btn_update_protein_scene.setText("")
        # self.ui.btn_update_protein_scene.setIconSize(
        #     self.ui.btn_update_protein_scene.icon().actualSize(QtCore.QSize(30, 30)))
        # </editor-fold>

        # <editor-fold desc="Protein Pair Session">
        # open
        open_protein_pair_session_icon = QtGui.QIcon(QtGui.QPixmap(":/icons/open_in_new_w200.svg"))
        open_protein_pair_session_icon.addPixmap(QtGui.QPixmap(":/icons/open_in_new_disabled_w200.svg"),
                                            mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_open_protein_pair_session.setIcon(open_protein_pair_session_icon)
        self.ui.btn_open_protein_pair_session.setText("")
        self.ui.btn_open_protein_pair_session.setIconSize(open_protein_pair_session_icon.actualSize(QtCore.QSize(30, 30)))
        # self.ui.btn_open_protein_pair_session.setIcon(
        #     QtGui.QIcon(":/icons/open_in_new_w200.svg"))
        # self.ui.btn_open_protein_pair_session.setText("")
        # self.ui.btn_open_protein_pair_session.setIconSize(
        #     self.ui.btn_open_protein_pair_session.icon().actualSize(QtCore.QSize(30, 30)))

        # create
        create_protein_pair_session_icon = QtGui.QIcon(QtGui.QPixmap(":/icons/add_circle_w200.svg"))
        create_protein_pair_session_icon.addPixmap(QtGui.QPixmap(":/icons/add_circle_disabled_w200.svg"),
                                              mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_create_protein_pair_scene.setIcon(create_protein_pair_session_icon)
        self.ui.btn_create_protein_pair_scene.setText("")
        self.ui.btn_create_protein_pair_scene.setIconSize(create_protein_pair_session_icon.actualSize(QtCore.QSize(30, 30)))
        # self.ui.btn_create_protein_pair_scene.setIcon(
        #     QtGui.QIcon(":/icons/add_circle_w200.svg"))
        # self.ui.btn_create_protein_pair_scene.setText("")
        # self.ui.btn_create_protein_pair_scene.setIconSize(
        #     self.ui.btn_create_protein_pair_scene.icon().actualSize(QtCore.QSize(30, 30)))

        # update
        update_protein_pair_session_icon = QtGui.QIcon(QtGui.QPixmap(":/icons/change_circle_w200.svg"))
        update_protein_pair_session_icon.addPixmap(QtGui.QPixmap(":/icons/change_circle_disabled_w200.svg"),
                                              mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_update_protein_pair_scene.setIcon(update_protein_pair_session_icon)
        self.ui.btn_update_protein_pair_scene.setText("")
        self.ui.btn_update_protein_pair_scene.setIconSize(update_protein_pair_session_icon.actualSize(QtCore.QSize(30, 30)))
        # self.ui.btn_update_protein_pair_scene.setIcon(
        #     QtGui.QIcon(":/icons/change_circle_w200.svg"))
        # self.ui.btn_update_protein_pair_scene.setText("")
        # self.ui.btn_update_protein_pair_scene.setIconSize(
        #     self.ui.btn_update_protein_pair_scene.icon().actualSize(QtCore.QSize(30, 30)))

        # delete
        delete_protein_pair_icon = QtGui.QIcon(QtGui.QPixmap(":icons/scan_delete_w200.svg"))
        delete_protein_pair_icon.addPixmap(QtGui.QPixmap(":icons/scan_delete_disabled_w200.svg"), mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_delete_protein_pair.setIcon(delete_protein_pair_icon)
        self.ui.btn_delete_protein_pair.setText("")
        self.ui.btn_delete_protein_pair.setIconSize(delete_protein_pair_icon.actualSize(QtCore.QSize(40, 40)))
        # self.ui.btn_delete_protein_pair.setIcon(
        #     QtGui.QIcon(":/icons/scan_deletew200.svg"))
        # self.ui.btn_delete_protein_pair.setText("")
        # self.ui.btn_delete_protein_pair.setIconSize(
        #     self.ui.btn_delete_protein_pair.icon().actualSize(QtCore.QSize(30, 30)))
        # </editor-fold>
        # </editor-fold>

        self._create_all_tooltips()
        pixmap = QtGui.QPixmap(str(constants.PLUGIN_LOGO_WITH_FONT_FILEPATH))
        # Resize the pixmap
        pixmap = QtGui.QPixmap(
            r"C:\ProgramData\pyssa\mambaforge_pyssa\pyssa-mamba-env\Lib\site-packages\pymol\pymol_path\data\startup\PySSA\assets\images\splash_screen_logo_002.png")
        scaled_pixmap = pixmap.scaled(700, 700, aspectRatioMode=Qt.KeepAspectRatio,
                                      transformMode=Qt.SmoothTransformation)
        # Set the scaled pixmap to the QLabel
        self.ui.lbl_logo.setPixmap(scaled_pixmap)
        self.ui.lbl_logo.setAlignment(Qt.AlignCenter)
        styles.set_stylesheet(self)
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        self.setWindowTitle("PySSA")
        constants.PYSSA_LOGGER.info(f"PySSA started with version {constants.VERSION_NUMBER}.")
        constants.PYSSA_LOGGER.info("Successful initialization of basic UI.")

    def build_sequence_table(self):
        self.line_edit_seq_name = custom_line_edit.CustomLineEdit()

        self.ui.seqs_table_widget.verticalHeader().setVisible(False)
        self.ui.seqs_table_widget.setColumnCount(2)
        self.ui.seqs_table_widget.setHorizontalHeaderLabels(["Name", "Value"])

    def build_proteins_table(self):
        self.cb_chain_color = QtWidgets.QComboBox()
        self.cb_chain_representation = QtWidgets.QComboBox()

        self.ui.proteins_table_widget.verticalHeader().setVisible(False)
        self.ui.proteins_table_widget.setColumnCount(2)
        self.ui.proteins_table_widget.setHorizontalHeaderLabels(["Name", "Value"])
        gui_utils.fill_combo_box(self.cb_chain_color, constants.PYMOL_COLORS)
        self.cb_chain_color.adjustSize()
        gui_utils.fill_combo_box(self.cb_chain_representation, constants.PYMOL_REPRESENTATIONS)
        self.cb_chain_representation.adjustSize()

    def build_protein_pairs_table(self):
        self.cb_chain_color_protein_pair = QtWidgets.QComboBox()
        self.cb_chain_representation_protein_pair = QtWidgets.QComboBox()

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
        # self.ui.seqs_list_view.setToolTip("A list of all sequences in the project")
        self.ui.seqs_table_widget.setToolTip("A table with additional information about the selected sequence")
        self.ui.btn_import_seq.setToolTip("Import an existing .fasta file")
        self.ui.btn_add_sequence.setToolTip("Add a sequence by pasting the sequence string")
        self.ui.btn_save_sequence.setToolTip("Save the selected sequence as .fasta file")
        self.ui.btn_delete_sequence.setToolTip("Delete the selected sequence from the project")

        # self.ui.proteins_tree_view.setToolTip("A tree of all proteins in the project")
        # self.ui.proteins_table_widget.setToolTip(
        #     "A table with changeable PyMOL parameters for the currently active session"
        # )
        self.ui.btn_import_protein.setToolTip("Import an existing .pdb file")
        self.ui.btn_save_protein.setToolTip("Save the selected protein as .pdb file")
        self.ui.btn_delete_protein.setToolTip("Delete the selected protein from the project")
        self.ui.btn_open_protein_session.setToolTip("Open protein PyMOL session")
        self.ui.btn_create_protein_scene.setToolTip("Create a new PyMOL scene")
        self.ui.btn_update_protein_scene.setToolTip("Update the current scene in PyMOL")

        # self.ui.protein_pairs_tree_view.setToolTip("A tree of all protein pairs in the project")
        # self.ui.protein_pairs_table_widget.setToolTip(
        #     "A table with changeable PyMOL parameters for the currently active session"
        # )
        self.ui.btn_delete_protein_pair.setToolTip("Delete the selected protein pair from the project")
        self.ui.btn_open_protein_pair_session.setToolTip("Open protein pair PyMOL session")
        self.ui.btn_create_protein_pair_scene.setToolTip("Create a new PyMOL scene")
        self.ui.btn_update_protein_pair_scene.setToolTip("Update the current scene in PyMOL")

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
