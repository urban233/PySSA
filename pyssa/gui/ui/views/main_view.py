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
"""Module for the main view of the PySSA plugin."""
import pathlib

from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QAbstractItemView
from pyssa.gui.ui.custom_widgets import custom_line_edit, toggle_button, color_grid, custom_label
from pyssa.gui.ui.forms.auto_generated import auto_main_view
from pyssa.gui.ui.styles import styles
from pyssa.gui.ui import icon_resources  # this import is used for the icons! DO NOT DELETE THIS
from pyssa.util import constants, gui_utils


class MainView(QtWidgets.QMainWindow):
    """Class representing the main view of PySSA."""
    dialogClosed = QtCore.pyqtSignal(tuple)
    """
    The spinner that shows up if something runs as Task
    """
    # wait_spinner: spinner.WaitingSpinner

    def __init__(self) -> None:
        """Constructor."""
        super().__init__()
        # build ui object
        self.ui = auto_main_view.Ui_MainWindow()
        self.ui.setupUi(self)
        self.status_bar = QtWidgets.QStatusBar()
        self.progress_bar = QtWidgets.QProgressBar()
        self.status_bar.addWidget(self.progress_bar)
        self.progress_bar.hide()
        # self.wait_spinner = spinner.WaitingSpinner(
        #     parent=self,
        #     center_on_parent=True,
        #     disable_parent_when_spinning=True,
        #     modality=Qt.ApplicationModal,
        #     roundness=100.0,
        #     fade=45.0,
        #     radius=14,
        #     lines=8,
        #     line_length=17,
        #     line_width=10,
        #     speed=1.25,
        #     color=QtGui.QColor(75, 145, 247),
        # )
        self.add_custom_widgets()
        self.initialize_ui()
        self.add_custom_job_panels()
        gui_utils.fill_combo_box(self.ui.box_protein_color, constants.PYMOL_COLORS)
        gui_utils.fill_combo_box(self.ui.box_protein_pair_color, constants.PYMOL_COLORS)

    def closeEvent(self, event):
        # Emit the custom signal when the window is closed
        self.dialogClosed.emit(("", event))

    def add_custom_job_panels(self):
        # job panel related gui elements
        self.btn_open_job_overview = QtWidgets.QPushButton()
        self.lbl_job_overview = QtWidgets.QLabel("No running jobs.")
        self.btn_open_job_notification = QtWidgets.QPushButton()
        self.lbl_job_notification = QtWidgets.QLabel("No new notifications.")
        self.icon_jobs = QtGui.QIcon(QtGui.QPixmap(":icons/play_circle_w200.png"))
        self.icon_jobs_running = QtGui.QIcon(QtGui.QPixmap(":icons/play_circle_run_w200_blue.png"))
        self.icon_notify = QtGui.QIcon(QtGui.QPixmap(":icons/notifications_w200.png"))
        self.icon_notify_unread = QtGui.QIcon(QtGui.QPixmap(":icons/notifications_unread_w200.png"))
        self.btn_open_job_overview.setIcon(self.icon_jobs)
        self.btn_open_job_overview.setText("")
        self.btn_open_job_overview.setIconSize(self.icon_jobs.actualSize(QtCore.QSize(30, 30)))
        self.btn_open_job_overview.setStyleSheet("""
            QPushButton {
                background-color: rgba(220, 219, 227, 0.01);
                border: none;
                border-radius: 4px;
                min-width: 20px;
                max-width: 20px;
                min-height: 20px;
                max-height: 20px;
            }
            QPushButton::hover {
                background-color: rgba(220, 219, 227, 0.5);
                border: none;
                min-width: 24px;
                max-width: 24px;
                min-height: 24px;
                max-height: 24px;
            }
        """)
        self.lbl_job_overview.setStyleSheet("""padding-top: 30px;""")
        self.btn_open_job_notification.setIcon(self.icon_notify)
        self.btn_open_job_notification.setText("")
        self.btn_open_job_notification.setIconSize(self.icon_notify.actualSize(QtCore.QSize(24, 24)))
        self.btn_open_job_notification.setStyleSheet("""
            QPushButton {
                background-color: rgba(220, 219, 227, 0.01);
                border: none;
                border-radius: 4px;
                min-width: 20px;
                max-width: 20px;
                min-height: 20px;
                max-height: 20px;
            }
            QPushButton::hover {
                background-color: rgba(220, 219, 227, 0.5);
                border: none;
                min-width: 24px;
                max-width: 24px;
                min-height: 24px;
                max-height: 24px;
            }
        """)
        """QPushButton {
                background-color: white;
                border: none;
                border-width: 2px;
                border-radius: 10px;
                padding: 2px;
                min-width: 20px;
                max-width: 20px;
                min-height: 20px;
                max-height: 20px
            }
            QPushButton::hover {
                background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
                                            stop: 0 #4B91F7, stop: 0.4 #367AF6,
                                            stop: 0.5 #367AF6, stop: 1.0 #4B91F7);
                background: white;
                color: white;
                color: #4B91F7;
                border: 2px solid #DCDBE3;
            }"""
        self.lbl_job_notification.setStyleSheet("""padding-top: 30px;""")
        self.ui.job_overview_layout.insertWidget(0, self.lbl_job_overview)  # After inserting the widget count is 2
        self.ui.job_overview_layout.setAlignment(self.lbl_job_overview, QtCore.Qt.AlignHCenter)
        self.ui.job_notification_layout.insertWidget(0,
                                                     self.lbl_job_notification)  # After inserting the widget count is 2
        self.ui.job_notification_layout.setAlignment(self.lbl_job_notification, QtCore.Qt.AlignHCenter)

    def add_custom_widgets(self):
        # Protein
        self.tg_protein_white_bg = toggle_button.ToggleWidget()
        self.ui.protein_white_bg_layout.addWidget(self.tg_protein_white_bg)
        self.tg_protein_color_atoms = toggle_button.ToggleWidget()
        self.ui.protein_atoms_layout.addWidget(self.tg_protein_color_atoms)
        # fixme: this code could become useful
        # self.tg_protein_hydrogen_atoms = toggle_button.ToggleWidget()
        # self.ui.protein_hydrogen_layout.addWidget(self.tg_protein_hydrogen_atoms)
        # representations
        self.tg_protein_cartoon = toggle_button.ToggleWidget()
        self.ui.protein_cartoon_layout.addWidget(self.tg_protein_cartoon)
        self.tg_protein_sticks = toggle_button.ToggleWidget()
        self.ui.protein_sticks_layout.addWidget(self.tg_protein_sticks)
        self.tg_protein_ribbon = toggle_button.ToggleWidget()
        self.ui.protein_ribbon_layout.addWidget(self.tg_protein_ribbon)
        self.tg_protein_lines = toggle_button.ToggleWidget()
        self.ui.protein_lines_layout.addWidget(self.tg_protein_lines)
        self.tg_protein_spheres = toggle_button.ToggleWidget()
        self.ui.protein_spheres_layout.addWidget(self.tg_protein_spheres)
        self.tg_protein_dots = toggle_button.ToggleWidget()
        self.ui.protein_dots_layout.addWidget(self.tg_protein_dots)
        self.tg_protein_mesh = toggle_button.ToggleWidget()
        self.ui.protein_mesh_layout.addWidget(self.tg_protein_mesh)
        self.tg_protein_surface = toggle_button.ToggleWidget()
        self.ui.protein_surface_layout.addWidget(self.tg_protein_surface)

        self.color_grid_proteins = color_grid.PyMOLColorGrid()
        self.ui.verticalLayout_6.insertWidget(2, self.color_grid_proteins)

        # Protein Pair
        self.tg_protein_pair_white_bg = toggle_button.ToggleWidget()
        self.ui.protein_pair_white_bg_layout.addWidget(self.tg_protein_pair_white_bg)
        self.tg_protein_pair_color_atoms = toggle_button.ToggleWidget()
        self.ui.protein_pair_atoms_layout.addWidget(self.tg_protein_pair_color_atoms)
        # fixme: this code could become useful
        # self.tg_protein_pair_hydrogen_atoms = toggle_button.ToggleWidget()
        # self.ui.protein_pair_hydrogen_layout.addWidget(self.tg_protein_pair_hydrogen_atoms)
        # representations
        self.tg_protein_pair_cartoon = toggle_button.ToggleWidget()
        self.ui.protein_pair_cartoon_layout.addWidget(self.tg_protein_pair_cartoon)
        self.tg_protein_pair_sticks = toggle_button.ToggleWidget()
        self.ui.protein_pair_sticks_layout.addWidget(self.tg_protein_pair_sticks)
        self.tg_protein_pair_ribbon = toggle_button.ToggleWidget()
        self.ui.protein_pair_ribbon_layout.addWidget(self.tg_protein_pair_ribbon)
        self.tg_protein_pair_lines = toggle_button.ToggleWidget()
        self.ui.protein_pair_lines_layout.addWidget(self.tg_protein_pair_lines)
        self.tg_protein_pair_spheres = toggle_button.ToggleWidget()
        self.ui.protein_pair_spheres_layout.addWidget(self.tg_protein_pair_spheres)
        self.tg_protein_pair_dots = toggle_button.ToggleWidget()
        self.ui.protein_pair_dots_layout.addWidget(self.tg_protein_pair_dots)
        self.tg_protein_pair_mesh = toggle_button.ToggleWidget()
        self.ui.protein_pair_mesh_layout.addWidget(self.tg_protein_pair_mesh)
        self.tg_protein_pair_surface = toggle_button.ToggleWidget()
        self.ui.protein_pair_surface_layout.addWidget(self.tg_protein_pair_surface)

        self.color_grid_protein_pairs = color_grid.PyMOLColorGrid()
        self.ui.verticalLayout_10.insertWidget(2, self.color_grid_protein_pairs)

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
        self.ui.action_abort_prediction.setVisible(False)

        self.ui.frame_job_overview.hide()
        self.ui.frame_job_notification.hide()

        # Project and Session name
        self.ui.lbl_project_name.setText("Project Name: No Project Opened")
        self.ui.lbl_session_name.setText("Session Name: No Session Loaded")

        # Sequences tab
        self.ui.btn_save_sequence.setEnabled(False)
        self.ui.btn_delete_sequence.setEnabled(False)
        self.ui.seqs_list_view.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)

        # Proteins tab
        self.ui.btn_protein_tree_view_expand.setEnabled(False)
        self.ui.btn_protein_tree_view_collapse.setEnabled(False)
        self.ui.btn_save_protein.setEnabled(False)
        self.ui.btn_delete_protein.setEnabled(False)
        self.ui.btn_open_protein_session.setEnabled(False)
        self.ui.btn_create_protein_scene.setEnabled(False)
        self.ui.btn_update_protein_scene.setEnabled(False)
        self.ui.btn_delete_protein_scene.setEnabled(False)
        self.ui.proteins_tree_view.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.ui.lbl_pymol_protein_scene.setText("PyMOL Scene: No Scene Loaded")

        # Hides all ui elements for the scene modifications
        self.ui.lbl_protein_color.hide()
        self.ui.lbl_protein_atoms.hide()
        self.ui.lbl_protein_cartoon.hide()
        self.ui.lbl_protein_sticks.hide()
        self.ui.lbl_protein_ribbon.hide()
        self.ui.lbl_protein_all_representations.hide()

        self.ui.frame_protein_color.hide()
        self.ui.frame_protein_repr.hide()
        self.ui.box_protein_color.hide()
        self.ui.lbl_protein_current_color.hide()
        self.ui.lbl_protein_current_color.setText("")
        # old
        self.ui.btn_protein_color_atoms.hide()
        self.ui.btn_protein_reset_atoms.hide()
        self.ui.btn_protein_show_cartoon.hide()
        self.ui.btn_protein_hide_cartoon.hide()
        self.ui.btn_protein_show_sticks.hide()
        self.ui.btn_protein_hide_sticks.hide()
        self.ui.btn_protein_show_ribbon.hide()
        self.ui.btn_protein_hide_ribbon.hide()
        # checkboxes
        self.ui.cb_protein_cartoon.hide()
        self.ui.cb_protein_sticks.hide()
        self.ui.cb_protein_ribbon.hide()
        self.ui.cb_protein_lines.hide()
        self.ui.cb_protein_spheres.hide()
        self.ui.cb_protein_dots.hide()
        self.ui.cb_protein_mesh.hide()
        self.ui.cb_protein_surface.hide()
        # toggle buttons
        self.tg_protein_cartoon.hide()
        self.tg_protein_sticks.hide()
        self.tg_protein_ribbon.hide()
        self.tg_protein_lines.hide()
        self.tg_protein_spheres.hide()
        self.tg_protein_dots.hide()
        self.tg_protein_mesh.hide()
        self.tg_protein_surface.hide()

        self.ui.btn_protein_hide_all_representations.hide()

        self.ui.lbl_info.setText("Please select a protein.")
        self.ui.lbl_info_2.hide()

        # Protein Pairs tab
        self.ui.btn_protein_pair_tree_view_expand.setEnabled(False)
        self.ui.btn_protein_pair_tree_view_collapse.setEnabled(False)
        self.ui.btn_delete_protein_pair.setEnabled(False)
        self.ui.btn_open_protein_pair_session.setEnabled(False)
        self.ui.btn_create_protein_pair_scene.setEnabled(False)
        self.ui.btn_update_protein_pair_scene.setEnabled(False)
        self.ui.btn_delete_protein_pair_scene.setEnabled(False)
        self.ui.protein_pairs_tree_view.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.ui.lbl_pymol_protein_pair_scene.setText("PyMOL Scene: No Scene Loaded")

        # Hides all ui elements for the scene modifications
        self.ui.lbl_protein_pair_current_color.hide()
        self.ui.lbl_protein_pair_current_color.setText("")
        self.ui.lbl_protein_pair_pymol_colors.hide()
        self.ui.lbl_protein_pair_color.hide()
        self.ui.lbl_protein_pair_atoms.hide()
        self.ui.lbl_protein_pair_cartoon.hide()
        self.ui.lbl_protein_pair_sticks.hide()
        self.ui.lbl_protein_pair_ribbon.hide()
        self.ui.lbl_protein_pair_lines.hide()
        self.ui.lbl_protein_pair_spheres.hide()
        self.ui.lbl_protein_pair_dots.hide()
        self.ui.lbl_protein_pair_mesh.hide()
        self.ui.lbl_protein_pair_surface.hide()
        self.ui.lbl_protein_pair_all_representations.hide()

        self.ui.box_protein_pair_color.hide()
        self.ui.btn_protein_pair_color_atoms.hide()
        self.ui.btn_protein_pair_reset_atoms.hide()
        self.ui.btn_protein_pair_show_cartoon.hide()
        self.ui.btn_protein_pair_hide_cartoon.hide()
        self.ui.btn_protein_pair_show_sticks.hide()
        self.ui.btn_protein_pair_hide_sticks.hide()
        self.ui.btn_protein_pair_show_ribbon.hide()
        self.ui.btn_protein_pair_hide_ribbon.hide()
        self.ui.btn_protein_pair_hide_all_representations.hide()
        # checkboxes
        self.ui.cb_protein_pair_cartoon.hide()
        self.ui.cb_protein_pair_sticks.hide()
        self.ui.cb_protein_pair_ribbon.hide()
        self.ui.cb_protein_pair_lines.hide()
        self.ui.cb_protein_pair_spheres.hide()
        self.ui.cb_protein_pair_dots.hide()
        self.ui.cb_protein_pair_mesh.hide()
        self.ui.cb_protein_pair_surface.hide()
        # toggle buttons
        #self.tg_protein_pair_color_atoms.hide()
        self.tg_protein_pair_cartoon.hide()
        self.tg_protein_pair_sticks.hide()
        self.tg_protein_pair_ribbon.hide()
        self.tg_protein_pair_lines.hide()
        self.tg_protein_pair_spheres.hide()
        self.tg_protein_pair_dots.hide()
        self.tg_protein_pair_mesh.hide()
        self.tg_protein_pair_surface.hide()

        self.ui.lbl_info_3.setText("Please select a protein pair.")
        self.ui.lbl_info_4.hide()

        # Extra UI elements
        self.cb_chain_color = QtWidgets.QComboBox()
        self.cb_chain_representation = QtWidgets.QComboBox()
        self.cb_chain_color_protein_pair = QtWidgets.QComboBox()
        self.cb_chain_representation_protein_pair = QtWidgets.QComboBox()
        self.line_edit_seq_name = custom_line_edit.CustomLineEdit()
        self.ui.lbl_pymol_protein_pair_scene.setWordWrap(False)

        self.build_sequence_table()
        #self.build_proteins_table()
        #self.build_protein_pairs_table()

        # Main tree views
        self.ui.seqs_list_view.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.ui.seqs_list_view.setContextMenuPolicy(Qt.CustomContextMenu)
        self.ui.proteins_tree_view.setContextMenuPolicy(Qt.CustomContextMenu)
        self.ui.protein_pairs_tree_view.setContextMenuPolicy(Qt.CustomContextMenu)
        self.setMinimumWidth(700)
        self.setMinimumHeight(900)

        # <editor-fold desc="Set icons">
        pixmapi = QtWidgets.QStyle.SP_MessageBoxQuestion
        # <editor-fold desc="Help">
        icon = self.style().standardIcon(pixmapi)
        self.ui.btn_help.setIcon(QtGui.QIcon(":/icons/help_w200.png"))
        self.ui.btn_help.setIconSize(self.ui.btn_help.icon().actualSize(QtCore.QSize(30, 30)))
        self.ui.btn_help.setText("")
        self.ui.btn_help_2.setIcon(QtGui.QIcon(":/icons/help_w200.png"))
        self.ui.btn_help_2.setIconSize(self.ui.btn_help_2.icon().actualSize(QtCore.QSize(30, 30)))
        self.ui.btn_help_2.setText("")
        self.ui.btn_help_3.setIcon(QtGui.QIcon(":/icons/help_w200.png"))
        self.ui.btn_help_3.setIconSize(self.ui.btn_help_3.icon().actualSize(QtCore.QSize(30, 30)))
        self.ui.btn_help_3.setText("")
        # </editor-fold>

        # <editor-fold desc="Sequence">
        # add
        add_sequence_icon = QtGui.QIcon(QtGui.QPixmap(":icons/note_add_w200.png"))
        add_sequence_icon.addPixmap(QtGui.QPixmap(":icons/note_add_disabled_w200.png"),
                                    mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_add_sequence.setIcon(add_sequence_icon)
        self.ui.btn_add_sequence.setText("")
        self.ui.btn_add_sequence.setIconSize(add_sequence_icon.actualSize(QtCore.QSize(30, 30)))
        # self.ui.btn_add_sequence.setIcon(
        #     QtGui.QIcon(":/icons/note_add_w200.png")
        # )
        # self.ui.btn_add_sequence.setText("")
        # self.ui.btn_add_sequence.setIconSize(self.ui.btn_add_sequence.icon().actualSize(QtCore.QSize(30, 30)))

        # import
        import_sequence_icon = QtGui.QIcon(QtGui.QPixmap(":icons/upload_file_w200.png"))
        import_sequence_icon.addPixmap(QtGui.QPixmap(":icons/upload_file_disabled_w200.png"),
                                       mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_import_seq.setIcon(import_sequence_icon)
        self.ui.btn_import_seq.setText("")
        self.ui.btn_import_seq.setIconSize(import_sequence_icon.actualSize(QtCore.QSize(30, 30)))
        # self.ui.btn_import_seq.setIcon(
        #     QtGui.QIcon(":/icons/upload_file_w200.png")
        # )
        # self.ui.btn_import_seq.setText("")
        # self.ui.btn_import_seq.setIconSize(self.ui.btn_import_seq.icon().actualSize(QtCore.QSize(30, 30)))

        # save
        save_seq_icon = QtGui.QIcon(QtGui.QPixmap(":icons/file_save_w200.png"))
        save_seq_icon.addPixmap(QtGui.QPixmap(":icons/file_save_disabled_w200.png"), mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_save_sequence.setIcon(save_seq_icon)
        self.ui.btn_save_sequence.setText("")
        self.ui.btn_save_sequence.setIconSize(save_seq_icon.actualSize(QtCore.QSize(30, 30)))
        # save_icon = QtGui.QIcon(":/icons/file_save_w200.png")
        # save_icon.addPixmap(QtGui.QPixmap(":icons/file_save_disabled_w200.png"), mode=QtGui.QIcon.Mode.Disabled)
        # self.ui.btn_save_sequence.setIcon(save_icon)
        # self.ui.btn_save_sequence.setText("")
        # self.ui.btn_save_sequence.setIconSize(self.ui.btn_save_sequence.icon().actualSize(QtCore.QSize(30, 30)))
        # self.ui.btn_save_sequence.setIcon(QtGui.QIcon(":/icons/file_save_w200.png"))

        # delete
        delete_seq_icon = QtGui.QIcon(QtGui.QPixmap(":icons/scan_delete_w200.png"))
        delete_seq_icon.addPixmap(QtGui.QPixmap(":icons/scan_delete_disabled_w200.png"), mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_delete_sequence.setIcon(delete_seq_icon)
        self.ui.btn_delete_sequence.setText("")
        self.ui.btn_delete_sequence.setIconSize(delete_seq_icon.actualSize(QtCore.QSize(30, 30)))
        # self.ui.btn_delete_sequence.setIcon(
        #     QtGui.QIcon(":/icons/scan_deletew200.png"))
        # self.ui.btn_delete_sequence.setText("")
        # self.ui.btn_delete_sequence.setIconSize(self.ui.btn_delete_sequence.icon().actualSize(QtCore.QSize(30, 30)))
        # </editor-fold>

        # <editor-fold desc="Proteins Tab">
        # expand all
        expand_all_protein_icon = QtGui.QIcon(QtGui.QPixmap(":icons/expand_all_w200.png"))
        self.ui.btn_protein_tree_view_expand.setIcon(expand_all_protein_icon)
        self.ui.btn_protein_tree_view_expand.setText("")
        self.ui.btn_protein_tree_view_expand.setIconSize(expand_all_protein_icon.actualSize(QtCore.QSize(14, 14)))

        # collapse all
        collapse_all_protein_icon = QtGui.QIcon(QtGui.QPixmap(":icons/collapse_all_w200.png"))
        self.ui.btn_protein_tree_view_collapse.setIcon(collapse_all_protein_icon)
        self.ui.btn_protein_tree_view_collapse.setText("")
        self.ui.btn_protein_tree_view_collapse.setIconSize(collapse_all_protein_icon.actualSize(QtCore.QSize(18, 18)))

        # import
        import_protein_icon = QtGui.QIcon(QtGui.QPixmap(":icons/upload_file_w200.png"))
        import_protein_icon.addPixmap(QtGui.QPixmap(":icons/upload_file_disabled_w200.png"),
                                      mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_import_protein.setIcon(import_protein_icon)
        self.ui.btn_import_protein.setText("")
        self.ui.btn_import_protein.setIconSize(import_protein_icon.actualSize(QtCore.QSize(30, 30)))

        # save
        save_protein_icon = QtGui.QIcon(QtGui.QPixmap(":icons/file_save_w200.png"))
        save_protein_icon.addPixmap(QtGui.QPixmap(":icons/file_save_disabled_w200.png"), mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_save_protein.setIcon(save_protein_icon)
        self.ui.btn_save_protein.setText("")
        self.ui.btn_save_protein.setIconSize(save_protein_icon.actualSize(QtCore.QSize(30, 30)))

        # delete
        delete_protein_icon = QtGui.QIcon(QtGui.QPixmap(":icons/scan_delete_w200.png"))
        delete_protein_icon.addPixmap(QtGui.QPixmap(":icons/scan_delete_disabled_w200.png"), mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_delete_protein.setIcon(delete_protein_icon)
        self.ui.btn_delete_protein.setText("")
        self.ui.btn_delete_protein.setIconSize(delete_protein_icon.actualSize(QtCore.QSize(30, 30)))

        # </editor-fold>

        # <editor-fold desc="Protein Session">
        # open
        open_protein_session_icon = QtGui.QIcon(QtGui.QPixmap(":/icons/open_in_new_w200.png"))
        open_protein_session_icon.addPixmap(QtGui.QPixmap(":/icons/open_in_new_disabled_w200.png"),
                                            mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_open_protein_session.setIcon(open_protein_session_icon)
        self.ui.btn_open_protein_session.setText("")
        self.ui.btn_open_protein_session.setIconSize(open_protein_session_icon.actualSize(QtCore.QSize(30, 30)))

        # create
        create_protein_session_icon = QtGui.QIcon(QtGui.QPixmap(":/icons/add_circle_w200.png"))
        create_protein_session_icon.addPixmap(QtGui.QPixmap(":/icons/add_circle_disabled_w200.png"),
                                            mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_create_protein_scene.setIcon(create_protein_session_icon)
        self.ui.btn_create_protein_scene.setText("")
        self.ui.btn_create_protein_scene.setIconSize(create_protein_session_icon.actualSize(QtCore.QSize(30, 30)))

        # refresh
        update_protein_session_icon = QtGui.QIcon(QtGui.QPixmap(":/icons/change_circle_w200.png"))
        update_protein_session_icon.addPixmap(QtGui.QPixmap(":/icons/change_circle_disabled_w200.png"),
                                              mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_update_protein_scene.setIcon(update_protein_session_icon)
        self.ui.btn_update_protein_scene.setText("")
        self.ui.btn_update_protein_scene.setIconSize(update_protein_session_icon.actualSize(QtCore.QSize(30, 30)))

        # delete scene
        delete_protein_session_icon = QtGui.QIcon(QtGui.QPixmap(":icons/cancel_w200.png"))
        delete_protein_session_icon.addPixmap(QtGui.QPixmap(":icons/cancel_disabled_w200.png"),
                                      mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_delete_protein_scene.setIcon(delete_protein_session_icon)
        self.ui.btn_delete_protein_scene.setText("")
        self.ui.btn_delete_protein_scene.setIconSize(delete_protein_session_icon.actualSize(QtCore.QSize(30, 30)))
        # </editor-fold>

        # <editor-fold desc="Protein Pairs Tab">
        # expand all
        expand_all_protein_pair_icon = QtGui.QIcon(QtGui.QPixmap(":icons/expand_all_w200.png"))
        self.ui.btn_protein_pair_tree_view_expand.setIcon(expand_all_protein_pair_icon)
        self.ui.btn_protein_pair_tree_view_expand.setText("")
        self.ui.btn_protein_pair_tree_view_expand.setIconSize(
            expand_all_protein_pair_icon.actualSize(QtCore.QSize(14, 14))
        )

        # collapse all
        collapse_all_protein_pair_icon = QtGui.QIcon(QtGui.QPixmap(":icons/collapse_all_w200.png"))
        self.ui.btn_protein_pair_tree_view_collapse.setIcon(collapse_all_protein_pair_icon)
        self.ui.btn_protein_pair_tree_view_collapse.setText("")
        self.ui.btn_protein_pair_tree_view_collapse.setIconSize(
            collapse_all_protein_pair_icon.actualSize(QtCore.QSize(18, 18))
        )
        # </editor-fold>

        # <editor-fold desc="Protein Pair Session">
        # open
        open_protein_pair_session_icon = QtGui.QIcon(QtGui.QPixmap(":/icons/open_in_new_w200.png"))
        open_protein_pair_session_icon.addPixmap(QtGui.QPixmap(":/icons/open_in_new_disabled_w200.png"),
                                            mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_open_protein_pair_session.setIcon(open_protein_pair_session_icon)
        self.ui.btn_open_protein_pair_session.setText("")
        self.ui.btn_open_protein_pair_session.setIconSize(open_protein_pair_session_icon.actualSize(QtCore.QSize(30, 30)))

        # create
        create_protein_pair_session_icon = QtGui.QIcon(QtGui.QPixmap(":/icons/add_circle_w200.png"))
        create_protein_pair_session_icon.addPixmap(QtGui.QPixmap(":/icons/add_circle_disabled_w200.png"),
                                              mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_create_protein_pair_scene.setIcon(create_protein_pair_session_icon)
        self.ui.btn_create_protein_pair_scene.setText("")
        self.ui.btn_create_protein_pair_scene.setIconSize(create_protein_pair_session_icon.actualSize(QtCore.QSize(30, 30)))

        # refresh
        update_protein_pair_session_icon = QtGui.QIcon(QtGui.QPixmap(":/icons/change_circle_w200.png"))
        update_protein_pair_session_icon.addPixmap(QtGui.QPixmap(":/icons/change_circle_disabled_w200.png"),
                                              mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_update_protein_pair_scene.setIcon(update_protein_pair_session_icon)
        self.ui.btn_update_protein_pair_scene.setText("")
        self.ui.btn_update_protein_pair_scene.setIconSize(update_protein_pair_session_icon.actualSize(QtCore.QSize(30, 30)))

        # delete
        delete_protein_pair_icon = QtGui.QIcon(QtGui.QPixmap(":icons/scan_delete_w200.png"))
        delete_protein_pair_icon.addPixmap(QtGui.QPixmap(":icons/scan_delete_disabled_w200.png"), mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_delete_protein_pair.setIcon(delete_protein_pair_icon)
        self.ui.btn_delete_protein_pair.setText("")
        self.ui.btn_delete_protein_pair.setIconSize(delete_protein_pair_icon.actualSize(QtCore.QSize(40, 40)))

        # delete scene
        delete_protein_pair_session_icon = QtGui.QIcon(QtGui.QPixmap(":icons/cancel_w200.png"))
        delete_protein_pair_session_icon.addPixmap(QtGui.QPixmap(":icons/cancel_disabled_w200.png"
                                                                 ""),
                                      mode=QtGui.QIcon.Mode.Disabled)
        self.ui.btn_delete_protein_pair_scene.setIcon(delete_protein_pair_session_icon)
        self.ui.btn_delete_protein_pair_scene.setText("")
        self.ui.btn_delete_protein_pair_scene.setIconSize(delete_protein_pair_session_icon.actualSize(QtCore.QSize(30, 30)))
        # </editor-fold>
        # </editor-fold>

        self._create_all_tooltips()
        pixmap = QtGui.QPixmap(str(constants.PLUGIN_LOGO_WITH_FONT_FILEPATH))
        # Resize the pixmap
        pixmap = QtGui.QPixmap(r"C:\ProgramData\pyssa\mambaforge_pyssa\pyssa-mamba-env\Lib\site-packages\pymol\pymol_path\data\startup\PySSA\assets\images\logo_type_2.tiff")
        scaled_pixmap = pixmap.scaled(700, 700, aspectRatioMode=Qt.KeepAspectRatio, transformMode=Qt.SmoothTransformation)
        # Set the scaled pixmap to the QLabel
        self.ui.lbl_logo.setPixmap(scaled_pixmap)
        self.ui.lbl_logo.setAlignment(Qt.AlignCenter)
        styles.set_stylesheet_homepage(self)
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        self.setWindowTitle("PySSA")
        constants.PYSSA_LOGGER.info(f"PySSA started with version {constants.VERSION_NUMBER}.")
        constants.PYSSA_LOGGER.info("Successful initialization of basic UI.")

    def disable_menu_bar_without_exit_application(self):
        self.ui.menuProject.setEnabled(True)
        self.ui.action_new_project.setEnabled(False)
        self.ui.action_open_project.setEnabled(False)
        self.ui.action_close_project.setEnabled(False)
        self.ui.action_use_project.setEnabled(False)
        self.ui.action_delete_project.setEnabled(False)
        self.ui.action_export_project.setEnabled(False)
        self.ui.action_import_project.setEnabled(False)
        self.ui.menuPrediction.setEnabled(False)
        self.ui.menuAnalysis.setEnabled(False)
        self.ui.menuResults.setEnabled(False)
        self.ui.menuImage.setEnabled(False)
        self.ui.menuHotspots.setEnabled(False)
        self.ui.menuSettings.setEnabled(False)
        self.ui.menuAbout.setEnabled(False)

    def disable_tab_widget(self):
        self.ui.project_tab_widget.setEnabled(False)

    def disable_job_panels(self):
        self.ui.frame_job_overview.setEnabled(False)
        self.ui.frame_job_notification.setEnabled(False)
        self.btn_open_job_overview.setEnabled(False)
        self.btn_open_job_notification.setEnabled(False)

    def enable_job_panels(self):
        self.ui.frame_job_overview.setEnabled(True)
        self.ui.frame_job_notification.setEnabled(True)
        self.btn_open_job_overview.setEnabled(True)
        self.btn_open_job_notification.setEnabled(True)

    def build_sequence_table(self):
        #self.line_edit_seq_name = custom_line_edit.CustomLineEdit()
        self.ui.seqs_table_widget.verticalHeader().setVisible(False)
        self.ui.seqs_table_widget.setColumnCount(2)
        self.ui.seqs_table_widget.setHorizontalHeaderLabels(["Name", "Value"])

    # def build_proteins_table(self):
    #     self.cb_chain_color = QtWidgets.QComboBox()
    #     self.cb_chain_representation = QtWidgets.QComboBox()
    #
    #     self.ui.proteins_table_widget.verticalHeader().setVisible(False)
    #     self.ui.proteins_table_widget.setColumnCount(2)
    #     self.ui.proteins_table_widget.setHorizontalHeaderLabels(["Name", "Value"])
    #     gui_utils.fill_combo_box(self.cb_chain_color, constants.PYMOL_COLORS)
    #     self.cb_chain_color.adjustSize()
    #     gui_utils.fill_combo_box(self.cb_chain_representation, constants.PYMOL_REPRESENTATIONS)
    #     self.cb_chain_representation.adjustSize()
    #
    # def build_protein_pairs_table(self):
    #     self.cb_chain_color_protein_pair = QtWidgets.QComboBox()
    #     self.cb_chain_representation_protein_pair = QtWidgets.QComboBox()
    #
    #     self.ui.protein_pairs_table_widget.setColumnCount(2)
    #     self.ui.protein_pairs_table_widget.verticalHeader().setVisible(False)
    #     self.ui.protein_pairs_table_widget.setHorizontalHeaderLabels(["Name", "Value"])
    #     gui_utils.fill_combo_box(self.cb_chain_color_protein_pair, constants.PYMOL_COLORS)
    #     self.cb_chain_color_protein_pair.adjustSize()
    #     gui_utils.fill_combo_box(self.cb_chain_representation_protein_pair, constants.PYMOL_REPRESENTATIONS)
    #     self.cb_chain_representation_protein_pair.adjustSize()

    def setup_sequences_table(self, row_count):
        #self.line_edit_seq_name.setStyleSheet("QLineEdit { background-color: white; border-radius: 0; }")
        self.ui.seqs_table_widget.setRowCount(row_count)
        #self.ui.seqs_table_widget.setCellWidget(0, 1, self.line_edit_seq_name)

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
        # Sequence Tab
        # self.ui.seqs_list_view.setToolTip("A list of all sequences in the project")
        self.ui.seqs_table_widget.setToolTip("A table with additional information about the selected sequence")
        self.ui.btn_import_seq.setToolTip("Import an existing .fasta file")
        self.ui.btn_add_sequence.setToolTip("Add a sequence by pasting the sequence string")
        self.ui.btn_save_sequence.setToolTip("Export the selected sequence as .fasta file")
        self.ui.btn_delete_sequence.setToolTip("Delete the selected sequence from the project")

        # Proteins Tab
        # self.ui.proteins_tree_view.setToolTip("A tree of all proteins in the project")
        # self.ui.proteins_table_widget.setToolTip(
        #     "A table with changeable PyMOL parameters for the currently active session"
        # )
        self.ui.btn_protein_tree_view_expand.setToolTip("Expand All")
        self.ui.btn_protein_tree_view_collapse.setToolTip("Collapse All")
        self.ui.btn_import_protein.setToolTip("Import an existing .pdb file")
        self.ui.btn_save_protein.setToolTip("Export the selected protein as .pdb file")
        self.ui.btn_delete_protein.setToolTip("Delete the selected protein from the project")
        self.ui.btn_open_protein_session.setToolTip("Open protein PyMOL session")
        self.ui.btn_create_protein_scene.setToolTip("Create a new PyMOL scene")
        self.ui.btn_update_protein_scene.setToolTip("Update the current scene in PyMOL")
        self.ui.btn_delete_protein_scene.setToolTip("Delete the current scene in PyMOL")

        # Protein Pairs Tab
        # self.ui.protein_pairs_tree_view.setToolTip("A tree of all protein pairs in the project")
        # self.ui.protein_pairs_table_widget.setToolTip(
        #     "A table with changeable PyMOL parameters for the currently active session"
        # )
        self.ui.btn_protein_pair_tree_view_expand.setToolTip("Expand All")
        self.ui.btn_protein_pair_tree_view_collapse.setToolTip("Collapse All")
        self.ui.btn_delete_protein_pair.setToolTip("Delete the selected protein pair from the project")
        self.ui.btn_open_protein_pair_session.setToolTip("Open protein pair PyMOL session")
        self.ui.btn_create_protein_pair_scene.setToolTip("Create a new PyMOL scene")
        self.ui.btn_update_protein_pair_scene.setToolTip("Update the current scene in PyMOL")
        self.ui.btn_delete_protein_pair_scene.setToolTip("Delete the current scene in PyMOL")

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
