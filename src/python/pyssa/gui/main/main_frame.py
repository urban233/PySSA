#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
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
"""Module for the main view of the PySSA plugin."""
import logging

from PyQt6 import QtWidgets
from PyQt6 import QtGui
from PyQt6 import QtCore
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import QAbstractItemView
# from pyssa.gui.custom_widgets import custom_line_edit
# from pyssa.gui.custom_widgets import toggle_button
from pyssa.gui.custom_widgets import color_grid
from pyssa.gui.main.forms.auto import auto_main_view
from pyssa.model.preference import model_definitions
from pyssa.model.util.gui_style import styles_utils, icons

from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)

__docformat__ = "google"


class MainFrame(QtWidgets.QMainWindow):
  """This class represents the main view of the application."""

  dialogClosed = QtCore.pyqtSignal(tuple)
  """A signal indicating that the dialog is closed."""

  def __init__(self) -> None:
    """Constructor."""
    super().__init__()
    # build ui object
    self.ui = auto_main_view.Ui_MainWindow()
    self.ui.setupUi(self)
    self.progress_bar = QtWidgets.QProgressBar()
    self.status_bar = QtWidgets.QStatusBar()
    self.status_bar.addWidget(self.progress_bar)
    self.progress_bar.hide()
    self.add_custom_widgets()
    self.initialize_ui()
    self.add_custom_job_panels()

  def closeEvent(self, event) -> None:  # noqa: ANN001
    """Overrides the closeEvent of the QMainWindow class.

    Args:
        event: The event object representing the close event.
    """
    # <editor-fold desc="Checks">
    if event is None:
      logger.error("event is None.")
      raise exception.IllegalArgumentError("event is None.")

    # </editor-fold>
    # Emit the custom signal when the window is closed
    self.dialogClosed.emit(("", event))

  def add_custom_job_panels(self) -> None:
    """Add custom job panels.

    Creates and configures the GUI elements for the job panels, including buttons, labels, and icons.
    """
    # job panel related gui elements
    self.btn_open_job_overview = QtWidgets.QPushButton()
    self.lbl_job_overview = QtWidgets.QLabel("No running jobs.")
    self.btn_open_job_notification = QtWidgets.QPushButton()
    self.lbl_job_notification = QtWidgets.QLabel("No new notifications.")
    icons.set_icon(self.btn_open_job_overview, model_definitions.IconsEnum.JOBS)
    icons.set_icon(self.btn_open_job_notification, model_definitions.IconsEnum.NOTIFY)
    # <editor-fold desc="Custom stylesheets">
    self.btn_open_job_overview.setStyleSheet(
        """
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
        """
    )
    self.lbl_job_overview.setStyleSheet("""padding-top: 30px;""")
    self.btn_open_job_notification.setStyleSheet(
        """
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
        """
    )
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
    # </editor-fold>
    self.ui.job_overview_layout.insertWidget(
        0, self.lbl_job_overview
    )  # After inserting the widget count is 2
    self.ui.job_overview_layout.setAlignment(
        self.lbl_job_overview, QtCore.Qt.AlignmentFlag.AlignHCenter
    )
    self.ui.job_notification_layout.insertWidget(
        0, self.lbl_job_notification
    )  # After inserting the widget count is 2
    self.ui.job_notification_layout.setAlignment(
        self.lbl_job_notification, QtCore.Qt.AlignmentFlag.AlignHCenter
    )

  def add_custom_widgets(self) -> None:
    """Add custom widgets to the UI for controlling protein and protein pair representations."""
    self.color_grid_proteins = color_grid.ColorGrid()
    self.ui.verticalLayout_6.insertWidget(2, self.color_grid_proteins)
    self.color_grid_protein_pairs = color_grid.ColorGrid()
    self.ui.verticalLayout_10.insertWidget(2, self.color_grid_protein_pairs)

  def initialize_ui(self) -> None:
    """Initialize the UI elements."""
    self.ui.lbl_project_name.hide()
    self.ui.lbl_session_name.hide()
    self.ui.project_tab_widget.hide()

    # fixme: the tutorial action of the menubar IS HIDDEN!!!
    self.ui.action_tutorials.setVisible(False)
    self.ui.action_abort_prediction.setVisible(False)
    self.ui.action_preview_image.setVisible(False)

    self.ui.frame_job_overview.hide()
    self.ui.frame_job_notification.hide()

    # Project and Session name
    self.ui.lbl_project_name.setText("Project Name: No Project Opened")
    self.ui.lbl_session_name.setText("Session Name: No Session Loaded")

    # Sequences tab
    self.ui.seqs_list_view.setEditTriggers(
        QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
    )

    # Proteins tab
    self.ui.proteins_tree_view.setEditTriggers(
        QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
    )
    self.ui.lbl_pymol_protein_scene.setText("PyMOL Scene: No Scene Loaded")

    self.ui.box_protein_color.hide()
    self.ui.lbl_protein_current_color.setText("")

    self.ui.lbl_info.setText("Please select a protein.")
    self.ui.lbl_info_2.hide()

    # Protein Pairs tab
    self.ui.protein_pairs_tree_view.setEditTriggers(
        QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
    )
    self.ui.lbl_pymol_protein_pair_scene.setText("PyMOL Scene: No Scene Loaded")

    # Hides all ui elements for the scene modifications
    self.ui.lbl_protein_pair_current_color.setText("")
    self.ui.box_protein_pair_color.hide()

    self.ui.lbl_info_3.setText("Please select a protein pair.")
    self.ui.lbl_info_4.hide()

    # Extra UI elements
    #self.line_edit_seq_name = custom_line_edit.CustomLineEdit()
    self.ui.lbl_pymol_protein_pair_scene.setWordWrap(False)

    self.build_sequence_table()
    # self.build_proteins_table()
    # self.build_protein_pairs_table()

    # Main tree views
    self.ui.seqs_list_view.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
    self.ui.seqs_list_view.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
    self.ui.proteins_tree_view.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
    self.ui.protein_pairs_tree_view.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
    self.setMinimumWidth(700)
    self.setMinimumHeight(900)

    self._create_all_tooltips()
    self.set_icons()
    pixmap = QtGui.QPixmap(str(model_definitions.ModelDefinitions.PLUGIN_LOGO_WITH_CAPTION_FILEPATH))
    scaled_pixmap = pixmap.scaled(
        450,
        450,
        aspectRatioMode=Qt.AspectRatioMode.KeepAspectRatio,
        transformMode=Qt.TransformationMode.SmoothTransformation,
    )
    # Set the scaled pixmap to the QLabel
    self.ui.lbl_logo.setPixmap(scaled_pixmap)
    self.ui.lbl_logo.setAlignment(Qt.AlignmentFlag.AlignCenter)
    self.setWindowIcon(QtGui.QIcon(model_definitions.ModelDefinitions.PLUGIN_LOGO_FILEPATH))
    self.setWindowTitle("PySSA")
    # constants.PYSSA_LOGGER.info(
    #     f"PySSA started with version {constants.VERSION_NUMBER}."
    # )
    # constants.PYSSA_LOGGER.info("Successful initialization of basic UI.")

  def set_icons(self) -> None:
    """Sets all icons for the main frame."""
    # <editor-fold desc="Help">
    icons.set_icon(self.ui.btn_help, model_definitions.IconsEnum.HELP)
    icons.set_icon(self.ui.btn_help_2, model_definitions.IconsEnum.HELP)
    icons.set_icon(self.ui.btn_help_3, model_definitions.IconsEnum.HELP)
    # </editor-fold>
    # <editor-fold desc="Sequence tab">
    # import
    icons.set_icon(
      self.ui.btn_import_seq,
      model_definitions.IconsEnum.IMPORT_SEQUENCE,
      model_definitions.IconsEnum.IMPORT_SEQUENCE_DISABLED,
    )
    # add
    icons.set_icon(
      self.ui.btn_add_sequence,
      model_definitions.IconsEnum.ADD_SEQUENCE,
      model_definitions.IconsEnum.ADD_SEQUENCE_DISABLED,
    )
    # save sequence
    icons.set_icon(
      self.ui.btn_save_sequence,
      model_definitions.IconsEnum.SAVE_SEQUENCE,
      model_definitions.IconsEnum.SAVE_SEQUENCE_DISABLED,
    )

    # delete sequence
    icons.set_icon(
      self.ui.btn_delete_sequence,
      model_definitions.IconsEnum.DELETE_SEQUENCE,
      model_definitions.IconsEnum.DELETE_SEQUENCE_DISABLED,
    )
    # </editor-fold>
    # <editor-fold desc="Proteins tab">
    # expand all
    icons.set_icon(
      self.ui.btn_protein_tree_view_expand,
      model_definitions.IconsEnum.EXPAND_ALL,
      size=QtCore.QSize(20, 20)
    )
    # collapse all
    icons.set_icon(
      self.ui.btn_protein_tree_view_collapse,
      model_definitions.IconsEnum.COLLAPSE_ALL,
      size=QtCore.QSize(20, 20)
    )
    # import
    icons.set_icon(
      self.ui.btn_import_protein,
      model_definitions.IconsEnum.IMPORT_PROTEIN,
      model_definitions.IconsEnum.IMPORT_PROTEIN_DISABLED,
    )
    # save
    icons.set_icon(
      self.ui.btn_save_protein,
      model_definitions.IconsEnum.SAVE_PROTEIN,
      model_definitions.IconsEnum.SAVE_PROTEIN_DISABLED,
    )
    # delete
    icons.set_icon(
      self.ui.btn_delete_protein,
      model_definitions.IconsEnum.DELETE_PROTEIN,
      model_definitions.IconsEnum.DELETE_PROTEIN_DISABLED,
    )
    # open
    icons.set_icon(
      self.ui.btn_open_protein_session,
      model_definitions.IconsEnum.OPEN_SESSION,
      model_definitions.IconsEnum.OPEN_SESSION_DISABLED,
    )
    # create
    icons.set_icon(
      self.ui.btn_create_protein_scene,
      model_definitions.IconsEnum.CREATE_SESSION_SCENE,
      model_definitions.IconsEnum.CREATE_SESSION_SCENE_DISABLED,
    )
    # update
    icons.set_icon(
      self.ui.btn_update_protein_scene,
      model_definitions.IconsEnum.UPDATE_SESSION_SCENE,
      model_definitions.IconsEnum.UPDATE_SESSION_SCENE_DISABLED,
    )
    # delete scene
    icons.set_icon(
      self.ui.btn_delete_protein_scene,
      model_definitions.IconsEnum.DELETE_SESSION_SCENE,
      model_definitions.IconsEnum.DELETE_SESSION_SCENE_DISABLED,
    )
    # </editor-fold>
    # <editor-fold desc="Protein-Protein complex tab">
    # expand all
    icons.set_icon(
      self.ui.btn_protein_pair_tree_view_expand,
      model_definitions.IconsEnum.EXPAND_ALL,
      size=QtCore.QSize(20, 20)
    )
    # collapse all
    icons.set_icon(
      self.ui.btn_protein_pair_tree_view_collapse,
      model_definitions.IconsEnum.COLLAPSE_ALL,
      size=QtCore.QSize(20, 20)
    )
    # delete
    icons.set_icon(
      self.ui.btn_delete_protein_pair,
      model_definitions.IconsEnum.DELETE_PROTEIN_PAIR,
      model_definitions.IconsEnum.DELETE_PROTEIN_PAIR_DISABLED,
    )
    # open
    icons.set_icon(
      self.ui.btn_open_protein_pair_session,
      model_definitions.IconsEnum.OPEN_SESSION,
      model_definitions.IconsEnum.OPEN_SESSION_DISABLED,
    )
    # create
    icons.set_icon(
      self.ui.btn_create_protein_pair_scene,
      model_definitions.IconsEnum.CREATE_SESSION_SCENE,
      model_definitions.IconsEnum.CREATE_SESSION_SCENE_DISABLED,
    )
    # update
    icons.set_icon(
      self.ui.btn_update_protein_pair_scene,
      model_definitions.IconsEnum.UPDATE_SESSION_SCENE,
      model_definitions.IconsEnum.UPDATE_SESSION_SCENE_DISABLED,
    )
    # delete scene
    icons.set_icon(
      self.ui.btn_delete_protein_pair_scene,
      model_definitions.IconsEnum.DELETE_SESSION_SCENE,
      model_definitions.IconsEnum.DELETE_SESSION_SCENE_DISABLED,
    )
    # </editor-fold>

  def disable_menu_bar_without_exit_application(self) -> None:
    """Disables the menu entries but not 'Exit Application'."""
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

  def disable_job_panels(self) -> None:
    """Disables the job overview panel, job notification panel, and corresponding buttons."""
    self.ui.frame_job_overview.setEnabled(False)
    self.ui.frame_job_notification.setEnabled(False)
    self.btn_open_job_overview.setEnabled(False)
    self.btn_open_job_notification.setEnabled(False)

  def enable_job_panels(self) -> None:
    """Enables the job overview panel, job notification panel,and the corresponding buttons."""
    self.ui.frame_job_overview.setEnabled(True)
    self.ui.frame_job_notification.setEnabled(True)
    self.btn_open_job_overview.setEnabled(True)
    self.btn_open_job_notification.setEnabled(True)

  def build_sequence_table(self) -> None:
    """Builds the sequence table.

    This method initializes and populates the sequence table widget with column labels.
    """
    # self.line_edit_seq_name = custom_line_edit.CustomLineEdit()
    self.ui.seqs_table_widget.verticalHeader().setVisible(False)
    self.ui.seqs_table_widget.setColumnCount(2)
    self.ui.seqs_table_widget.setHorizontalHeaderLabels(["Name", "Value"])

  def setup_sequences_table(self, row_count: int) -> None:
    """Sets up the sequences table widget by setting the desired number of rows using the given row_count.

    Args:
        row_count: The number of rows to set for the sequences table widget.

    Raises:
        exception.IllegalArgumentError: If `row_count` is either None or has a value less than 0.
    """
    # <editor-fold desc="Checks">
    if row_count is None or row_count < 0:
      logger.error("row_count is either None or has a value less than 0.")
      raise exception.IllegalArgumentError(
          "row_count is either None or has a value less than 0."
      )

    # </editor-fold>

    # self.line_edit_seq_name.setStyleSheet("QLineEdit { background-color: white; border-radius: 0; }")
    self.ui.seqs_table_widget.setRowCount(row_count)
    # self.ui.seqs_table_widget.setCellWidget(0, 1, self.line_edit_seq_name)

  def _create_all_tooltips(self) -> None:
    """Creates all tooltips for the gui elements."""
    # <editor-fold desc="Sequence tab">
    self.ui.seqs_list_view.setToolTip("Import or add sequences to do structure predictions")  # fixme: is this a useful tooltip?
    self.ui.seqs_table_widget.setToolTip(
        "A table with additional information about the selected sequence"
    )
    self.ui.btn_import_seq.setToolTip("Import an existing .fasta file")
    self.ui.btn_add_sequence.setToolTip(
        "Add a sequence by pasting the sequence string"
    )
    self.ui.btn_save_sequence.setToolTip(
        "Export the selected sequence as .fasta file"
    )
    self.ui.btn_delete_sequence.setToolTip(
        "Delete the selected sequence from the project"
    )
    # </editor-fold>
    # <editor-fold desc="Proteins tab">
    self.ui.proteins_tree_view.setToolTip("Get proteins by running a structure prediction or importing an existing pdb file")  # fixme: is this a useful tooltip?
    self.ui.btn_protein_tree_view_expand.setToolTip("Expand All")
    self.ui.btn_protein_tree_view_collapse.setToolTip("Collapse All")
    self.ui.btn_import_protein.setToolTip("Import an existing .pdb file")
    self.ui.btn_save_protein.setToolTip(
        "Export the selected protein as .pdb file"
    )
    self.ui.btn_delete_protein.setToolTip(
        "Delete the selected protein from the project"
    )
    self.ui.btn_open_protein_session.setToolTip("Open protein PyMOL session")
    self.ui.btn_create_protein_scene.setToolTip("Create a new PyMOL scene")
    self.ui.btn_update_protein_scene.setToolTip(
        "Update the current scene in PyMOL"
    )
    self.ui.btn_delete_protein_scene.setToolTip(
        "Delete the current scene in PyMOL"
    )
    # </editor-fold>
    # <editor-fold desc="Protein pairs tab">
    self.ui.protein_pairs_tree_view.setToolTip("Get protein pairs by running a distance analysis")  # fixme: is this a useful tooltip?
    self.ui.btn_protein_pair_tree_view_expand.setToolTip("Expand All")
    self.ui.btn_protein_pair_tree_view_collapse.setToolTip("Collapse All")
    self.ui.btn_delete_protein_pair.setToolTip(
        "Delete the selected protein pair from the project"
    )
    self.ui.btn_open_protein_pair_session.setToolTip(
        "Open protein pair PyMOL session"
    )
    self.ui.btn_create_protein_pair_scene.setToolTip("Create a new PyMOL scene")
    self.ui.btn_update_protein_pair_scene.setToolTip(
        "Update the current scene in PyMOL"
    )
    self.ui.btn_delete_protein_pair_scene.setToolTip(
        "Delete the current scene in PyMOL"
    )
    # </editor-fold>
