import logging
import os
import pathlib
import shutil
import subprocess
import platform
from io import BytesIO

import pygetwindow
import requests
import zmq
import pymol
from Bio.Seq import Seq
from pymol import cmd
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from PyQt5 import QtGui
from Bio import SeqRecord
from Bio import SeqIO
from xml import sax

from pyssa.gui.ui.custom_dialogs import custom_message_box
from pyssa.internal.thread.async_pyssa import util_async
from pyssa.controller import results_view_controller, rename_protein_view_controller, use_project_view_controller, \
    pymol_session_manager, hotspots_protein_regions_view_controller, predict_multimer_view_controller, \
    add_sequence_view_controller, add_scene_view_controller, add_protein_view_controller, settings_view_controller
from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.gui.ui.styles import styles
from pyssa.gui.ui.views import predict_monomer_view, delete_project_view, rename_protein_view
from pyssa.gui.ui.dialogs import dialog_startup, dialog_settings_global, dialog_tutorial_videos, dialog_about, \
    dialog_rename_protein, dialog_help
from pyssa.internal.data_structures import project, settings, protein, protein_pair, chain, selection
from pyssa.internal.data_structures.data_classes import prediction_protein_info, database_operation
from pyssa.internal.portal import graphic_operations, pymol_io
from pyssa.gui.ui.dialogs import dialog_settings_global, dialog_tutorial_videos, dialog_about
from pyssa.internal.data_structures import project, settings, protein, protein_pair
from pyssa.internal.data_structures.data_classes import database_operation
from pyssa.internal.portal import graphic_operations
from pyssa.internal.thread import tasks, task_workers, database_thread
from pyssa.io_pyssa import safeguard, filesystem_io
from pyssa.logging_pyssa import log_handlers
from pyssa.presenter import main_presenter_async
from pyssa.util import constants, enums, exit_codes, gui_utils, tools
from pyssa.gui.ui.views import main_view
from pyssa.model import application_model
from pyssa.controller import interface_manager, distance_analysis_view_controller, predict_monomer_view_controller, \
    delete_project_view_controller, create_project_view_controller, open_project_view_controller, database_manager
from pyssa.util import globals

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class MainViewController:
    """Class for main presenter of the pyssa plugin."""

    """
    The main view of the pyssa plugin.
    """
    _view: "main_view.MainView"

    """
    
    """
    _app_model: "application_model.ApplicationModel"

    """
    The path to the active workspace.
    """
    _workspace_path: pathlib.Path

    """
    The application settings object.
    """
    _application_settings: "settings.Settings"

    """
    A watcher for the project state.
    """
    _interface_manager: "interface_manager.InterfaceManager"

    """
    A manager for the pymol session.
    """
    _pymol_session_manager: "pymol_session_manager.PymolSessionManager"

    """
    The active task of the application.
    """
    _active_task: tasks.Task

    """
    a thread for database related processes
    """
    _database_thread: "database_thread.DatabaseThread"

    def __init__(self,
                 the_interface_manager: "interface_manager.InterfaceManager",
                 the_pymol_session_manager: "pymol_session_manager.PymolSessionManager") -> None:
        """Constructor.

        Args:
            a_view: the main view of the pyssa plugin.

        Raises:
            IllegalArgumentException: If an argument is illegal.
        """
        # <editor-fold desc="Checks">
        safeguard.Safeguard.check_if_value_is_not_none(the_interface_manager, logger)

        # </editor-fold>

        self._view: "main_view.MainView" = the_interface_manager.get_main_view()
        self._interface_manager: "interface_manager.InterfaceManager" = the_interface_manager
        self._pymol_session_manager = the_pymol_session_manager
        self._database_manager = database_manager.DatabaseManager("")
        self._database_manager.set_application_settings(self._interface_manager.get_application_settings())
        self._database_thread: "database_thread.DatabaseThread" = database_thread.DatabaseThread("")

        self._app_model: "application_model.ApplicationModel" = application_model.ApplicationModel(project.Project())
        self._external_view = None
        self._protein_model = QtGui.QStandardItemModel()
        self._workspace_path = constants.DEFAULT_WORKSPACE_PATH
        self._workspace_status = f"Current workspace: {str(self._workspace_path)}"
        self._workspace_label = QtWidgets.QLabel(f"Current Workspace: {self._workspace_path}")

        self._setup_statusbar()
        self._connect_all_ui_elements_with_slot_functions()
        self._init_context_menus()
        self._interface_manager.refresh_main_view()
        if self._interface_manager.get_application_settings().start_help_at_startup == 1:
            self._start_documentation_server()

    def _connect_all_ui_elements_with_slot_functions(self):
        self._view.dialogClosed.connect(self._close_main_window)
        # <editor-fold desc="Menu">
        self._view.ui.action_new_project.triggered.connect(self._create_project)
        self._interface_manager.get_create_view().dialogClosed.connect(self.__await_create_project)
        self._view.ui.action_open_project.triggered.connect(self._open_project)
        self._interface_manager.get_open_view().dialogClosed.connect(self._post_open_project)
        self._view.ui.action_use_project.triggered.connect(self._use_project)
        self._view.ui.action_delete_project.triggered.connect(self._delete_project)
        self._interface_manager.get_delete_view().dialogClosed.connect(self._post_delete_project)
        self._view.ui.action_import_project.triggered.connect(self.import_project)
        self._view.ui.action_export_project.triggered.connect(self.export_current_project)
        self._view.ui.action_close_project.triggered.connect(self._close_project)

        self._view.ui.action_results_summary.triggered.connect(self._results_summary)
        self._view.ui.action_preview_image.triggered.connect(self._preview_image)
        self._view.ui.action_ray_tracing_image.triggered.connect(self._create_ray_traced_image)
        self._view.ui.action_simple_image.triggered.connect(self._create_drawn_image)
        self._view.ui.action_protein_regions.triggered.connect(self._hotspots_protein_regions)
        self._interface_manager.get_hotspots_protein_regions_view().dialogClosed.connect(self.post_hotspots_protein_regions)

        self._view.ui.action_edit_settings.triggered.connect(self.open_settings_global)
        self._view.ui.action_restore_settings.triggered.connect(self.restore_settings)
        self._view.ui.action_show_log_in_explorer.triggered.connect(self.open_logs)
        self._view.ui.action_clear_logs.triggered.connect(self.clear_all_log_files)
        self._view.ui.action_documentation.triggered.connect(self._open_help_center)
        self._view.ui.action_tutorials.triggered.connect(self.open_tutorial)
        self._view.ui.action_get_demo_projects.triggered.connect(self.get_demo_projects)
        self._view.ui.action_about.triggered.connect(self.open_about)
        self._view.ui.action_predict_monomer.triggered.connect(self._predict_monomer)
        self._view.ui.action_predict_multimer.triggered.connect(self._predict_multimer)
        self._view.ui.action_distance_analysis.triggered.connect(self._distance_analysis)
        self._view.ui.action_arrange_windows.triggered.connect(self.arrange_windows)

        self._view.ui.project_tab_widget.currentChanged.connect(self._update_tab)
        # </editor-fold>

        # <editor-fold desc="Sequence Tab">
        self._view.ui.seqs_list_view.clicked.connect(self._show_sequence_information)
        self._view.ui.btn_add_sequence.clicked.connect(self._add_sequence)
        self._view.ui.btn_import_seq.clicked.connect(self._import_sequence)
        self._view.ui.btn_save_sequence.clicked.connect(self._save_selected_sequence_as_fasta_file)
        self._view.ui.btn_delete_sequence.clicked.connect(self._delete_selected_sequence)
        self._view.ui.seqs_table_widget.cellClicked.connect(self._open_text_editor_for_seq)
        self._view.line_edit_seq_name.textChanged.connect(self._set_new_sequence_name_in_table_item)
        self._view.ui.seqs_table_widget.cellChanged.connect(self._rename_sequence)
        self._view.ui.btn_help.clicked.connect(self._open_sequences_tab_help)
        # </editor-fold>

        # <editor-fold desc="Proteins Tab">
        self._view.ui.proteins_tree_view.customContextMenuRequested.connect(self.open_context_menu_for_proteins)
        self._view.ui.proteins_tree_view.clicked.connect(self.__slot_get_information_about_selected_object_in_protein_branch)
        self._view.ui.btn_save_protein.clicked.connect(self._save_selected_protein_structure_as_pdb_file)
        # import
        self._view.ui.btn_import_protein.clicked.connect(self._import_protein_structure)
        self._interface_manager.get_add_protein_view().return_value.connect(self._post_import_protein_structure)
        self._view.ui.btn_open_protein_session.clicked.connect(self._open_protein_pymol_session)
        self._view.ui.btn_create_protein_scene.clicked.connect(self.save_scene)
        self._view.ui.btn_delete_protein.clicked.connect(self._delete_protein)
        self._view.ui.btn_update_protein_scene.clicked.connect(self._update_scene)
        self._view.ui.btn_delete_protein_scene.clicked.connect(self.delete_current_scene)
        self._view.ui.box_protein_color.currentIndexChanged.connect(self._change_chain_color_proteins)
        self._view.ui.btn_protein_color_atoms.clicked.connect(self._change_chain_color_proteins_atoms)
        self._view.ui.btn_protein_reset_atoms.clicked.connect(self._change_chain_reset_proteins_atoms)
        # self._view.ui.btn_protein_show_cartoon.clicked.connect(self.__slot_show_protein_chain_as_cartoon)
        # self._view.ui.btn_protein_hide_cartoon.clicked.connect(self.__slot_hide_protein_chain_as_cartoon)
        # self._view.ui.btn_protein_show_sticks.clicked.connect(self.__slot_show_protein_chain_as_sticks)
        # self._view.ui.btn_protein_hide_sticks.clicked.connect(self.__slot_hide_protein_chain_as_sticks)
        # self._view.ui.btn_protein_show_ribbon.clicked.connect(self.__slot_show_protein_chain_as_ribbon)
        # self._view.ui.btn_protein_hide_ribbon.clicked.connect(self.__slot_hide_protein_chain_as_ribbon)
        # self._view.ui.cb_protein_atoms.stateChanged.connect(self.__slot_color_protein_chain_atoms_by_element)
        self._view.ui.cb_protein_cartoon.stateChanged.connect(self.__slot_protein_chain_as_cartoon)
        self._view.ui.cb_protein_sticks.stateChanged.connect(self.__slot_protein_chain_as_sticks)
        self._view.ui.cb_protein_ribbon.stateChanged.connect(self.__slot_protein_chain_as_ribbon)
        self._view.ui.cb_protein_lines.stateChanged.connect(self.__slot_protein_chain_as_lines)
        self._view.ui.cb_protein_spheres.stateChanged.connect(self.__slot_protein_chain_as_spheres)
        self._view.ui.cb_protein_dots.stateChanged.connect(self.__slot_protein_chain_as_dots)
        self._view.ui.cb_protein_mesh.stateChanged.connect(self.__slot_protein_chain_as_mesh)
        self._view.ui.cb_protein_surface.stateChanged.connect(self.__slot_protein_chain_as_surface)
        # self._view.tg_protein_color_atoms.toggleChanged.connect(self.__slot_color_protein_chain_atoms_by_element)
        self._view.tg_protein_cartoon.toggleChanged.connect(self.__slot_protein_chain_as_cartoon)
        self._view.tg_protein_sticks.toggleChanged.connect(self.__slot_protein_chain_as_sticks)
        self._view.tg_protein_ribbon.toggleChanged.connect(self.__slot_protein_chain_as_ribbon)
        self._view.tg_protein_lines.toggleChanged.connect(self.__slot_protein_chain_as_lines)
        self._view.tg_protein_spheres.toggleChanged.connect(self.__slot_protein_chain_as_spheres)
        self._view.tg_protein_dots.toggleChanged.connect(self.__slot_protein_chain_as_dots)
        self._view.tg_protein_mesh.toggleChanged.connect(self.__slot_protein_chain_as_mesh)
        self._view.tg_protein_surface.toggleChanged.connect(self.__slot_protein_chain_as_surface)

        self._view.ui.btn_protein_hide_all_representations.clicked.connect(self.__slot_hide_protein_chain_all)

        # <editor-fold desc="Color Grid">
        self._view.color_grid.c_red.clicked.connect(self.set_color_name_in_label_red)
        self._view.color_grid.c_tv_red.clicked.connect(self.set_color_name_in_label_tv_red)
        self._view.color_grid.c_salomon.clicked.connect(self.set_color_name_in_label_salmon)
        self._view.color_grid.c_raspberry.clicked.connect(self.set_color_name_in_label_raspberry)

        self._view.color_grid.c_green.clicked.connect(self.set_color_name_in_label_green)
        self._view.color_grid.c_tv_green.clicked.connect(self.set_color_name_in_label_tv_green)
        self._view.color_grid.c_palegreen.clicked.connect(self.set_color_name_in_label_palegreen)
        self._view.color_grid.c_forest.clicked.connect(self.set_color_name_in_label_forest)

        self._view.color_grid.c_blue.clicked.connect(self.set_color_name_in_label_blue)
        self._view.color_grid.c_tv_blue.clicked.connect(self.set_color_name_in_label_tv_blue)
        self._view.color_grid.c_lightblue.clicked.connect(self.set_color_name_in_label_lightblue)
        self._view.color_grid.c_skyblue.clicked.connect(self.set_color_name_in_label_skyblue)

        self._view.color_grid.c_yellow.clicked.connect(self.set_color_name_in_label_yellow)
        self._view.color_grid.c_tv_yellow.clicked.connect(self.set_color_name_in_label_tv_yellow)
        self._view.color_grid.c_paleyellow.clicked.connect(self.set_color_name_in_label_paleyellow)
        self._view.color_grid.c_sand.clicked.connect(self.set_color_name_in_label_sand)

        self._view.color_grid.c_magenta.clicked.connect(self.set_color_name_in_label_magenta)
        self._view.color_grid.c_purple.clicked.connect(self.set_color_name_in_label_purple)
        self._view.color_grid.c_pink.clicked.connect(self.set_color_name_in_label_pink)
        self._view.color_grid.c_hotpink.clicked.connect(self.set_color_name_in_label_hotpink)

        self._view.color_grid.c_cyan.clicked.connect(self.set_color_name_in_label_cyan)
        self._view.color_grid.c_aquamarine.clicked.connect(self.set_color_name_in_label_aquamarine)
        self._view.color_grid.c_palecyan.clicked.connect(self.set_color_name_in_label_palecyan)
        self._view.color_grid.c_teal.clicked.connect(self.set_color_name_in_label_teal)

        self._view.color_grid.c_orange.clicked.connect(self.set_color_name_in_label_orange)
        self._view.color_grid.c_tv_orange.clicked.connect(self.set_color_name_in_label_tv_orange)
        self._view.color_grid.c_lightorange.clicked.connect(self.set_color_name_in_label_lightorange)
        self._view.color_grid.c_olive.clicked.connect(self.set_color_name_in_label_olive)

        self._view.color_grid.c_white.clicked.connect(self.set_color_name_in_label_white)
        self._view.color_grid.c_grey_70.clicked.connect(self.set_color_name_in_label_grey_70)
        self._view.color_grid.c_grey_30.clicked.connect(self.set_color_name_in_label_grey_30)
        self._view.color_grid.c_black.clicked.connect(self.set_color_name_in_label_black)
        # </editor-fold>

        # </editor-fold>

        # <editor-fold desc="Context Menu">
        self._interface_manager.get_rename_protein_view().dialogClosed.connect(
            self.post_rename_selected_protein_structure)

        self._view.ui.btn_help_2.clicked.connect(self._open_proteins_tab_help)
        # </editor-fold>

        # <editor-fold desc="Proteins Pair Tab">
        self._view.ui.protein_pairs_tree_view.customContextMenuRequested.connect(self.open_context_menu_for_protein_pairs)
        self._view.ui.protein_pairs_tree_view.clicked.connect(self.__slot_get_information_about_selected_object_in_protein_pair_branch)
        self._view.ui.btn_delete_protein_pair.clicked.connect(self._delete_protein_pair_from_project)
        self._view.ui.btn_open_protein_pair_session.clicked.connect(self._open_protein_pair_pymol_session)
        self._view.ui.btn_create_protein_pair_scene.clicked.connect(self.save_scene)
        self._view.ui.btn_update_protein_pair_scene.clicked.connect(self._update_scene)
        self._view.ui.btn_delete_protein_pair_scene.clicked.connect(self.delete_current_scene)
        self._view.ui.protein_pairs_tree_view.clicked.connect(self._check_for_results)
        self._view.ui.box_protein_pair_color.currentIndexChanged.connect(self._change_chain_color_protein_pairs)
        self._view.ui.btn_protein_pair_color_atoms.clicked.connect(self._change_chain_color_proteins_pair_atoms)
        self._view.ui.btn_protein_pair_reset_atoms.clicked.connect(self._change_chain_reset_proteins_pair_atoms)
        # self._view.ui.btn_protein_pair_show_cartoon.clicked.connect(self.__slot_show_protein_pair_chain_as_cartoon)
        # self._view.ui.btn_protein_pair_hide_cartoon.clicked.connect(self.__slot_hide_protein_pair_chain_as_cartoon)
        # self._view.ui.btn_protein_pair_show_sticks.clicked.connect(self.__slot_show_protein_pair_chain_as_sticks)
        # self._view.ui.btn_protein_pair_hide_sticks.clicked.connect(self.__slot_hide_protein_pair_chain_as_sticks)
        # self._view.ui.btn_protein_pair_show_ribbon.clicked.connect(self.__slot_show_protein_pair_chain_as_ribbon)
        # self._view.ui.btn_protein_pair_hide_ribbon.clicked.connect(self.__slot_hide_protein_pair_chain_as_ribbon)
        self._view.ui.cb_protein_pair_cartoon.stateChanged.connect(self.__slot_protein_pair_chain_as_cartoon)
        self._view.ui.cb_protein_pair_sticks.stateChanged.connect(self.__slot_protein_pair_chain_as_sticks)
        self._view.ui.cb_protein_pair_ribbon.stateChanged.connect(self.__slot_protein_pair_chain_as_ribbon)
        self._view.ui.cb_protein_pair_lines.stateChanged.connect(self.__slot_protein_pair_chain_as_lines)
        self._view.ui.cb_protein_pair_spheres.stateChanged.connect(self.__slot_protein_pair_chain_as_spheres)
        self._view.ui.cb_protein_pair_dots.stateChanged.connect(self.__slot_protein_pair_chain_as_dots)
        self._view.ui.cb_protein_pair_mesh.stateChanged.connect(self.__slot_protein_pair_chain_as_mesh)
        self._view.ui.cb_protein_pair_surface.stateChanged.connect(self.__slot_protein_pair_chain_as_surface)
        self._view.ui.btn_protein_pair_hide_all_representations.clicked.connect(self.__slot_hide_protein_pair_chain_all)
        self._view.ui.btn_help_3.clicked.connect(self._open_protein_pairs_tab_help)
        # </editor-fold>

    @staticmethod
    def _close_main_window():
        """Cleans after the main window closes."""
        # Closes the documentation browser if it is still open
        if len(pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)) == 1:
            pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)[0].close()

    # <editor-fold desc="Util methods">
    def update_status(self, message: str) -> None:
        """Updates the status bar of the main view with a custom message."""
        self._view.status_bar.showMessage(message)

    def _update_tab(self):
        self._interface_manager.current_tab_index = self._view.ui.project_tab_widget.currentIndex()
        if self._pymol_session_manager.session_object_type == "protein" and self._interface_manager.current_tab_index == 2:
            self._interface_manager.hide_protein_pair_pymol_scene_configuration()
            self._view.ui.lbl_info_protein_pair.setText(
                "Please load the PyMOL session of the \nselected protein pair.")
        elif self._pymol_session_manager.session_object_type == "protein_pair" and self._interface_manager.current_tab_index == 1:
            self._interface_manager.hide_protein_pymol_scene_configuration()
            self._view.ui.lbl_info.setText("Please load the PyMOL session of the selected protein.")
        elif self._pymol_session_manager.is_the_current_session_empty():
            self._interface_manager.hide_protein_pymol_scene_configuration()
            self._interface_manager.hide_protein_pair_pymol_scene_configuration()
            self._view.ui.lbl_info.setText("Please load the PyMOL session of the selected protein.")
            self._view.ui.lbl_info_protein_pair.setText(
                "Please load the PyMOL session of the \nselected protein pair.")

    def _setup_statusbar(self) -> None:
        """Sets up the status bar and fills it with the current workspace."""
        self._interface_manager.get_main_view().setStatusBar(self._interface_manager.get_main_view().status_bar)

    # <editor-fold desc="Help related methods">
    def _start_documentation_server(self):
        self._active_task = tasks.Task(
            target=util_async.start_documentation_server,
            args=(0, 0),
            post_func=self.__await_start_documentation_server,
        )
        self._active_task.start()

    def __await_start_documentation_server(self):
        self._interface_manager.update_status_bar("Opening help center finished.")

    def open_help(self, a_page_name: str):
        """Opens the pyssa documentation window if it's not already open.

        Args:
            a_page_name (str): a name of a documentation page to display
        """
        self._interface_manager.start_wait_spinner()
        self._interface_manager.update_status_bar("Opening help center ...")
        self._active_task = tasks.Task(
            target=util_async.open_documentation_on_certain_page,
            args=(a_page_name, 0),
            post_func=self.__await_open_help,
        )
        self._active_task.start()

    def __await_open_help(self):
        subprocess.run([constants.HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH])
        self._interface_manager.stop_wait_spinner()
        self._interface_manager.update_status_bar("Opening help center finished.")

    def _init_context_menus(self):
        # <editor-fold desc="General context menu setup">
        self.context_menu = QtWidgets.QMenu()
        self.help_context_action = self.context_menu.addAction(self._view.tr("Get Help"))
        self.help_context_action.triggered.connect(self._open_sequences_tab_help)

        # </editor-fold>

        # Set the context menu for the buttons
        self._view.ui.seqs_list_view.setContextMenuPolicy(3)
        self._view.ui.seqs_list_view.customContextMenuRequested.connect(self._show_context_menu_for_seq_list)
        self._view.ui.seqs_table_widget.setContextMenuPolicy(3)
        self._view.ui.seqs_table_widget.customContextMenuRequested.connect(self._show_context_menu_for_seq_table)
        self._view.ui.btn_import_seq.setContextMenuPolicy(3)  # 3 corresponds to Qt.CustomContextMenu
        self._view.ui.btn_import_seq.customContextMenuRequested.connect(self._show_context_menu_for_seq_import)
        self._view.ui.btn_add_sequence.setContextMenuPolicy(3)  # 3 corresponds to Qt.CustomContextMenu
        self._view.ui.btn_add_sequence.customContextMenuRequested.connect(self._show_context_menu_for_seq_add)
        self._view.ui.btn_save_sequence.setContextMenuPolicy(3)  # 3 corresponds to Qt.CustomContextMenu
        self._view.ui.btn_save_sequence.customContextMenuRequested.connect(self._show_context_menu_for_seq_save)
        self._view.ui.btn_delete_sequence.setContextMenuPolicy(3)  # 3 corresponds to Qt.CustomContextMenu
        self._view.ui.btn_delete_sequence.customContextMenuRequested.connect(self._show_context_menu_for_seq_delete)
        self._view.ui.proteins_tree_view.setContextMenuPolicy(3)
        self._view.ui.proteins_tree_view.customContextMenuRequested.connect(
            self._show_context_menu_for_proteins_tree_view
        )
        self._view.ui.btn_import_protein.setContextMenuPolicy(3)
        self._view.ui.btn_import_protein.customContextMenuRequested.connect(
            self._show_context_menu_for_protein_import
        )
        self._view.ui.btn_save_protein.setContextMenuPolicy(3)
        self._view.ui.btn_save_protein.customContextMenuRequested.connect(
            self._show_context_menu_for_protein_save
        )
        self._view.ui.btn_delete_protein.setContextMenuPolicy(3)
        self._view.ui.btn_delete_protein.customContextMenuRequested.connect(
            self._show_context_menu_for_protein_delete
        )
        self._view.ui.btn_open_protein_session.setContextMenuPolicy(3)
        self._view.ui.btn_open_protein_session.customContextMenuRequested.connect(
            self._show_context_menu_for_protein_load_session
        )
        self._view.ui.btn_create_protein_scene.setContextMenuPolicy(3)
        self._view.ui.btn_create_protein_scene.customContextMenuRequested.connect(
            self._show_context_menu_for_protein_add_scene
        )
        self._view.ui.btn_update_protein_scene.setContextMenuPolicy(3)
        self._view.ui.btn_update_protein_scene.customContextMenuRequested.connect(
            self._show_context_menu_for_protein_update_scene
        )
        self._view.ui.btn_delete_protein_scene.setContextMenuPolicy(3)
        self._view.ui.btn_delete_protein_scene.customContextMenuRequested.connect(
            self._show_context_menu_for_protein_delete_scene
        )
        self._view.ui.frame_protein_pymol_scene.setContextMenuPolicy(3)
        self._view.ui.frame_protein_pymol_scene.customContextMenuRequested.connect(
            self._show_context_menu_for_protein_pymol_scene_config
        )
        # add more buttons here ...
        self._view.ui.protein_pairs_tree_view.setContextMenuPolicy(3)
        self._view.ui.protein_pairs_tree_view.customContextMenuRequested.connect(
            self._show_context_menu_for_protein_pairs_tree_view
        )
        self._view.ui.btn_delete_protein_pair.setContextMenuPolicy(3)
        self._view.ui.btn_delete_protein_pair.customContextMenuRequested.connect(
            self._show_context_menu_for_protein_pair_delete
        )
        self._view.ui.btn_open_protein_pair_session.setContextMenuPolicy(3)
        self._view.ui.btn_open_protein_pair_session.customContextMenuRequested.connect(
            self._show_context_menu_for_protein_pair_load_session
        )
        self._view.ui.btn_create_protein_pair_scene.setContextMenuPolicy(3)
        self._view.ui.btn_create_protein_pair_scene.customContextMenuRequested.connect(
            self._show_context_menu_for_protein_pair_add_scene
        )
        self._view.ui.btn_update_protein_pair_scene.setContextMenuPolicy(3)
        self._view.ui.btn_update_protein_pair_scene.customContextMenuRequested.connect(
            self._show_context_menu_for_protein_pair_update_scene
        )
        self._view.ui.btn_delete_protein_pair_scene.setContextMenuPolicy(3)
        self._view.ui.btn_delete_protein_pair_scene.customContextMenuRequested.connect(
            self._show_context_menu_for_protein_pair_delete_scene
        )
        self._view.ui.frame_protein_pair_pymol_scene.setContextMenuPolicy(3)
        self._view.ui.frame_protein_pair_pymol_scene.customContextMenuRequested.connect(
            self._show_context_menu_for_protein_pair_pymol_scene_config
        )

    # <editor-fold desc="Help pages">
    def _open_help_center(self):
        self.open_help("help/")

    def _open_sequences_tab_help(self):
        self.open_help("help/sequences/sequences_tab/")

    def _open_additional_information_table_help(self):
        self.open_help("help/sequences/additional_sequence_information/")

    def _open_sequence_import_help(self):
        self.open_help("help/sequences/sequence_import/")

    def _open_sequence_add_help(self):
        self.open_help("help/sequences/sequence_add/")

    def _open_sequence_save_help(self):
        self.open_help("help/sequences/sequence_save/")

    def _open_sequence_delete_help(self):
        self.open_help("help/sequences/sequence_delete/")

    def _open_proteins_tab_help(self):
        self.open_help("help/proteins/proteins_tab/")

    def _open_protein_import_help(self):
        self.open_help("help/proteins/protein_import/")

    def _open_protein_save_help(self):
        self.open_help("help/proteins/protein_save/")

    def _open_protein_delete_help(self):
        self.open_help("help/proteins/protein_delete/")

    def _open_protein_pymol_scene_config_help(self):
        self.open_help("help/proteins/protein_pymol_scene_configuration/")

    def _open_protein_load_session_help(self):
        self.open_help("help/proteins/protein_load_session/")

    def _open_protein_add_scene_help(self):
        self.open_help("help/proteins/protein_add_scene/")

    def _open_protein_update_scene_help(self):
        self.open_help("help/proteins/protein_update_scene/")

    def _open_protein_delete_scene_help(self):
        self.open_help("help/proteins/protein_delete_scene/")

    def _open_protein_pairs_tab_help(self):
        self.open_help("help/protein_pairs/protein_pairs_tab/")

    def _open_protein_pair_delete_help(self):
        self.open_help("help/protein_pairs/protein_pair_delete/")

    def _open_protein_pair_pymol_scene_config_help(self):
        self.open_help("help/protein_pairs/protein_pair_pymol_scene_configuration/")

    def _open_protein_pair_load_session_help(self):
        self.open_help("help/protein_pairs/protein_pair_load_session/")

    def _open_protein_pair_add_scene_help(self):
        self.open_help("help/protein_pairs/protein_pair_add_scene/")

    def _open_protein_pair_update_scene_help(self):
        self.open_help("help/protein_pairs/protein_pair_update_scene/")

    def _open_protein_pair_delete_scene_help(self):
        self.open_help("help/protein_pairs/protein_pair_delete_scene/")

    # </editor-fold>

    # <editor-fold desc="Context menu connections">
    def _show_context_menu_for_seq_list(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_sequences_tab_help)
        self.context_menu.exec_(self._view.ui.seqs_list_view.mapToGlobal(a_point))

    def _show_context_menu_for_seq_table(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_additional_information_table_help)
        self.context_menu.exec_(self._view.ui.seqs_table_widget.mapToGlobal(a_point))

    def _show_context_menu_for_seq_import(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_sequence_import_help)
        self.context_menu.exec_(self._view.ui.btn_import_seq.mapToGlobal(a_point))

    def _show_context_menu_for_seq_add(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_sequence_add_help)
        self.context_menu.exec_(self._view.ui.btn_add_sequence.mapToGlobal(a_point))

    def _show_context_menu_for_seq_save(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_sequence_save_help)
        self.context_menu.exec_(self._view.ui.btn_save_sequence.mapToGlobal(a_point))

    def _show_context_menu_for_seq_delete(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_sequence_delete_help)
        self.context_menu.exec_(self._view.ui.btn_delete_sequence.mapToGlobal(a_point))

    def _show_context_menu_for_proteins_tree_view(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_proteins_tab_help)
        self.context_menu.exec_(self._view.ui.proteins_tree_view.mapToGlobal(a_point))

    def _show_context_menu_for_protein_import(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_protein_import_help)
        self.context_menu.exec_(self._view.ui.btn_import_protein.mapToGlobal(a_point))

    def _show_context_menu_for_protein_save(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_protein_save_help)
        self.context_menu.exec_(self._view.ui.btn_save_protein.mapToGlobal(a_point))

    def _show_context_menu_for_protein_delete(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_protein_delete_help)
        self.context_menu.exec_(self._view.ui.btn_delete_protein.mapToGlobal(a_point))

    def _show_context_menu_for_protein_load_session(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_protein_load_session_help)
        self.context_menu.exec_(self._view.ui.btn_open_protein_session.mapToGlobal(a_point))

    def _show_context_menu_for_protein_add_scene(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_protein_add_scene_help)
        self.context_menu.exec_(self._view.ui.btn_create_protein_scene.mapToGlobal(a_point))

    def _show_context_menu_for_protein_update_scene(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_protein_update_scene_help)
        self.context_menu.exec_(self._view.ui.btn_update_protein_scene.mapToGlobal(a_point))

    def _show_context_menu_for_protein_delete_scene(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_protein_delete_scene_help)
        self.context_menu.exec_(self._view.ui.btn_delete_protein_scene.mapToGlobal(a_point))

    def _show_context_menu_for_protein_pymol_scene_config(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_protein_pymol_scene_config_help)
        self.context_menu.exec_(self._view.ui.frame_protein_pymol_scene.mapToGlobal(a_point))

    def _show_context_menu_for_protein_pairs_tree_view(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_protein_pairs_tab_help)
        self.context_menu.exec_(self._view.ui.protein_pairs_tree_view.mapToGlobal(a_point))

    def _show_context_menu_for_protein_pair_delete(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_protein_pair_delete_help)
        self.context_menu.exec_(self._view.ui.btn_delete_protein_pair.mapToGlobal(a_point))

    def _show_context_menu_for_protein_pair_load_session(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_protein_pair_load_session_help)
        self.context_menu.exec_(self._view.ui.btn_open_protein_pair_session.mapToGlobal(a_point))

    def _show_context_menu_for_protein_pair_add_scene(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_protein_pair_add_scene_help)
        self.context_menu.exec_(self._view.ui.btn_create_protein_pair_scene.mapToGlobal(a_point))

    def _show_context_menu_for_protein_pair_update_scene(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_protein_pair_update_scene_help)
        self.context_menu.exec_(self._view.ui.btn_update_protein_pair_scene.mapToGlobal(a_point))

    def _show_context_menu_for_protein_pair_delete_scene(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_protein_pair_delete_scene_help)
        self.context_menu.exec_(self._view.ui.btn_delete_protein_pair_scene.mapToGlobal(a_point))

    def _show_context_menu_for_protein_pair_pymol_scene_config(self, a_point):
        self.help_context_action.triggered.disconnect()
        self.help_context_action.triggered.connect(self._open_protein_pair_pymol_scene_config_help)
        self.context_menu.exec_(self._view.ui.frame_protein_pair_pymol_scene.mapToGlobal(a_point))

    # </editor-fold>

    # </editor-fold>

    # </editor-fold>

    # <editor-fold desc="Menu bar methods">

    # <editor-fold desc="Project menu">
    def _close_project(self):
        """Closes the current project"""
        self._active_task = tasks.Task(
            target=main_presenter_async.close_project,
            args=(self._database_thread, self._pymol_session_manager),
            post_func=self.__await_close_project,
        )
        self._active_task.start()
        # self.msg_box = basic_boxes.no_buttons("Saving Project",
        #                                       "Please wait the program is saving your project.",
        #                                       QtWidgets.QMessageBox.Information)
        # # self.msg_box.show()
        self._interface_manager.restore_default_main_view()
        self.update_status("Saving current project ...")
        self._view.wait_spinner.start()

    def __await_close_project(self):
        """Await the async closing process."""
        self._interface_manager.set_new_project(project.Project())
        self._interface_manager.refresh_main_view()
        # self.msg_box.hide()
        self.update_status("Closing project finished.")
        self._view.wait_spinner.stop()

    def _create_project(self) -> None:
        self._external_controller = create_project_view_controller.CreateProjectViewController(self._interface_manager)
        self._external_controller.user_input.connect(self._post_create_project)
        self._interface_manager.get_create_view().show()

    def _post_create_project(self, user_input: tuple) -> None:
        self._interface_manager.start_wait_spinner()
        tmp_project_name, tmp_protein_name = user_input
        tmp_project_database_filepath = str(
            pathlib.Path(f"{self._interface_manager.get_application_settings().workspace_path}/{tmp_project_name}.db"))
        with database_manager.DatabaseManager(tmp_project_database_filepath) as db_manager:
            db_manager.build_new_database()
            db_manager.close_project_database()
        self._database_thread = database_thread.DatabaseThread(tmp_project_database_filepath)
        self._database_thread.start()
        self._active_task = tasks.Task(
            target=main_presenter_async.create_new_project,
            args=(
                tmp_project_name,
                self._interface_manager.get_application_settings().get_workspace_path(),
                tmp_protein_name
            ),
            post_func=self.__await_create_project,
        )
        self._active_task.start()

        #
        # tmp_project_name, tmp_protein_name = user_input
        # tmp_project_database_filepath = str(pathlib.Path(f"{self._interface_manager.get_application_settings().workspace_path}/{tmp_project_name}.db"))
        # self._database_manager.set_database_filepath(tmp_project_database_filepath)
        # self._database_manager.build_new_database()
        # self._database_thread = database_thread.DatabaseThread(tmp_project_database_filepath)
        # self._database_thread.start()
        # self._database_manager.open_project_database()
        # tmp_project = project.Project(tmp_project_name,
        #                               self._interface_manager.get_application_settings().workspace_path)
        # tmp_project.set_id(self._database_manager.insert_new_project(tmp_project.get_project_name(),
        #                                                              platform.system()))
        # if len(tmp_protein_name) == 4:
        #     tmp_ref_protein = protein.Protein(tmp_protein_name.upper())
        #     tmp_ref_protein.db_project_id = tmp_project.get_id()
        #     tmp_ref_protein.add_protein_structure_data_from_pdb_db(tmp_protein_name.upper())
        #     tmp_ref_protein.create_new_pymol_session()
        #     tmp_ref_protein.save_pymol_session_as_base64_string()
        #     tmp_project.add_existing_protein(tmp_ref_protein)
        #     tmp_ref_protein.set_id(self._database_manager.insert_new_protein(tmp_ref_protein))
        #     constants.PYSSA_LOGGER.info("Create project finished with protein from the PDB.")
        # elif len(tmp_protein_name) > 0:
        #     # local pdb file as input
        #     pdb_filepath = pathlib.Path(tmp_protein_name)
        #     graphic_operations.setup_default_session_graphic_settings()
        #     tmp_ref_protein = protein.Protein(
        #         pdb_filepath.name.replace(".pdb","")
        #     )
        #     tmp_ref_protein.db_project_id = tmp_project.get_id()
        #     tmp_ref_protein.add_protein_structure_data_from_local_pdb_file(pathlib.Path(tmp_protein_name))
        #     tmp_ref_protein.create_new_pymol_session()
        #     tmp_ref_protein.save_pymol_session_as_base64_string()
        #     tmp_project.add_existing_protein(tmp_ref_protein)
        #     tmp_ref_protein.db_project_id = self._database_manager.insert_new_protein(tmp_ref_protein)
        #     constants.PYSSA_LOGGER.info("Create project finished with protein from local filesystem.")
        # else:
        #     constants.PYSSA_LOGGER.info("Create empty project finished.")

    def __await_create_project(self, return_value: tuple):
        if return_value[1] is False:
            self._interface_manager.stop_wait_spinner()
            return

        _, tmp_project = return_value
        self._interface_manager.set_new_project(tmp_project)
        self._interface_manager.refresh_workspace_model()
        self._interface_manager.refresh_main_view()
        self._pymol_session_manager.reinitialize_session()
        self._interface_manager.stop_wait_spinner()

    def _open_project(self) -> None:
        self._external_controller = open_project_view_controller.OpenProjectViewController(self._interface_manager)
        self._external_controller.return_value.connect(self._post_open_project)
        self._interface_manager.get_open_view().show()

    def _post_open_project(self, return_value: str):
        if return_value[1] is False:
            return

        self._interface_manager.update_status_bar("Opening existing project ...")
        #self._interface_manager.update_progress_bar(10, "Opening existing project ...")
        self._interface_manager.start_wait_spinner()
        tmp_project_name = return_value
        tmp_project_database_filepath = str(
            pathlib.Path(
                f"{self._interface_manager.get_application_settings().workspace_path}/{tmp_project_name}.db"
            )
        )
        self._database_thread = database_thread.DatabaseThread(tmp_project_database_filepath)
        self._database_thread.start()
        self._interface_manager.update_progress_of_progress_bar(20)
        self._database_manager.set_database_filepath(tmp_project_database_filepath)
        self._active_task = tasks.Task(
            target=main_presenter_async.open_project,
            args=(
                tmp_project_name,
                tmp_project_database_filepath,
                self._interface_manager,
                self._pymol_session_manager,
            ),
            post_func=self.__await_open_project,
        )
        self._active_task.start()
        self._interface_manager.update_progress_of_progress_bar(50)

    def __await_open_project(self, return_value: tuple):
        self._interface_manager.update_progress_of_progress_bar(60)
        exit_code, tmp_project, tmp_interface_manager = return_value
        if exit_code == 0:
            self._interface_manager = tmp_interface_manager
            self._interface_manager.refresh_main_view()
            self._interface_manager.hide_progress_bar()
            self._interface_manager.update_status_bar("Opening existing project finished.")
        else:
            self._interface_manager.update_status_bar("Opening existing project failed!")
        self._interface_manager.stop_wait_spinner()

    def _use_project(self) -> None:
        self._external_controller = use_project_view_controller.UseProjectViewController(self._interface_manager)
        self._external_controller.user_input.connect(self._post_use_project)
        self._interface_manager.get_use_project_view().show()

    def _post_use_project(self, user_input: tuple) -> None:
        self._interface_manager.start_wait_spinner()
        tmp_project_database_filepath = str(pathlib.Path(f"{self._interface_manager.get_application_settings().get_workspace_path()}/{user_input[0]}.db"))
        with database_manager.DatabaseManager(tmp_project_database_filepath) as db_manager:
            db_manager.build_new_database()
            db_manager.close_project_database()

        self._active_task = tasks.Task(
            target=main_presenter_async.create_use_project,
            args=(
                user_input[0],
                self._interface_manager.get_application_settings().get_workspace_path(),
                user_input[1]
            ),
            post_func=self.__await_use_project,
        )
        self._active_task.start()

    def __await_use_project(self, return_value: tuple):
        _, tmp_project = return_value
        self._interface_manager.set_new_project(tmp_project)
        self._interface_manager.refresh_main_view()
        self._pymol_session_manager.reinitialize_session()
        self._interface_manager.stop_wait_spinner()
        self._interface_manager.update_status_bar("Use process finished.")

    # def __await_open_project(self, a_result: tuple) -> None:
    #     self._interface_manager.set_new_project(a_result[1])
    #     self._interface_manager.update_status_bar(self._workspace_status)
    #     self._interface_manager.refresh_main_view()
    #     self._interface_manager.stop_wait_spinner()

    def _delete_project(self) -> None:
        self._external_controller = delete_project_view_controller.DeleteProjectViewController(self._interface_manager)
        self._interface_manager.get_delete_view().show()

    def _post_delete_project(self) -> None:
        self._interface_manager.refresh_main_view()

    def import_project(self) -> None:
        """Imports a project.xml into the current workspace."""
        file_dialog = QtWidgets.QFileDialog()
        desktop_path = QtCore.QStandardPaths.standardLocations(QtCore.QStandardPaths.DesktopLocation)[0]
        file_dialog.setDirectory(desktop_path)
        file_path, _ = file_dialog.getOpenFileName(
            self._view,
            "Select a project file to import",
            "",
            "Project Database File (*.db)",
        )
        if file_path:
            file = QtCore.QFile(file_path)
            if not file.open(QtCore.QFile.ReadOnly | QtCore.QFile.Text):
                print("Error: Cannot open file for reading")
                return
            tmp_import_filepath = pathlib.Path(file_path)
            tmp_project_name_input_dialog = QtWidgets.QInputDialog()
            tmp_new_project_name = tmp_project_name_input_dialog.getText(
                self._view,
                "Project Name",
                "Enter A Project Name:",
                text=tmp_import_filepath.name.replace(".db", "")
            )[0]
            # Copy db file into new workspace
            tmp_project_database_filepath = str(
                pathlib.Path(
                    f"{self._interface_manager.get_application_settings().workspace_path}/{tmp_new_project_name}.db"
                )
            )
            shutil.copyfile(file_path, tmp_project_database_filepath)
            # Open db and create project
            self._database_thread = database_thread.DatabaseThread(tmp_project_database_filepath)
            self._database_thread.start()
            self._database_manager.set_database_filepath(tmp_project_database_filepath)
            self._database_manager.open_project_database()
            self._database_manager.update_project_name(tmp_new_project_name)
            tmp_project = self._database_manager.get_project_as_object(
                tmp_new_project_name,
                self._interface_manager.get_application_settings().workspace_path,
                self._interface_manager.get_application_settings()
            )
            self._interface_manager.set_new_project(tmp_project)

            self._interface_manager.refresh_main_view()
            self._pymol_session_manager.reinitialize_session()
            self._interface_manager.stop_wait_spinner()
            self._interface_manager.update_status_bar("Importing project finished.")

            # tmp_project = project.Project()
            # handler = filesystem_io.ProjectParserHandler(tmp_project,
            #                                              self._interface_manager.get_application_settings())
            # parser = sax.make_parser()
            # parser.setContentHandler(handler)
            # parser.parse(file_path)
            # file.close()
            # tmp_project = handler.get_project()
            #
            # tmp_project.set_workspace_path(self._workspace_path)
            # if len(tmp_project.proteins) <= 1:
            #     if self._interface_manager.get_application_settings().wsl_install == 0:
            #         basic_boxes.ok(
            #             "Create new project",
            #             "Please install local colabfold to import this project!",
            #             QtWidgets.QMessageBox.Warning,
            #         )
            #         return
            #     elif self._interface_manager.get_application_settings().local_colabfold == 0:  # noqa: RET505
            #         basic_boxes.ok(
            #             "Create new project",
            #             "Please install local colabfold to import this project!",
            #             QtWidgets.QMessageBox.Warning,
            #         )
            #         return
            # new_filepath = pathlib.Path(f"{self._workspace_path}/{tmp_project.get_project_name()}.xml")
            # tmp_project.serialize_project(new_filepath)
            # self._interface_manager.set_new_project(
            #     self._interface_manager.get_current_project().deserialize_project(
            #         new_filepath, self._interface_manager.get_application_settings()
            #     )
            # )
            # constants.PYSSA_LOGGER.info(
            #     f"Opening the project {self._interface_manager.get_current_project().get_project_name()}."
            # )
            # self._view.ui.lbl_project_name.setText(self._interface_manager.get_current_project().get_project_name())
            # self._interface_manager.refresh_main_view()
            # basic_boxes.ok(
            #     "Import Project",
            #     "The project was successfully imported.",
            #     QtWidgets.QMessageBox.Information,
            # )

    def export_current_project(self) -> None:
        """Exports the current project to an importable format."""
        file_dialog = QtWidgets.QFileDialog()
        desktop_path = QtCore.QStandardPaths.standardLocations(QtCore.QStandardPaths.DesktopLocation)[0]
        file_dialog.setDirectory(desktop_path)
        file_path, _ = file_dialog.getSaveFileName(self._view, "Export current project", "", "Project Database File (*.db)")
        if file_path:
            shutil.copyfile(self._interface_manager.get_current_project().get_database_filepath(),
                            file_path)
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "The project was successfully exported.", "Export Project",
                custom_message_box.CustomMessageBoxIcons.INFORMATION.value
            )
            tmp_dialog.exec_()

    # </editor-fold>

    # <editor-fold desc="Analysis menu">
    def _distance_analysis(self):
        self._external_controller = distance_analysis_view_controller.DistanceAnalysisViewController(
            self._interface_manager
        )
        self._external_controller.job_input.connect(self.start_process_batch)
        self._interface_manager.get_distance_analysis_view().show()

    def start_process_batch(self, job_input: tuple) -> None:
        """Sets up the worker for the analysis task."""
        constants.PYSSA_LOGGER.info("Begin analysis process.")

        _, tmp_raw_analysis_run_names, tmp_checkbox_state = job_input

        self._active_task = tasks.Task(
            target=main_presenter_async.run_distance_analysis,
            args=(
                tmp_raw_analysis_run_names,
                self._interface_manager.get_current_project(),
                self._interface_manager.get_application_settings(),
                tmp_checkbox_state,
            ),
            post_func=self.post_analysis_process,
        )
        self._active_task.start()

        if not os.path.exists(constants.SCRATCH_DIR_ANALYSIS):
            os.mkdir(constants.SCRATCH_DIR_ANALYSIS)
        #self.block_box_analysis.exec_()

    # </editor-fold>

    # <editor-fold desc="Prediction menu">

    # <editor-fold desc="Monomer">
    def _predict_monomer(self):
        self._external_controller = predict_monomer_view_controller.PredictMonomerViewController(
            self._interface_manager
        )
        self._external_controller.job_input.connect(self._post_predict_monomer)
        self._interface_manager.get_predict_monomer_view().show()

    def _post_predict_monomer(self, result: tuple) -> None:
        """Sets up the worker for the prediction of the proteins."""
        self._view.wait_spinner.start()

        # <editor-fold desc="Check if WSL2 and ColabFold are installed">
        if globals.g_os == "win32":
            constants.PYSSA_LOGGER.info("Checking if WSL2 is installed ...")
            if not dialog_settings_global.is_wsl2_installed():
                constants.PYSSA_LOGGER.warning("WSL2 is NOT installed.")
                self._interface_manager.get_application_settings().wsl_install = 0
                tmp_dialog = custom_message_box.CustomMessageBoxOk(
                    "Prediction failed because the WSL2 environment is not installed!",
                    "Structure Prediction",
                    custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
                )
                tmp_dialog.exec_()
                return
            constants.PYSSA_LOGGER.info("Checking if Local Colabfold is installed ...")
            if not dialog_settings_global.is_local_colabfold_installed():
                constants.PYSSA_LOGGER.warning("Local Colabfold is NOT installed.")
                self._interface_manager.get_application_settings().local_colabfold = 0
                tmp_dialog = custom_message_box.CustomMessageBoxOk(
                    "Prediction failed because the ColabFold is not installed!",
                    "Structure Prediction",
                    custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
                )
                tmp_dialog.exec_()
                return

        # </editor-fold>

        self.prediction_type = constants.PREDICTION_TYPE_PRED_MONO_ANALYSIS
        constants.PYSSA_LOGGER.info("Begin prediction process.")
        self._interface_manager.update_status_bar("Begin prediction process ...")
        if result[3] is True:
            constants.PYSSA_LOGGER.info("Running prediction with subsequent analysis.")
            # Analysis should be run after the prediction
            self._active_task = tasks.Task(
                target=main_presenter_async.predict_protein_with_colabfold,
                args=(
                    result[1],
                    result[2],
                    self._interface_manager.get_current_project(),
                ),
                post_func=self.__await_monomer_prediction_for_subsequent_analysis,
            )
            self._active_task.start()
        else:
            constants.PYSSA_LOGGER.info("Running only a prediction.")
            # No analysis after prediction
            self._active_task = tasks.Task(
                target=main_presenter_async.predict_protein_with_colabfold,
                args=(
                    result[1],
                    result[2],
                    self._interface_manager.get_current_project(),
                ),
                post_func=self.__await_predict_protein_with_colabfold,
            )
            self._active_task.start()

        self._view.status_bar.showMessage("A prediction is currently running ...")
        self.block_box_prediction = QtWidgets.QMessageBox()
        self.block_box_prediction.setIcon(QtWidgets.QMessageBox.Information)
        self.block_box_prediction.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        styles.set_stylesheet(self.block_box_prediction)
        self.block_box_prediction.setWindowTitle("Structure Prediction")
        self.block_box_prediction.setText("A prediction is currently running.")
        btn_abort = self.block_box_prediction.addButton("Abort", QtWidgets.QMessageBox.ActionRole)
        self.block_box_prediction.exec_()
        if self.block_box_prediction.clickedButton() == btn_abort:
            self.abort_prediction()
            self.block_box_prediction.close()
            self._view.wait_spinner.stop()
        else:
            self.block_box_prediction.close()
            self._view.wait_spinner.stop()
        # self._interface_manager.update_status_bar("")

    def abort_prediction(self) -> None:
        """Aborts the running prediction."""
        constants.PYSSA_LOGGER.info("Structure prediction process was aborted manually.")
        subprocess.run(["wsl", "--shutdown"])
        constants.PYSSA_LOGGER.info("Shutdown of wsl environment.")
        filesystem_io.FilesystemCleaner.clean_prediction_scratch_folder()
        constants.PYSSA_LOGGER.info("Cleaned scratch directory.")
        tmp_dialog = custom_message_box.CustomMessageBoxOk(
            "The structure prediction was aborted.",
            "Abort Structure Prediction",
            custom_message_box.CustomMessageBoxIcons.INFORMATION.value
        )
        tmp_dialog.exec_()
        self._interface_manager.refresh_main_view()

    def __await_monomer_prediction_for_subsequent_analysis(self, result: tuple) -> None:
        tmp_exit_code = result[0]
        tmp_exit_code_description = [1]
        if tmp_exit_code == exit_codes.EXIT_CODE_ZERO[0]:
            # Prediction was successful
            self.block_box_prediction.destroy(True)
            constants.PYSSA_LOGGER.info("All structure predictions are done.")
            self._interface_manager.update_status_bar("All structure predictions are done.")
            constants.PYSSA_LOGGER.info("Begin analysis process.")
            self._interface_manager.update_status_bar("Running analysis process ...")

            tmp_raw_analysis_run_names: list = []
            for row_no in range(self._interface_manager.get_predict_monomer_view().ui.list_pred_analysis_mono_overview.count()):
                tmp_raw_analysis_run_names.append(
                    self._interface_manager.get_predict_monomer_view().ui.list_pred_analysis_mono_overview.item(row_no).text())

            self._database_manager.close_project_database()
            self._interface_manager.refresh_protein_model()
            self._active_task = tasks.Task(
                target=main_presenter_async.run_distance_analysis,
                args=(
                    tmp_raw_analysis_run_names,
                    self._interface_manager.get_current_project(),
                    self._interface_manager.get_application_settings(),
                    self._interface_manager.get_predict_monomer_view().ui.cb_pred_analysis_mono_images.isChecked(),
                ),
                post_func=self.post_analysis_process,
            )
            self._active_task.start()
            if not os.path.exists(constants.SCRATCH_DIR_ANALYSIS):
                os.mkdir(constants.SCRATCH_DIR_ANALYSIS)

        elif tmp_exit_code == exit_codes.ERROR_WRITING_FASTA_FILES[0]:
            self.block_box_prediction.destroy(True)
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Prediction failed because there was an error writing the fasta file(s)!",
                "Structure Prediction",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
        elif tmp_exit_code == exit_codes.ERROR_FASTA_FILES_NOT_FOUND[0]:
            self.block_box_prediction.destroy(True)
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Prediction failed because the fasta file(s) could not be found!",
                "Structure Prediction",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
        elif tmp_exit_code == exit_codes.ERROR_PREDICTION_FAILED[0]:
            self.block_box_prediction.destroy(True)
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Prediction failed because a subprocess failed!",
                "Structure Prediction",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
            self._view.wait_spinner.stop()
        elif tmp_exit_code == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
            self.block_box_prediction.destroy(True)
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Prediction failed because of an unknown error!",
                "Structure Prediction",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
        self._view.wait_spinner.stop()
        self._interface_manager.refresh_main_view()

    def post_analysis_process(self, an_exit_code: tuple[int, str, list]) -> None:
        """Post process after the analysis thread finished."""
        constants.PYSSA_LOGGER.debug("post_analysis_process() started ...")
        if an_exit_code[0] == exit_codes.ERROR_DISTANCE_ANALYSIS_FAILED[0]:
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Distance analysis failed because there was an error during the analysis!",
                "Distance Analysis",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
            constants.PYSSA_LOGGER.error(
                f"Distance analysis ended with exit code {an_exit_code[0]}: {an_exit_code[1]}",
            )
            self._interface_manager.update_status_bar(
                f"Distance analysis ended with exit code {an_exit_code[0]}: {an_exit_code[1]}")
        elif an_exit_code[0] == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Distance analysis failed because of an unknown error!",
                "Distance Analysis",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
            constants.PYSSA_LOGGER.error(
                f"Distance analysis ended with exit code {an_exit_code[0]}: {an_exit_code[1]}",
            )
            self._interface_manager.update_status_bar(
                f"Distance analysis ended with exit code {an_exit_code[0]}: {an_exit_code[1]}")
        elif an_exit_code[0] == exit_codes.EXIT_CODE_ZERO[0]:
            constants.PYSSA_LOGGER.info("Project has been saved to project database.")
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "All structure analysis' are done. Go to results to check the new results.",
                "Distance Analysis",
                custom_message_box.CustomMessageBoxIcons.INFORMATION.value
            )
            tmp_dialog.exec_()
            constants.PYSSA_LOGGER.info("All structure analysis' are done.")
            self._interface_manager.update_status_bar("All structure analysis' are done.")
        self._database_manager.open_project_database()
        self._interface_manager.refresh_protein_pair_model()
        self._interface_manager.refresh_main_view()
        self._interface_manager.stop_wait_spinner()

    def __await_predict_protein_with_colabfold(self, result: tuple) -> None:
        """Process which runs after each prediction job."""
        print(result)
        tmp_exit_code = result[0]
        tmp_exit_code_description = result[1]
        if tmp_exit_code == exit_codes.ERROR_WRITING_FASTA_FILES[0]:
            self.block_box_prediction.destroy(True)
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Prediction failed because there was an error writing the fasta file(s)!",
                "Structure Prediction",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
            self._interface_manager.update_status_bar(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}")
        elif tmp_exit_code == exit_codes.ERROR_FASTA_FILES_NOT_FOUND[0]:
            self.block_box_prediction.destroy(True)
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Prediction failed because the fasta file(s) could not be found!",
                "Structure Prediction",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
            self._interface_manager.update_status_bar(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}")
        elif tmp_exit_code == exit_codes.ERROR_PREDICTION_FAILED[0]:
            self.block_box_prediction.destroy(True)
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Prediction failed because a subprocess failed!",
                "Structure Prediction",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
            self._interface_manager.update_status_bar(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}")
        elif tmp_exit_code == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
            self.block_box_prediction.destroy(True)
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Prediction failed because of an unknown error!",
                "Structure Prediction",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
            self._interface_manager.update_status_bar(f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}")
        elif tmp_exit_code == exit_codes.EXIT_CODE_ZERO[0]:
            # Prediction was successful
            self._interface_manager.refresh_protein_model()
            self._interface_manager.refresh_main_view()
            self.block_box_prediction.destroy(True)
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "All structure predictions are done. Go to View to the Proteins tab to see the new protein(s).",
                "Structure Prediction",
                custom_message_box.CustomMessageBoxIcons.INFORMATION.value
            )
            tmp_dialog.exec_()
            constants.PYSSA_LOGGER.info("All structure predictions are done.")
            self._interface_manager.update_status_bar("All structure predictions are done.")
        else:
            self.block_box_prediction.destroy(True)
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Prediction failed because of an unknown case!",
                "Structure Prediction",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
            self._interface_manager.update_status_bar("Prediction failed because of an unknown case!")
        self._interface_manager.stop_wait_spinner()

    # </editor-fold>

    # <editor-fold desc="Multimer">
    def _predict_multimer(self):
        self._external_controller = predict_multimer_view_controller.PredictMultimerViewController(
            self._interface_manager
        )
        self._external_controller.job_input.connect(self._post_predict_monomer)
        self._interface_manager.get_predict_multimer_view().show()

    def _post_predict_multimer(self, result: tuple):
        self._view.wait_spinner.start()

        # <editor-fold desc="Check if WSL2 and ColabFold are installed">
        if globals.g_os == "win32":
            constants.PYSSA_LOGGER.info("Checking if WSL2 is installed ...")
            if not dialog_settings_global.is_wsl2_installed():
                constants.PYSSA_LOGGER.warning("WSL2 is NOT installed.")
                self._interface_manager.get_application_settings().wsl_install = 0
                tmp_dialog = custom_message_box.CustomMessageBoxOk(
                    "Prediction failed because the WSL2 environment is not installed!",
                    "Structure Prediction",
                    custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
                )
                tmp_dialog.exec_()
                return
            constants.PYSSA_LOGGER.info("Checking if Local Colabfold is installed ...")
            if not dialog_settings_global.is_local_colabfold_installed():
                constants.PYSSA_LOGGER.warning("Local Colabfold is NOT installed.")
                self._interface_manager.get_application_settings().local_colabfold = 0
                tmp_dialog = custom_message_box.CustomMessageBoxOk(
                    "Prediction failed because the ColabFold is not installed!",
                    "Structure Prediction",
                    custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
                )
                tmp_dialog.exec_()
                return

        # </editor-fold>

        self.prediction_type = constants.PREDICTION_TYPE_PRED_MULTI_ANALYSIS
        constants.PYSSA_LOGGER.info("Begin prediction process.")
        self._interface_manager.update_status_bar("Begin prediction process ...")
        if result[3] is True:
            constants.PYSSA_LOGGER.info("Running prediction with subsequent analysis.")
            # Analysis should be run after the prediction
            self._active_task = tasks.Task(
                target=main_presenter_async.predict_protein_with_colabfold,
                args=(
                    result[1],
                    result[2],
                    self._interface_manager.get_current_project(),
                ),
                post_func=self.__await_multimer_prediction_for_subsequent_analysis,
            )
            self._active_task.start()
        else:
            constants.PYSSA_LOGGER.info("Running only a prediction.")
            # No analysis after prediction
            self._active_task = tasks.Task(
                target=main_presenter_async.predict_protein_with_colabfold,
                args=(
                    result[1],
                    result[2],
                    self._interface_manager.get_current_project(),
                ),
                post_func=self.__await_predict_protein_with_colabfold,
            )
            self._active_task.start()

        self._view.status_bar.showMessage("A prediction is currently running ...")
        self.block_box_prediction = QtWidgets.QMessageBox()
        self.block_box_prediction.setIcon(QtWidgets.QMessageBox.Information)
        self.block_box_prediction.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        styles.set_stylesheet(self.block_box_prediction)
        self.block_box_prediction.setWindowTitle("Structure Prediction")
        self.block_box_prediction.setText("A prediction is currently running.")
        btn_abort = self.block_box_prediction.addButton("Abort", QtWidgets.QMessageBox.ActionRole)
        self.block_box_prediction.exec_()
        if self.block_box_prediction.clickedButton() == btn_abort:
            self.abort_prediction()
            self.block_box_prediction.close()
            self._view.wait_spinner.stop()
        else:
            self.block_box_prediction.close()
            self._view.wait_spinner.stop()

    def __await_multimer_prediction_for_subsequent_analysis(self, result: tuple) -> None:
        tmp_exit_code = result[0]
        tmp_exit_code_description = [1]
        if tmp_exit_code == exit_codes.EXIT_CODE_ZERO[0]:
            # Prediction was successful
            self.block_box_prediction.destroy(True)
            constants.PYSSA_LOGGER.info("All structure predictions are done.")
            self.update_status("All structure predictions are done.")
            constants.PYSSA_LOGGER.info("Begin analysis process.")
            self.update_status("Begin analysis process ...")
            tmp_raw_analysis_run_names: list = []
            for row_no in range(self._view.ui.list_pred_analysis_multi_overview.count()):
                tmp_raw_analysis_run_names.append(self._view.ui.list_pred_analysis_multi_overview.item(row_no).text())

            self._active_task = tasks.Task(
                target=main_presenter_async.run_distance_analysis,
                args=(
                    tmp_raw_analysis_run_names,
                    self._interface_manager.get_current_project(),
                    self._interface_manager.get_application_settings(),
                    self._interface_manager.get_predict_multimer_view().cb_pred_analysis_multi_images.isChecked(),
                ),
                post_func=self.post_analysis_process,
            )
            self._active_task.start()

            if not os.path.exists(constants.SCRATCH_DIR_ANALYSIS):
                os.mkdir(constants.SCRATCH_DIR_ANALYSIS)

        elif tmp_exit_code == exit_codes.ERROR_WRITING_FASTA_FILES[0]:
            self.block_box_prediction.destroy(True)
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Prediction failed because there was an error writing the fasta file(s)!",
                "Structure Prediction",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
            self._interface_manager.update_status_bar(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}")
            self._view.wait_spinner.stop()
        elif tmp_exit_code == exit_codes.ERROR_FASTA_FILES_NOT_FOUND[0]:
            self.block_box_prediction.destroy(True)
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Prediction failed because the fasta file(s) could not be found!",
                "Structure Prediction",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
            self._interface_manager.update_status_bar(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}")
            self._view.wait_spinner.stop()
        elif tmp_exit_code == exit_codes.ERROR_PREDICTION_FAILED[0]:
            self.block_box_prediction.destroy(True)
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Prediction failed because a subprocess failed!",
                "Structure Prediction",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
            self._interface_manager.update_status_bar(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}")
            self._view.wait_spinner.stop()
        elif tmp_exit_code == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
            self.block_box_prediction.destroy(True)
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Prediction failed because of an unknown error!",
                "Structure Prediction",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
            constants.PYSSA_LOGGER.error(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
            )
            self._interface_manager.update_status_bar(
                f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}")
            self._view.wait_spinner.stop()


    # </editor-fold>

    # </editor-fold>

    # <editor-fold desc="Hotspots">
    def _hotspots_protein_regions(self) -> None:
        self._external_controller = hotspots_protein_regions_view_controller.HotspotsProteinRegionsViewController(
            self._interface_manager
        )
        self._interface_manager.get_hotspots_protein_regions_view().show()
        self._pymol_session_manager.show_sequence_view()

    def post_hotspots_protein_regions(self) -> None:
        self._pymol_session_manager.hide_sequence_view()
        cmd.select(name="", selection="none")

    # </editor-fold>

    # <editor-fold desc="Settings menu methods">
    def open_settings_global(self) -> None:
        """Opens the dialog for the global settings."""
        self._external_controller = settings_view_controller.SettingsViewController(self._interface_manager)
        self._external_controller.user_input.connect(self.post_open_settings_global)
        self._external_controller.restore_ui()
        self._interface_manager.get_settings_view().show()
        # dialog = dialog_settings_global.DialogSettingsGlobal(self._interface_manager)
        # dialog.exec_()
        # self._interface_manager.update_settings()
        # self._workspace_label = QtWidgets.QLabel(f"Current Workspace: {self._workspace_path}")

    def post_open_settings_global(self, return_value):
        try:
            tmp_type = self._interface_manager.get_current_protein_tree_index_type()
        except AttributeError:
            return
        if self._view.ui.cb_protein_atoms.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
            self._view.tg_protein_color_atoms.toggle_button.setChecked(True)
            self._view.ui.cb_protein_atoms.setChecked(False)
        elif self._view.tg_protein_color_atoms.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
            self._view.tg_protein_color_atoms.toggle_button.setChecked(False)
            self._view.ui.cb_protein_atoms.setChecked(True)
        else:
            self._view.tg_protein_color_atoms.toggle_button.setChecked(False)
            self._view.ui.cb_protein_atoms.setChecked(False)

        if tmp_type == "chain":
            self._interface_manager.set_index_of_protein_color_combo_box(self._pymol_session_manager)
            self._interface_manager.set_repr_state_in_ui_for_protein_chain(self._pymol_session_manager)
            self._interface_manager.show_protein_pymol_scene_configuration()

    def restore_settings(self) -> None:
        """Restores the settings.xml file to the default values."""
        tmp_dialog = custom_message_box.CustomMessageBoxYesNo(
            "Are you sure you want to restore all settings?", "Restore Settings",
            custom_message_box.CustomMessageBoxIcons.INFORMATION.value
        )
        tmp_dialog.exec_()
        if tmp_dialog.dialogClosed is not None:
            return
        if tmp_dialog.response:
            tools.restore_default_settings(self._interface_manager.get_application_settings())
            self._view.status_bar.showMessage("Settings were successfully restored.")
            logging.info("Settings were successfully restored.")
        else:
            self._view.status_bar.showMessage("Settings were not modified.")
            logging.info("Settings were not modified.")

    # </editor-fold>

    # <editor-fold desc="Help menu methods">
    def arrange_windows(self):
        if not os.path.exists(constants.ARRANGE_WINDOWS_EXE_FILEPATH):
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "The script for arranging the windows could not be found!", "Arrange Windows",
                custom_message_box.CustomMessageBoxIcons.ERROR.value
            )
            tmp_dialog.exec_()
        else:
            logger.debug("Started script to arrange window ...")
            subprocess.Popen([constants.ARRANGE_WINDOWS_EXE_FILEPATH])
            logger.debug("Script to arrange windows finished.")

    def open_logs(self) -> None:
        """Opens a file explorer with all log files and can open a log file in the default application."""
        file_dialog = QtWidgets.QFileDialog()
        log_path = str(constants.LOG_PATH)
        file_dialog.setDirectory(log_path)
        file_path, _ = file_dialog.getOpenFileName(self._view, "Select a log file to open", "", "LOG File (*.log)")
        if file_path:
            os.startfile(file_path)

    @staticmethod
    def clear_all_log_files() -> None:
        """Clears all log files generated under .pyssa/logs."""
        response = basic_boxes.yes_or_no(
            "Clear log files",
            "Are you sure you want to delete all log files?",
            QtWidgets.QMessageBox.Information,
        )
        if response:
            try:
                shutil.rmtree(str(constants.LOG_PATH))
            except PermissionError:
                print("The active log file was not deleted.")
            if len(os.listdir(str(constants.LOG_PATH))) == 1:
                tmp_dialog = custom_message_box.CustomMessageBoxOk(
                    "All log files could be deleted.", "Clear Log Files",
                    custom_message_box.CustomMessageBoxIcons.INFORMATION.value
                )
                tmp_dialog.exec_()
                constants.PYSSA_LOGGER.info("All log files were deleted.")
            else:
                tmp_dialog = custom_message_box.CustomMessageBoxOk(
                    "Not all log files could be deleted.",
                    "Clear Log Files",
                    custom_message_box.CustomMessageBoxIcons.WARNING.value
                )
                tmp_dialog.exec_()
                constants.PYSSA_LOGGER.warning("Not all log files were deleted!")

    @staticmethod
    def open_tutorial() -> None:
        """Opens the official tutorial pdf file."""
        tmp_dialog = dialog_tutorial_videos.TutorialVideosDialog()
        tmp_dialog.exec_()

    @staticmethod
    def open_documentation() -> None:
        """Opens the official plugin documentation as PDF."""
        os.startfile(constants.DOCS_PATH)

    @staticmethod
    def open_about() -> None:
        """Opens the About dialog."""
        dialog = dialog_about.DialogAbout()
        dialog.exec_()

    def get_demo_projects(self):
        self._interface_manager.update_status_bar("Getting demo projects ...")
        import zipfile
        download_dest = pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects.zip")
        if not os.path.exists(download_dest):
            # download demo projects
            url = f'https://w-hs.sciebo.de/s/ZHJa6XB9SKWtqGi/download'
            tmp_error_flag = False
            try:
                response = requests.get(url)
                response.raise_for_status()  # Check for errors
                zipfile = zipfile.ZipFile(BytesIO(response.content))
                zipfile.extractall(pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects"))
            except requests.exceptions.HTTPError as errh:
                constants.PYSSA_LOGGER.error(f"HTTP Error: {errh}")
                tmp_error_flag = True
            except requests.exceptions.ConnectionError as errc:
                constants.PYSSA_LOGGER.error(f"Error Connecting: {errc}")
                tmp_error_flag = True
            except requests.exceptions.Timeout as errt:
                constants.PYSSA_LOGGER.error(f"Timeout Error: {errt}")
                tmp_error_flag = True
            except requests.exceptions.RequestException as err:
                constants.PYSSA_LOGGER.error(f"Error: {err}")
                tmp_error_flag = True
            else:
                constants.PYSSA_LOGGER.info(f"Demo projects downloaded and extracted successfully.")

            if tmp_error_flag:
                tmp_dialog = custom_message_box.CustomMessageBoxOk(
                    "The download of the demo projects failed. Please try again later.",
                    "Get Demo Projects",
                    custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
                )
                tmp_dialog.exec_()
                self._interface_manager.update_status_bar("The download of the demo projects failed.")
                return
        else:
            constants.PYSSA_LOGGER.info("Demo projects are getting extracted ...")
            try:
                with zipfile.ZipFile(pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects.zip"), "r") as zip_ref:
                    zip_ref.extractall(pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects"))
                constants.PYSSA_LOGGER.info(
                    "Demo projects are downloaded and extracted.\n Import of demo projects started ...",
                )
            except Exception as e:
                constants.PYSSA_LOGGER.error(f"Extraction process of demo projects finished with the error: {e}.")
                tmp_dialog = custom_message_box.CustomMessageBoxOk(
                    "Extraction process of demo projects finished with an error. Check the logs to get more information.",
                    "Get Demo Projects",
                    custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
                )
                tmp_dialog.exec_()
                self._interface_manager.update_status_bar("Extraction process of demo projects finished with an error.")
                return
        try:
            path_of_demo_projects = pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects")
            for tmp_filename in os.listdir(path_of_demo_projects):
                # Copy db file into new workspace
                tmp_project_database_filepath = str(
                    pathlib.Path(
                        f"{self._interface_manager.get_application_settings().workspace_path}/{tmp_filename}"
                    )
                )
                tmp_src_filepath = str(pathlib.Path(f"{path_of_demo_projects}/{tmp_filename}"))
                shutil.copyfile(tmp_src_filepath, tmp_project_database_filepath)
            constants.PYSSA_LOGGER.info("Import process of demo projects finished.")
        except Exception as e:
            constants.PYSSA_LOGGER.error(f"Import process of demo projects finished with the error: {e}.")
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Import process of demo projects finished with an error. Check the logs to get more information.",
                "Get Demo Projects",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
            self._interface_manager.update_status_bar("Import process of demo projects finished with an error.")
        else:
            self._interface_manager.refresh_workspace_model()
            self._interface_manager.refresh_main_view()
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Getting demo projects finished successfully.", "Get Demo Projects",
                custom_message_box.CustomMessageBoxIcons.INFORMATION.value
            )
            tmp_dialog.exec_()
            self._interface_manager.update_status_bar("Getting demo projects finished successfully.")

    # </editor-fold>

    # <editor-fold desc="Image menu methods">
    # TODO: images need to be reimplemented
    # def post_preview_image(self) -> None:
    #     """Hides the block box of the preview process."""
    #     # self.block_box_uni.hide()
    #     # self.block_box_uni.destroy(True)
    #     self._view.status_bar.showMessage("Finished preview of ray-traced image.")
    #     QtWidgets.QApplication.restoreOverrideCursor()
    #
    # def preview_image(self) -> None:
    #     """Previews the image."""
    #     QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
    #     if self._view.ui.cb_ray_tracing.isChecked():
    #         self._view.status_bar.showMessage("Preview ray-traced image ...")
    #         # <editor-fold desc="Worker setup">
    #         # --Begin: worker setup
    #         self.tmp_thread = QtCore.QThread()
    #         self.tmp_worker = task_workers.PreviewRayImageWorker(self.renderer)
    #         self.tmp_thread = task_workers.setup_worker_for_work(
    #             self.tmp_thread,
    #             self.tmp_worker,
    #             self.display_view_page,
    #         )
    #         self.tmp_worker.finished.connect(self.post_preview_image)
    #         self.tmp_thread.start()
    #         # --End: worker setup
    #
    #         # </editor-fold>
    #         gui_utils.setup_standard_block_box(
    #             self.block_box_uni,
    #             "Preview ray-trace image",
    #             "Creating preview for the ray-traced image ...",
    #         )
    #         # self.block_box_uni.exec_()
    #     else:
    #         self._view.status_bar.showMessage("Preview draw image ...")
    #         cmd.draw(2400, 2400)
    #         self._view.status_bar.showMessage("Finished preview of drawn image.")
    #         QtWidgets.QApplication.restoreOverrideCursor()

    def _preview_image(self):
        self._active_task = tasks.Task(
            target=main_presenter_async.preview_image,
            args=(0, 0),
            post_func=self.__await_preview_image,
        )
        self._active_task.start()
        self.update_status("Creating preview of image ...")
        self._view.wait_spinner.start()

    def __await_preview_image(self, return_value: tuple):
        self._view.wait_spinner.stop()
        self.update_status("Preview finished.")

    def _create_ray_traced_image(self) -> None:
        save_dialog = QtWidgets.QFileDialog()
        full_file_name = save_dialog.getSaveFileName(caption="Save Image", filter="Image (*.png)")
        if full_file_name == ("", ""):
            tools.quick_log_and_display(
                "info",
                "No file has been selected.",
                self._view.status_bar,
                "No file has been selected.",
            )
            return

        self._active_task = tasks.Task(
            target=main_presenter_async.create_ray_traced_image,
            args=(full_file_name[0], self._interface_manager.get_application_settings()),
            post_func=self.__await_create_ray_traced_image,
        )
        self._active_task.start()
        self.update_status("Creating ray-traced image ...")
        self._view.wait_spinner.start()

    def __await_create_ray_traced_image(self, return_value: tuple) -> None:
        self._view.wait_spinner.stop()
        self.update_status("Image creation finished.")

    def _create_drawn_image(self) -> None:
        save_dialog = QtWidgets.QFileDialog()
        full_file_name = save_dialog.getSaveFileName(caption="Save Image", filter="Image (*.png)")
        if full_file_name == ("", ""):
            tools.quick_log_and_display(
                "info",
                "No file has been selected.",
                self._view.status_bar,
                "No file has been selected.",
            )
            return

        self._active_task = tasks.Task(
            target=main_presenter_async.create_drawn_image,
            args=(full_file_name[0], self._interface_manager.get_application_settings()),
            post_func=self.__await_create_drawn_image,
        )
        self._active_task.start()
        self.update_status("Creating simple image ...")
        self._view.wait_spinner.start()

    def __await_create_drawn_image(self, return_value: tuple) -> None:
        self._view.wait_spinner.stop()
        self.update_status("Image creation finished.")


    def post_save_image(self) -> None:
        """Displays a message box which informs that the process has finished."""
        self.block_box_uni.hide()
        self.block_box_uni.destroy(True)
        self._view.status_bar.showMessage("Finished image creation.")
        QtWidgets.QApplication.restoreOverrideCursor()
        tmp_dialog = custom_message_box.CustomMessageBoxOk(
            "The image has been created.", "Image Creation",
            custom_message_box.CustomMessageBoxIcons.INFORMATION.value
        )
        tmp_dialog.exec_()

    def save_image(self) -> None:
        """Saves the image as a png file."""
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        if self._view.ui.cb_ray_tracing.isChecked():
            save_dialog = QtWidgets.QFileDialog()
            try:
                full_file_name = save_dialog.getSaveFileName(caption="Save Image", filter="Image (*.png)")
                if full_file_name == ("", ""):
                    tools.quick_log_and_display(
                        "info",
                        "No file has been selected.",
                        self._view.status_bar,
                        "No file has been selected.",
                    )
                    return
                self._view.status_bar.showMessage("Creating ray-traced image ...")

                # <editor-fold desc="Worker setup">
                # --Begin: worker setup
                self.tmp_thread = QtCore.QThread()
                self.tmp_worker = task_workers.SaveRayImageWorker(self.renderer, full_file_name[0])
                self.tmp_thread = task_workers.setup_worker_for_work(
                    self.tmp_thread,
                    self.tmp_worker,
                    self.display_view_page,
                )
                self.tmp_worker.finished.connect(self.post_save_image)
                self.tmp_thread.start()
                # --End: worker setup

                # </editor-fold>
                gui_utils.setup_standard_block_box(
                    self.block_box_uni,
                    "Save ray-trace image",
                    "Creating the ray-traced image ...",
                )
                self.block_box_uni.exec_()

                # cmd.ray(2400, 2400, renderer=int(self.renderer))
                # cmd.png(full_file_name[0], dpi=300)

            except FileExistsError:
                tools.quick_log_and_display(
                    "error",
                    "File exists already.",
                    self._view.status_bar,
                    "File exists already.",
                )
            except pymol.CmdException:
                tools.quick_log_and_display(
                    "error",
                    "Unexpected Error from PyMOL while saving the " "an image",
                    self._view.status_bar,
                    "Unexpected Error from PyMOL",
                )
        else:
            save_dialog = QtWidgets.QFileDialog()
            try:
                full_file_name = save_dialog.getSaveFileName(caption="Save Image", filter="Image (*.png)")
                if full_file_name == ("", ""):
                    tools.quick_log_and_display(
                        "info",
                        "No file has been selected.",
                        self._view.status_bar,
                        "No file has been selected.",
                    )
                    return
                self._view.status_bar.showMessage("Creating draw image ...")
                cmd.draw(2400, 2400)
                cmd.png(full_file_name[0], dpi=300)
                self._view.status_bar.showMessage("Finished image creation.")
                tmp_dialog = custom_message_box.CustomMessageBoxOk(
                    "The image has been created.", "Image Creation",
                    custom_message_box.CustomMessageBoxIcons.INFORMATION.value
                )
                tmp_dialog.exec_()
            except FileExistsError:
                tools.quick_log_and_display(
                    "error",
                    "File exists already.",
                    self._view.status_bar,
                    "File exists already.",
                )
            except pymol.CmdException:
                tools.quick_log_and_display(
                    "error",
                    "Unexpected Error from PyMOL while saving the " "an image",
                    self._view.status_bar,
                    "Unexpected Error from PyMOL",
                )
            finally:
                QtWidgets.QApplication.restoreOverrideCursor()

    # </editor-fold>

    # <editor-fold desc="Sequences tab methods">
    def _open_text_editor_for_seq(self):
        self.tmp_txt_browser = QtWidgets.QTextBrowser()
        try:
            self.tmp_txt_browser.setText(
                self._view.ui.seqs_table_widget.currentItem().data(enums.ModelEnum.OBJECT_ROLE).seq
            )
        except AttributeError:
            return
        else:
            self.tmp_txt_browser.setWindowTitle("View Sequence")
            self.tmp_txt_browser.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
            self.tmp_txt_browser.resize(500, 150)
            self.tmp_txt_browser.show()

    def _set_new_sequence_name_in_table_item(self):
        try:
            tmp_new_seq = self._view.line_edit_seq_name.text()
            self._view.ui.seqs_table_widget.currentItem().setText(tmp_new_seq)
        except AttributeError:
            return

    def _rename_sequence(self):
        tmp_old_name = self._view.ui.seqs_list_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE).name
        tmp_new_name = self._view.ui.seqs_table_widget.currentItem().text()
        tmp_seq = self._view.ui.seqs_list_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE).seq
        self._view.ui.seqs_list_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE).name = tmp_new_name
        self._view.ui.seqs_list_view.model().setData(self._view.ui.seqs_list_view.currentIndex(), tmp_new_name, Qt.DisplayRole)
        tmp_database_operation = database_operation.DatabaseOperation(
            enums.SQLQueryType.UPDATE_SEQUENCE_NAME, (0, tmp_new_name, tmp_old_name, tmp_seq)
        )
        self._database_thread.put_database_operation_into_queue(tmp_database_operation)

    def _show_sequence_information(self):
        self._view.ui.seqs_table_widget.cellChanged.disconnect(self._rename_sequence)
        self._interface_manager.show_sequence_parameters(
            self._view.ui.seqs_list_view.currentIndex()
        )
        self._view.ui.seqs_table_widget.cellChanged.connect(self._rename_sequence)
        self._view.ui.btn_save_sequence.setEnabled(True)
        self._view.ui.btn_delete_sequence.setEnabled(True)

    def _import_sequence(self) -> None:
        self._interface_manager.get_import_sequence_view().return_value.connect(self._post_import_sequence)
        self._interface_manager.get_import_sequence_view().show()

    def _post_import_sequence(self, return_value: tuple):
        tmp_fasta_filepath, _ = return_value
        with open(tmp_fasta_filepath, "r") as handle:
            for tmp_record in SeqIO.parse(handle, "fasta"):
                # Append each SeqRecord object to the list
                self._interface_manager.get_current_project().sequences.append(tmp_record)
                self._database_manager.insert_new_sequence(tmp_record)
        self._interface_manager.refresh_sequence_model()
        self._interface_manager.refresh_main_view()

    def _add_sequence(self):
        self._external_controller = add_sequence_view_controller.AddSequenceViewController(self._interface_manager)
        self._external_controller.return_value.connect(self._post_add_sequence)
        self._external_controller.restore_default_view()
        self._interface_manager.get_add_sequence_view().show()

    def _post_add_sequence(self, return_value: tuple):
        logger.info(f"Adding new sequence {return_value[0]} with {return_value[1]} to the current project.")
        tmp_seq_name = return_value[0]
        tmp_sequence = return_value[1]
        tmp_seq_record = SeqRecord.SeqRecord(tmp_sequence, name=tmp_seq_name)
        self._interface_manager.get_current_project().sequences.append(tmp_seq_record)
        tmp_database_operation = database_operation.DatabaseOperation(
            enums.SQLQueryType.INSERT_NEW_SEQUENCE,
            (0, tmp_seq_record)
        )
        self._database_thread.put_database_operation_into_queue(tmp_database_operation)
        self._interface_manager.refresh_sequence_model()
        self._interface_manager.refresh_main_view()

    def _save_selected_sequence_as_fasta_file(self):
        self._view.wait_spinner.start()
        file_dialog = QtWidgets.QFileDialog()
        desktop_path = QtCore.QStandardPaths.standardLocations(QtCore.QStandardPaths.DesktopLocation)[0]
        file_dialog.setDirectory(desktop_path)
        file_path, _ = file_dialog.getSaveFileName(
            self._view,
            "Save Protein Sequence",
            "",
            "FASTA File (*.fasta)",
        )
        if file_path:
            tmp_seq_record = self._interface_manager.get_current_sequence_list_index_object()
            # pre-process seq record object for the SeqIO module
            if tmp_seq_record.id == "<unknown id>":
                tmp_seq_record.id = tmp_seq_record.name
            tmp_seq_record.seq = Seq(tmp_seq_record.seq)
            # defines the task to save the sequence as .fasta file
            self._active_task = tasks.Task(
                target=main_presenter_async.save_selected_protein_sequence_as_fasta_file,
                args=(
                    tmp_seq_record,
                    file_path,
                    self._interface_manager.get_current_project().get_database_filepath()
                ),
                post_func=self.__await_save_selected_sequence_as_fasta_file,
            )
            self._active_task.start()
        else:
            self._view.wait_spinner.stop()

    def __await_save_selected_sequence_as_fasta_file(self, result: tuple):
        self._view.wait_spinner.stop()
        if result[0] == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Saving the sequence as .fasta file failed!",
                "Save Protein Sequence",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
        elif result[0] == exit_codes.EXIT_CODE_ZERO[0]:
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "The sequence was successfully saved as .fasta file.",
                "Save Protein Sequence",
                custom_message_box.CustomMessageBoxIcons.INFORMATION.value
            )
            tmp_dialog.exec_()
        else:
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Saving the sequence as .fasta file failed with an unexpected error!",
                "Save Protein Sequence",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
        self._interface_manager.refresh_sequence_model()
        self._interface_manager.refresh_main_view()

    def _delete_selected_sequence(self):
        response: bool = gui_utils.warning_message_sequence_gets_deleted()
        if response:
            tmp_seq_record: "SeqRecord.SeqRecord" = self._interface_manager.get_current_sequence_list_index_object()
            tmp_database_operation = database_operation.DatabaseOperation(enums.SQLQueryType.DELETE_EXISTING_SEQUENCE,
                                                                          (0, tmp_seq_record))
            self._database_thread.put_database_operation_into_queue(tmp_database_operation)
            self._interface_manager.get_current_project().delete_specific_sequence(tmp_seq_record.name)
            self._interface_manager.refresh_sequence_model()
            self._interface_manager.refresh_main_view()
            # extra ui changes
            self._view.ui.seqs_table_widget.setRowCount(0)
            self._view.build_sequence_table()

    # </editor-fold>

    # <editor-fold desc="Proteins tab methods">
    def open_context_menu_for_proteins(self, position):
        indexes = self._view.ui.proteins_tree_view.selectedIndexes()
        if len(indexes) > 0:
            level = 0
            index = indexes[0]
            while index.parent().isValid():
                index = index.parent()
                level += 1
        else:
            return

        tmp_type = self._interface_manager.get_current_protein_tree_index_type()

        if tmp_type == "protein" or tmp_type == "chain":
            tmp_protein = self._interface_manager.get_current_active_protein_object()
        elif tmp_type == "header" or tmp_type == "scene":
            return
        else:
            logger.warning("Unknown object type occurred in Protein tab.")
            return
        tmp_is_protein_in_any_pair: bool = self._interface_manager.get_current_project().check_if_protein_is_in_any_protein_pair(
            tmp_protein.get_molecule_object()
        )
        self.protein_context_menu = QtWidgets.QMenu()
        if level == 0:
            self.proteins_context_menu_clean_action = self.protein_context_menu.addAction(self._view.tr("Clean selected protein"))
            self.proteins_context_menu_clean_action.triggered.connect(self.clean_protein_update)
            self.proteins_context_menu_rename_action = self.protein_context_menu.addAction(self._view.tr("Rename selected protein"))
            self.proteins_context_menu_rename_action.triggered.connect(self.rename_selected_protein_structure)

            # <editor-fold desc="Check if protein is in any protein pair">
            if tmp_is_protein_in_any_pair:
                self.proteins_context_menu_rename_action.setEnabled(False)
            else:
                self.proteins_context_menu_rename_action.setEnabled(True)

            # </editor-fold>

        elif level == 1:
            # header level
            pass
        elif level == 2:
            tmp_show_sequence_action = self.protein_context_menu.addAction(self._view.tr("Show sequence"))
            tmp_show_sequence_action.triggered.connect(self._show_protein_chain_sequence)

        self.protein_context_menu.exec_(self._view.ui.proteins_tree_view.viewport().mapToGlobal(position))

    def __slot_get_information_about_selected_object_in_protein_branch(self) -> None:
        tmp_type = self._interface_manager.get_current_protein_tree_index_type()

        if tmp_type == "protein":
            pass
        elif tmp_type == "scene":
            if self._pymol_session_manager.is_the_current_protein_in_session():
                tmp_scene_name = self._interface_manager.get_current_active_scene_name()
                self._pymol_session_manager.current_scene_name = tmp_scene_name
                self._view.ui.lbl_pymol_protein_scene.setText(f"PyMOL Scene: {tmp_scene_name}")
                self._pymol_session_manager.load_scene(tmp_scene_name)

        elif tmp_type == "chain":
            if self._pymol_session_manager.current_scene_name != "":
                self.set_icon_for_current_color()
                self._interface_manager.set_index_of_protein_color_combo_box(self._pymol_session_manager)
                self._interface_manager.set_repr_state_in_ui_for_protein_chain(self._pymol_session_manager)
        elif tmp_type == "header":
            pass

        else:
            logger.warning("Unknown object type occurred in Protein tab.")
            return

        self._interface_manager.manage_ui_of_protein_tab(
            tmp_type,
            self._interface_manager.get_current_project().check_if_protein_is_in_any_protein_pair(
                self._interface_manager.get_current_active_protein_object().get_molecule_object()
            ),
            self._pymol_session_manager
        )

    def _open_protein_pymol_session(self):
        tmp_protein: "protein.Protein" = self._interface_manager.get_current_active_protein_object()
        # fixme: I am no sure if the code below is needed
        # if not self._pymol_session_manager.is_the_current_session_empty():
        #     tmp_flag = True  # Session is NOT empty and needs reinitialization
        # else:
        #     tmp_flag = False  # Session is empty

        tmp_flag = False
        self._active_task = tasks.Task(
            target=main_presenter_async.load_protein_pymol_session,
            args=(
                tmp_protein,
                self._pymol_session_manager,
                tmp_flag
            ),
            post_func=self.__await_open_protein_pymol_session,
        )
        self._active_task.start()
        self._pymol_session_manager.current_scene_name = "base"
        self._interface_manager.start_wait_spinner()
        self._interface_manager.update_status_bar(f"Loading PyMOL session of {tmp_protein.get_molecule_object()} ...")

    def __await_open_protein_pymol_session(self, return_value: tuple):
        _, exit_boolean = return_value
        self._view.ui.action_protein_regions.setEnabled(False)
        if exit_boolean:
            self._view.cb_chain_color.setEnabled(True)
            self._view.cb_chain_representation.setEnabled(True)
            self._view.ui.action_protein_regions.setEnabled(True)
            self._view.ui.btn_create_protein_scene.setEnabled(True)
            self._view.ui.btn_update_protein_scene.setEnabled(True)
            self._view.ui.lbl_session_name.setText(f"Session Name: {self._pymol_session_manager.session_name}")
            self._view.ui.lbl_pymol_protein_scene.setText(f"PyMOL Scene: base")
            self._view.ui.lbl_info.setText("Please select a chain.")
            self._pymol_session_manager.get_all_scenes_in_current_session()
            logger.info("Successfully opened protein session.")
            self._interface_manager.update_status_bar("Loading the PyMOL session was successful.")
        else:
            logger.error("The protein name could not be found in the object list in PyMOL!")
            self._view.cb_chain_color.setEnabled(False)
            self._view.cb_chain_representation.setEnabled(False)
            self._view.ui.btn_create_protein_scene.setEnabled(False)
            self._view.ui.btn_update_protein_scene.setEnabled(False)
            self._interface_manager.update_status_bar("Loading the PyMOL session failed! Check out the log file to get more information.")
            self._view.ui.lbl_info.setText("Please load the PyMOL session of the selected protein.")
        self._interface_manager.stop_wait_spinner()

    def _change_chain_color_proteins(self) -> None:
        tmp_protein = self._interface_manager.get_current_active_protein_object()
        tmp_chain = self._interface_manager.get_current_active_chain_object()
        if self._interface_manager._settings_manager.settings.proteins_tab_use_combobox_for_colors == 1:
            tmp_color = self._view.ui.box_protein_color.currentText()
        else:
            tmp_color = self._view.ui.lbl_protein_current_color.text()

        if self._pymol_session_manager.session_object_type == "protein" and self._pymol_session_manager.session_name == tmp_protein.get_molecule_object():
            # Update pymol parameter in PyMOL
            tmp_protein.pymol_selection.set_selection_for_a_single_chain(tmp_chain.chain_letter)
            try:
                tmp_protein.pymol_selection.color_selection(tmp_color)
            except pymol.CmdException:
                logger.warning("No protein in session found. This can lead to more serious problems.")
            else:
                # Update pymol parameter in memory
                tmp_chain.pymol_parameters["chain_color"] = tmp_color
                # Update pymol parameter in database
                with database_manager.DatabaseManager(str(self._interface_manager.get_current_project().get_database_filepath())) as db_manager:
                    db_manager.open_project_database()
                    db_manager.update_protein_chain_color(tmp_chain.get_id(), tmp_color)
                    db_manager.close_project_database()
                self._update_scene()
                self._save_protein_pymol_session()
        else:
            logger.warning("The color of a protein chain could not be changed. This can be due to UI setup reasons.")

    def _change_chain_color_proteins_atoms(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object().chain_letter)
        cmd.color(color='atomic', selection=f"{tmp_selection.selection_string} and not elem C")

    def _change_chain_reset_proteins_atoms(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object().chain_letter)
        if self._interface_manager.get_settings_manager().settings.proteins_tab_use_combobox_for_colors == 1:
            tmp_color_name = self._view.ui.box_protein_color.currentText()
        else:
            tmp_color_name = self._view.ui.lbl_protein_current_color.text().strip()
        cmd.color(color=tmp_color_name, selection=f"{tmp_selection.selection_string}")

    # <editor-fold desc="Color Grid slot methods">
    def set_icon_for_current_color(self):
        color_index_functions = {
            "red": self.set_color_name_in_label_red,
            "tv_red": self.set_color_name_in_label_tv_red,
            "salmon": self.set_color_name_in_label_salmon,
            "raspberry": self.set_color_name_in_label_raspberry,
            "green": self.set_color_name_in_label_green,
            "tv_green": self.set_color_name_in_label_tv_green,
            "palegreen": self.set_color_name_in_label_palegreen,
            "forest": self.set_color_name_in_label_forest,
            "blue": self.set_color_name_in_label_blue,
            "tv_blue": self.set_color_name_in_label_tv_blue,
            "lightblue": self.set_color_name_in_label_lightblue,
            "skyblue": self.set_color_name_in_label_skyblue,
            "yellow": self.set_color_name_in_label_yellow,
            "tv_yellow": self.set_color_name_in_label_tv_yellow,
            "paleyellow": self.set_color_name_in_label_paleyellow,
            "sand": self.set_color_name_in_label_sand,
            "magenta": self.set_color_name_in_label_magenta,
            "purple": self.set_color_name_in_label_purple,
            "pink": self.set_color_name_in_label_pink,
            "hotpink": self.set_color_name_in_label_hotpink,
            "cyan": self.set_color_name_in_label_cyan,
            "aquamarine": self.set_color_name_in_label_aquamarine,
            "palecyan": self.set_color_name_in_label_palecyan,
            "teal": self.set_color_name_in_label_teal,
            "orange": self.set_color_name_in_label_orange,
            "tv_orange": self.set_color_name_in_label_tv_orange,
            "lightorange": self.set_color_name_in_label_lightorange,
            "olive": self.set_color_name_in_label_olive,
            "white": self.set_color_name_in_label_white,
            "grey70": self.set_color_name_in_label_grey_70,
            "grey30": self.set_color_name_in_label_grey_30,
            "black": self.set_color_name_in_label_black
        }
        tmp_protein = self._interface_manager.get_current_active_protein_object()
        tmp_chain = self._interface_manager.get_current_active_chain_object()
        tmp_protein.pymol_selection.set_selection_for_a_single_chain(tmp_chain.chain_letter)
        color_index_functions[
            tmp_chain.get_color(tmp_protein.pymol_selection.selection_string)
        ]()

    # Fixme: write function for uncheck colorbutton
    

    def set_color_name_in_label_red(self):
        self._view.ui.lbl_protein_current_color.setText("red")
        self._change_chain_color_proteins()
        self._view.color_grid.c_red.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_red.setIconSize(self._view.color_grid.c_red.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_tv_red(self):
        self._view.ui.lbl_protein_current_color.setText("tv_red")
        self._change_chain_color_proteins()
        self._view.color_grid.c_tv_red.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_tv_red.setIconSize(self._view.color_grid.c_tv_red.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_salmon(self):
        self._view.ui.lbl_protein_current_color.setText("salmon")
        self._change_chain_color_proteins()
        self._view.color_grid.c_salomon.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_salomon.setIconSize(self._view.color_grid.c_salomon.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_raspberry(self):
        self._view.ui.lbl_protein_current_color.setText("raspberry")
        self._change_chain_color_proteins()
        self._view.color_grid.c_raspberry.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_raspberry.setIconSize(self._view.color_grid.c_raspberry.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_green(self):
        self._view.ui.lbl_protein_current_color.setText("green")
        self._change_chain_color_proteins()
        self._view.color_grid.c_green.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_green.setIconSize(
            self._view.color_grid.c_green.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_tv_green(self):
        self._view.ui.lbl_protein_current_color.setText("tv_green")
        self._change_chain_color_proteins()
        self._view.color_grid.c_tv_green.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_tv_green.setIconSize(
            self._view.color_grid.c_tv_green.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_palegreen(self):
        self._view.ui.lbl_protein_current_color.setText("palegreen")
        self._change_chain_color_proteins()
        self._view.color_grid.c_palegreen.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_palegreen.setIconSize(
            self._view.color_grid.c_palegreen.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_forest(self):
        self._view.ui.lbl_protein_current_color.setText("forest")
        self._change_chain_color_proteins()
        self._view.color_grid.c_forest.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_forest.setIconSize(
            self._view.color_grid.c_forest.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_blue(self):
        self._view.ui.lbl_protein_current_color.setText("blue")
        self._change_chain_color_proteins()
        self._view.color_grid.c_blue.setIcon(QtGui.QIcon(
            ":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_blue.setIconSize(self._view.color_grid.c_blue.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_tv_blue(self):
        self._view.ui.lbl_protein_current_color.setText("tv_blue")
        self._change_chain_color_proteins()
        self._view.color_grid.c_tv_blue.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_tv_blue.setIconSize(
            self._view.color_grid.c_tv_blue.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_lightblue(self):
        self._view.ui.lbl_protein_current_color.setText("lightblue")
        self._change_chain_color_proteins()
        self._view.color_grid.c_lightblue.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_lightblue.setIconSize(
            self._view.color_grid.c_lightblue.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_skyblue(self):
        self._view.ui.lbl_protein_current_color.setText("skyblue")
        self._change_chain_color_proteins()
        self._view.color_grid.c_skyblue.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_skyblue.setIconSize(
            self._view.color_grid.c_skyblue.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_yellow(self):
        self._view.ui.lbl_protein_current_color.setText("yellow")
        self._change_chain_color_proteins()
        self._view.color_grid.c_yellow.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_yellow.setIconSize(
            self._view.color_grid.c_yellow.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_tv_yellow(self):
        self._view.ui.lbl_protein_current_color.setText("tv_yellow")
        self._change_chain_color_proteins()
        self._view.color_grid.c_tv_yellow.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_tv_yellow.setIconSize(
            self._view.color_grid.c_tv_yellow.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_paleyellow(self):
        self._view.ui.lbl_protein_current_color.setText("paleyellow")
        self._change_chain_color_proteins()
        self._view.color_grid.c_paleyellow.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_paleyellow.setIconSize(
            self._view.color_grid.c_paleyellow.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_sand(self):
        self._view.ui.lbl_protein_current_color.setText("sand")
        self._change_chain_color_proteins()
        self._view.color_grid.c_sand.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_sand.setIconSize(
            self._view.color_grid.c_sand.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_magenta(self):
        self._view.ui.lbl_protein_current_color.setText("magenta")
        self._change_chain_color_proteins()
        self._view.color_grid.c_magenta.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_magenta.setIconSize(
            self._view.color_grid.c_magenta.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_purple(self):
        self._view.ui.lbl_protein_current_color.setText("purple")
        self._change_chain_color_proteins()
        self._view.color_grid.c_purple.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_purple.setIconSize(
            self._view.color_grid.c_purple.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_pink(self):
        self._view.ui.lbl_protein_current_color.setText("pink")
        self._change_chain_color_proteins()
        self._view.color_grid.c_pink.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_pink.setIconSize(
            self._view.color_grid.c_pink.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_hotpink(self):
        self._view.ui.lbl_protein_current_color.setText("hotpink")
        self._change_chain_color_proteins()
        self._view.color_grid.c_hotpink.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_hotpink.setIconSize(
            self._view.color_grid.c_hotpink.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_cyan(self):
        self._view.ui.lbl_protein_current_color.setText("cyan")
        self._change_chain_color_proteins()
        self._view.color_grid.c_cyan.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_cyan.setIconSize(
            self._view.color_grid.c_cyan.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_aquamarine(self):
        self._view.ui.lbl_protein_current_color.setText("aquamarine")
        self._change_chain_color_proteins()
        self._view.color_grid.c_aquamarine.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_aquamarine.setIconSize(
            self._view.color_grid.c_aquamarine.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_palecyan(self):
        self._view.ui.lbl_protein_current_color.setText("palecyan")
        self._change_chain_color_proteins()
        self._view.color_grid.c_palecyan.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_palecyan.setIconSize(
            self._view.color_grid.c_palecyan.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_teal(self):
        self._view.ui.lbl_protein_current_color.setText("teal")
        self._change_chain_color_proteins()
        self._view.color_grid.c_teal.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_teal.setIconSize(
            self._view.color_grid.c_teal.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_orange(self):
        self._view.ui.lbl_protein_current_color.setText("orange")
        self._change_chain_color_proteins()
        self._view.color_grid.c_orange.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_orange.setIconSize(
            self._view.color_grid.c_orange.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_tv_orange(self):
        self._view.ui.lbl_protein_current_color.setText("tv_orange")
        self._change_chain_color_proteins()
        self._view.color_grid.c_tv_orange.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_tv_orange.setIconSize(
            self._view.color_grid.c_tv_orange.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_lightorange(self):
        self._view.ui.lbl_protein_current_color.setText("lightorange")
        self._change_chain_color_proteins()
        self._view.color_grid.c_lightorange.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_lightorange.setIconSize(
            self._view.color_grid.c_lightorange.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_olive(self):
        self._view.ui.lbl_protein_current_color.setText("olive")
        self._change_chain_color_proteins()
        self._view.color_grid.c_olive.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_olive.setIconSize(
            self._view.color_grid.c_olive.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_white(self):
        self._view.ui.lbl_protein_current_color.setText("white")
        self._change_chain_color_proteins()
        self._view.color_grid.c_white.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_white.setIconSize(
            self._view.color_grid.c_white.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_grey_70(self):
        self._view.ui.lbl_protein_current_color.setText("grey70")
        self._change_chain_color_proteins()
        self._view.color_grid.c_grey_70.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_grey_70.setIconSize(
            self._view.color_grid.c_grey_70.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_grey_30(self):
        self._view.ui.lbl_protein_current_color.setText("grey30")
        self._change_chain_color_proteins()
        self._view.color_grid.c_grey_30.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_grey_30.setIconSize(
            self._view.color_grid.c_grey_30.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_black(self):
        self._view.ui.lbl_protein_current_color.setText("black")
        self._change_chain_color_proteins()
        self._view.color_grid.c_black.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid.c_black.setIconSize(
            self._view.color_grid.c_black.icon().actualSize(QtCore.QSize(14, 14)))
    # </editor-fold>

    # <editor-fold desc="Representations">
    # cartoon
    def __slot_protein_chain_as_cartoon(self):
        if self._view.ui.cb_protein_cartoon.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
            self.__slot_show_protein_chain_as_cartoon()
        elif self._view.tg_protein_cartoon.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
            self.__slot_show_protein_chain_as_cartoon()
        else:
            self.__slot_hide_protein_chain_as_cartoon()
        self._update_scene()
        self._save_protein_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_chain()

    def __slot_show_protein_chain_as_cartoon(self):
        # TODO: should be refactored after the checkbox approach is the way to go
        tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object().chain_letter)
        tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.CARTOON.value)

    def __slot_hide_protein_chain_as_cartoon(self):
        # TODO: should be refactored after the checkbox approach is the way to go
        tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object().chain_letter)
        tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.CARTOON.value)

    # sticks
    def __slot_protein_chain_as_sticks(self):
        if self._view.ui.cb_protein_sticks.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
            self.__slot_show_protein_chain_as_sticks()
        elif self._view.tg_protein_sticks.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
            self.__slot_show_protein_chain_as_sticks()
        else:
            self.__slot_hide_protein_chain_as_sticks()
        self._update_scene()
        self._save_protein_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_chain()

    def __slot_show_protein_chain_as_sticks(self):
        # TODO: should be refactored after the checkbox approach is the way to go
        tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object().chain_letter)
        tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.STICKS.value)

    def __slot_hide_protein_chain_as_sticks(self):
        # TODO: should be refactored after the checkbox approach is the way to go
        tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object().chain_letter)
        tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.STICKS.value)

    # ribbon
    def __slot_protein_chain_as_ribbon(self):
        if self._view.ui.cb_protein_ribbon.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
            self.__slot_show_protein_chain_as_ribbon()
        elif self._view.tg_protein_ribbon.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
            self.__slot_show_protein_chain_as_ribbon()
        else:
            self.__slot_hide_protein_chain_as_ribbon()
        self._update_scene()
        self._save_protein_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_chain()

    def __slot_show_protein_chain_as_ribbon(self):
        # TODO: should be refactored after the checkbox approach is the way to go
        tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object().chain_letter)
        tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.RIBBON.value)

    def __slot_hide_protein_chain_as_ribbon(self):
        # TODO: should be refactored after the checkbox approach is the way to go
        tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object().chain_letter)
        tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.RIBBON.value)

    # others
    def __slot_protein_chain_as_lines(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object().chain_letter)
        if self._view.ui.cb_protein_lines.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
            tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.LINES.value)
        elif self._view.tg_protein_lines.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
            tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.LINES.value)
        else:
            tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.LINES.value)
        self._update_scene()
        self._save_protein_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_chain()

    def __slot_protein_chain_as_spheres(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object().chain_letter)
        if self._view.ui.cb_protein_spheres.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
            tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SPHERES.value)
        elif self._view.tg_protein_spheres.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
            tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SPHERES.value)
        else:
            tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.SPHERES.value)
        self._update_scene()
        self._save_protein_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_chain()

    def __slot_protein_chain_as_dots(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object().chain_letter)
        if self._view.ui.cb_protein_dots.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
            tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.DOTS.value)
        elif self._view.tg_protein_dots.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
            tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.DOTS.value)
        else:
            tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.DOTS.value)
        self._update_scene()
        self._save_protein_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_chain()

    def __slot_protein_chain_as_mesh(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object().chain_letter)
        if self._view.ui.cb_protein_mesh.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
            tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.MESH.value)
        elif self._view.tg_protein_mesh.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
            tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.MESH.value)
        else:
            tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.MESH.value)
        self._update_scene()
        self._save_protein_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_chain()

    def __slot_protein_chain_as_surface(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object().chain_letter)
        if self._view.ui.cb_protein_surface.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
            tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SURFACE.value)
        elif self._view.tg_protein_surface.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
            tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SURFACE.value)
        else:
            tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.SURFACE.value)
        self._update_scene()
        self._save_protein_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_chain()

    # all
    def __slot_hide_protein_chain_all(self):
        self._view.ui.cb_protein_cartoon.setChecked(False)
        self._view.ui.cb_protein_sticks.setChecked(False)
        self._view.ui.cb_protein_ribbon.setChecked(False)
        self._view.ui.cb_protein_lines.setChecked(False)
        self._view.ui.cb_protein_spheres.setChecked(False)
        self._view.ui.cb_protein_dots.setChecked(False)
        self._view.ui.cb_protein_mesh.setChecked(False)
        self._view.ui.cb_protein_surface.setChecked(False)
        self._view.tg_protein_cartoon.toggle_button.setChecked(False)
        self._view.tg_protein_sticks.toggle_button.setChecked(False)
        self._view.tg_protein_ribbon.toggle_button.setChecked(False)
        self._view.tg_protein_lines.toggle_button.setChecked(False)
        self._view.tg_protein_spheres.toggle_button.setChecked(False)
        self._view.tg_protein_dots.toggle_button.setChecked(False)
        self._view.tg_protein_mesh.toggle_button.setChecked(False)
        self._view.tg_protein_surface.toggle_button.setChecked(False)
        self.__slot_hide_protein_chain_as_cartoon()
        self.__slot_hide_protein_chain_as_sticks()
        self.__slot_hide_protein_chain_as_ribbon()
        self.__slot_protein_chain_as_lines()
        self.__slot_protein_chain_as_spheres()
        self.__slot_protein_chain_as_dots()
        self.__slot_protein_chain_as_mesh()
        self.__slot_protein_chain_as_surface()
        self._update_scene()
        self._save_protein_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_chain()

    # </editor-fold>

    # def _change_chain_representation_proteins(self) -> None:
    #     tmp_type = self._interface_manager.get_current_protein_tree_index_type()
    #     tmp_representation = self._view.cb_chain_representation.currentText()
    #     if tmp_type == "chain":
    #         tmp_protein: "protein.Protein" = self._interface_manager.get_parent_index_object_of_current_protein_tree_index()
    #         tmp_raw_chain = self._interface_manager.get_current_protein_tree_index_object()
    #         tmp_chain = tmp_protein.get_chain_by_letter(tmp_raw_chain.chain_letter)
    #     elif tmp_type == "protein":
    #         tmp_protein = self._interface_manager.get_current_protein_tree_index_object()
    #         tmp_chain = tmp_protein.chains[0]
    #     else:
    #         return
    #
    #     if self._pymol_session_manager.session_object_type == "protein" and self._pymol_session_manager.session_name == tmp_protein.get_molecule_object():
    #         # Update pymol parameter in PyMOL
    #         tmp_protein.pymol_selection.set_selection_for_a_single_chain(tmp_chain.chain_letter)
    #         try:
    #             tmp_protein.pymol_selection.change_representaion_of_selection(tmp_representation)
    #         except pymol.CmdException:
    #             logger.warning("No protein in session found. This can lead to more serious problems.")
    #         else:
    #             # Update pymol parameter in memory
    #             tmp_chain.pymol_parameters[enums.PymolParameterEnum.REPRESENTATION.value] = tmp_representation
    #             # Update pymol parameter in database
    #             with database_manager.DatabaseManager(
    #                     str(self._interface_manager.get_current_project().get_database_filepath())) as db_manager:
    #                 db_manager.open_project_database()
    #                 db_manager.update_protein_chain_representation(tmp_chain.get_id(), tmp_representation)
    #                 db_manager.close_project_database()
    #             self._save_protein_pymol_session()
    #     else:
    #         logger.warning("The representation of a protein chain could not be changed. This can be due to UI setup reasons.")

    def _import_protein_structure(self):
        self._external_controller = add_protein_view_controller.AddProteinViewController(self._interface_manager)
        self._external_controller.user_input.connect(self._post_import_protein_structure)
        self._external_controller.restore_ui()
        self._interface_manager.get_add_protein_view().show()

    def _post_import_protein_structure(self, return_value: tuple):
        # TODO: this function needs an async/await part
        tmp_protein_name, tmp_name_len = return_value
        if tmp_name_len == 4:
            self._active_task = tasks.Task(
                target=main_presenter_async.add_protein_from_pdb_to_project,
                args=(
                    tmp_protein_name,
                    self._database_manager,
                    self._interface_manager
                ),
                post_func=self.__await_post_import_protein_structure,
            )
            self._active_task.start()
            constants.PYSSA_LOGGER.info("Create project finished with protein from the PDB.")
        elif tmp_name_len > 0:
            # local pdb file as input
            self._active_task = tasks.Task(
                target=main_presenter_async.add_protein_from_local_filesystem_to_project,
                args=(
                    tmp_protein_name,
                    self._database_manager,
                    self._interface_manager
                ),
                post_func=self.__await_post_import_protein_structure,
            )
            self._active_task.start()
            constants.PYSSA_LOGGER.info("Create project finished with protein from local filesystem.")
        else:
            logger.warning("No protein object was created.")
            return
        self._interface_manager.update_status_bar("Importing protein structure ...")
        self._interface_manager.start_wait_spinner()

    def __await_post_import_protein_structure(self, return_value: tuple):
        tmp_protein: "protein.Protein" = return_value[1]
        self._interface_manager.get_current_project().add_existing_protein(tmp_protein)
        self._database_thread.put_database_operation_into_queue(
            database_operation.DatabaseOperation(enums.SQLQueryType.INSERT_NEW_PROTEIN,
                                                 (0, tmp_protein)))
        self._interface_manager.refresh_main_view()
        self._interface_manager.update_status_bar("Importing protein structure finished.")
        self._interface_manager.stop_wait_spinner()

    def _delete_protein(self):
        tmp_dialog = custom_message_box.CustomMessageBoxDelete(
            "Are you sure you want to delete this protein?", "Delete Protein",
            custom_message_box.CustomMessageBoxIcons.WARNING.value
        )
        tmp_dialog.exec_()
        response: bool = tmp_dialog.response
        if response:
            tmp_protein: "protein.Protein" = self._view.ui.proteins_tree_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE)
            tmp_database_operation = database_operation.DatabaseOperation(enums.SQLQueryType.DELETE_EXISTING_PROTEIN,
                                                                          (0, tmp_protein.get_id()))
            self._database_thread.put_database_operation_into_queue(tmp_database_operation)
            self._interface_manager.get_current_project().delete_specific_protein(tmp_protein.get_molecule_object())
            self._interface_manager.remove_protein_from_proteins_model()
            self._interface_manager.refresh_main_view()

    def _save_selected_protein_structure_as_pdb_file(self) -> None:
        """Saves selected protein as pdb file."""
        file_dialog = QtWidgets.QFileDialog()
        desktop_path = QtCore.QStandardPaths.standardLocations(QtCore.QStandardPaths.DesktopLocation)[0]
        file_dialog.setDirectory(desktop_path)
        file_path, _ = file_dialog.getSaveFileName(
            self._view,
            "Save Protein Structure",
            "",
            "Protein Data Bank File (*.pdb)",
        )
        if file_path:
            tmp_protein: "protein.Protein" = self._interface_manager.get_current_protein_tree_index_object()
            self._active_task = tasks.Task(
                target=main_presenter_async.save_selected_protein_structure_as_pdb_file,
                args=(
                    tmp_protein,
                    file_path,
                    self._interface_manager.get_current_project().get_database_filepath()
                ),
                post_func=self.__await_save_selected_protein_structure_as_pdb_file,
            )
            self._active_task.start()
            self._interface_manager.start_wait_spinner()
        else:
            self._interface_manager.stop_wait_spinner()

    def __await_save_selected_protein_structure_as_pdb_file(self, result: tuple) -> None:
        if result[0] == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Saving the protein as .pdb file failed!",
                "Save Protein Structure",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
        elif result[0] == exit_codes.EXIT_CODE_ZERO[0]:
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "The protein was successfully saved as .pdb file.",
                "Save Protein Structure",
                custom_message_box.CustomMessageBoxIcons.INFORMATION.value
            )
            tmp_dialog.exec_()
        else:
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Saving the protein as .pdb file failed with an unexpected error!",
                "Save Protein Structure",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
        self._interface_manager.refresh_protein_model()
        self._interface_manager.refresh_main_view()
        self._interface_manager.stop_wait_spinner()

    def clean_protein_update(self) -> None:
        """Cleans the selected protein structure."""
        tmp_dialog = custom_message_box.CustomMessageBoxYesNo(
            "Are you sure you want to clean this protein?\n" "This will remove all organic and solvent components!", "Clean Protein",
            custom_message_box.CustomMessageBoxIcons.WARNING.value
        )
        tmp_dialog.exec_()
        if tmp_dialog.response:
            self._active_task = tasks.Task(
                target=main_presenter_async.clean_protein_update,
                args=(
                    self._interface_manager.get_current_protein_tree_index_object(),
                    self._interface_manager.get_current_project().get_database_filepath(),
                ),
                post_func=self.__await_clean_protein_update,
            )
            self._active_task.start()
            self._interface_manager.start_wait_spinner()
            self.update_status("Cleaning protein ...")
        else:
            constants.PYSSA_LOGGER.info("No protein has been cleaned.")
            self._interface_manager.stop_wait_spinner()

    def __await_clean_protein_update(self) -> None:
        self.update_status("Cleaning protein finished.")
        self._interface_manager.stop_wait_spinner()

    def rename_selected_protein_structure(self) -> None:
        """Opens a new view to rename the selected protein."""
        self._external_controller = rename_protein_view_controller.RenameProteinViewController(self._interface_manager)
        self._external_controller.user_input.connect(self.post_rename_selected_protein_structure)
        self._interface_manager.get_rename_protein_view().show()

    def post_rename_selected_protein_structure(self, return_value: tuple) -> None:
        """Renames a selected protein structure."""
        if return_value[1] is True:
            self._active_task = tasks.Task(
                target=main_presenter_async.rename_selected_protein_structure,
                args=(
                    self._interface_manager.get_current_protein_tree_index_object(),
                    return_value[0],
                    self._interface_manager.get_current_project().get_database_filepath()
                ),
                post_func=self.__await_post_rename_selected_protein_structure,
            )
            self._active_task.start()
            self._interface_manager.start_wait_spinner()
        else:
            self._interface_manager.start_wait_spinner()

    def __await_post_rename_selected_protein_structure(self, result: tuple) -> None:
        self._view.ui.proteins_tree_view.model().setData(
            self._interface_manager.get_current_protein_tree_index(), result[1], enums.ModelEnum.OBJECT_ROLE
        )
        tmp_database_operation = database_operation.DatabaseOperation(
            enums.SQLQueryType.UPDATE_PYMOL_SESSION_PROTEIN,
            (
                0,
                self._view.ui.proteins_tree_view.model().data(
                    self._interface_manager.get_current_protein_tree_index(),
                    enums.ModelEnum.OBJECT_ROLE
                )
            )
        )
        self._database_thread.put_database_operation_into_queue(tmp_database_operation)
        self._interface_manager.refresh_protein_model()
        self._interface_manager.refresh_main_view()
        self._interface_manager.start_wait_spinner()

    def _show_protein_chain_sequence(self) -> None:
        self.tmp_txt_browser = QtWidgets.QTextBrowser()
        try:
            tmp_chain: "chain.Chain" = self._interface_manager.get_current_protein_tree_index_object()
            if tmp_chain.chain_sequence.sequence == "":
                self.tmp_txt_browser.setText(
                    "This chain is a non-protein chain."
                )
            else:
                self.tmp_txt_browser.setText(
                    tmp_chain.chain_sequence.sequence
                )
        except AttributeError:
            return
        else:
            self.tmp_txt_browser.setWindowTitle("View Protein Sequence")
            self.tmp_txt_browser.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
            self.tmp_txt_browser.resize(500, 150)
            self.tmp_txt_browser.show()

    @staticmethod
    def _update_scene() -> None:
        """Updates the current selected PyMOL scene."""
        cmd.scene(key="auto", action="update")

    def save_scene(self) -> None:
        """Saves the current view as a new PyMOL scene."""
        self._external_controller = add_scene_view_controller.AddSceneViewController(self._interface_manager)
        self._external_controller.user_input.connect(self.post_save_scene)
        self._interface_manager.get_add_scene_view().show()

    def post_save_scene(self, return_value: tuple):
        tmp_scene_name, _ = return_value
        cmd.scene(key=tmp_scene_name, action="append")

        if self._interface_manager.current_tab_index == 1:
            self._active_task = tasks.Task(
                target=main_presenter_async.save_protein_pymol_session_to_database,
                args=(
                    self._interface_manager,
                    0
                ),
                post_func=self.__await_save_scene_protein,
            )
            self._active_task.start()
            self._interface_manager.update_status_bar("Adding new scene to protein ...")
            self._interface_manager.start_wait_spinner()
            self._interface_manager.add_scene_to_proteins_model(tmp_scene_name)
        elif self._interface_manager.current_tab_index == 2:
            # The database thread cannot be used here because the session gets loaded again
            # before the new data is in the db
            self._active_task = tasks.Task(
                target=main_presenter_async.save_protein_pair_pymol_session_to_database,
                args=(
                    self._interface_manager,
                    0
                ),
                post_func=self.__await_save_scene_protein_pair,
            )
            self._active_task.start()
            self._interface_manager.update_status_bar("Adding new scene to protein pair ...")
            self._interface_manager.start_wait_spinner()
            self._interface_manager.add_scene_to_protein_pairs_model(tmp_scene_name)
        else:
            logger.warning("The current tab index is not for the proteins nor for the protein pairs tab?!")
            return

    def __await_save_scene_protein(self, return_value: tuple):
        _, exit_flag = return_value
        if exit_flag:
            self._interface_manager.refresh_main_view()
            self._interface_manager.update_status_bar("Adding new scene to protein finished.")
        else:
            self._interface_manager.update_status_bar("Adding new scene to protein failed!")
        self._interface_manager.stop_wait_spinner()

    def __await_save_scene_protein_pair(self, return_value: tuple):
        _, exit_flag = return_value
        if exit_flag:
            self._interface_manager.refresh_main_view()
            self._interface_manager.update_status_bar("Adding new scene to protein pair finished.")
        else:
            self._interface_manager.update_status_bar("Adding new scene to protein pair failed!")
        self._interface_manager.stop_wait_spinner()

    def delete_current_scene(self):
        tmp_dialog = custom_message_box.CustomMessageBoxDelete(
            "Are you sure you want to delete this scene?", "Delete PyMOL Scene",
            custom_message_box.CustomMessageBoxIcons.WARNING.value
        )
        tmp_dialog.exec_()
        response: bool = tmp_dialog.response
        if response:
            cmd.scene(key=self._pymol_session_manager.current_scene_name, action="clear")  # TODO: Does not work as expected!

            if self._interface_manager.current_tab_index == 1:
                self._save_protein_pymol_session()
                self._interface_manager.remove_scene_from_proteins_model(self._interface_manager.get_current_protein_tree_index())
                self._interface_manager.refresh_main_view()
            elif self._interface_manager.current_tab_index == 2:
                # The database thread cannot be used here because the session gets loaded again
                # before the new data is in the db
                self._active_task = tasks.Task(
                    target=main_presenter_async.save_protein_pair_pymol_session_to_database,
                    args=(
                        self._interface_manager,
                        0
                    ),
                    post_func=self.__await_delete_current_scene,
                )
                self._active_task.start()
                self._interface_manager.update_status_bar("Deleting selected scene ...")
                self._interface_manager.start_wait_spinner()
                self._interface_manager.remove_scene_from_protein_pairs_model(
                    self._interface_manager.get_current_protein_pair_tree_index()
                )
            else:
                logger.warning("The current tab index is not for the proteins nor for the protein pairs tab?!")
                return

    def __await_delete_current_scene(self, return_value: tuple):
        _, exit_flag = return_value
        if exit_flag:
            self._interface_manager.refresh_main_view()
            self._interface_manager.update_status_bar("Deleted the scene successfully.")
        else:
            self._interface_manager.update_status_bar("Deleting the scene failed!")
        self._interface_manager.stop_wait_spinner()

    def _save_protein_pymol_session(self):
        """Saves the session as base64 string and updates the database"""
        tmp_database_operation = database_operation.DatabaseOperation(
            enums.SQLQueryType.UPDATE_PYMOL_SESSION_PROTEIN,
            (0, self._interface_manager.get_current_active_protein_object())
        )
        self._database_thread.put_database_operation_into_queue(tmp_database_operation)

    # </editor-fold>

    # <editor-fold desc="Protein Pairs tab methods">
    def open_context_menu_for_protein_pairs(self, position):
        indexes = self._view.ui.protein_pairs_tree_view.selectedIndexes()
        if len(indexes) > 0:
            level = 0
            index = indexes[0]
            while index.parent().isValid():
                index = index.parent()
                level += 1
        else:
            return

        # elif tmp_type == "chain":
        #     tmp_protein = self._interface_manager.get_parent_index_object_of_current_protein_tree_index()
        # else:
        #     logger.warning("Unknown object type occurred in Protein tab.")
        #     return
        # tmp_is_protein_in_any_pair: bool = self._interface_manager.get_current_project().check_if_protein_is_in_any_protein_pair(
        #     tmp_protein.get_molecule_object()
        # )

        self.protein_pair_context_menu = QtWidgets.QMenu()
        if level == 0:
            # protein pair level
            self.protein_pair_context_open_results_summary_action = self.protein_pair_context_menu.addAction(self._view.tr("Open Results Summary"))
            self.protein_pair_context_open_results_summary_action.triggered.connect(self._results_summary)
            self.protein_pair_context_color_based_on_rmsd_action = self.protein_pair_context_menu.addAction(self._view.tr("Color By RMSD"))
            self.protein_pair_context_color_based_on_rmsd_action.triggered.connect(self._color_protein_pair_by_rmsd)

            if not self._pymol_session_manager.is_the_current_protein_pair_in_session():
                self.protein_pair_context_color_based_on_rmsd_action.setEnabled(False)
            else:
                self.protein_pair_context_color_based_on_rmsd_action.setEnabled(True)
        else:
            pass

        self.protein_pair_context_menu.exec_(self._view.ui.protein_pairs_tree_view.viewport().mapToGlobal(position))

    def __slot_get_information_about_selected_object_in_protein_pair_branch(self):
        tmp_type = self._interface_manager.get_current_protein_pair_tree_index_type()

        if tmp_type == "protein_pair":
            pass
        elif tmp_type == "protein":
            pass
        elif tmp_type == "scene":
            if self._pymol_session_manager.is_the_current_protein_pair_in_session():
                tmp_scene_name = self._interface_manager.get_current_active_scene_name_of_protein_pair()
                self._pymol_session_manager.current_scene_name = tmp_scene_name
                self._view.ui.lbl_pymol_protein_pair_scene.setText(f"PyMOL Scene: {tmp_scene_name}")
                self._pymol_session_manager.load_scene(tmp_scene_name)
        elif tmp_type == "chain":
            if self._pymol_session_manager.current_scene_name != "":
                self._interface_manager.show_chain_pymol_parameter_for_protein_pairs(self._pymol_session_manager)
                self._interface_manager.set_repr_state_in_ui_for_protein_pair_chain(self._pymol_session_manager)
        elif tmp_type == "header":
            pass
        else:
            logger.warning("Unknown object type occurred in Protein Pairs tab.")
            return
        self._interface_manager.manage_ui_of_protein_pairs_tab(tmp_type, self._pymol_session_manager)

    def _color_protein_pair_by_rmsd(self) -> None:
        """Colors the residues in 5 colors depending on their distance to the reference."""
        self._active_task = tasks.Task(
            target=main_presenter_async.color_protein_pair_by_rmsd_value,
            args=(
                self._interface_manager.get_current_protein_pair_tree_index_object(),
                0
            ),
            post_func=self.__await_color_protein_pair_by_rmsd,
        )
        self._active_task.start()
        self._interface_manager.start_wait_spinner()

    def __await_color_protein_pair_by_rmsd(self, result: tuple) -> None:
        self._interface_manager.stop_wait_spinner()

    def _open_protein_pair_pymol_session(self):
        tmp_protein_pair: "protein_pair.ProteinPair" = self._interface_manager.get_current_active_protein_pair_object()
        # fixme: I am no sure if the code below is needed
        # if not self._pymol_session_manager.is_the_current_session_empty():
        #     tmp_flag = True  # Session is NOT empty and needs reinitialization
        # else:
        #     tmp_flag = False  # Session is empty

        tmp_flag = False
        self._active_task = tasks.Task(
            target=main_presenter_async.load_protein_pair_pymol_session,
            args=(
                tmp_protein_pair,
                self._pymol_session_manager,
                tmp_flag
            ),
            post_func=self.__await_open_protein_pair_pymol_session,
        )
        self._active_task.start()
        self._interface_manager.start_wait_spinner()
        self._interface_manager.update_status_bar(f"Loading PyMOL session of {tmp_protein_pair.name} ...")

    def __await_open_protein_pair_pymol_session(self, return_value: tuple):
        _, exit_boolean = return_value
        self._view.ui.action_protein_regions.setEnabled(False)
        if exit_boolean:
            self._view.ui.action_protein_regions.setEnabled(True)
            self._view.ui.btn_create_protein_pair_scene.setEnabled(True)
            self._view.ui.btn_update_protein_pair_scene.setEnabled(True)
            self._view.ui.lbl_session_name.setText(f"Session Name: {self._pymol_session_manager.session_name}")
            tmp_protein_pair = self._interface_manager.get_current_active_protein_pair_object()
            self._pymol_session_manager.current_scene_name = tmp_protein_pair.name
            self._view.ui.lbl_pymol_protein_pair_scene.setText(f"PyMOL Scene: {tmp_protein_pair.name}")
            logger.info("Successfully opened protein pair session.")
            self._interface_manager.update_status_bar("Loading the PyMOL session was successful.")
            self._view.ui.lbl_info_protein_pair.setText("Please select a chain.")
        else:
            logger.error("The protein name could not be found in the object list in PyMOL!")
            self._view.ui.btn_create_protein_pair_scene.setEnabled(False)
            self._view.ui.btn_update_protein_pair_scene.setEnabled(False)
            self._interface_manager.update_status_bar(
                "Loading the PyMOL session failed! Check out the log file to get more information.")
            self._view.ui.lbl_info_protein_pair.setText("Please load the PyMOL session of the selected protein.")
        self._interface_manager.stop_wait_spinner()

    def _change_chain_color_protein_pairs(self) -> None:
        tmp_protein_pair = self._interface_manager.get_current_active_protein_pair_object()
        tmp_protein = self._interface_manager.get_current_active_protein_object_of_protein_pair()
        tmp_chain = self._interface_manager.get_current_active_chain_object_of_protein_pair()
        tmp_color: str = self._view.ui.box_protein_pair_color.currentText()
        if self._pymol_session_manager.session_object_type == "protein_pair" and self._pymol_session_manager.session_name == tmp_protein_pair.name:
            # Update pymol parameter in PyMOL
            tmp_protein.pymol_selection.set_selection_for_a_single_chain(tmp_chain.chain_letter)
            try:
                tmp_protein.pymol_selection.color_selection(tmp_color)
            except pymol.CmdException:
                logger.warning("No protein in session found. This can lead to more serious problems.")
            else:
                self._update_scene()
                self._save_protein_pair_pymol_session()
        else:
            logger.warning("The color of a protein chain could not be changed. This can be due to UI setup reasons.")

    def _change_chain_color_proteins_pair_atoms(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
        cmd.color(color='atomic',selection=f"{tmp_selection.selection_string} and not elem C")

    def _change_chain_reset_proteins_pair_atoms(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
        cmd.color(color=self._view.ui.box_protein_pair_color.currentText(), selection=f"{tmp_selection.selection_string}")

    # <editor-fold desc="Representations">
    # cartoon
    def __slot_protein_pair_chain_as_cartoon(self):
        if self._view.ui.cb_protein_pair_cartoon.isChecked():
            self.__slot_show_protein_pair_chain_as_cartoon()
        else:
            self.__slot_hide_protein_pair_chain_as_cartoon()
        self._update_scene()
        self._save_protein_pair_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()

    def __slot_show_protein_pair_chain_as_cartoon(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
        tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.CARTOON.value)

    def __slot_hide_protein_pair_chain_as_cartoon(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
        tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.CARTOON.value)

    # sticks
    def __slot_protein_pair_chain_as_sticks(self):
        if self._view.ui.cb_protein_pair_sticks.isChecked():
            self.__slot_show_protein_pair_chain_as_sticks()
        else:
            self.__slot_hide_protein_pair_chain_as_sticks()
        self._update_scene()
        self._save_protein_pair_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()

    def __slot_show_protein_pair_chain_as_sticks(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
        tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.STICKS.value)

    def __slot_hide_protein_pair_chain_as_sticks(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
        tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.STICKS.value)

    # ribbon
    def __slot_protein_pair_chain_as_ribbon(self):
        if self._view.ui.cb_protein_pair_ribbon.isChecked():
            self.__slot_show_protein_pair_chain_as_ribbon()
        else:
            self.__slot_hide_protein_pair_chain_as_ribbon()
        self._update_scene()
        self._save_protein_pair_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()

    def __slot_show_protein_pair_chain_as_ribbon(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
        tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.RIBBON.value)

    def __slot_hide_protein_pair_chain_as_ribbon(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
        tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.RIBBON.value)

    # others
    def __slot_protein_pair_chain_as_lines(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
        if self._view.ui.cb_protein_pair_lines.isChecked():
            tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.LINES.value)
        else:
            tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.LINES.value)
        self._update_scene()
        self._save_protein_pair_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()

    def __slot_protein_pair_chain_as_spheres(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
        if self._view.ui.cb_protein_pair_spheres.isChecked():
            tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SPHERES.value)
        else:
            tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.SPHERES.value)
        self._update_scene()
        self._save_protein_pair_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()

    def __slot_protein_pair_chain_as_dots(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
        if self._view.ui.cb_protein_pair_dots.isChecked():
            tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.DOTS.value)
        else:
            tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.DOTS.value)
        self._update_scene()
        self._save_protein_pair_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()

    def __slot_protein_pair_chain_as_mesh(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
        if self._view.ui.cb_protein_pair_mesh.isChecked():
            tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.MESH.value)
        else:
            tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.MESH.value)
        self._update_scene()
        self._save_protein_pair_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()

    def __slot_protein_pair_chain_as_surface(self):
        tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
        if self._view.ui.cb_protein_pair_surface.isChecked():
            tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SURFACE.value)
        else:
            tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.SURFACE.value)
        self._update_scene()
        self._save_protein_pair_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()

    # all
    def __slot_hide_protein_pair_chain_all(self):
        self._view.ui.cb_protein_pair_cartoon.setChecked(False)
        self._view.ui.cb_protein_pair_sticks.setChecked(False)
        self._view.ui.cb_protein_pair_ribbon.setChecked(False)
        self._view.ui.cb_protein_pair_lines.setChecked(False)
        self._view.ui.cb_protein_pair_spheres.setChecked(False)
        self._view.ui.cb_protein_pair_dots.setChecked(False)
        self._view.ui.cb_protein_pair_mesh.setChecked(False)
        self._view.ui.cb_protein_pair_surface.setChecked(False)
        self.__slot_hide_protein_pair_chain_as_cartoon()
        self.__slot_hide_protein_pair_chain_as_sticks()
        self.__slot_hide_protein_pair_chain_as_ribbon()
        self.__slot_protein_pair_chain_as_lines()
        self.__slot_protein_pair_chain_as_spheres()
        self.__slot_protein_pair_chain_as_dots()
        self.__slot_protein_pair_chain_as_mesh()
        self.__slot_protein_pair_chain_as_surface()
        self._update_scene()
        self._save_protein_pair_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()

    # </editor-fold>

    def _delete_protein_pair_from_project(self):
        tmp_dialog = custom_message_box.CustomMessageBoxDelete(
            "Are you sure you want to delete this protein pair?","Delete Protein Pair",
            custom_message_box.CustomMessageBoxIcons.WARNING.value
        )
        tmp_dialog.exec_()
        response: bool = tmp_dialog.response
        if response:
            tmp_protein_pair: "protein_pair.ProteinPair" = self._interface_manager.get_current_protein_pair_tree_index_object()
            tmp_database_operation = database_operation.DatabaseOperation(
                enums.SQLQueryType.DELETE_EXISTING_PROTEIN_PAIR, (0, tmp_protein_pair.get_id())
            )
            self._database_thread.put_database_operation_into_queue(tmp_database_operation)
            self._interface_manager.get_current_project().delete_specific_protein_pair(tmp_protein_pair.name)
            self._interface_manager.remove_protein_pair_from_protein_pairs_model()
            self._interface_manager.refresh_main_view()

    def _check_for_results(self) -> None:
        if self._view.ui.protein_pairs_tree_view.model().data(self._view.ui.protein_pairs_tree_view.currentIndex(), Qt.DisplayRole).find("_vs_") != -1:
            self._view.ui.action_results_summary.setEnabled(True)
        else:
            self._view.ui.action_results_summary.setEnabled(False)

    def _results_summary(self) -> None:
        tmp_protein_pair = self._view.ui.protein_pairs_tree_view.model().data(self._view.ui.protein_pairs_tree_view.currentIndex(),
                                                                              enums.ModelEnum.OBJECT_ROLE)
        self._external_controller = results_view_controller.ResultsViewController(self._interface_manager,
                                                                                  tmp_protein_pair,
                                                                                  self._pymol_session_manager)
        self._interface_manager.get_results_view().show()

    def _save_protein_pair_pymol_session(self):
        tmp_protein_pair = self._interface_manager.get_current_active_protein_pair_object()
        tmp_database_operation = database_operation.DatabaseOperation(
            enums.SQLQueryType.UPDATE_PYMOL_SESSION_PROTEIN_PAIR,
            (0, tmp_protein_pair.get_id(), tmp_protein_pair)
        )
        self._database_thread.put_database_operation_into_queue(tmp_database_operation)

    # </editor-fold>
