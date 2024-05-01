import copy
import logging
import os
import pathlib
import shutil
import subprocess
import pygetwindow
import requests
import pymol
from io import BytesIO
from Bio import SeqRecord
from Bio.Seq import Seq

from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from PyQt5 import QtGui

from pyssa.gui.ui.custom_context_menus import protein_tree_context_menu, protein_pair_tree_context_menu, \
    sequence_list_context_menu
from pyssa.gui.ui.custom_dialogs import custom_message_box
from pyssa.internal.thread.async_pyssa import util_async, custom_signals, project_async, image_async, \
    pymol_session_async, protein_async, sequence_async, protein_pair_async
from pyssa.controller import results_view_controller, rename_protein_view_controller, use_project_view_controller, \
    pymol_session_manager, predict_multimer_view_controller, \
    add_sequence_view_controller, add_scene_view_controller, add_protein_view_controller, settings_view_controller, \
    predict_protein_view_controller, import_sequence_view_controller, rename_sequence_view_controller, watcher
from pyssa.internal.data_structures import chain
from pyssa.internal.data_structures.data_classes import main_view_state
from pyssa.gui.ui.dialogs import dialog_settings_global, dialog_tutorial_videos, dialog_about
from pyssa.internal.data_structures import project, settings, protein, protein_pair
from pyssa.internal.data_structures.data_classes import database_operation
from pyssa.internal.thread import tasks, task_workers, database_thread
from pyssa.io_pyssa import safeguard, filesystem_io
from pyssa.logging_pyssa import log_handlers, log_levels
from pyssa.util import constants, enums, exit_codes, tools, ui_util
from pyssa.gui.ui.views import main_view
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
    The active task which controls the help window of the application.
    """
    _help_task: tasks.Task

    """
    a thread for database related processes
    """
    _database_thread: "database_thread.DatabaseThread"

    def __init__(self,
                 the_interface_manager: "interface_manager.InterfaceManager") -> None:
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
        self._database_manager = database_manager.DatabaseManager("")
        self._database_manager.set_application_settings(self._interface_manager.get_application_settings())
        self._database_thread: "database_thread.DatabaseThread" = database_thread.DatabaseThread("")
        self._external_view = None
        self.active_custom_message_box: "custom_message_box.CustomMessageBoxOk" = None
        self._main_view_state = main_view_state.MainViewState(
            self._view.ui.seqs_list_view,
            self.__slot_show_sequence_information,
            self._view.ui.proteins_tree_view,
            self.__slot_get_information_about_selected_object_in_protein_branch,
            self._view.ui.protein_pairs_tree_view,
            self.__slot_get_information_about_selected_object_in_protein_pair_branch
        )
        self.custom_progress_signal = custom_signals.ProgressSignal()
        self.abort_signal = custom_signals.AbortSignal()
        self.disable_pymol_signal = custom_signals.DisablePyMOLSignal()
        self._sequence_list_context_menu = sequence_list_context_menu.SequenceListContextMenu()
        self._protein_tree_context_menu = protein_tree_context_menu.ProteinTreeContextMenu()
        self._protein_pair_tree_context_menu = protein_pair_tree_context_menu.ProteinPairTreeContextMenu()

        self._setup_statusbar()
        self._init_generic_help_context_menus()
        self._interface_manager.refresh_main_view()
        self._connect_all_ui_elements_with_slot_functions()
        if self._interface_manager.get_application_settings().start_help_at_startup == 1:
            self._start_documentation_server()

    def _connect_all_ui_elements_with_slot_functions(self):
        self._view.dialogClosed.connect(self._close_main_window)
        self.custom_progress_signal.progress.connect(self._update_progress_bar)
        self.abort_signal.abort.connect(self._abort_task)
        self._interface_manager.refresh_after_job_finished_signal.refresh.connect(self._update_main_view_ui)
        self._view.btn_open_job_overview.clicked.connect(self.__slot_open_job_overview_panel)
        self._view.btn_open_job_notification.clicked.connect(self.__slot_open_notification_panel)

        # <editor-fold desc="Menu">
        self._view.ui.action_new_project.triggered.connect(self.__slot_create_project)
        self._interface_manager.get_create_view().dialogClosed.connect(self.__await_create_project)
        self._view.ui.action_open_project.triggered.connect(self.__slot_open_project)
        self._interface_manager.get_open_view().dialogClosed.connect(self._post_open_project)
        self._view.ui.action_use_project.triggered.connect(self.__slot_use_project)
        self._view.ui.action_delete_project.triggered.connect(self.__slot_delete_project)
        self._interface_manager.get_delete_view().dialogClosed.connect(self._post_delete_project)
        self._view.ui.action_import_project.triggered.connect(self.__slot_import_project)
        self._view.ui.action_export_project.triggered.connect(self.__slot_export_current_project)
        self._view.ui.action_close_project.triggered.connect(self.__slot_close_project)
        self._view.ui.action_exit_application.triggered.connect(self.__slot_close_all)

        self._view.ui.action_results_summary.triggered.connect(self.__slot_results_summary)
        self._view.ui.action_preview_image.triggered.connect(self.__slot_preview_image)
        self._view.ui.action_ray_tracing_image.triggered.connect(self.__slot_create_ray_traced_image)
        self._view.ui.action_simple_image.triggered.connect(self.__slot_create_drawn_image)
        self._view.ui.action_protein_regions.triggered.connect(self.__slot_hotspots_protein_regions)
        self._interface_manager.get_hotspots_protein_regions_view().dialogClosed.connect(self.post_hotspots_protein_regions)

        self._view.ui.action_edit_settings.triggered.connect(self.__slot_open_settings_global)
        self._view.ui.action_restore_settings.triggered.connect(self.__slot_restore_settings)
        self._view.ui.action_show_log_in_explorer.triggered.connect(self.__slot_open_logs)
        self._view.ui.action_clear_logs.triggered.connect(self.__slot_clear_all_log_files)
        self._view.ui.action_documentation.triggered.connect(self.__slot_open_help_center)
        self._view.ui.action_tutorials.triggered.connect(self.__slot_open_tutorial)
        self._view.ui.action_get_demo_projects.triggered.connect(self.__slot_get_demo_projects)
        self._view.ui.action_about.triggered.connect(self.__slot_open_about)
        self._view.ui.action_predict_monomer.triggered.connect(self.__slot_predict_monomer)
        self._view.ui.action_predict_multimer.triggered.connect(self.__slot_predict_multimer)
        self._view.ui.action_abort_prediction.triggered.connect(self.__slot_abort_prediction)
        self._view.ui.action_distance_analysis.triggered.connect(self.__slot_distance_analysis)
        self._view.ui.action_arrange_windows.triggered.connect(self.__slot_arrange_windows)

        self._view.ui.project_tab_widget.currentChanged.connect(self._update_tab)
        # </editor-fold>

        # <editor-fold desc="Sequence Tab">
        self._view.ui.seqs_list_view.customContextMenuRequested.connect(self.open_context_menu_for_sequences)
        self._view.ui.seqs_list_view.clicked.connect(self.__slot_show_sequence_information)
        self._view.ui.btn_add_sequence.clicked.connect(self.__slot_add_sequence)
        self._view.ui.btn_import_seq.clicked.connect(self.__slot_import_sequence)
        self._view.ui.btn_save_sequence.clicked.connect(self.__slot_save_selected_sequence_as_fasta_file)
        self._view.ui.btn_delete_sequence.clicked.connect(self.__slot_delete_selected_sequence)
        self._view.ui.seqs_table_widget.cellClicked.connect(self.__slot_open_text_editor_for_seq)
        self._view.line_edit_seq_name.textChanged.connect(self._set_new_sequence_name_in_table_item)
        #self._view.ui.seqs_table_widget.cellChanged.connect(self._rename_sequence)
        self._view.ui.btn_help.clicked.connect(self.__slot_open_sequences_tab_help)

        # <editor-fold desc="Context menu">
        self._sequence_list_context_menu.connect_rename_sequence_action(self.__slot_rename_selected_sequence)
        self._sequence_list_context_menu.connect_help_action(self.__slot_open_sequences_tab_help)
        # </editor-fold>

        # </editor-fold>

        # <editor-fold desc="Proteins Tab">
        self._view.ui.proteins_tree_view.customContextMenuRequested.connect(self.open_context_menu_for_proteins)
        self._view.ui.proteins_tree_view.clicked.connect(self.__slot_get_information_about_selected_object_in_protein_branch)
        self._view.ui.btn_save_protein.clicked.connect(self.__slot_save_selected_protein_structure_as_pdb_file)
        # import
        self._view.ui.btn_import_protein.clicked.connect(self.__slot_import_protein_structure)
        self._interface_manager.get_add_protein_view().return_value.connect(self._post_import_protein_structure)
        self._view.ui.btn_open_protein_session.clicked.connect(self.__slot_open_protein_pymol_session)
        self._view.ui.btn_create_protein_scene.clicked.connect(self.__slot_save_scene)
        self._view.ui.btn_delete_protein.clicked.connect(self.__slot_delete_protein)
        self._view.ui.btn_update_protein_scene.clicked.connect(self.__slot_update_protein_scene)
        self._view.ui.btn_delete_protein_scene.clicked.connect(self.__slot_delete_current_scene)
        self._view.ui.box_protein_color.currentIndexChanged.connect(self.__slot_change_chain_color_proteins)
        self._view.ui.btn_protein_color_atoms.clicked.connect(self.__slot_change_chain_color_proteins_atoms)
        self._view.ui.btn_protein_reset_atoms.clicked.connect(self.__slot_change_chain_reset_proteins_atoms)
        self._view.tg_protein_white_bg.toggleChanged.connect(self.__slot_protein_change_background_color)
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
        # representation
        # self._view.ui.btn_protein_show_hydrogens.clicked.connect(self.__slot_show_protein_chain_with_hydrogens)
        # self._view.ui.btn_protein_hide_hydrogens.clicked.connect(self.__slot_hide_protein_chain_with_hydrogens)
        self._view.tg_protein_hydrogen_atoms.toggleChanged.connect(self.__slot_protein_chain_with_hydrogens)
        self._view.tg_protein_cartoon.toggleChanged.connect(self.__slot_protein_chain_as_cartoon)
        self._view.tg_protein_sticks.toggleChanged.connect(self.__slot_protein_chain_as_sticks)
        self._view.tg_protein_ribbon.toggleChanged.connect(self.__slot_protein_chain_as_ribbon)
        self._view.tg_protein_lines.toggleChanged.connect(self.__slot_protein_chain_as_lines)
        self._view.tg_protein_spheres.toggleChanged.connect(self.__slot_protein_chain_as_spheres)
        self._view.tg_protein_dots.toggleChanged.connect(self.__slot_protein_chain_as_dots)
        self._view.tg_protein_mesh.toggleChanged.connect(self.__slot_protein_chain_as_mesh)
        self._view.tg_protein_surface.toggleChanged.connect(self.__slot_protein_chain_as_surface)

        self._view.ui.btn_protein_hide_all_representations.clicked.connect(self.__slot_hide_protein_chain_all)
        self._view.ui.btn_help_2.clicked.connect(self.__slot_open_proteins_tab_help)

        # <editor-fold desc="Color Grid">
        self._view.color_grid_proteins.c_red.clicked.connect(self.set_color_name_in_label_red_in_proteins_tab)
        self._view.color_grid_proteins.c_tv_red.clicked.connect(self.set_color_name_in_label_tv_red_in_proteins_tab)
        self._view.color_grid_proteins.c_salomon.clicked.connect(self.set_color_name_in_label_salmon_in_proteins_tab)
        self._view.color_grid_proteins.c_raspberry.clicked.connect(self.set_color_name_in_label_raspberry_in_proteins_tab)

        self._view.color_grid_proteins.c_green.clicked.connect(self.set_color_name_in_label_green_in_proteins_tab)
        self._view.color_grid_proteins.c_tv_green.clicked.connect(self.set_color_name_in_label_tv_green_in_proteins_tab)
        self._view.color_grid_proteins.c_palegreen.clicked.connect(self.set_color_name_in_label_palegreen_in_proteins_tab)
        self._view.color_grid_proteins.c_forest.clicked.connect(self.set_color_name_in_label_forest_in_proteins_tab)

        self._view.color_grid_proteins.c_blue.clicked.connect(self.set_color_name_in_label_blue_in_proteins_tab)
        self._view.color_grid_proteins.c_tv_blue.clicked.connect(self.set_color_name_in_label_tv_blue_in_proteins_tab)
        self._view.color_grid_proteins.c_lightblue.clicked.connect(self.set_color_name_in_label_lightblue_in_proteins_tab)
        self._view.color_grid_proteins.c_skyblue.clicked.connect(self.set_color_name_in_label_skyblue_in_proteins_tab)

        self._view.color_grid_proteins.c_yellow.clicked.connect(self.set_color_name_in_label_yellow_in_proteins_tab)
        self._view.color_grid_proteins.c_tv_yellow.clicked.connect(self.set_color_name_in_label_tv_yellow_in_proteins_tab)
        self._view.color_grid_proteins.c_paleyellow.clicked.connect(self.set_color_name_in_label_paleyellow_in_proteins_tab)
        self._view.color_grid_proteins.c_sand.clicked.connect(self.set_color_name_in_label_sand_in_proteins_tab)

        self._view.color_grid_proteins.c_magenta.clicked.connect(self.set_color_name_in_label_magenta_in_proteins_tab)
        self._view.color_grid_proteins.c_purple.clicked.connect(self.set_color_name_in_label_purple_in_proteins_tab)
        self._view.color_grid_proteins.c_pink.clicked.connect(self.set_color_name_in_label_pink_in_proteins_tab)
        self._view.color_grid_proteins.c_hotpink.clicked.connect(self.set_color_name_in_label_hotpink_in_proteins_tab)

        self._view.color_grid_proteins.c_cyan.clicked.connect(self.set_color_name_in_label_cyan_in_proteins_tab)
        self._view.color_grid_proteins.c_aquamarine.clicked.connect(self.set_color_name_in_label_aquamarine_in_proteins_tab)
        self._view.color_grid_proteins.c_palecyan.clicked.connect(self.set_color_name_in_label_palecyan_in_proteins_tab)
        self._view.color_grid_proteins.c_teal.clicked.connect(self.set_color_name_in_label_teal_in_proteins_tab)

        self._view.color_grid_proteins.c_orange.clicked.connect(self.set_color_name_in_label_orange_in_proteins_tab)
        self._view.color_grid_proteins.c_tv_orange.clicked.connect(self.set_color_name_in_label_tv_orange_in_proteins_tab)
        self._view.color_grid_proteins.c_lightorange.clicked.connect(self.set_color_name_in_label_lightorange_in_proteins_tab)
        self._view.color_grid_proteins.c_olive.clicked.connect(self.set_color_name_in_label_olive_in_proteins_tab)

        self._view.color_grid_proteins.c_white.clicked.connect(self.set_color_name_in_label_white_in_proteins_tab)
        self._view.color_grid_proteins.c_grey_70.clicked.connect(self.set_color_name_in_label_grey_70_in_proteins_tab)
        self._view.color_grid_proteins.c_grey_30.clicked.connect(self.set_color_name_in_label_grey_30_in_proteins_tab)
        self._view.color_grid_proteins.c_black.clicked.connect(self.set_color_name_in_label_black_in_proteins_tab)
        # </editor-fold>

        # <editor-fold desc="Protein tree context menu">
        self._protein_tree_context_menu.connect_clean_protein_action(self.__slot_clean_protein_update)
        self._protein_tree_context_menu.connect_rename_protein_action(self.__slot_rename_selected_protein_structure)
        self._protein_tree_context_menu.connect_show_sequence_action(self.__slot_show_protein_chain_sequence)
        self._protein_tree_context_menu.connect_help_action(self.__slot_open_proteins_tab_help)
        # </editor-fold>

        self._view.ui.btn_protein_sticks_show.clicked.connect(self.__slot_show_protein_regions_resi_sticks)
        self._view.ui.btn_protein_sticks_hide.clicked.connect(self.__slot_hide_protein_regions_resi_sticks)
        self._view.ui.btn_protein_disulfide_bonds_show.clicked.connect(
            self.__slot_show_protein_regions_disulfide_bonds)
        self._view.ui.btn_protein_disulfide_bonds_hide.clicked.connect(
            self.__slot_hide_protein_regions_disulfide_bonds)
        self._view.ui.btn_protein_position_zoom.clicked.connect(self.__slot_zoom_protein_regions_resi_position)

        # </editor-fold>

        # <editor-fold desc="Proteins Pair Tab">
        self._view.ui.protein_pairs_tree_view.customContextMenuRequested.connect(self.open_context_menu_for_protein_pairs)
        self._view.ui.protein_pairs_tree_view.clicked.connect(self.__slot_get_information_about_selected_object_in_protein_pair_branch)
        self._view.ui.btn_delete_protein_pair.clicked.connect(self.__slot_delete_protein_pair_from_project)
        self._view.ui.btn_open_protein_pair_session.clicked.connect(self.__slot_open_protein_pair_pymol_session)
        self._view.ui.btn_create_protein_pair_scene.clicked.connect(self.__slot_save_scene)
        self._view.ui.btn_update_protein_pair_scene.clicked.connect(self.__slot_update_protein_pair_scene)
        self._view.ui.btn_delete_protein_pair_scene.clicked.connect(self.__slot_delete_current_scene)
        self._view.ui.protein_pairs_tree_view.clicked.connect(self._check_for_results)
        self._view.ui.box_protein_pair_color.currentIndexChanged.connect(self.__slot_change_chain_color_protein_pairs)
        self._view.ui.btn_protein_pair_color_atoms.clicked.connect(self.__slot_change_chain_color_protein_pairs_atoms)
        self._view.ui.btn_protein_pair_reset_atoms.clicked.connect(self.__slot_change_chain_reset_protein_pairs_atoms)
        self._view.tg_protein_pair_white_bg.toggleChanged.connect(self.__slot_protein_pair_change_background_color)
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
        # toggle representation
        # self._view.ui.btn_protein_pair_show_hydrogens.clicked.connect(self.__slot_show_protein_pair_chain_with_hydrogens)
        # self._view.ui.btn_protein_pair_hide_hydrogens.clicked.connect(self.__slot_hide_protein_pair_chain_with_hydrogens)
        self._view.tg_protein_pair_hydrogen_atoms.toggleChanged.connect(self.__slot_protein_pair_chain_with_hydrogens)
        self._view.tg_protein_pair_cartoon.toggleChanged.connect(self.__slot_protein_pair_chain_as_cartoon)
        self._view.tg_protein_pair_sticks.toggleChanged.connect(self.__slot_protein_pair_chain_as_sticks)
        self._view.tg_protein_pair_ribbon.toggleChanged.connect(self.__slot_protein_pair_chain_as_ribbon)
        self._view.tg_protein_pair_lines.toggleChanged.connect(self.__slot_protein_pair_chain_as_lines)
        self._view.tg_protein_pair_spheres.toggleChanged.connect(self.__slot_protein_pair_chain_as_spheres)
        self._view.tg_protein_pair_dots.toggleChanged.connect(self.__slot_protein_pair_chain_as_dots)
        self._view.tg_protein_pair_mesh.toggleChanged.connect(self.__slot_protein_pair_chain_as_mesh)
        self._view.tg_protein_pair_surface.toggleChanged.connect(self.__slot_protein_pair_chain_as_surface)

        self._view.ui.btn_protein_pair_hide_all_representations.clicked.connect(self.__slot_hide_protein_pair_chain_all)
        self._view.ui.btn_help_3.clicked.connect(self._open_protein_pairs_tab_help)

        # <editor-fold desc="Color Grid">
        self._view.color_grid_protein_pairs.c_red.clicked.connect(self.set_color_name_in_label_red_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_tv_red.clicked.connect(self.set_color_name_in_label_tv_red_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_salomon.clicked.connect(self.set_color_name_in_label_salmon_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_raspberry.clicked.connect(self.set_color_name_in_label_raspberry_in_protein_pairs_tab)

        self._view.color_grid_protein_pairs.c_green.clicked.connect(self.set_color_name_in_label_green_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_tv_green.clicked.connect(self.set_color_name_in_label_tv_green_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_palegreen.clicked.connect(self.set_color_name_in_label_palegreen_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_forest.clicked.connect(self.set_color_name_in_label_forest_in_protein_pairs_tab)

        self._view.color_grid_protein_pairs.c_blue.clicked.connect(self.set_color_name_in_label_blue_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_tv_blue.clicked.connect(self.set_color_name_in_label_tv_blue_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_lightblue.clicked.connect(self.set_color_name_in_label_lightblue_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_skyblue.clicked.connect(self.set_color_name_in_label_skyblue_in_protein_pairs_tab)

        self._view.color_grid_protein_pairs.c_yellow.clicked.connect(self.set_color_name_in_label_yellow_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_tv_yellow.clicked.connect(self.set_color_name_in_label_tv_yellow_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_paleyellow.clicked.connect(self.set_color_name_in_label_paleyellow_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_sand.clicked.connect(self.set_color_name_in_label_sand_in_protein_pairs_tab)

        self._view.color_grid_protein_pairs.c_magenta.clicked.connect(self.set_color_name_in_label_magenta_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_purple.clicked.connect(self.set_color_name_in_label_purple_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_pink.clicked.connect(self.set_color_name_in_label_pink_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_hotpink.clicked.connect(self.set_color_name_in_label_hotpink_in_protein_pairs_tab)

        self._view.color_grid_protein_pairs.c_cyan.clicked.connect(self.set_color_name_in_label_cyan_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_aquamarine.clicked.connect(self.set_color_name_in_label_aquamarine_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_palecyan.clicked.connect(self.set_color_name_in_label_palecyan_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_teal.clicked.connect(self.set_color_name_in_label_teal_in_protein_pairs_tab)

        self._view.color_grid_protein_pairs.c_orange.clicked.connect(self.set_color_name_in_label_orange_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_tv_orange.clicked.connect(self.set_color_name_in_label_tv_orange_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_lightorange.clicked.connect(self.set_color_name_in_label_lightorange_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_olive.clicked.connect(self.set_color_name_in_label_olive_in_protein_pairs_tab)

        self._view.color_grid_protein_pairs.c_white.clicked.connect(self.set_color_name_in_label_white_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_grey_70.clicked.connect(self.set_color_name_in_label_grey_70_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_grey_30.clicked.connect(self.set_color_name_in_label_grey_30_in_protein_pairs_tab)
        self._view.color_grid_protein_pairs.c_black.clicked.connect(self.set_color_name_in_label_black_in_protein_pairs_tab)

        # </editor-fold>

        # <editor-fold desc="Context menu">
        self._protein_pair_tree_context_menu.connect_open_results_summary_action(self.__slot_results_summary)
        self._protein_pair_tree_context_menu.connect_color_based_on_rmsd_action(self.__slot_color_protein_pair_by_rmsd)
        self._protein_pair_tree_context_menu.connect_help_action(self._open_protein_pairs_tab_help)
        # </editor-fold>

        self._view.ui.btn_protein_pair_sticks_show.clicked.connect(self.__slot_show_protein_regions_resi_sticks)
        self._view.ui.btn_protein_pair_sticks_hide.clicked.connect(self.__slot_hide_protein_regions_resi_sticks)
        self._view.ui.btn_protein_pair_disulfide_bonds_show.clicked.connect(
            self.__slot_show_protein_regions_disulfide_bonds)
        self._view.ui.btn_protein_pair_disulfide_bonds_hide.clicked.connect(
            self.__slot_hide_protein_regions_disulfide_bonds)
        self._view.ui.btn_protein_pair_position_zoom.clicked.connect(self.__slot_zoom_protein_regions_resi_position)

        # </editor-fold>

    def _close_main_window(self, return_value):
        """Cleans after the main window closes."""
        _, tmp_event = return_value
        logger.info("Check if any jobs are running before closing PySSA.")
        if self._interface_manager.job_manager.there_are_jobs_running():
            logger.info("Running jobs found!")
            tmp_dialog = custom_message_box.CustomMessageBoxYesNo(
                "There are still jobs running!\n\nAre you sure you want to close PySSA?\nThis could lead to data loss and a damaged project!", "Close PySSA",
                custom_message_box.CustomMessageBoxIcons.WARNING.value
            )
            tmp_dialog.exec_()
            if not tmp_dialog.response:
                logger.info("PySSA should not be closed.")
                tmp_event.ignore()
                return
            logger.info("There are jobs running but PySSA will closed! (Requested by the user)")
        # Closes the documentation browser if it is still open
        if len(pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)) == 1:
            pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)[0].close()
        self._interface_manager.job_manager.stop_auxiliary_pymol()
        tmp_event.accept()

    def __slot_close_all(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Project/Exit Application' clicked.")
        tmp_dialog = custom_message_box.CustomMessageBoxYesNo(
            "Are you sure you want to close PySSA?", "Close PySSA",
            custom_message_box.CustomMessageBoxIcons.WARNING.value
        )
        tmp_dialog.exec_()
        if tmp_dialog.response:
            tmp_number_of_help_windows = len(pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER))
            tmp_number_of_pymol_windows = len(pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYMOL_PART))
            # PySSA should be closed
            if not self._view.ui.lbl_logo.isVisible():
                logger.info("A project is currently opend. It will now be saved and the application exists afterwards.")
                self.__slot_close_project()
            # Help windows
            if tmp_number_of_help_windows == 1:
                logger.info("The documentation window is open. It will be closed now.")
                pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)[0].close()
            elif tmp_number_of_help_windows > 1:
                for tmp_window_index in range(tmp_number_of_help_windows):
                    pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)[tmp_window_index].close()
            # PyMOL windows
            if tmp_number_of_pymol_windows == 1:
                logger.info("PyMOL will be closed now.")
                pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYMOL_PART)[0].close()
            elif tmp_number_of_pymol_windows > 1:
                tmp_dialog = custom_message_box.CustomMessageBoxYesNo(
                    "There are multiple windows open which contain PyMOL as window title.\nDo you want to close all?", "Close PySSA",
                    custom_message_box.CustomMessageBoxIcons.WARNING.value
                )
                tmp_dialog.exec_()
                if tmp_dialog.response:
                    for tmp_window_index in range(tmp_number_of_help_windows + 1):
                        pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYMOL_PART)[tmp_window_index].close()

    def _abort_task(self, return_value):
        if return_value[0] is True and return_value[1] == "ColabFold Prediction":
            self.__slot_abort_prediction()

    def _update_main_view_ui(self, refresh_after_job_finished_signal_values):
        tmp_job_is_for_current_project_flag, tmp_job_base_information, tmp_job_notification_widget = refresh_after_job_finished_signal_values
        self._interface_manager.remove_job_notification_widget(tmp_job_notification_widget)
        if tmp_job_base_information.job_progress == enums.JobProgress.FAILED:
            return
        if tmp_job_is_for_current_project_flag:
            if tmp_job_base_information.job_type == enums.JobType.PREDICTION:
                # refresh protein model
                self._active_task = tasks.Task(
                    target=util_async.add_proteins_to_project_and_model,
                    args=(self._interface_manager, tmp_job_base_information.protein_names),
                    post_func=self._post_update_project_and_model,
                )
                self._active_task.start()
            elif tmp_job_base_information.job_type == enums.JobType.DISTANCE_ANALYSIS:
                # refresh protein pair model
                self._active_task = tasks.Task(
                    target=util_async.add_protein_pairs_to_project_and_model,
                    args=(self._interface_manager, tmp_job_base_information.protein_pair_names),
                    post_func=self._post_update_project_and_model,
                )
                self._active_task.start()
            elif tmp_job_base_information.job_type == enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS:
                # refresh protein and protein pair model
                self._active_task = tasks.Task(
                    target=util_async.add_proteins_and_protein_pairs_to_project_and_model,
                    args=(
                        self._interface_manager,
                        tmp_job_base_information.protein_names,
                        tmp_job_base_information.protein_pair_names
                    ),
                    post_func=self._post_update_project_and_model,
                )
                self._active_task.start()
            elif tmp_job_base_information.job_type == enums.JobType.RAY_TRACING:
                if not os.path.exists(tmp_job_base_information.get_image_filepath()):
                    self._interface_manager.status_bar_manager.show_error_message("Could not find image file!")
                os.startfile(tmp_job_base_information.get_image_filepath())
                return
        else:
            # open other project
            tmp_a_project_is_open = self._view.ui.lbl_project_name.isVisible()
            if tmp_a_project_is_open:
                self._disconnect_sequence_selection_model()
            self._active_task = tasks.Task(
                target=util_async.close_project_automatically,
                args=(
                    tmp_a_project_is_open,
                    self._database_thread,
                    self._interface_manager.pymol_session_manager,
                    tmp_job_base_information.project_name
                ),
                post_func=self._post_close_project_automatically,
            )
            self._active_task.start()
            self.update_status("Saving current project ...")
            self._interface_manager.restore_default_main_view()
            self._interface_manager.close_job_notification_panel()
        self._interface_manager.start_wait_cursor()

    def _post_update_project_and_model(self, return_value):
        self._interface_manager = return_value[0]
        self._interface_manager.refresh_main_view()
        self._interface_manager.stop_wait_cursor()

    def _post_close_project_automatically(self, return_value):
        tmp_project_name = return_value[0]
        self._interface_manager.status_bar_manager.show_temporary_message(
            enums.StatusMessages.OPENING_PROJECT.value, False
        )
        tmp_project_database_filepath = str(
            pathlib.Path(
                f"{self._interface_manager.get_application_settings().workspace_path}/{tmp_project_name}.db"
            )
        )
        self._database_thread = database_thread.DatabaseThread(tmp_project_database_filepath)
        #self._database_thread.start()
        self._database_manager.set_database_filepath(tmp_project_database_filepath)
        self._active_task = tasks.Task(
            target=project_async.open_project,
            args=(
                tmp_project_name,
                tmp_project_database_filepath,
                self._interface_manager,
                self._interface_manager.pymol_session_manager,
                self.custom_progress_signal,
                self._interface_manager.watcher
            ),
            post_func=self.__await_open_project,
        )
        self._active_task.start()

    def _update_progress_bar(self, return_value):
        self._interface_manager.status_bar_manager.update_progress_bar(return_value)

    # <editor-fold desc="Util methods">
    def update_status(self, message: str) -> None:
        """Updates the status bar of the main view with a custom message."""
        self._interface_manager.status_bar_manager.show_temporary_message(message)

    def _update_tab(self):
        self._interface_manager.current_tab_index = self._view.ui.project_tab_widget.currentIndex()
        if self._interface_manager.pymol_session_manager.session_object_type == "protein" and self._interface_manager.current_tab_index == 2:
            self._interface_manager.hide_protein_pair_pymol_scene_configuration()
            self._view.ui.lbl_info_3.setText(
                "Please load the PyMOL session of the \nselected protein pair.")
        elif self._interface_manager.pymol_session_manager.session_object_type == "protein_pair" and self._interface_manager.current_tab_index == 1:
            self._interface_manager.hide_protein_pymol_scene_configuration()
            self._view.ui.lbl_info.setText("Please load the PyMOL session of the selected protein.")
        elif self._interface_manager.pymol_session_manager.is_the_current_session_empty():
            self._interface_manager.hide_protein_pymol_scene_configuration()
            self._interface_manager.hide_protein_pair_pymol_scene_configuration()
            self._view.ui.lbl_info.setText("Please load the PyMOL session of the selected protein.")
            self._view.ui.lbl_info_3.setText(
                "Please load the PyMOL session of the \nselected protein pair.")

    def _setup_statusbar(self) -> None:
        """Sets up the status bar and fills it with the current workspace."""
        self._interface_manager.get_main_view().setStatusBar(self._interface_manager.get_main_view().status_bar)

    # <editor-fold desc="Help related methods">
    def _start_documentation_server(self):
        self._help_task = tasks.Task(
            target=util_async.start_documentation_server,
            args=(0, 0),
            post_func=self.__await_start_documentation_server,
        )
        self._help_task.start()

    def __await_start_documentation_server(self, return_value: tuple):
        self._interface_manager.documentation_window = return_value[1]
        self._interface_manager.status_bar_manager.show_temporary_message("Opening help center finished.")

    def open_help(self, a_page_name: str):
        """Opens the pyssa documentation window if it's not already open.

        Args:
            a_page_name (str): a name of a documentation page to display
        """
        self._interface_manager.start_wait_cursor()
        self._interface_manager.status_bar_manager.show_temporary_message(
            "Opening help center ...", False)

        if len(pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)) != 1:
            self._interface_manager.documentation_window = None
        self._help_task = tasks.Task(
            target=util_async.open_documentation_on_certain_page,
            args=(a_page_name, self._interface_manager.documentation_window),
            post_func=self.__await_open_help,
        )
        self._help_task.start()

    def __await_open_help(self, return_value):
        self._interface_manager.documentation_window = return_value[2]
        if not os.path.exists(constants.HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH):
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "The script for bringing the documentation window in front could not be found!", "Documentation",
                custom_message_box.CustomMessageBoxIcons.ERROR.value
            )
            tmp_dialog.exec_()
        else:
            self._interface_manager.documentation_window.restore()
            subprocess.run([constants.HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH])
            self._interface_manager.status_bar_manager.show_temporary_message("Opening help center finished.")
        self._interface_manager.stop_wait_cursor()
        self._interface_manager.refresh_main_view()

    def _init_generic_help_context_menus(self):
        # <editor-fold desc="General context menu setup">
        context_menu = QtWidgets.QMenu()
        self.help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        self.help_context_action.triggered.connect(self.__slot_open_sequences_tab_help)

        # </editor-fold>

        # Set the context menu for the buttons
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
    def __slot_open_help_center(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Help/Documentation' clicked.")
        self.open_help("help/")

    def __slot_open_sequences_tab_help(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Help' button on the 'Sequence Tab' was clicked.")
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

    def __slot_open_proteins_tab_help(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,"'Help' button on the 'Proteins Tab' was clicked.")
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
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Help' button on the 'Protein Pairs Tab' was clicked.")
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
    def _show_context_menu_for_seq_table(self, a_point):
        context_menu = QtWidgets.QMenu()
        help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        help_context_action.triggered.connect(self._open_additional_information_table_help)
        context_menu.exec_(self._view.ui.seqs_table_widget.mapToGlobal(a_point))

    def _show_context_menu_for_seq_import(self, a_point):
        context_menu = QtWidgets.QMenu()
        help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        help_context_action.triggered.connect(self._open_sequence_import_help)
        context_menu.exec_(self._view.ui.btn_import_seq.mapToGlobal(a_point))

    def _show_context_menu_for_seq_add(self, a_point):
        context_menu = QtWidgets.QMenu()
        help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        help_context_action.triggered.connect(self._open_sequence_add_help)
        context_menu.exec_(self._view.ui.btn_add_sequence.mapToGlobal(a_point))

    def _show_context_menu_for_seq_save(self, a_point):
        context_menu = QtWidgets.QMenu()
        help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        help_context_action.triggered.connect(self._open_sequence_save_help)
        context_menu.exec_(self._view.ui.btn_save_sequence.mapToGlobal(a_point))

    def _show_context_menu_for_seq_delete(self, a_point):
        context_menu = QtWidgets.QMenu()
        help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        help_context_action.triggered.connect(self._open_sequence_delete_help)
        context_menu.exec_(self._view.ui.btn_delete_sequence.mapToGlobal(a_point))

    def _show_context_menu_for_protein_import(self, a_point):
        context_menu = QtWidgets.QMenu()
        help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        help_context_action.triggered.connect(self._open_protein_import_help)
        context_menu.exec_(self._view.ui.btn_import_protein.mapToGlobal(a_point))

    def _show_context_menu_for_protein_save(self, a_point):
        context_menu = QtWidgets.QMenu()
        help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        help_context_action.triggered.connect(self._open_protein_save_help)
        context_menu.exec_(self._view.ui.btn_save_protein.mapToGlobal(a_point))

    def _show_context_menu_for_protein_delete(self, a_point):
        context_menu = QtWidgets.QMenu()
        help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        help_context_action.triggered.connect(self._open_protein_delete_help)
        context_menu.exec_(self._view.ui.btn_delete_protein.mapToGlobal(a_point))

    def _show_context_menu_for_protein_load_session(self, a_point):
        context_menu = QtWidgets.QMenu()
        help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        help_context_action.triggered.connect(self._open_protein_load_session_help)
        context_menu.exec_(self._view.ui.btn_open_protein_session.mapToGlobal(a_point))

    def _show_context_menu_for_protein_add_scene(self, a_point):
        context_menu = QtWidgets.QMenu()
        help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        help_context_action.triggered.connect(self._open_protein_add_scene_help)
        context_menu.exec_(self._view.ui.btn_create_protein_scene.mapToGlobal(a_point))

    def _show_context_menu_for_protein_update_scene(self, a_point):
        context_menu = QtWidgets.QMenu()
        help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        help_context_action.triggered.connect(self._open_protein_update_scene_help)
        context_menu.exec_(self._view.ui.btn_update_protein_scene.mapToGlobal(a_point))

    def _show_context_menu_for_protein_delete_scene(self, a_point):
        context_menu = QtWidgets.QMenu()
        help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        help_context_action.triggered.connect(self._open_protein_delete_scene_help)
        context_menu.exec_(self._view.ui.btn_delete_protein_scene.mapToGlobal(a_point))

    def _show_context_menu_for_protein_pymol_scene_config(self, a_point):
        context_menu = QtWidgets.QMenu()
        help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        help_context_action.triggered.connect(self._open_protein_pymol_scene_config_help)
        context_menu.exec_(self._view.ui.frame_protein_pymol_scene.mapToGlobal(a_point))

    def _show_context_menu_for_protein_pair_delete(self, a_point):
        context_menu = QtWidgets.QMenu()
        help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        help_context_action.triggered.connect(self._open_protein_pair_delete_help)
        context_menu.exec_(self._view.ui.btn_delete_protein_pair.mapToGlobal(a_point))

    def _show_context_menu_for_protein_pair_load_session(self, a_point):
        context_menu = QtWidgets.QMenu()
        help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        help_context_action.triggered.connect(self._open_protein_pair_load_session_help)
        context_menu.exec_(self._view.ui.btn_open_protein_pair_session.mapToGlobal(a_point))

    def _show_context_menu_for_protein_pair_add_scene(self, a_point):
        context_menu = QtWidgets.QMenu()
        help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        help_context_action.triggered.connect(self._open_protein_pair_add_scene_help)
        context_menu.exec_(self._view.ui.btn_create_protein_pair_scene.mapToGlobal(a_point))

    def _show_context_menu_for_protein_pair_update_scene(self, a_point):
        context_menu = QtWidgets.QMenu()
        help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        help_context_action.triggered.connect(self._open_protein_pair_update_scene_help)
        context_menu.exec_(self._view.ui.btn_update_protein_pair_scene.mapToGlobal(a_point))

    def _show_context_menu_for_protein_pair_delete_scene(self, a_point):
        context_menu = QtWidgets.QMenu()
        help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        help_context_action.triggered.connect(self._open_protein_pair_delete_scene_help)
        context_menu.exec_(self._view.ui.btn_delete_protein_pair_scene.mapToGlobal(a_point))

    def _show_context_menu_for_protein_pair_pymol_scene_config(self, a_point):
        context_menu = QtWidgets.QMenu()
        help_context_action = context_menu.addAction(self._view.tr("Get Help"))
        help_context_action.triggered.connect(self._open_protein_pair_pymol_scene_config_help)
        context_menu.exec_(self._view.ui.frame_protein_pair_pymol_scene.mapToGlobal(a_point))

    # </editor-fold>

    # </editor-fold>

    # </editor-fold>

    # <editor-fold desc="Job panels">
    def __slot_open_job_overview_panel(self):
        if self._view.ui.frame_job_overview.isVisible():
            self._view.ui.frame_job_overview.hide()
        elif not self._view.ui.frame_job_overview.isVisible():
            self._view.ui.frame_job_notification.hide()
            self._view.ui.frame_job_overview.show()

    def __slot_open_notification_panel(self):
        if self._view.ui.frame_job_notification.isVisible():
            self._view.ui.frame_job_notification.hide()
            self._view.btn_open_job_notification.setIcon(self._view.icon_notify)
        elif not self._view.ui.frame_job_notification.isVisible():
            self._view.ui.frame_job_overview.hide()
            self._view.ui.frame_job_notification.show()
            self._view.btn_open_job_notification.setIcon(self._view.icon_notify)
    # </editor-fold>

    # <editor-fold desc="Project menu">
    def __slot_close_project(self):
        """Closes the current project"""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Project/Close' clicked.")
        self._active_task = tasks.Task(
            target=project_async.close_project,
            args=(self._database_thread, self._interface_manager.pymol_session_manager),
            post_func=self.__await_close_project,
        )
        self._active_task.start()
        # self.msg_box = basic_boxes.no_buttons("Saving Project",
        #                                       "Please wait the program is saving your project.",
        #                                       QtWidgets.QMessageBox.Information)
        # # self.msg_box.show()
        self._interface_manager.restore_default_main_view()
        self._interface_manager.close_job_notification_panel()
        self._interface_manager.close_job_overview_panel()
        self._disconnect_sequence_selection_model()
        self.update_status("Saving current project ...")
        self._interface_manager.start_wait_cursor()

    def __await_close_project(self):
        """Await the async closing process."""
        self._interface_manager.set_new_project(project.Project())
        self._view.ui.project_tab_widget.setCurrentIndex(0)
        self._interface_manager.refresh_main_view()
        self.update_status("Closing project finished.")
        self._interface_manager.stop_wait_cursor()

    def __slot_create_project(self) -> None:
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Project/Create' clicked.")
        self._external_controller = create_project_view_controller.CreateProjectViewController(self._interface_manager)
        self._external_controller.user_input.connect(self._post_create_project)
        self._interface_manager.get_create_view().show()

    def _post_create_project(self, user_input: tuple) -> None:
        self._interface_manager.start_wait_cursor()
        tmp_project_name, tmp_protein_name = user_input
        tmp_project_database_filepath = str(
            pathlib.Path(f"{self._interface_manager.get_application_settings().workspace_path}/{tmp_project_name}.db"))
        with database_manager.DatabaseManager(tmp_project_database_filepath) as db_manager:
            db_manager.build_new_database()
        self._database_thread = database_thread.DatabaseThread(tmp_project_database_filepath)
        #self._database_thread.start()
        self._active_task = tasks.Task(
            target=project_async.create_new_project,
            args=(
                tmp_project_name,
                self._interface_manager.get_application_settings().get_workspace_path(),
                self._interface_manager.watcher,
                self._interface_manager
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
        # #self._database_thread.start()
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
            self._interface_manager.refresh_main_view()
            return

        _, tmp_project, self._interface_manager.watcher, self._interface_manager = return_value
        self._interface_manager.set_new_project(tmp_project)
        self._interface_manager.refresh_workspace_model()
        self._interface_manager.refresh_main_view()
        self._interface_manager.pymol_session_manager.reinitialize_session()
        self._connect_sequence_selection_model()
        self._interface_manager.stop_wait_cursor()

    def __slot_open_project(self) -> None:
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Project/Open' clicked.")
        self._external_controller = open_project_view_controller.OpenProjectViewController(self._interface_manager)
        self._external_controller.return_value.connect(self._post_open_project)
        self._interface_manager.get_open_view().show()

    def _post_open_project(self, return_value: str):
        if return_value[1] is False:
            self._interface_manager.refresh_main_view()
            return

        self._interface_manager.status_bar_manager.show_temporary_message(
            enums.StatusMessages.OPENING_PROJECT.value, False
        )
        self._interface_manager.start_wait_cursor()
        tmp_project_name = return_value
        tmp_project_database_filepath = str(
            pathlib.Path(
                f"{self._interface_manager.get_application_settings().workspace_path}/{tmp_project_name}.db"
            )
        )
        self._database_thread = database_thread.DatabaseThread(tmp_project_database_filepath)
        #self._database_thread.start()
        self._database_manager.set_database_filepath(tmp_project_database_filepath)
        self._active_task = tasks.Task(
            target=project_async.open_project,
            args=(
                tmp_project_name,
                tmp_project_database_filepath,
                self._interface_manager,
                self._interface_manager.pymol_session_manager,
                self.custom_progress_signal,
                self._interface_manager.watcher
            ),
            post_func=self.__await_open_project,
        )
        self._active_task.start()

    def __await_open_project(self, return_value: tuple):
        self._interface_manager.status_bar_manager.hide_progress_bar()
        exit_code, tmp_project, tmp_interface_manager, tmp_watcher = return_value
        if exit_code == 0:
            self._interface_manager = tmp_interface_manager
            self._interface_manager.watcher = tmp_watcher
            self._interface_manager.refresh_main_view()
            self._interface_manager.hide_progress_bar()
            self._interface_manager.status_bar_manager.show_temporary_message(
                enums.StatusMessages.OPENING_PROJECT_FINISHED.value
            )
            self._connect_sequence_selection_model()
        else:
            self._interface_manager.status_bar_manager.show_error_message(
                enums.StatusMessages.OPENING_PROJECT_FAILED.value,
            )
            self._interface_manager.refresh_main_view()
        self._interface_manager.stop_wait_cursor()

    def __slot_use_project(self) -> None:
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Project/Use' clicked.")
        self._external_controller = use_project_view_controller.UseProjectViewController(self._interface_manager)
        self._external_controller.user_input.connect(self._post_use_project)
        self._interface_manager.get_use_project_view().show()

    def _post_use_project(self, user_input: tuple) -> None:
        self._interface_manager.start_wait_cursor()
        tmp_project_database_filepath = str(pathlib.Path(f"{self._interface_manager.get_application_settings().get_workspace_path()}/{user_input[0]}.db"))
        with database_manager.DatabaseManager(tmp_project_database_filepath) as db_manager:
            db_manager.build_new_database()
            db_manager.close_project_database()

        self._active_task = tasks.Task(
            target=project_async.create_use_project,
            args=(
                user_input[0],
                self._interface_manager.get_application_settings().get_workspace_path(),
                user_input[1]
            ),
            post_func=self.__await_use_project,
        )
        self._active_task.start()

    def __await_use_project(self, return_value: tuple):
        _, tmp_project, self._interface_manager.watcher, self._interface_manager = return_value
        self._interface_manager.set_new_project(tmp_project)
        self._interface_manager.add_project_to_workspace_model(tmp_project.get_project_name())
        self._interface_manager.refresh_main_view()
        self._interface_manager.pymol_session_manager.reinitialize_session()
        self._connect_sequence_selection_model()
        self._interface_manager.stop_wait_cursor()
        self._interface_manager.status_bar_manager.show_temporary_message("Use process finished.")

    def __slot_delete_project(self) -> None:
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Project/Delete' clicked.")
        self._external_controller = delete_project_view_controller.DeleteProjectViewController(self._interface_manager)
        self._interface_manager.get_delete_view().show()

    def _post_delete_project(self) -> None:
        self._interface_manager.refresh_main_view()

    def __slot_import_project(self) -> None:
        """Imports a project.xml into the current workspace."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Project/Import' clicked.")
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
            #self._database_thread.start()
            self._database_manager.set_database_filepath(tmp_project_database_filepath)
            #self._database_manager.open_project_database()
            self._database_manager.update_project_name(tmp_new_project_name)
            tmp_project = self._database_manager.get_project_as_object(
                tmp_new_project_name,
                self._interface_manager.get_application_settings().workspace_path,
                self._interface_manager.get_application_settings()
            )
            self._interface_manager.set_new_project(tmp_project)

            self._interface_manager.refresh_main_view()
            self._interface_manager.pymol_session_manager.reinitialize_session()
            self._interface_manager.stop_wait_cursor()
            self._interface_manager.status_bar_manager.show_temporary_message("Importing project finished.")

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

    def __slot_export_current_project(self) -> None:
        """Exports the current project to an importable format."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Project/Export' clicked.")
        file_dialog = QtWidgets.QFileDialog()
        desktop_path = QtCore.QStandardPaths.standardLocations(QtCore.QStandardPaths.DesktopLocation)[0]
        file_dialog.setDirectory(desktop_path)
        file_path, _ = file_dialog.getSaveFileName(self._view, "Export current project", "", "Project Database File (*.db)")
        if file_path:
            shutil.copyfile(self._interface_manager.get_current_project().get_database_filepath(),
                            file_path)
            # tmp_dialog = custom_message_box.CustomMessageBoxOk(
            #     "The project was successfully exported.", "Export Project",
            #     custom_message_box.CustomMessageBoxIcons.INFORMATION.value
            # )
            # tmp_dialog.exec_()
            self._interface_manager.status_bar_manager.show_temporary_message("The project was successfully exported.")

    # </editor-fold>

    # <editor-fold desc="Analysis menu">
    def __slot_distance_analysis(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Analysis/Distance' clicked.")
        self._external_controller = distance_analysis_view_controller.DistanceAnalysisViewController(
            self._interface_manager, self._interface_manager.watcher
        )
        self._external_controller.job_input.connect(self._post_distance_analysis)
        self._interface_manager.get_distance_analysis_view().show()

    def _post_distance_analysis(self, job_input: tuple) -> None:
        """Sets up the worker for the analysis task."""
        constants.PYSSA_LOGGER.info("Begin analysis process.")

        _, tmp_raw_analysis_run_names, tmp_checkbox_state = job_input

        # --- New job approach
        #self._interface_manager.get_current_project()
        self._interface_manager.watcher.add_protein_pairs_from_new_job(tmp_raw_analysis_run_names)
        tmp_distance_analysis_job, tmp_distance_analysis_entry_widget = self._interface_manager.job_manager.create_distance_analysis_job(
            self._interface_manager.get_current_project(),
            self._interface_manager.project_lock,
            self._interface_manager,
            tmp_raw_analysis_run_names,
            self._interface_manager.get_settings_manager().settings.cutoff,
            self._interface_manager.get_settings_manager().settings.cycles
        )
        self._interface_manager.job_manager.put_job_into_queue(tmp_distance_analysis_job)
        self._interface_manager.add_job_entry_to_job_overview_layout(tmp_distance_analysis_entry_widget)

        # tmp_distance_analysis_task = tasks.Task(
        #     target=main_tasks_async.run_distance_analysis,
        #     args=(
        #         tmp_raw_analysis_run_names,
        #         self._interface_manager.get_current_project(),
        #         self._interface_manager.get_application_settings(),
        #         tmp_checkbox_state,
        #         self.custom_progress_signal,
        #         self._interface_manager.pymol_lock,
        #         self.disable_pymol_signal
        #     ),
        #     post_func=self.__await_run_distance_analysis,
        # )
        # self._interface_manager.main_tasks_manager.start_distance_analysis_task(tmp_distance_analysis_task)
        #
        # if not os.path.exists(constants.SCRATCH_DIR_ANALYSIS):
        #     os.mkdir(constants.SCRATCH_DIR_ANALYSIS)
        # self._interface_manager.status_bar_manager.update_progress_bar(
        #     (enums.StatusMessages.DISTANCE_ANALYSIS_IS_RUNNING.value, 10)
        # )
        # self._main_view_state.set_protein_pairs_list(self._interface_manager.get_current_project().protein_pairs)
        # self._interface_manager.refresh_main_view()

    def _add_new_protein_pairs_to_protein_pair_model(self):
        """Adds the new protein pairs to the interface manager's protein pair model."""
        tmp_protein_pairs_to_add = self._main_view_state.get_not_matching_protein_pairs(
            self._interface_manager.get_current_project().protein_pairs
        )
        for tmp_protein_pair in tmp_protein_pairs_to_add:
            self._interface_manager.add_protein_pair_to_protein_pairs_model(tmp_protein_pair)

    def __await_unfreeze_pymol_session_after_analysis(self):
        self._interface_manager.main_tasks_manager.distance_analysis_task = None
        self._interface_manager.refresh_main_view()
        # self.active_custom_message_box = custom_message_box.CustomMessageBoxOk(
        #     "All structure analysis' are done. \nGo to the Protein Pairs tab to view the new results.",
        #     "Distance Analysis",
        #     custom_message_box.CustomMessageBoxIcons.INFORMATION.value
        # )
        # self.active_custom_message_box.exec_()
        constants.PYSSA_LOGGER.info("All structure analysis' are done.")
        self._interface_manager.status_bar_manager.show_temporary_message("All structure analysis' are done.")

    # </editor-fold>

    # <editor-fold desc="Prediction menu">
    def _connect_sequence_selection_model(self):
        self._view.ui.seqs_list_view.selectionModel().selectionChanged.connect(
            self._check_options_for_sequence_selection)

    def _disconnect_sequence_selection_model(self):
        try:
            self._view.ui.seqs_list_view.selectionModel().selectionChanged.disconnect(
                self._check_options_for_sequence_selection)
        except TypeError:
            logger.warning("Catching error because no sequences exist in the project. Therefore the selectionChanged signal does not need to be disconnected.")

    def _check_options_for_sequence_selection(self):
        if len(self._view.ui.seqs_list_view.selectedIndexes()) > 0:
            tmp_enable_monomer_flag = False
            tmp_enable_multimer_flag = False
            for tmp_model_index in self._view.ui.seqs_list_view.selectedIndexes():
                tmp_type = tmp_model_index.data(enums.ModelEnum.TYPE_ROLE)
                tmp_sequence_name: SeqRecord.SeqRecord = tmp_model_index.data(enums.ModelEnum.OBJECT_ROLE).name
                if tmp_type == enums.ModelTypeEnum.MONOMER_SEQ and tmp_enable_monomer_flag is False and not self._interface_manager.get_current_project().is_sequence_as_protein_in_project(tmp_sequence_name):
                    tmp_enable_monomer_flag = True
                elif tmp_model_index.data(enums.ModelEnum.TYPE_ROLE) == enums.ModelTypeEnum.MULTIMER_SEQ and tmp_enable_multimer_flag is False and not self._interface_manager.get_current_project().is_sequence_as_protein_in_project(tmp_sequence_name):
                    tmp_enable_multimer_flag = True
            self._view.ui.action_predict_monomer.setEnabled(tmp_enable_monomer_flag)
            self._view.ui.action_predict_multimer.setEnabled(tmp_enable_multimer_flag)
            if tmp_enable_monomer_flag or tmp_enable_multimer_flag:
                self._view.ui.menuPrediction.setEnabled(True)
            elif not tmp_enable_monomer_flag and not tmp_enable_multimer_flag:
                self._view.ui.menuPrediction.setEnabled(False)
        else:
            self._view.ui.btn_save_sequence.setEnabled(False)
            self._view.ui.btn_delete_sequence.setEnabled(False)
            self._interface_manager.refresh_main_view()

    # <editor-fold desc="Protein structure prediction">
    def __slot_predict_monomer(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Prediction/Monomer' clicked.")
        tmp_indexes = []
        if len(self._view.ui.seqs_list_view.selectedIndexes()) == 0:
            tmp_model = self._interface_manager.get_main_view().ui.seqs_list_view.model()
            for tmp_row_no in range(tmp_model.rowCount()):
                tmp_index = tmp_model.index(tmp_row_no, 0)
                if tmp_index.data(enums.ModelEnum.TYPE_ROLE) == enums.ModelTypeEnum.MONOMER_SEQ:
                    tmp_indexes.append(tmp_index)
        else:
            tmp_indexes = self._view.ui.seqs_list_view.selectedIndexes()
        self._external_controller = predict_protein_view_controller.PredictProteinViewController(
            self._interface_manager, self._interface_manager.watcher, tmp_indexes, "monomer"
        )
        self._external_controller.job_input.connect(self._post_predict_protein)
        self._interface_manager.get_predict_protein_view().show()

    def __slot_predict_multimer(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Prediction/Multimer' clicked.")
        tmp_indexes = []
        if len(self._view.ui.seqs_list_view.selectedIndexes()) == 0:
            tmp_model = self._interface_manager.get_main_view().ui.seqs_list_view.model()
            for tmp_row_no in range(tmp_model.rowCount()):
                tmp_index = tmp_model.index(tmp_row_no, 0)
                if tmp_index.data(enums.ModelEnum.TYPE_ROLE) == enums.ModelTypeEnum.MULTIMER_SEQ:
                    tmp_indexes.append(tmp_index)
        else:
            tmp_indexes = self._view.ui.seqs_list_view.selectedIndexes()
        self._external_controller = predict_protein_view_controller.PredictProteinViewController(
            self._interface_manager, self._interface_manager.watcher, tmp_indexes, "multimer"
        )
        self._external_controller.job_input.connect(self._post_predict_protein)
        self._interface_manager.get_predict_protein_view().show()

    def _post_predict_protein(self, result: tuple):

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
        if result[3] is True:
            constants.PYSSA_LOGGER.info("Running prediction with subsequent analysis.")
            # --- New job approach
            _, tmp_prediction_protein_infos, tmp_prediction_configuration, _ = result
            self._interface_manager.watcher.add_proteins_from_new_job(tmp_prediction_protein_infos)
            tmp_prediction_job, _ = self._interface_manager.job_manager.create_prediction_job(
                self._interface_manager.get_current_project(),
                tmp_prediction_protein_infos,
                tmp_prediction_configuration,
                self._interface_manager.project_lock,
                self._interface_manager
            )
            tmp_raw_analysis_run_names: list = []
            for row_no in range(self._interface_manager.get_predict_protein_view().ui.list_analysis_overview.count()):
                tmp_raw_analysis_run_names.append(
                    self._interface_manager.get_predict_protein_view().ui.list_analysis_overview.item(row_no).text())

            self._interface_manager.watcher.add_protein_pairs_from_new_job(tmp_raw_analysis_run_names)
            tmp_distance_analysis_job, _ = self._interface_manager.job_manager.create_distance_analysis_job(
                self._interface_manager.get_current_project(),
                self._interface_manager.project_lock,
                self._interface_manager,
                tmp_raw_analysis_run_names,
                self._interface_manager.get_settings_manager().settings.cutoff,
                self._interface_manager.get_settings_manager().settings.cycles
            )
            tmp_job = self._interface_manager.job_manager.create_prediction_and_distance_analysis_job(
                tmp_prediction_job,
                tmp_distance_analysis_job,
                self._interface_manager
            )
            tmp_prediction_and_distance_analysis_job = tmp_job[0]
            tmp_prediction_and_distance_analysis_job_entry_widget = tmp_job[1]
            self._interface_manager.job_manager.put_job_into_queue(tmp_prediction_and_distance_analysis_job)
            self._interface_manager.add_job_entry_to_job_overview_layout(tmp_prediction_and_distance_analysis_job_entry_widget)
        else:
            constants.PYSSA_LOGGER.info("Running prediction without subsequent analysis.")
            # --- New job approach
            _, tmp_prediction_protein_infos, tmp_prediction_configuration, _ = result
            self._interface_manager.watcher.add_proteins_from_new_job(tmp_prediction_protein_infos)
            tmp_prediction_job, tmp_prediction_entry_widget = self._interface_manager.job_manager.create_prediction_job(
                self._interface_manager.get_current_project(),
                tmp_prediction_protein_infos,
                tmp_prediction_configuration,
                self._interface_manager.project_lock,
                self._interface_manager
            )
            self._interface_manager.job_manager.put_job_into_queue(tmp_prediction_job)
            self._interface_manager.add_job_entry_to_job_overview_layout(tmp_prediction_entry_widget)
        #     tmp_prediction_task = tasks.Task(
        #         target=main_tasks_async.predict_protein_with_colabfold,
        #         args=(
        #             result[1],
        #             result[2],
        #             self._interface_manager.get_current_project(),
        #             self.custom_progress_signal,
        #             self._interface_manager.pymol_lock,
        #             self.disable_pymol_signal
        #         ),
        #         post_func=self.__await_predict_protein_with_colabfold,
        #     )
        # self._interface_manager.main_tasks_manager.start_prediction_task(tmp_prediction_task)
        self._main_view_state.set_proteins_list(self._interface_manager.get_current_project().proteins)
        self._main_view_state.set_protein_pairs_list(self._interface_manager.get_current_project().protein_pairs)
        self._interface_manager.refresh_main_view()

    def __slot_abort_prediction(self) -> None:
        """Aborts the running prediction."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Prediction/Abort' clicked.")
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
        self._interface_manager.status_bar_manager.hide_progress_bar()
        self._interface_manager.refresh_main_view()

    # </editor-fold>

    # <editor-fold desc="Old code">
    # def _post_predict_monomer(self, result: tuple) -> None:
    #     """Sets up the worker for the prediction of the proteins."""
    #
    #     # <editor-fold desc="Check if WSL2 and ColabFold are installed">
    #     if globals.g_os == "win32":
    #         constants.PYSSA_LOGGER.info("Checking if WSL2 is installed ...")
    #         if not dialog_settings_global.is_wsl2_installed():
    #             constants.PYSSA_LOGGER.warning("WSL2 is NOT installed.")
    #             self._interface_manager.get_application_settings().wsl_install = 0
    #             tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #                 "Prediction failed because the WSL2 environment is not installed!",
    #                 "Structure Prediction",
    #                 custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #             )
    #             tmp_dialog.exec_()
    #             return
    #         constants.PYSSA_LOGGER.info("Checking if Local Colabfold is installed ...")
    #         if not dialog_settings_global.is_local_colabfold_installed():
    #             constants.PYSSA_LOGGER.warning("Local Colabfold is NOT installed.")
    #             self._interface_manager.get_application_settings().local_colabfold = 0
    #             tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #                 "Prediction failed because the ColabFold is not installed!",
    #                 "Structure Prediction",
    #                 custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #             )
    #             tmp_dialog.exec_()
    #             return
    #
    #     # </editor-fold>
    #
    #     self.prediction_type = constants.PREDICTION_TYPE_PRED_MONO_ANALYSIS
    #     constants.PYSSA_LOGGER.info("Begin prediction process.")
    #     if result[3] is True:
    #         constants.PYSSA_LOGGER.info("Running prediction with subsequent analysis.")
    #         # --- New job approach
    #         _, tmp_prediction_protein_infos, tmp_prediction_configuration, _ = result
    #         self._interface_manager.watcher.add_proteins_from_new_job(tmp_prediction_protein_infos)
    #         tmp_prediction_job, _ = self._interface_manager.job_manager.create_prediction_job(
    #             self._interface_manager.get_current_project(),
    #             tmp_prediction_protein_infos,
    #             tmp_prediction_configuration,
    #             self._interface_manager.project_lock,
    #             self._interface_manager
    #         )
    #         tmp_raw_analysis_run_names: list = []
    #         for row_no in range(self._interface_manager.get_predict_protein_view().ui.list_analysis_overview.count()):
    #             tmp_raw_analysis_run_names.append(
    #                 self._interface_manager.get_predict_protein_view().ui.list_analysis_overview.item(row_no).text())
    #
    #         self._interface_manager.watcher.add_protein_pairs_from_new_job(tmp_raw_analysis_run_names)
    #         tmp_distance_analysis_job, _ = self._interface_manager.job_manager.create_distance_analysis_job(
    #             self._interface_manager.get_current_project(),
    #             self._interface_manager.project_lock,
    #             self._interface_manager,
    #             tmp_raw_analysis_run_names,
    #             self._interface_manager.get_settings_manager().settings.cutoff,
    #             self._interface_manager.get_settings_manager().settings.cycles
    #         )
    #         tmp_job = self._interface_manager.job_manager.create_prediction_and_distance_analysis_job(
    #             tmp_prediction_job,
    #             tmp_distance_analysis_job,
    #             self._interface_manager
    #         )
    #         tmp_prediction_and_distance_analysis_job = tmp_job[0]
    #         tmp_prediction_and_distance_analysis_job_entry_widget = tmp_job[1]
    #         self._interface_manager.job_manager.put_job_into_queue(tmp_prediction_and_distance_analysis_job)
    #         self._interface_manager.add_job_entry_to_job_overview_layout(tmp_prediction_and_distance_analysis_job_entry_widget)
    #     else:
    #         constants.PYSSA_LOGGER.info("Running prediction without subsequent analysis.")
    #         # --- New job approach
    #         _, tmp_prediction_protein_infos, tmp_prediction_configuration, _ = result
    #         self._interface_manager.watcher.add_proteins_from_new_job(tmp_prediction_protein_infos)
    #         tmp_prediction_job, tmp_prediction_entry_widget = self._interface_manager.job_manager.create_prediction_job(
    #             self._interface_manager.get_current_project(),
    #             tmp_prediction_protein_infos,
    #             tmp_prediction_configuration,
    #             self._interface_manager.project_lock,
    #             self._interface_manager
    #         )
    #         self._interface_manager.job_manager.put_job_into_queue(tmp_prediction_job)
    #         self._interface_manager.add_job_entry_to_job_overview_layout(tmp_prediction_entry_widget)
    #     #     tmp_prediction_task = tasks.Task(
    #     #         target=main_tasks_async.predict_protein_with_colabfold,
    #     #         args=(
    #     #             result[1],
    #     #             result[2],
    #     #             self._interface_manager.get_current_project(),
    #     #             self.custom_progress_signal,
    #     #             self._interface_manager.pymol_lock,
    #     #             self.disable_pymol_signal
    #     #         ),
    #     #         post_func=self.__await_predict_protein_with_colabfold,
    #     #     )
    #     # self._interface_manager.main_tasks_manager.start_prediction_task(tmp_prediction_task)
    #     self._main_view_state.set_proteins_list(self._interface_manager.get_current_project().proteins)
    #     self._main_view_state.set_protein_pairs_list(self._interface_manager.get_current_project().protein_pairs)
    #     self._interface_manager.refresh_main_view()


    # def __await_monomer_prediction_for_subsequent_analysis(self, result: tuple) -> None:
    #     print(result)
    #     tmp_exit_code, tmp_exit_code_description = result
    #     self._interface_manager.stop_wait_cursor()
    #     if tmp_exit_code == exit_codes.EXIT_CODE_ZERO[0]:
    #         if self.active_custom_message_box is not None:
    #             self.active_custom_message_box.close()
    #         # Prediction was successful
    #         constants.PYSSA_LOGGER.info("All structure predictions are done.")
    #         self._interface_manager.status_bar_manager.show_temporary_message("All structure predictions are done.")
    #         constants.PYSSA_LOGGER.info("Begin analysis process.")
    #
    #         # <editor-fold desc="Analysis preperations">
    #         tmp_raw_analysis_run_names: list = []
    #         for row_no in range(self._interface_manager.get_predict_protein_view().ui.list_analysis_overview.count()):
    #             tmp_raw_analysis_run_names.append(
    #                 self._interface_manager.get_predict_protein_view().ui.list_analysis_overview.item(row_no).text())
    #
    #         if os.path.exists(constants.SCRATCH_DIR_ANALYSIS):
    #             shutil.rmtree(constants.SCRATCH_DIR_ANALYSIS)
    #             os.mkdir(constants.SCRATCH_DIR_ANALYSIS)
    #         else:
    #             os.mkdir(constants.SCRATCH_DIR_ANALYSIS)
    #         # </editor-fold>
    #
    #         self._add_new_proteins_to_protein_model()
    #
    #         tmp_distance_analysis_task = tasks.Task(
    #             target=main_tasks_async.run_distance_analysis,
    #             args=(
    #                 tmp_raw_analysis_run_names,
    #                 self._interface_manager.get_current_project(),
    #                 self._interface_manager.get_application_settings(),
    #                 False,
    #                                     self.custom_progress_signal,
    #                 self._interface_manager.pymol_lock,
    #                 self.disable_pymol_signal
    #             ),
    #             post_func=self.__await_run_distance_analysis_after_prediction,
    #         )
    #         self._interface_manager.main_tasks_manager.start_distance_analysis_task(tmp_distance_analysis_task)
    #         self._interface_manager.status_bar_manager.show_permanent_message(
    #             enums.StatusMessages.DISTANCE_ANALYSIS_IS_RUNNING.value
    #         )
    #     elif tmp_exit_code == exit_codes.ERROR_WRITING_FASTA_FILES[0]:
    #         # tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #         #     "Prediction failed because there was an error writing the fasta file(s)!",
    #         #     "Structure Prediction",
    #         #     custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #         # )
    #         # tmp_dialog.exec_()
    #         constants.PYSSA_LOGGER.error(
    #             f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
    #         )
    #         self._interface_manager.status_bar_manager.show_error_message("Prediction failed because there was an error writing the fasta file(s)!")
    #     elif tmp_exit_code == exit_codes.ERROR_FASTA_FILES_NOT_FOUND[0]:
    #         # tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #         #     "Prediction failed because the fasta file(s) could not be found!",
    #         #     "Structure Prediction",
    #         #     custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #         # )
    #         # tmp_dialog.exec_()
    #         constants.PYSSA_LOGGER.error(
    #             f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
    #         )
    #         self._interface_manager.status_bar_manager.show_error_message(
    #             "Prediction failed because the fasta file(s) could not be found!")
    #     elif tmp_exit_code == exit_codes.ERROR_PREDICTION_FAILED[0]:
    #         # tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #         #     "Prediction failed because a subprocess failed!",
    #         #     "Structure Prediction",
    #         #     custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #         # )
    #         # tmp_dialog.exec_()
    #         constants.PYSSA_LOGGER.error(
    #             f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
    #         )
    #         self._interface_manager.status_bar_manager.show_error_message(
    #             "Prediction failed because a subprocess failed!")
    #     elif tmp_exit_code == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
    #         # tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #         #     "Prediction failed because of an unknown error!",
    #         #     "Structure Prediction",
    #         #     custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #         # )
    #         # tmp_dialog.exec_()
    #         constants.PYSSA_LOGGER.error(
    #             f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
    #         )
    #         self._interface_manager.status_bar_manager.show_error_message(
    #             "Prediction failed because of an unknown error!")
    #     self._interface_manager.refresh_main_view()
    #
    # def __await_run_distance_analysis_after_prediction(self, an_exit_code: tuple[int, str, list]) -> None:
    #     """Post process after the analysis thread finished."""
    #     self._interface_manager.status_bar_manager.hide_progress_bar()
    #     constants.PYSSA_LOGGER.debug("__await_run_distance_analysis() started ...")
    #     if an_exit_code[0] == exit_codes.EXIT_CODE_ZERO[0]:
    #         # Analysis was successful
    #         self._add_new_protein_pairs_to_protein_pair_model()
    #         self._active_task = tasks.Task(
    #             target=util_async.unfreeze_pymol_session,
    #             args=(
    #                 self._interface_manager.pymol_session_manager, 0
    #             ),
    #             post_func=self.__await_unfreeze_pymol_session_after_prediction_and_analysis,
    #         )
    #         self._active_task.start()
    #     elif an_exit_code[0] == exit_codes.ERROR_DISTANCE_ANALYSIS_FAILED[0]:
    #         # tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #         #     "Distance analysis failed because there was an error during the analysis!",
    #         #     "Distance Analysis",
    #         #     custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #         # )
    #         # tmp_dialog.exec_()
    #         constants.PYSSA_LOGGER.error(
    #             f"Distance analysis ended with exit code {an_exit_code[0]}: {an_exit_code[1]}",
    #         )
    #         self._interface_manager.status_bar_manager.show_error_message(
    #             "Distance analysis failed!")
    #     elif an_exit_code[0] == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
    #         # tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #         #     "Distance analysis failed because of an unknown error!",
    #         #     "Distance Analysis",
    #         #     custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #         # )
    #         # tmp_dialog.exec_()
    #         constants.PYSSA_LOGGER.error(
    #             f"Distance analysis ended with exit code {an_exit_code[0]}: {an_exit_code[1]}",
    #         )
    #         self._interface_manager.status_bar_manager.show_error_message(
    #             "Distance analysis failed because of an unknown error!")
    #
    # def __await_unfreeze_pymol_session_after_prediction_and_analysis(self):
    #     self._interface_manager.main_tasks_manager.prediction_task = None
    #     self._interface_manager.main_tasks_manager.distance_analysis_task = None
    #     self._interface_manager.refresh_main_view()
    #     self._main_view_state.restore_main_view_state()
    #
    #     # self.active_custom_message_box = custom_message_box.CustomMessageBoxOk(
    #     #     "All structure analysis' are done. \nGo to the Protein Pairs tab to view the new results.",
    #     #     "Distance Analysis",
    #     #     custom_message_box.CustomMessageBoxIcons.INFORMATION.value
    #     # )
    #     # self.active_custom_message_box.exec_()
    #     constants.PYSSA_LOGGER.info("All structure analysis' are done.")
    #     self._interface_manager.status_bar_manager.show_temporary_message("All structure analysis' are done.")
    #     self._interface_manager.stop_wait_cursor()
    #
    # def __await_predict_protein_with_colabfold(self, result: tuple) -> None:
    #     """Process which runs after each prediction job."""
    #     print(result)
    #     self._interface_manager.status_bar_manager.hide_progress_bar()
    #     tmp_exit_code = result[0]
    #     tmp_exit_code_description = result[1]
    #     if tmp_exit_code == exit_codes.ERROR_WRITING_FASTA_FILES[0]:
    #
    #         # tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #         #     "Prediction failed because there was an error writing the fasta file(s)!",
    #         #     "Structure Prediction",
    #         #     custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #         # )
    #         # tmp_dialog.exec_()
    #         constants.PYSSA_LOGGER.error(
    #             f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
    #         )
    #         self._interface_manager.status_bar_manager.show_error_message(
    #             "Prediction failed because there was an error writing the fasta file(s)!")
    #     elif tmp_exit_code == exit_codes.ERROR_FASTA_FILES_NOT_FOUND[0]:
    #         # tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #         #     "Prediction failed because the fasta file(s) could not be found!",
    #         #     "Structure Prediction",
    #         #     custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #         # )
    #         # tmp_dialog.exec_()
    #         constants.PYSSA_LOGGER.error(
    #             f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
    #         )
    #         self._interface_manager.status_bar_manager.show_error_message(
    #             "Prediction failed because the fasta file(s) could not be found!")
    #     elif tmp_exit_code == exit_codes.ERROR_PREDICTION_FAILED[0]:
    #         # tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #         #     "Prediction failed because a subprocess failed!",
    #         #     "Structure Prediction",
    #         #     custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #         # )
    #         # tmp_dialog.exec_()
    #         constants.PYSSA_LOGGER.error(
    #             f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
    #         )
    #         self._interface_manager.status_bar_manager.show_error_message(
    #             "Prediction failed because a subprocess failed!")
    #     elif tmp_exit_code == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
    #         # tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #         #     "Prediction failed because of an unknown error!",
    #         #     "Structure Prediction",
    #         #     custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #         # )
    #         # tmp_dialog.exec_()
    #         constants.PYSSA_LOGGER.error(
    #             f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
    #         )
    #         self._interface_manager.status_bar_manager.show_error_message("Prediction failed because of an unknown error!")
    #     elif tmp_exit_code == exit_codes.EXIT_CODE_ZERO[0]:
    #         # Prediction was successful
    #         self._add_new_proteins_to_protein_model()
    #         self._active_task = tasks.Task(
    #             target=util_async.unfreeze_pymol_session,
    #             args=(
    #                 self._interface_manager.pymol_session_manager, 0
    #             ),
    #             post_func=self.__await_unfreeze_pymol_session_after_prediction,
    #         )
    #         self._active_task.start()
    #     else:
    #         # tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #         #     "Prediction failed because of an unknown case!",
    #         #     "Structure Prediction",
    #         #     custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #         # )
    #         # tmp_dialog.exec_()
    #         self._interface_manager.status_bar_manager.show_error_message("Prediction failed because of an unknown case!")
    #
    # def __await_unfreeze_pymol_session_after_prediction(self):
    #     self._interface_manager.main_tasks_manager.prediction_task = None
    #     self._interface_manager.refresh_main_view()
    #
    #     if self.active_custom_message_box is not None:
    #         self.active_custom_message_box.close()
    #
    #     self._main_view_state.restore_main_view_state()
    #
    #     # self.active_custom_message_box = custom_message_box.CustomMessageBoxOk(
    #     #     "All structure predictions are done.\nGo to the Proteins tab to see the new protein(s).",
    #     #     "Structure Prediction",
    #     #     custom_message_box.CustomMessageBoxIcons.INFORMATION.value
    #     # )
    #     # self.active_custom_message_box.exec_()
    #     constants.PYSSA_LOGGER.info("All structure predictions are done.")
    #     self._interface_manager.status_bar_manager.show_temporary_message("All structure predictions are done.")
    #     self._interface_manager.stop_wait_cursor()

    # def _add_new_proteins_to_protein_model(self):
    #     """Adds the new predicted proteins to the interface manager's protein model."""
    #     tmp_proteins_to_add = self._main_view_state.get_not_matching_proteins(
    #         self._interface_manager.get_current_project().proteins
    #     )
    #     for tmp_protein in tmp_proteins_to_add:
    #         self._interface_manager.add_protein_to_proteins_model(tmp_protein)

    # <editor-fold desc="Multimer">


    # def _post_predict_multimer(self, result: tuple):
    #
    #     # <editor-fold desc="Check if WSL2 and ColabFold are installed">
    #     if globals.g_os == "win32":
    #         constants.PYSSA_LOGGER.info("Checking if WSL2 is installed ...")
    #         if not dialog_settings_global.is_wsl2_installed():
    #             constants.PYSSA_LOGGER.warning("WSL2 is NOT installed.")
    #             self._interface_manager.get_application_settings().wsl_install = 0
    #             tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #                 "Prediction failed because the WSL2 environment is not installed!",
    #                 "Structure Prediction",
    #                 custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #             )
    #             tmp_dialog.exec_()
    #             return
    #         constants.PYSSA_LOGGER.info("Checking if Local Colabfold is installed ...")
    #         if not dialog_settings_global.is_local_colabfold_installed():
    #             constants.PYSSA_LOGGER.warning("Local Colabfold is NOT installed.")
    #             self._interface_manager.get_application_settings().local_colabfold = 0
    #             tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #                 "Prediction failed because the ColabFold is not installed!",
    #                 "Structure Prediction",
    #                 custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #             )
    #             tmp_dialog.exec_()
    #             return
    #
    #     # </editor-fold>
    #
    #     self.prediction_type = constants.PREDICTION_TYPE_PRED_MULTI_ANALYSIS
    #     constants.PYSSA_LOGGER.info("Begin prediction process.")
    #     if result[3] is True:
    #         constants.PYSSA_LOGGER.info("Running prediction with subsequent analysis.")
    #         # Analysis should be run after the prediction
    #         self._active_task = tasks.Task(
    #             target=main_tasks_async.predict_protein_with_colabfold,
    #             args=(
    #                 result[1],
    #                 result[2],
    #                 self._interface_manager.get_current_project(),
    #                 self.custom_progress_signal,
    #                 self._interface_manager.pymol_lock,
    #                 self.disable_pymol_signal
    #             ),
    #             post_func=self.__await_monomer_prediction_for_subsequent_analysis,
    #         )
    #         self._active_task.start()
    #     else:
    #         constants.PYSSA_LOGGER.info("Running only a prediction.")
    #         # No analysis after prediction
    #         self._active_task = tasks.Task(
    #             target=main_tasks_async.predict_protein_with_colabfold,
    #             args=(
    #                 result[1],
    #                 result[2],
    #                 self._interface_manager.get_current_project(),
    #                 self.custom_progress_signal,
    #                 self._interface_manager.pymol_lock,
    #                 self.disable_pymol_signal
    #             ),
    #             post_func=self.__await_predict_protein_with_colabfold,
    #         )
    #         self._active_task.start()
    #
    # def __await_multimer_prediction_for_subsequent_analysis(self, result: tuple) -> None:
    #     tmp_exit_code = result[0]
    #     tmp_exit_code_description = [1]
    #     self._interface_manager.stop_wait_cursor()
    #     if tmp_exit_code == exit_codes.EXIT_CODE_ZERO[0]:
    #         # Prediction was successful
    #         constants.PYSSA_LOGGER.info("All structure predictions are done.")
    #         self.update_status("All structure predictions are done.")
    #         constants.PYSSA_LOGGER.info("Begin analysis process.")
    #         self.update_status("Begin analysis process ...")
    #         tmp_raw_analysis_run_names: list = []
    #         for row_no in range(self._view.ui.list_pred_analysis_multi_overview.count()):
    #             tmp_raw_analysis_run_names.append(self._view.ui.list_pred_analysis_multi_overview.item(row_no).text())
    #
    #         self._active_task = tasks.Task(
    #             target=main_tasks_async.run_distance_analysis,
    #             args=(
    #                 tmp_raw_analysis_run_names,
    #                 self._interface_manager.get_current_project(),
    #                 self._interface_manager.get_application_settings(),
    #                 False,
    #                 self.custom_progress_signal,
    #                 self._interface_manager.pymol_lock,
    #                 self.disable_pymol_signal
    #             ),
    #             post_func=self.__await_run_distance_analysis,
    #         )
    #         self._active_task.start()
    #
    #         if not os.path.exists(constants.SCRATCH_DIR_ANALYSIS):
    #             os.mkdir(constants.SCRATCH_DIR_ANALYSIS)
    #
    #     elif tmp_exit_code == exit_codes.ERROR_WRITING_FASTA_FILES[0]:
    #         # tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #         #     "Prediction failed because there was an error writing the fasta file(s)!",
    #         #     "Structure Prediction",
    #         #     custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #         # )
    #         # tmp_dialog.exec_()
    #         constants.PYSSA_LOGGER.error(
    #             f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
    #         )
    #         self._interface_manager.status_bar_manager.show_error_message(
    #             "Prediction failed because there was an error writing the fasta file(s)!")
    #     elif tmp_exit_code == exit_codes.ERROR_FASTA_FILES_NOT_FOUND[0]:
    #         # tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #         #     "Prediction failed because the fasta file(s) could not be found!",
    #         #     "Structure Prediction",
    #         #     custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #         # )
    #         # tmp_dialog.exec_()
    #         constants.PYSSA_LOGGER.error(
    #             f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
    #         )
    #         self._interface_manager.status_bar_manager.show_error_message(
    #             "Prediction failed because the fasta file(s) could not be found!")
    #     elif tmp_exit_code == exit_codes.ERROR_PREDICTION_FAILED[0]:
    #         # tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #         #     "Prediction failed because a subprocess failed!",
    #         #     "Structure Prediction",
    #         #     custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #         # )
    #         # tmp_dialog.exec_()
    #         constants.PYSSA_LOGGER.error(
    #             f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
    #         )
    #         self._interface_manager.status_bar_manager.show_error_message(
    #             "Prediction failed because a subprocess failed!")
    #     elif tmp_exit_code == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
    #         # tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #         #     "Prediction failed because of an unknown error!",
    #         #     "Structure Prediction",
    #         #     custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #         # )
    #         # tmp_dialog.exec_()
    #         constants.PYSSA_LOGGER.error(
    #             f"Prediction ended with exit code {tmp_exit_code}: {tmp_exit_code_description}",
    #         )
    #         self._interface_manager.status_bar_manager.show_error_message(
    #             "Prediction failed because of an unknown error!")
    #     self._interface_manager.refresh_main_view()

    # </editor-fold>
    # </editor-fold>
    # </editor-fold>

    # <editor-fold desc="Hotspots">
    def __slot_hotspots_protein_regions(self) -> None:
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Hotspots/Protein Regions' clicked.")
        if self._interface_manager.current_tab_index == 1 and self._interface_manager.pymol_session_manager.session_object_type == "protein" and self._interface_manager.pymol_session_manager.is_the_current_protein_in_session(self._interface_manager.get_current_active_protein_object().get_molecule_object()):
            # Proteins tab
            if self._view.ui.lbl_protein_protein_regions.isVisible():
                self._view.ui.lbl_protein_protein_regions.hide()
                self._view.ui.frame_protein_protein_regions.hide()
                self._interface_manager.pymol_session_manager.hide_sequence_view()
            else:
                self._view.ui.lbl_protein_protein_regions.show()
                self._view.ui.frame_protein_protein_regions.show()
                self._interface_manager.pymol_session_manager.show_sequence_view()
        elif self._interface_manager.current_tab_index == 2 and self._interface_manager.pymol_session_manager.session_object_type == "protein_pair" and self._interface_manager.pymol_session_manager.is_the_current_protein_pair_in_session(self._interface_manager.get_current_active_protein_pair_object().name):
            if self._view.ui.lbl_protein_pair_protein_regions.isVisible():
                # Protein Pairs tab
                self._view.ui.lbl_protein_pair_protein_regions.hide()
                self._view.ui.frame_protein_pair_protein_regions.hide()
                self._interface_manager.pymol_session_manager.hide_sequence_view()
            else:
                self._view.ui.lbl_protein_pair_protein_regions.show()
                self._view.ui.frame_protein_pair_protein_regions.show()
                self._interface_manager.pymol_session_manager.show_sequence_view()



    def post_hotspots_protein_regions(self) -> None:
        self._interface_manager.pymol_session_manager.hide_sequence_view()
        self._interface_manager.pymol_session_manager.pymol_interface.select(
            "", "none"
        )

    def __slot_show_protein_regions_resi_sticks(self) -> None:
        """Shows the pymol selection as sticks."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Show' sticks button was clicked.")
        if self._interface_manager.pymol_session_manager.check_if_sele_is_empty():
            return
        self._interface_manager.pymol_session_manager.hide_specific_representation(
            "sticks", "sele and not hydrogens"
        )
        self._interface_manager.pymol_session_manager.pymol_interface.select(
            "sele", "sele and not hydrogens"
        )
        self._interface_manager.pymol_session_manager.color_protein(
            "atomic", "disulfides and not elem C"
        )
        self._interface_manager.pymol_session_manager.pymol_interface.set_custom_setting(
            "valence", 0
        )


    def __slot_hide_protein_regions_resi_sticks(self) -> None:
        """Hides the balls and sticks representation of the pymol selection."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Hide' sticks button was clicked.")
        if self._interface_manager.pymol_session_manager.check_if_sele_is_empty():
            return
        self._interface_manager.pymol_session_manager.hide_specific_representation(
            "sticks", "sele"
        )

    def __slot_show_protein_regions_disulfide_bonds(self) -> None:
        """Shows all disulfid bonds within the pymol session."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Show' disulfide bonds button was clicked.")
        if self._interface_manager.current_tab_index == 1:
            tmp_protein_names: tuple = self._interface_manager.get_current_active_protein_object().get_molecule_object(), 0
        elif self._interface_manager.current_tab_index == 2:
            tmp_protein_pair: "protein_pair.ProteinPair" = self._interface_manager.get_current_active_protein_pair_object()
            tmp_protein_names = (tmp_protein_pair.protein_1.get_molecule_object(), tmp_protein_pair.protein_2.get_molecule_object())
        else:
            return
        tmp_pymol_selection_option: str = "byres (resn CYS and name SG) within 2 of (resn CYS and name SG)"
        for tmp_protein_name in tmp_protein_names:
            if tmp_protein_name != 0:
                self._interface_manager.pymol_session_manager.pymol_interface.select(
                    "disulfides", f"{tmp_protein_name} & {tmp_pymol_selection_option}"
                )
                self._interface_manager.pymol_session_manager.color_protein(
                    "atomic", "disulfides and not elem C"
                )
                self._interface_manager.pymol_session_manager.pymol_interface.set_custom_setting(
                    "valence", 0
                )
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    "sticks", "disulfides"
                )
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    "sticks", "elem H"
                )

    def __slot_hide_protein_regions_disulfide_bonds(self) -> None:
        """Hides all disulfid bonds within the pymol session."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Hide' disulfide bonds button was clicked.")

        if self._interface_manager.current_tab_index == 1:
            tmp_protein_names: tuple = self._interface_manager.get_current_active_protein_object().get_molecule_object(), 0
        elif self._interface_manager.current_tab_index == 2:
            tmp_protein_pair: "protein_pair.ProteinPair" = self._interface_manager.get_current_active_protein_pair_object()
            tmp_protein_names = (tmp_protein_pair.protein_1.get_molecule_object(), tmp_protein_pair.protein_2.get_molecule_object())
        else:
            return
        tmp_pymol_selection_option: str = "byres (resn CYS and name SG) within 2 of (resn CYS and name SG)"
        for tmp_protein_name in tmp_protein_names:
            if tmp_protein_name != 0:
                self._interface_manager.pymol_session_manager.pymol_interface.select(
                    "disulfides", f"{tmp_protein_name} & {tmp_pymol_selection_option}"
                )
                self._interface_manager.pymol_session_manager.hide_specific_representation("sticks", "disulfides")

    def __slot_zoom_protein_regions_resi_position(self) -> None:
        """Zooms to the pymol selection."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Zoom' button was clicked.")
        if self._interface_manager.pymol_session_manager.check_if_sele_is_empty():
            return
        self._interface_manager.pymol_session_manager.zoom_to_residue_in_protein_position("sele")

    # </editor-fold>

    # <editor-fold desc="Settings menu methods">
    def __slot_open_settings_global(self) -> None:
        """Opens the dialog for the global settings."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Settings/Edit' clicked.")
        self._external_controller = settings_view_controller.SettingsViewController(self._interface_manager)
        self._external_controller.user_input.connect(self.post_open_settings_global)
        self._external_controller.restore_ui()
        self._interface_manager.get_settings_view().show()
        # dialog = dialog_settings_global.DialogSettingsGlobal(self._interface_manager)
        # dialog.exec_()
        # self._interface_manager.update_settings()
        # self._workspace_label = QtWidgets.QLabel(f"Current Workspace: {self._workspace_path}")

    def post_open_settings_global(self, return_value):
        self._interface_manager.refresh_workspace_model()
        self._interface_manager.refresh_main_view()
        try:
            tmp_type = self._interface_manager.get_current_protein_tree_index_type()
        except AttributeError:
            return
        # if self._view.ui.cb_protein_atoms.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
        #     self._view.tg_protein_color_atoms.toggle_button.setChecked(True)
        #     self._view.ui.cb_protein_atoms.setChecked(False)
        # elif self._view.tg_protein_color_atoms.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
        #     self._view.tg_protein_color_atoms.toggle_button.setChecked(False)
        #     self._view.ui.cb_protein_atoms.setChecked(True)
        # else:
        #     self._view.tg_protein_color_atoms.toggle_button.setChecked(False)
        #     self._view.ui.cb_protein_atoms.setChecked(False)

        if tmp_type == "chain":
            self._interface_manager.set_index_of_protein_color_combo_box(self._interface_manager.pymol_session_manager)
            self._interface_manager.set_repr_state_in_ui_for_protein_chain(self._interface_manager.pymol_session_manager)
            self._interface_manager.show_protein_pymol_scene_configuration()

    def __slot_restore_settings(self) -> None:
        """Restores the settings.xml file to the default values."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Settings/Restore' clicked.")
        tmp_dialog = custom_message_box.CustomMessageBoxYesNo(
            "Are you sure you want to restore all settings?", "Restore Settings",
            custom_message_box.CustomMessageBoxIcons.INFORMATION.value
        )
        tmp_dialog.exec_()
        if tmp_dialog.response:
            tools.restore_default_settings(self._interface_manager.get_application_settings())
            self._view.status_bar.showMessage("Settings were successfully restored.")
            logging.info("Settings were successfully restored.")
        else:
            self._view.status_bar.showMessage("Settings were not modified.")
            logging.info("Settings were not modified.")

    # </editor-fold>

    # <editor-fold desc="Help menu methods">
    def __slot_arrange_windows(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Help/Arrange Windows' clicked.")
        if not os.path.exists(constants.ARRANGE_WINDOWS_EXE_FILEPATH):
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "The script for arranging the windows could not be found!", "Arrange Windows",
                custom_message_box.CustomMessageBoxIcons.ERROR.value
            )
            tmp_dialog.exec_()
        else:
            logger.debug("Started script to arrange window ...")
            subprocess.Popen([constants.ARRANGE_WINDOWS_EXE_FILEPATH], creationflags=subprocess.CREATE_NO_WINDOW)
            logger.debug("Script to arrange windows finished.")

    def __slot_open_logs(self) -> None:
        """Opens a file explorer with all log files and can open a log file in the default application."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Help/Show Logs in Explorer' clicked.")
        file_dialog = QtWidgets.QFileDialog()
        log_path = str(constants.LOG_PATH)
        file_dialog.setDirectory(log_path)
        file_path, _ = file_dialog.getOpenFileName(self._view, "Select a log file to open", "", "LOG File (*.log)")
        if file_path:
            os.startfile(file_path)

    def __slot_clear_all_log_files(self) -> None:
        """Clears all log files generated under .pyssa/logs."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Help/Clear All Logs' clicked.")
        tmp_dialog = custom_message_box.CustomMessageBoxYesNo(
            "Are you sure you want to delete all log files?",
            "Clear Log Files",
            custom_message_box.CustomMessageBoxIcons.WARNING.value
        )
        tmp_dialog.exec_()
        if tmp_dialog.response:
            try:
                shutil.rmtree(str(constants.LOG_PATH))
            except PermissionError:
                print("The active log file was not deleted.")
            if len(os.listdir(str(constants.LOG_PATH))) == 1:
                # tmp_dialog = custom_message_box.CustomMessageBoxOk(
                #     "All log files could be deleted.", "Clear Log Files",
                #     custom_message_box.CustomMessageBoxIcons.INFORMATION.value
                # )
                # tmp_dialog.exec_()
                self._interface_manager.status_bar_manager.show_temporary_message("All log files could be deleted.")
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
    def __slot_open_tutorial() -> None:
        """Opens the official tutorial pdf file."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Help/Tutorials' clicked.")
        tmp_dialog = dialog_tutorial_videos.TutorialVideosDialog()
        tmp_dialog.exec_()

    @staticmethod
    def open_documentation() -> None:
        """Opens the official plugin documentation as PDF."""
        os.startfile(constants.DOCS_PATH)

    @staticmethod
    def __slot_open_about() -> None:
        """Opens the About dialog."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Help/About' clicked.")
        dialog = dialog_about.DialogAbout()
        dialog.exec_()

    def __slot_get_demo_projects(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Help/Get Demo Projects' clicked.")
        self._interface_manager.status_bar_manager.show_temporary_message(
            "Getting demo projects ...", False)
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
                self._interface_manager.status_bar_manager.show_error_message("The download of the demo projects failed.")
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
                self._interface_manager.status_bar_manager.show_error_message("Extraction process of demo projects finished with an error.")
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
            self._interface_manager.status_bar_manager.show_error_message("Import process of demo projects finished with an error.")
        else:
            self._interface_manager.refresh_workspace_model()
            self._interface_manager.refresh_main_view()
            # tmp_dialog = custom_message_box.CustomMessageBoxOk(
            #     "Getting demo projects finished successfully.", "Get Demo Projects",
            #     custom_message_box.CustomMessageBoxIcons.INFORMATION.value
            # )
            # tmp_dialog.exec_()
            self._interface_manager.status_bar_manager.show_temporary_message("Getting demo projects finished successfully.")

    # </editor-fold>

    # <editor-fold desc="Image menu methods">
    def __slot_preview_image(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Image/Preview' clicked.")
            self._active_task = tasks.Task(
                target=image_async.preview_image,
                args=(self._interface_manager.pymol_session_manager, 0),
                post_func=self.__await_preview_image,
            )
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
        else:
            self._active_task.start()
            self.update_status("Creating preview of image ...")
            self._interface_manager.start_wait_cursor()

    def __await_preview_image(self, return_value: tuple):
        self._interface_manager.stop_wait_cursor()
        self._interface_manager.refresh_main_view()
        self.update_status("Preview finished.")

    def __slot_create_ray_traced_image(self) -> None:
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Image/Ray-Traced' clicked.")
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

            # --- New job approach
            tmp_session_filepath = self._interface_manager.pymol_session_manager.save_current_pymol_session_as_pse_cache_file()
            tmp_ray_tracing_job, tmp_ray_tracing_entry_widget = self._interface_manager.job_manager.create_ray_tracing_job(
                full_file_name[0],
                tmp_session_filepath,
                self._interface_manager.get_application_settings().image_ray_trace_mode,
                self._interface_manager.get_application_settings().image_ray_texture,
                self._interface_manager.get_application_settings().image_renderer,
                self._interface_manager,
                self._interface_manager.get_current_project().get_project_name()
            )
            self._interface_manager.job_manager.put_job_into_queue(tmp_ray_tracing_job)
            self._interface_manager.add_job_entry_to_job_overview_layout(tmp_ray_tracing_entry_widget)
        #
        #
        # #self.add_job_entry_to_overview("Ray-tracing for scene 6OMNvsBMP2", "BMP2-validation")
        #
        # tmp_session_filepath = self._interface_manager.pymol_session_manager.save_current_pymol_session_as_pse_cache_file()
        # self._active_task = tasks.Task(
        #     target=main_presenter_async.create_ray_traced_image,
        #     args=(full_file_name[0], self._interface_manager.get_application_settings(), tmp_session_filepath),
        #     post_func=self.__await_create_ray_traced_image,
        # )
        # self._active_task.start()
        # self.update_status("Creating ray-traced image ...")
        # print(os.cpu_count())
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    def __await_create_ray_traced_image(self, return_value: tuple) -> None:
        self._interface_manager.stop_wait_cursor()
        self._interface_manager.refresh_main_view()
        self.update_status("Image creation finished.")

    def __slot_create_drawn_image(self) -> None:
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Image/Simple' clicked.")
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
                target=image_async.create_drawn_image,
                args=(full_file_name[0], self._interface_manager.pymol_session_manager),
                post_func=self.__await_create_drawn_image,
            )
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
        else:
            self._active_task.start()
            self.update_status("Creating simple image ...")
            self._interface_manager.start_wait_cursor()

    def __await_create_drawn_image(self, return_value: tuple) -> None:
        self._interface_manager.stop_wait_cursor()
        self._interface_manager.refresh_main_view()
        self.update_status("Image creation finished.")

    # </editor-fold>

    # <editor-fold desc="Sequences tab methods">
    def __slot_open_text_editor_for_seq(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "A sequence in the 'Addition Information' table was clicked.")
        if self._view.ui.seqs_table_widget.currentColumn() == 1 and self._view.ui.seqs_table_widget.currentRow() == 0:
            self.__slot_rename_selected_sequence()
        elif self._view.ui.seqs_table_widget.currentColumn() == 1 and self._view.ui.seqs_table_widget.currentRow() == 1:
            self.tmp_txt_browser = QtWidgets.QTextBrowser()
            try:
                tmp_seq = self._view.ui.seqs_table_widget.currentItem().data(enums.ModelEnum.OBJECT_ROLE).seq
                tmp_seqs = tmp_seq.split(",")
                tmp_seq = ",\n\n".join(tmp_seqs)
                self.tmp_txt_browser.setText(tmp_seq)
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
            self._view.ui.seqs_table_widget.item(0, 1).setText(tmp_new_seq)
        except AttributeError:
            return

    def _rename_sequence(self):
        tmp_old_name = self._view.ui.seqs_list_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE).name
        try:
            # this is needed because the signal is fired even if the current item is None
            tmp_new_name = self._view.ui.seqs_table_widget.item(0, 1).text()
        except AttributeError:
            return
        else:
            tmp_seq = self._view.ui.seqs_list_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE).seq
            self._view.ui.seqs_list_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE).name = tmp_new_name
            self._view.ui.seqs_list_view.model().setData(self._view.ui.seqs_list_view.currentIndex(), tmp_new_name, Qt.DisplayRole)
            tmp_database_operation = database_operation.DatabaseOperation(
                enums.SQLQueryType.UPDATE_SEQUENCE_NAME, (0, tmp_new_name, tmp_old_name, tmp_seq)
            )
            self._database_thread.put_database_operation_into_queue(tmp_database_operation)

    def __slot_show_sequence_information(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, f"The sequence '{self._view.ui.seqs_list_view.currentIndex().data(Qt.DisplayRole)}' on the 'Sequence Tab' was clicked.")
        self._interface_manager.show_sequence_parameters(
            self._view.ui.seqs_list_view.currentIndex()
        )
        self._view.ui.btn_save_sequence.setEnabled(True)
        self._view.ui.btn_delete_sequence.setEnabled(True)

    def __slot_import_sequence(self) -> None:
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Import sequence' button on the 'Sequence Tab' was clicked.")
        self._external_controller = import_sequence_view_controller.ImportSequenceViewController(self._interface_manager)
        self._external_controller.user_input.connect(self._post_import_sequence)
        self._external_controller.restore_ui()
        self._interface_manager.get_import_sequence_view().show()

    def _post_import_sequence(self, return_value: tuple):
        # tmp_fasta_filepath, _ = return_value
        # with open(tmp_fasta_filepath, "r") as handle:
        #     for tmp_record in SeqIO.parse(handle, "fasta"):
        #         # Append each SeqRecord object to the list
        #         self._interface_manager.get_current_project().sequences.append(tmp_record)
        #         self._database_manager.insert_new_sequence(tmp_record)
        # self._interface_manager.refresh_sequence_model()
        # self._interface_manager.refresh_main_view()
        for tmp_seq_record in return_value[1]:
            logger.info(f"Adding new sequence {tmp_seq_record.name} with {tmp_seq_record.seq} to the current project.")
            self._interface_manager.get_current_project().sequences.append(tmp_seq_record)
            tmp_database_operation = database_operation.DatabaseOperation(
                enums.SQLQueryType.INSERT_NEW_SEQUENCE,
                (0, tmp_seq_record)
            )
            self._database_thread.put_database_operation_into_queue(tmp_database_operation)
        self._interface_manager.refresh_sequence_model()
        self._interface_manager.show_menu_options_with_seq()
        self._interface_manager.refresh_main_view()

    def __slot_add_sequence(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Add sequence' button on the 'Sequence Tab' was clicked.")
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
        #self._interface_manager.show_menu_options_with_seq()
        self._interface_manager.refresh_main_view()
        #self._show_temporary_message("Adding a sequence was successful", "A prediction is currently running ...")

    def __slot_save_selected_sequence_as_fasta_file(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Export sequence' button on the 'Sequence Tab' was clicked.")
        self._interface_manager.start_wait_cursor()
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
                target=sequence_async.save_selected_protein_sequence_as_fasta_file,
                args=(
                    tmp_seq_record,
                    file_path,
                    self._interface_manager.get_current_project().get_database_filepath()
                ),
                post_func=self.__await_save_selected_sequence_as_fasta_file,
            )
            self._active_task.start()
        else:
            self._interface_manager.stop_wait_cursor()
            self._interface_manager.refresh_main_view()

    def __await_save_selected_sequence_as_fasta_file(self, result: tuple):
        self._interface_manager.stop_wait_cursor()
        if result[0] == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Saving the sequence as .fasta file failed!",
                "Save Protein Sequence",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
        elif result[0] == exit_codes.EXIT_CODE_ZERO[0]:
            # tmp_dialog = custom_message_box.CustomMessageBoxOk(
            #     "The sequence was successfully saved as .fasta file.",
            #     "Save Protein Sequence",
            #     custom_message_box.CustomMessageBoxIcons.INFORMATION.value
            # )
            # tmp_dialog.exec_()
            self._interface_manager.status_bar_manager.show_temporary_message("The sequence was successfully saved as .fasta file.")
        else:
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Saving the sequence as .fasta file failed with an unexpected error!",
                "Save Protein Sequence",
                custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
            )
            tmp_dialog.exec_()
        self._interface_manager.refresh_sequence_model()
        self._interface_manager.refresh_main_view()

    def __slot_delete_selected_sequence(self):
        # popup message which warns the user that the selected sequence gets deleted
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Delete sequence' button on the 'Sequence Tab' was clicked.")
        tmp_dialog = custom_message_box.CustomMessageBoxDelete(
            "Are you sure you want to delete this sequence?",
            "Delete Sequence",
            custom_message_box.CustomMessageBoxIcons.WARNING.value
        )
        tmp_dialog.exec_()
        response: bool = tmp_dialog.response

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
            # hide not available menu labels
        else:
            constants.PYSSA_LOGGER.info("No sequence has been deleted. No changes were made.")

    def __slot_rename_selected_sequence(self) -> None:
        """Opens a new view to rename the selected sequence."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Rename sequence' context menu action was clicked.")
        self._external_controller = rename_sequence_view_controller.RenameSequenceViewController(self._interface_manager)
        self._external_controller.user_input.connect(self.post_rename_selected_sequence_structure)
        self._external_controller.restore_ui()
        self._interface_manager.get_rename_sequence_view().show()

    def post_rename_selected_sequence_structure(self, return_value: tuple):
        tmp_new_name = return_value[0]
        tmp_old_name = self._view.ui.seqs_list_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE).name
        tmp_seq = self._view.ui.seqs_list_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE).seq
        self._view.ui.seqs_list_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE).name = tmp_new_name
        self._view.ui.seqs_list_view.model().setData(self._view.ui.seqs_list_view.currentIndex(), tmp_new_name,
                                                     Qt.DisplayRole)
        tmp_database_operation = database_operation.DatabaseOperation(
            enums.SQLQueryType.UPDATE_SEQUENCE_NAME, (0, tmp_new_name, tmp_old_name, tmp_seq)
        )
        self._database_thread.put_database_operation_into_queue(tmp_database_operation)
        self._view.ui.seqs_table_widget.item(0, 1).setText(tmp_new_name)
        self._interface_manager.refresh_main_view()

    def open_context_menu_for_sequences(self, position):
        tmp_context_menu = self._sequence_list_context_menu.get_context_menu(
            self._view.ui.seqs_list_view.selectedIndexes()
        )
        tmp_context_menu.exec_(self._view.ui.seqs_list_view.viewport().mapToGlobal(position))
        # sequence_context_menu = QtWidgets.QMenu()
        #
        # self.sequences_context_menu_rename_action = sequence_context_menu.addAction(
        #     self._view.tr("Rename selected sequence")
        # )
        # self.sequences_context_menu_rename_action.triggered.connect(self.rename_selected_sequence)
        #
        # sequence_context_menu.exec_(self._view.ui.seqs_list_view.viewport().mapToGlobal(position))

    # </editor-fold>

    # <editor-fold desc="Proteins tab methods">
    def open_context_menu_for_proteins(self, position):
        try:
            tmp_protein = self._interface_manager.get_current_active_protein_object()
        except ValueError:
            tmp_is_protein_in_any_pair_flag = True
        else:
            tmp_is_protein_in_any_pair_flag = self._interface_manager.get_current_project().check_if_protein_is_in_any_protein_pair(
                tmp_protein.get_molecule_object()
            )
        tmp_context_menu = self._protein_tree_context_menu.get_context_menu(
            self._view.ui.proteins_tree_view.selectedIndexes(),
            self._interface_manager.get_current_protein_tree_index_type(),
            tmp_is_protein_in_any_pair_flag,
            self._interface_manager.pymol_session_manager.is_the_current_protein_in_session(self._interface_manager.get_current_active_protein_object().get_molecule_object())
        )
        tmp_context_menu.exec_(self._view.ui.proteins_tree_view.viewport().mapToGlobal(position))

    def __slot_open_protein_pymol_session(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Open protein pymol session' button on the 'Proteins Tab' was clicked.")
            self._view.tg_protein_white_bg.toggle_button.setCheckState(False)
            tmp_protein: "protein.Protein" = self._interface_manager.get_current_active_protein_object()
            tmp_flag = False
            self._active_task = tasks.Task(
                target=pymol_session_async.load_protein_pymol_session,
                args=(
                    tmp_protein,
                    self._interface_manager.pymol_session_manager,
                    tmp_flag
                ),
                post_func=self.__await_open_protein_pymol_session,
            )
            self._interface_manager.status_bar_manager.show_temporary_message(
                f"Loading PyMOL session of {tmp_protein.get_molecule_object()} ...", False
            )
            self._interface_manager.start_wait_cursor()
        except Exception as e:
            logger.error(f"The error {e} occurred during the pymol session loading!")
            self._interface_manager.status_bar_manager.show_error_message("Loading the PyMOL session failed!")
            self._interface_manager.stop_wait_cursor()
        else:
            self._active_task.start()

    def __await_open_protein_pymol_session(self, return_value: tuple):
        try:
            logger.debug("Returning from async function.")
            tmp_pymol_session_manager, exit_boolean = return_value
            self._interface_manager.pymol_session_manager = tmp_pymol_session_manager
            self._view.ui.action_protein_regions.setEnabled(False)
            if exit_boolean:
                self._view.cb_chain_color.setEnabled(True)
                self._view.cb_chain_representation.setEnabled(True)
                self._view.ui.action_protein_regions.setEnabled(True)
                self._view.ui.btn_create_protein_scene.setEnabled(True)
                self._view.ui.btn_update_protein_scene.setEnabled(True)
                self._view.ui.lbl_session_name.setText(f"Session Name: {self._interface_manager.pymol_session_manager.session_name}")
                self._view.ui.lbl_pymol_protein_scene.setText(f"PyMOL Scene: base")
                self._view.ui.lbl_info.setText("Please select a chain.")
                logger.info("Successfully opened protein session.")
                self._interface_manager.status_bar_manager.show_temporary_message("Loading the PyMOL session was successful.")
            else:
                logger.error("The protein name could not be found in the object list in PyMOL!")
                self._view.cb_chain_color.setEnabled(False)
                self._view.cb_chain_representation.setEnabled(False)
                self._view.ui.btn_create_protein_scene.setEnabled(False)
                self._view.ui.btn_update_protein_scene.setEnabled(False)
                self._interface_manager.status_bar_manager.show_error_message("Loading the PyMOL session failed!")
                self._view.ui.lbl_info.setText("Please load the PyMOL session of the selected protein.")
        except Exception as e:
            logger.error(f"The error {e} occurred during the pymol session loading!")
            self._interface_manager.status_bar_manager.show_error_message("Loading the PyMOL session failed!")
        finally:
            self._interface_manager.refresh_main_view()
            self._interface_manager.stop_wait_cursor()

    def __slot_get_information_about_selected_object_in_protein_branch(self) -> None:
        try:
            tmp_type = self._interface_manager.get_current_protein_tree_index_type()

            if tmp_type == "protein":
                logger.log(
                    log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                    f"The protein object '{self._view.ui.proteins_tree_view.currentIndex().data(Qt.DisplayRole)}' on the 'Proteins Tab' was clicked."
                )

            elif tmp_type == "scene":
                tmp_scene_name = self._view.ui.proteins_tree_view.currentIndex().data(Qt.DisplayRole)
                tmp_protein_name = self._view.ui.proteins_tree_view.currentIndex().parent().parent().data(Qt.DisplayRole)
                logger.log(
                    log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                    f"The scene '{tmp_scene_name}' of the protein '{tmp_protein_name}' on the 'Proteins Tab' was clicked."
                )
                if self._interface_manager.pymol_session_manager.is_the_current_protein_in_session(self._interface_manager.get_current_active_protein_object().get_molecule_object()):
                    tmp_scene_name = self._interface_manager.get_current_active_scene_name()
                    self._interface_manager.pymol_session_manager.current_scene_name = tmp_scene_name
                    ui_util.set_pymol_scene_name_into_label(self._interface_manager.pymol_session_manager.current_scene_name,
                                                            self._view.ui.lbl_pymol_protein_scene)
                    self._interface_manager.pymol_session_manager.load_scene(tmp_scene_name)

            elif tmp_type == "chain":
                tmp_chain_letter = self._view.ui.proteins_tree_view.currentIndex().data(Qt.DisplayRole)
                tmp_protein_name = self._view.ui.proteins_tree_view.currentIndex().parent().parent().data(Qt.DisplayRole)
                logger.log(
                    log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                    f"The chain object '{tmp_chain_letter}' of the protein '{tmp_protein_name}' on the 'Proteins Tab' was clicked."
                )
                if self._interface_manager.pymol_session_manager.current_scene_name != "" and self._interface_manager.pymol_session_manager.is_the_current_protein_in_session(self._interface_manager.get_current_active_protein_object().get_molecule_object()):
                    self.set_icon_for_current_color_in_proteins_tab()
                    self._interface_manager.set_index_of_protein_color_combo_box(self._interface_manager.pymol_session_manager)
                    self._interface_manager.set_repr_state_in_ui_for_protein_chain(self._interface_manager.pymol_session_manager)
                    self._main_view_state.selected_chain_proteins = self._interface_manager.get_current_protein_tree_index()

            elif tmp_type == "header":
                logger.log(
                    log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                    f"The header '{self._view.ui.proteins_tree_view.currentIndex().data(Qt.DisplayRole)}' on the 'Proteins Tab' was clicked."
                )

            else:
                logger.warning("Unknown object type occurred in Protein tab.")
                return

            self._interface_manager.manage_ui_of_protein_tab(
                tmp_type,
                self._interface_manager.get_current_project().check_if_protein_is_in_any_protein_pair(
                    self._interface_manager.get_current_active_protein_object().get_molecule_object()
                ),
                self._interface_manager.pymol_session_manager
            )
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    def __slot_change_chain_color_proteins(self) -> None:
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                   "The 'Color' attribute of a protein chain on the 'Proteins Tab' changed.")
        tmp_protein = self._interface_manager.get_current_active_protein_object()
        tmp_chain = self._interface_manager.get_current_active_chain_object()
        tmp_color = self._view.ui.lbl_protein_current_color.text().strip()
        # fixme : below is legacy code!
        # if self._interface_manager.get_settings_manager().settings.proteins_tab_use_combobox_for_colors == 1:
        #     tmp_color = self._view.ui.box_protein_color.currentText()
        # else:

        if self._interface_manager.pymol_session_manager.session_object_type == "protein" and self._interface_manager.pymol_session_manager.session_name == tmp_protein.get_molecule_object():
            # Update pymol parameter in PyMOL
            tmp_protein.pymol_selection.set_selection_for_a_single_chain(tmp_chain.chain_letter)
            self._interface_manager.pymol_session_manager.color_protein(
                tmp_color, tmp_protein.pymol_selection.selection_string
            )
            # Update pymol parameter in memory
            tmp_chain.pymol_parameters["chain_color"] = tmp_color
            # Update pymol parameter in database
            with database_manager.DatabaseManager(str(self._interface_manager.get_current_project().get_database_filepath())) as db_manager:
                db_manager.update_protein_chain_color(tmp_chain.get_id(), tmp_color)
            self._update_scene()
            self._save_protein_pymol_session()
        else:
            logger.warning("The color of a protein chain could not be changed. This can be due to UI setup reasons.")

    def __slot_change_chain_color_proteins_atoms(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Color atoms by element' button on the 'Proteins Tab' was clicked.")
        tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object().chain_letter)
        self._interface_manager.pymol_session_manager.color_protein(
            "atomic", f"{tmp_selection.selection_string} and not elem C"
        )
        self._interface_manager.pymol_session_manager.color_protein(
            "grey70", f"{tmp_selection.selection_string} and elem C"
        )
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("grey70    ")
        self._view.color_grid_proteins.c_grey_70.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_grey_70.setIconSize(
            self._view.color_grid_proteins.c_grey_70.icon().actualSize(QtCore.QSize(14, 14)))

    def __slot_change_chain_reset_proteins_atoms(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                   "'Reset color atoms by element' button on the 'Proteins Tab' was clicked.")
        tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
        tmp_selection.set_selection_for_a_single_chain(
            self._interface_manager.get_current_active_chain_object().chain_letter)
        self._interface_manager.pymol_session_manager.color_protein(
            self._view.color_grid_proteins.last_clicked_color, f"{tmp_selection.selection_string}"
        )
        self.set_icon_for_current_color_in_proteins_tab()

    def __slot_protein_change_background_color(self):
        try:
            if self._view.tg_protein_white_bg.toggle_button.isChecked():
                self._interface_manager.pymol_session_manager.pymol_interface.set_background_color("white")
            else:
                self._interface_manager.pymol_session_manager.pymol_interface.set_background_color("black")
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # <editor-fold desc="Color Grid slot methods">
    def set_icon_for_current_color_in_proteins_tab(self):
        color_index_functions = {
            "red": self.set_color_name_in_label_red_in_proteins_tab,
            "tv_red": self.set_color_name_in_label_tv_red_in_proteins_tab,
            "salmon": self.set_color_name_in_label_salmon_in_proteins_tab,
            "raspberry": self.set_color_name_in_label_raspberry_in_proteins_tab,
            "green": self.set_color_name_in_label_green_in_proteins_tab,
            "tv_green": self.set_color_name_in_label_tv_green_in_proteins_tab,
            "palegreen": self.set_color_name_in_label_palegreen_in_proteins_tab,
            "forest": self.set_color_name_in_label_forest_in_proteins_tab,
            "blue": self.set_color_name_in_label_blue_in_proteins_tab,
            "tv_blue": self.set_color_name_in_label_tv_blue_in_proteins_tab,
            "lightblue": self.set_color_name_in_label_lightblue_in_proteins_tab,
            "skyblue": self.set_color_name_in_label_skyblue_in_proteins_tab,
            "yellow": self.set_color_name_in_label_yellow_in_proteins_tab,
            "tv_yellow": self.set_color_name_in_label_tv_yellow_in_proteins_tab,
            "paleyellow": self.set_color_name_in_label_paleyellow_in_proteins_tab,
            "sand": self.set_color_name_in_label_sand_in_proteins_tab,
            "magenta": self.set_color_name_in_label_magenta_in_proteins_tab,
            "purple": self.set_color_name_in_label_purple_in_proteins_tab,
            "pink": self.set_color_name_in_label_pink_in_proteins_tab,
            "hotpink": self.set_color_name_in_label_hotpink_in_proteins_tab,
            "cyan": self.set_color_name_in_label_cyan_in_proteins_tab,
            "aquamarine": self.set_color_name_in_label_aquamarine_in_proteins_tab,
            "palecyan": self.set_color_name_in_label_palecyan_in_proteins_tab,
            "teal": self.set_color_name_in_label_teal_in_proteins_tab,
            "orange": self.set_color_name_in_label_orange_in_proteins_tab,
            "tv_orange": self.set_color_name_in_label_tv_orange_in_proteins_tab,
            "lightorange": self.set_color_name_in_label_lightorange_in_proteins_tab,
            "olive": self.set_color_name_in_label_olive_in_proteins_tab,
            "white": self.set_color_name_in_label_white_in_proteins_tab,
            "grey70": self.set_color_name_in_label_grey_70_in_proteins_tab,
            "grey30": self.set_color_name_in_label_grey_30_in_proteins_tab,
            "black": self.set_color_name_in_label_black_in_proteins_tab,
        }
        tmp_protein = self._interface_manager.get_current_active_protein_object()
        tmp_chain = self._interface_manager.get_current_active_chain_object()
        tmp_protein.pymol_selection.set_selection_for_a_single_chain(tmp_chain.chain_letter)
        color_index_functions[
            tmp_chain.get_color(tmp_protein.pymol_selection.selection_string,
                                self._interface_manager.pymol_session_manager)[0]
        ]()

    def reset_icon_for_last_color_in_proteins_tab(self):
        tmp_color_name = self._view.ui.lbl_protein_current_color.text()
        if tmp_color_name == "":
            return
        color_index_functions = {
            "red": self.reset_icon_for_red_in_proteins_tab,
            "tv_red": self.reset_icon_for_tv_red_in_proteins_tab,
            "salmon": self.reset_icon_for_salmon_in_proteins_tab,
            "raspberry": self.reset_icon_for_raspberry_in_proteins_tab,
            "green": self.reset_icon_for_green_in_proteins_tab,
            "tv_green": self.reset_icon_for_tv_green_in_proteins_tab,
            "palegreen": self.reset_icon_for_palegreen_in_proteins_tab,
            "forest": self.reset_icon_for_forest_in_proteins_tab,
            "blue": self.reset_icon_for_blue_in_proteins_tab,
            "tv_blue": self.reset_icon_for_tv_blue_in_proteins_tab,
            "lightblue": self.reset_icon_for_lightblue_in_proteins_tab,
            "skyblue": self.reset_icon_for_skyblue_in_proteins_tab,
            "yellow": self.reset_icon_for_yellow_in_proteins_tab,
            "tv_yellow": self.reset_icon_for_tv_yellow_in_proteins_tab,
            "paleyellow": self.reset_icon_for_paleyellow_in_proteins_tab,
            "sand": self.reset_icon_for_sand_in_proteins_tab,
            "magenta": self.reset_icon_for_magenta_in_proteins_tab,
            "purple": self.reset_icon_for_purple_in_proteins_tab,
            "pink": self.reset_icon_for_pink_in_proteins_tab,
            "hotpink": self.reset_icon_for_hotpink_in_proteins_tab,
            "cyan": self.reset_icon_for_cyan_in_proteins_tab,
            "aquamarine": self.reset_icon_for_aquamarine_in_proteins_tab,
            "palecyan": self.reset_icon_for_palecyan_in_proteins_tab,
            "teal": self.reset_icon_for_teal_in_proteins_tab,
            "orange": self.reset_icon_for_orange_in_proteins_tab,
            "tv_orange": self.reset_icon_for_tv_orange_in_proteins_tab,
            "lightorange": self.reset_icon_for_lightorange_in_proteins_tab,
            "olive": self.reset_icon_for_olive_in_proteins_tab,
            "white": self.reset_icon_for_white_in_proteins_tab,
            "grey70": self.reset_icon_for_grey_70_in_proteins_tab,
            "grey30": self.reset_icon_for_grey_30_in_proteins_tab,
            "black": self.reset_icon_for_black_in_proteins_tab
        }
        color_index_functions[tmp_color_name.strip()]()

    # <editor-fold desc="Set color and icon">
    def set_color_name_in_label_red_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("red    ")
        self._view.color_grid_proteins.last_clicked_color = "red"
        self._view.color_grid_proteins.c_red.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_red.setIconSize(self._view.color_grid_proteins.c_red.icon().actualSize(QtCore.QSize(14, 14)))
        self.__slot_change_chain_color_proteins()

    def set_color_name_in_label_tv_red_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("tv_red    ")
        self._view.color_grid_proteins.last_clicked_color = "tv_red"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_tv_red.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_tv_red.setIconSize(self._view.color_grid_proteins.c_tv_red.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_salmon_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("salmon    ")
        self._view.color_grid_proteins.last_clicked_color = "salmon"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_salomon.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_salomon.setIconSize(self._view.color_grid_proteins.c_salomon.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_raspberry_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("raspberry    ")
        self._view.color_grid_proteins.last_clicked_color = "raspberry"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_raspberry.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_raspberry.setIconSize(self._view.color_grid_proteins.c_raspberry.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_green_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("green    ")
        self._view.color_grid_proteins.last_clicked_color = "green"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_green.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_green.setIconSize(
            self._view.color_grid_proteins.c_green.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_tv_green_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("tv_green    ")
        self._view.color_grid_proteins.last_clicked_color = "tv_green"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_tv_green.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_tv_green.setIconSize(
            self._view.color_grid_proteins.c_tv_green.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_palegreen_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("palegreen    ")
        self._view.color_grid_proteins.last_clicked_color = "palegreen"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_palegreen.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_palegreen.setIconSize(
            self._view.color_grid_proteins.c_palegreen.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_forest_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("forest    ")
        self._view.color_grid_proteins.last_clicked_color = "forest"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_forest.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_forest.setIconSize(
            self._view.color_grid_proteins.c_forest.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_blue_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("blue    ")
        self._view.color_grid_proteins.last_clicked_color = "blue"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_blue.setIcon(QtGui.QIcon(
            ":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_blue.setIconSize(self._view.color_grid_proteins.c_blue.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_tv_blue_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("tv_blue    ")
        self._view.color_grid_proteins.last_clicked_color = "tv_blue"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_tv_blue.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_tv_blue.setIconSize(
            self._view.color_grid_proteins.c_tv_blue.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_lightblue_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("lightblue    ")
        self._view.color_grid_proteins.last_clicked_color = "lightblue"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_lightblue.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_lightblue.setIconSize(
            self._view.color_grid_proteins.c_lightblue.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_skyblue_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("skyblue    ")
        self._view.color_grid_proteins.last_clicked_color = "skyblue"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_skyblue.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_skyblue.setIconSize(
            self._view.color_grid_proteins.c_skyblue.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_yellow_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("yellow    ")
        self._view.color_grid_proteins.last_clicked_color = "yellow"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_yellow.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_yellow.setIconSize(
            self._view.color_grid_proteins.c_yellow.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_tv_yellow_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("tv_yellow    ")
        self._view.color_grid_proteins.last_clicked_color = "tv_yellow"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_tv_yellow.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_tv_yellow.setIconSize(
            self._view.color_grid_proteins.c_tv_yellow.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_paleyellow_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("paleyellow    ")
        self._view.color_grid_proteins.last_clicked_color = "paleyellow"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_paleyellow.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_paleyellow.setIconSize(
            self._view.color_grid_proteins.c_paleyellow.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_sand_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("sand    ")
        self._view.color_grid_proteins.last_clicked_color = "sand"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_sand.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_sand.setIconSize(
            self._view.color_grid_proteins.c_sand.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_magenta_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("magenta    ")
        self._view.color_grid_proteins.last_clicked_color = "magenta"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_magenta.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_magenta.setIconSize(
            self._view.color_grid_proteins.c_magenta.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_purple_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("purple    ")
        self._view.color_grid_proteins.last_clicked_color = "purple"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_purple.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_purple.setIconSize(
            self._view.color_grid_proteins.c_purple.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_pink_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("pink    ")
        self._view.color_grid_proteins.last_clicked_color = "pink"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_pink.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_pink.setIconSize(
            self._view.color_grid_proteins.c_pink.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_hotpink_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("hotpink    ")
        self._view.color_grid_proteins.last_clicked_color = "hotpink"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_hotpink.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_hotpink.setIconSize(
            self._view.color_grid_proteins.c_hotpink.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_cyan_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("cyan    ")
        self._view.color_grid_proteins.last_clicked_color = "cyan"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_cyan.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_cyan.setIconSize(
            self._view.color_grid_proteins.c_cyan.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_aquamarine_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("aquamarine    ")
        self._view.color_grid_proteins.last_clicked_color = "aquamarine"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_aquamarine.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_aquamarine.setIconSize(
            self._view.color_grid_proteins.c_aquamarine.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_palecyan_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("palecyan    ")
        self._view.color_grid_proteins.last_clicked_color = "palecyan"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_palecyan.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_palecyan.setIconSize(
            self._view.color_grid_proteins.c_palecyan.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_teal_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("teal    ")
        self._view.color_grid_proteins.last_clicked_color = "teal"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_teal.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_teal.setIconSize(
            self._view.color_grid_proteins.c_teal.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_orange_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("orange    ")
        self._view.color_grid_proteins.last_clicked_color = "orange"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_orange.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_orange.setIconSize(
            self._view.color_grid_proteins.c_orange.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_tv_orange_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("tv_orange    ")
        self._view.color_grid_proteins.last_clicked_color = "tv_orange"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_tv_orange.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_tv_orange.setIconSize(
            self._view.color_grid_proteins.c_tv_orange.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_lightorange_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("lightorange    ")
        self._view.color_grid_proteins.last_clicked_color = "lightorange"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_lightorange.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_lightorange.setIconSize(
            self._view.color_grid_proteins.c_lightorange.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_olive_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("olive    ")
        self._view.color_grid_proteins.last_clicked_color = "olive"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_olive.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_olive.setIconSize(
            self._view.color_grid_proteins.c_olive.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_white_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("white    ")
        self._view.color_grid_proteins.last_clicked_color = "white"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_white.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_white.setIconSize(
            self._view.color_grid_proteins.c_white.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_grey_70_in_proteins_tab(self):
        tmp_protein = self._interface_manager.get_current_active_protein_object()
        tmp_chain = self._interface_manager.get_current_active_chain_object()
        tmp_protein.pymol_selection.set_selection_for_a_single_chain(tmp_chain.chain_letter)
        tmp_is_colored_by_elements: bool = tmp_chain.get_color(tmp_protein.pymol_selection.selection_string,
                                                               self._interface_manager.pymol_session_manager)[1]
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("grey70    ")
        self.__slot_change_chain_color_proteins()
        if tmp_is_colored_by_elements:
            self.__slot_change_chain_color_proteins_atoms()
        else:
            self._view.color_grid_proteins.last_clicked_color = "grey70"
        self._view.color_grid_proteins.c_grey_70.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_grey_70.setIconSize(
            self._view.color_grid_proteins.c_grey_70.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_grey_30_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("grey30    ")
        self._view.color_grid_proteins.last_clicked_color = "grey30"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_grey_30.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_grey_30.setIconSize(
            self._view.color_grid_proteins.c_grey_30.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_black_in_proteins_tab(self):
        self.reset_icon_for_last_color_in_proteins_tab()
        self._view.ui.lbl_protein_current_color.setText("black    ")
        self._view.color_grid_proteins.last_clicked_color = "black"
        self.__slot_change_chain_color_proteins()
        self._view.color_grid_proteins.c_black.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_proteins.c_black.setIconSize(
            self._view.color_grid_proteins.c_black.icon().actualSize(QtCore.QSize(14, 14)))
    # </editor-fold>

    # <editor-fold desc="Reset Icon">
    def reset_icon_for_red_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("red")
        self._view.color_grid_proteins.c_red.setIcon(QtGui.QIcon())

    def reset_icon_for_tv_red_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("tv_red")
        self._view.color_grid_proteins.c_tv_red.setIcon(QtGui.QIcon())

    def reset_icon_for_salmon_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("salmon")
        self._view.color_grid_proteins.c_salomon.setIcon(QtGui.QIcon())

    def reset_icon_for_raspberry_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("raspberry")
        self._view.color_grid_proteins.c_raspberry.setIcon(QtGui.QIcon())

    def reset_icon_for_green_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("green")
        self._view.color_grid_proteins.c_green.setIcon(QtGui.QIcon())

    def reset_icon_for_tv_green_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("tv_green")
        self._view.color_grid_proteins.c_tv_green.setIcon(QtGui.QIcon())

    def reset_icon_for_palegreen_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("palegreen")
        self._view.color_grid_proteins.c_palegreen.setIcon(QtGui.QIcon())

    def reset_icon_for_forest_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("forest")
        self._view.color_grid_proteins.c_forest.setIcon(QtGui.QIcon())

    def reset_icon_for_blue_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("blue")
        self._view.color_grid_proteins.c_blue.setIcon(QtGui.QIcon())

    def reset_icon_for_tv_blue_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("tv_blue")
        self._view.color_grid_proteins.c_tv_blue.setIcon(QtGui.QIcon())

    def reset_icon_for_lightblue_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("lightblue")
        self._view.color_grid_proteins.c_lightblue.setIcon(QtGui.QIcon())

    def reset_icon_for_skyblue_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("skyblue")
        self._view.color_grid_proteins.c_skyblue.setIcon(QtGui.QIcon())

    def reset_icon_for_yellow_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("yellow")
        self._view.color_grid_proteins.c_yellow.setIcon(QtGui.QIcon())

    def reset_icon_for_tv_yellow_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("tv_yellow")
        self._view.color_grid_proteins.c_tv_yellow.setIcon(QtGui.QIcon())

    def reset_icon_for_paleyellow_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("paleyellow")
        self._view.color_grid_proteins.c_paleyellow.setIcon(QtGui.QIcon())

    def reset_icon_for_sand_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("sand")
        self._view.color_grid_proteins.c_sand.setIcon(QtGui.QIcon())

    def reset_icon_for_magenta_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("magenta")
        self._view.color_grid_proteins.c_magenta.setIcon(QtGui.QIcon())

    def reset_icon_for_purple_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("purple")
        self._view.color_grid_proteins.c_purple.setIcon(QtGui.QIcon())

    def reset_icon_for_pink_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("pink")
        self._view.color_grid_proteins.c_pink.setIcon(QtGui.QIcon())

    def reset_icon_for_hotpink_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("hotpink")
        self._view.color_grid_proteins.c_hotpink.setIcon(QtGui.QIcon())

    def reset_icon_for_cyan_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("cyan")
        self._view.color_grid_proteins.c_cyan.setIcon(QtGui.QIcon())

    def reset_icon_for_aquamarine_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("aquamarine")
        self._view.color_grid_proteins.c_aquamarine.setIcon(QtGui.QIcon())

    def reset_icon_for_palecyan_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("palecyan")
        self._view.color_grid_proteins.c_palecyan.setIcon(QtGui.QIcon())

    def reset_icon_for_teal_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("teal")
        self._view.color_grid_proteins.c_teal.setIcon(QtGui.QIcon())

    def reset_icon_for_orange_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("orange")
        self._view.color_grid_proteins.c_orange.setIcon(QtGui.QIcon())

    def reset_icon_for_tv_orange_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("tv_orange")
        self._view.color_grid_proteins.c_tv_orange.setIcon(QtGui.QIcon())

    def reset_icon_for_lightorange_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("lightorange")
        self._view.color_grid_proteins.c_lightorange.setIcon(QtGui.QIcon())

    def reset_icon_for_olive_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("olive")
        self._view.color_grid_proteins.c_olive.setIcon(QtGui.QIcon())

    def reset_icon_for_white_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("white")
        self._view.color_grid_proteins.c_white.setIcon(QtGui.QIcon())

    def reset_icon_for_grey_70_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("grey70")
        self._view.color_grid_proteins.c_grey_70.setIcon(QtGui.QIcon())

    def reset_icon_for_grey_30_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("grey30")
        self._view.color_grid_proteins.c_grey_30.setIcon(QtGui.QIcon())

    def reset_icon_for_black_in_proteins_tab(self):
        self._view.ui.lbl_protein_current_color.setText("black")
        self._view.color_grid_proteins.c_black.setIcon(QtGui.QIcon())
    # </editor-fold>
    # </editor-fold>

    # <editor-fold desc="Representations">
    # hydrogens
    def __slot_protein_chain_with_hydrogens(self):
        pass



    def __slot_show_protein_chain_with_hydrogens(self):
        try:
            tmp_protein = self._interface_manager.get_current_active_protein_object()
            tmp_chain = self._interface_manager.get_current_active_chain_object()
            if self._view.tg_protein_sticks.toggle_button.isChecked():
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    "sticks", f"h. and chain {tmp_chain.chain_letter} and {tmp_protein.get_molecule_object()}"
                )
            if self._view.tg_protein_lines.toggle_button.isChecked():
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    "lines", f"h. and chain {tmp_chain.chain_letter} and {tmp_protein.get_molecule_object()}"
                )
            if self._view.tg_protein_spheres.toggle_button.isChecked():
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    "spheres", f"h. and chain {tmp_chain.chain_letter} and {tmp_protein.get_molecule_object()}"
                )
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    def __slot_hide_protein_chain_with_hydrogens(self):
        try:
            tmp_protein = self._interface_manager.get_current_active_protein_object()
            tmp_chain = self._interface_manager.get_current_active_chain_object()
            if self._view.tg_protein_sticks.toggle_button.isChecked():
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    "sticks", f"h. and chain {tmp_chain.chain_letter} and {tmp_protein.get_molecule_object()}"
                )
            if self._view.tg_protein_lines.toggle_button.isChecked():
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    "lines", f"h. and chain {tmp_chain.chain_letter} and {tmp_protein.get_molecule_object()}"
                )
            if self._view.tg_protein_spheres.toggle_button.isChecked():
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    "spheres", f"h. and chain {tmp_chain.chain_letter} and {tmp_protein.get_molecule_object()}"
                )
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # cartoon
    def __slot_protein_chain_as_cartoon(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Cartoon' toggle on the 'Proteins Tab' was clicked.")
            tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
            tmp_selection.set_selection_for_a_single_chain(
                self._interface_manager.get_current_active_chain_object().chain_letter)
            if self._view.ui.cb_protein_cartoon.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.CARTOON.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.CARTOON.value, tmp_selection.selection_string
                )
            elif self._view.tg_protein_cartoon.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.CARTOON.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.CARTOON.value, tmp_selection.selection_string
                )
            else:
                #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.CARTOON.value)
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    enums.PyMOLRepresentation.CARTOON.value, tmp_selection.selection_string
                )
            self._update_scene()
            self._save_protein_pymol_session()
            self._interface_manager.manage_coloring_by_element_option_for_protein_chain()
            self._interface_manager.manage_hydrogen_representation_for_protein_chain()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # sticks
    def __slot_protein_chain_as_sticks(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Sticks' toggle on the 'Proteins Tab' was clicked.")
            tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
            tmp_selection.set_selection_for_a_single_chain(
                self._interface_manager.get_current_active_chain_object().chain_letter)
            if self._view.ui.cb_protein_sticks.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.STICKS.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.STICKS.value, tmp_selection.selection_string
                )
            elif self._view.tg_protein_sticks.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.STICKS.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.STICKS.value, tmp_selection.selection_string
                )
            else:
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    enums.PyMOLRepresentation.STICKS.value, tmp_selection.selection_string
                )
            self._update_scene()
            self._save_protein_pymol_session()
            self._interface_manager.manage_coloring_by_element_option_for_protein_chain()
            self._interface_manager.manage_hydrogen_representation_for_protein_chain()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # ribbon
    def __slot_protein_chain_as_ribbon(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Ribbon' toggle on the 'Proteins Tab' was clicked.")
            tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
            tmp_selection.set_selection_for_a_single_chain(
                self._interface_manager.get_current_active_chain_object().chain_letter)
            if self._view.ui.cb_protein_ribbon.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.RIBBON.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.RIBBON.value, tmp_selection.selection_string
                )
            elif self._view.tg_protein_ribbon.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.RIBBON.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.RIBBON.value, tmp_selection.selection_string
                )
            else:
                #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.RIBBON.value)
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    enums.PyMOLRepresentation.RIBBON.value, tmp_selection.selection_string
                )
            self._update_scene()
            self._save_protein_pymol_session()
            self._interface_manager.manage_coloring_by_element_option_for_protein_chain()
            self._interface_manager.manage_hydrogen_representation_for_protein_chain()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # lines
    def __slot_protein_chain_as_lines(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Lines' toggle on the 'Proteins Tab' was clicked.")
            tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
            tmp_selection.set_selection_for_a_single_chain(
                self._interface_manager.get_current_active_chain_object().chain_letter)
            if self._view.ui.cb_protein_lines.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.LINES.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.LINES.value, tmp_selection.selection_string
                )
            elif self._view.tg_protein_lines.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.LINES.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.LINES.value, tmp_selection.selection_string
                )
            else:
                #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.LINES.value)
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    enums.PyMOLRepresentation.LINES.value, tmp_selection.selection_string
                )
            self._update_scene()
            self._save_protein_pymol_session()
            self._interface_manager.manage_coloring_by_element_option_for_protein_chain()
            self._interface_manager.manage_hydrogen_representation_for_protein_chain()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # spheres
    def __slot_protein_chain_as_spheres(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Spheres' toggle on the 'Proteins Tab' was clicked.")
            tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
            tmp_selection.set_selection_for_a_single_chain(
                self._interface_manager.get_current_active_chain_object().chain_letter)
            if self._view.ui.cb_protein_spheres.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SPHERES.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.SPHERES.value, tmp_selection.selection_string
                )
            elif self._view.tg_protein_spheres.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SPHERES.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.SPHERES.value, tmp_selection.selection_string
                )
            else:
                #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.SPHERES.value)
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    enums.PyMOLRepresentation.SPHERES.value, tmp_selection.selection_string
                )
            self._update_scene()
            self._save_protein_pymol_session()
            self._interface_manager.manage_coloring_by_element_option_for_protein_chain()
            self._interface_manager.manage_hydrogen_representation_for_protein_chain()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # dots
    def __slot_protein_chain_as_dots(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Dots' toggle on the 'Proteins Tab' was clicked.")
            tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
            tmp_selection.set_selection_for_a_single_chain(
                self._interface_manager.get_current_active_chain_object().chain_letter)
            if self._view.ui.cb_protein_dots.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.DOTS.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.DOTS.value, tmp_selection.selection_string
                )
            elif self._view.tg_protein_dots.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.DOTS.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.DOTS.value, tmp_selection.selection_string
                )
            else:
                #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.DOTS.value)
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    enums.PyMOLRepresentation.DOTS.value, tmp_selection.selection_string
                )
            self._update_scene()
            self._save_protein_pymol_session()
            self._interface_manager.manage_coloring_by_element_option_for_protein_chain()
            self._interface_manager.manage_hydrogen_representation_for_protein_chain()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # mesh
    def __slot_protein_chain_as_mesh(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Mesh' toggle on the 'Proteins Tab' was clicked.")
            tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
            tmp_selection.set_selection_for_a_single_chain(
                self._interface_manager.get_current_active_chain_object().chain_letter)
            if self._view.ui.cb_protein_mesh.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.MESH.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.MESH.value, tmp_selection.selection_string
                )
            elif self._view.tg_protein_mesh.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.MESH.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.MESH.value, tmp_selection.selection_string
                )
            else:
                #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.MESH.value)
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    enums.PyMOLRepresentation.MESH.value, tmp_selection.selection_string
                )
            self._update_scene()
            self._save_protein_pymol_session()
            self._interface_manager.manage_coloring_by_element_option_for_protein_chain()
            self._interface_manager.manage_hydrogen_representation_for_protein_chain()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # surface
    def __slot_protein_chain_as_surface(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Surface' toggle on the 'Proteins Tab' was clicked.")
            tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
            tmp_selection.set_selection_for_a_single_chain(
                self._interface_manager.get_current_active_chain_object().chain_letter)
            if self._view.ui.cb_protein_surface.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SURFACE.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.SURFACE.value, tmp_selection.selection_string
                )
            elif self._view.tg_protein_surface.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SURFACE.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.SURFACE.value, tmp_selection.selection_string
                )
            else:
                #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.SURFACE.value)
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    enums.PyMOLRepresentation.SURFACE.value, tmp_selection.selection_string
                )
            self._update_scene()
            self._save_protein_pymol_session()
            self._interface_manager.manage_coloring_by_element_option_for_protein_chain()
            self._interface_manager.manage_hydrogen_representation_for_protein_chain()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # all
    def __slot_hide_protein_chain_all(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Hide all' representations button on the 'Proteins Tab' was clicked.")
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
            self.__slot_protein_chain_as_cartoon()
            self.__slot_protein_chain_as_sticks()
            self.__slot_protein_chain_as_ribbon()
            self.__slot_protein_chain_as_lines()
            self.__slot_protein_chain_as_spheres()
            self.__slot_protein_chain_as_dots()
            self.__slot_protein_chain_as_mesh()
            self.__slot_protein_chain_as_surface()
            self._update_scene()
            self._save_protein_pymol_session()
            self._interface_manager.manage_coloring_by_element_option_for_protein_chain()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

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
    #     if self._interface_manager.pymol_session_manager.session_object_type == "protein" and self._interface_manager.pymol_session_manager.session_name == tmp_protein.get_molecule_object():
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

    def __slot_import_protein_structure(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Import protein' button on the 'Proteins Tab' was clicked.")
            self._external_controller = add_protein_view_controller.AddProteinViewController(self._interface_manager)
            self._external_controller.user_input.connect(self._post_import_protein_structure)
            self._external_controller.restore_ui()
            self._interface_manager.get_add_protein_view().show()
        except Exception as e:
            logger.error(f"An error occurred during the protein import: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    def _post_import_protein_structure(self, return_value: tuple):
        try:
            tmp_protein_name, tmp_name_len = return_value
            if tmp_name_len == 4:
                self._active_task = tasks.Task(
                    target=project_async.add_protein_from_pdb_to_project,
                    args=(
                        tmp_protein_name,
                        self._interface_manager
                    ),
                    post_func=self.__await_post_import_protein_structure,
                )
                self._active_task.start()
                constants.PYSSA_LOGGER.info("Create project finished with protein from the PDB.")
            elif tmp_name_len > 0:
                # local pdb file as input
                self._active_task = tasks.Task(
                    target=project_async.add_protein_from_local_filesystem_to_project,
                    args=(
                        tmp_protein_name,
                        self._interface_manager
                    ),
                    post_func=self.__await_post_import_protein_structure,
                )
                self._active_task.start()
                constants.PYSSA_LOGGER.info("Create project finished with protein from local filesystem.")
            else:
                logger.warning("No protein object was created.")
                return
            self._interface_manager.status_bar_manager.show_temporary_message(
                "Importing protein structure ...", False
            )
        except Exception as e:
            logger.error(f"An error occurred during the protein import: {e}")
            self._interface_manager.status_bar_manager.show_error_message("Protein import failed!")
        else:
            self._interface_manager.start_wait_cursor()

    def __await_post_import_protein_structure(self, return_value: tuple):
        try:
            tmp_protein: "protein.Protein" = return_value[1]
            self._interface_manager.get_current_project().add_existing_protein(tmp_protein)
            self._interface_manager.watcher.add_protein(tmp_protein.get_molecule_object())
            self._database_thread.put_database_operation_into_queue(
                database_operation.DatabaseOperation(enums.SQLQueryType.INSERT_NEW_PROTEIN,
                                                     (0, tmp_protein)))
            self._interface_manager.refresh_main_view()
            #self._interface_manager.pymol_session_manager.reinitialize_session()
            #self._interface_manager.pymol_session_manager.unfreeze_current_protein_pymol_session()
            #self._interface_manager.pymol_session_manager.unfreeze_current_protein_pair_pymol_session()
            #self._main_view_state.restore_main_view_state()
            self._interface_manager.status_bar_manager.show_temporary_message("Importing protein structure finished.")
        except Exception as e:
            logger.error(f"An error occurred during the protein import: {e}")
            self._interface_manager.status_bar_manager.show_error_message("Protein import failed!")
        finally:
            self._interface_manager.stop_wait_cursor()

    def __slot_delete_protein(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Delete protein' button on the 'Proteins Tab' was clicked.")
            tmp_dialog = custom_message_box.CustomMessageBoxDelete(
                "Are you sure you want to delete this protein?", "Delete Protein",
                custom_message_box.CustomMessageBoxIcons.WARNING.value
            )
            tmp_dialog.exec_()
            response: bool = tmp_dialog.response
            if response:
                tmp_protein: "protein.Protein" = self._interface_manager.get_current_active_protein_object()
                if self._interface_manager.pymol_session_manager.is_the_current_protein_in_session(tmp_protein.get_molecule_object()):
                    self._interface_manager.pymol_session_manager.reinitialize_session()
                tmp_database_operation = database_operation.DatabaseOperation(enums.SQLQueryType.DELETE_EXISTING_PROTEIN,
                                                                              (0, tmp_protein.get_id()))
                self._database_thread.put_database_operation_into_queue(tmp_database_operation)
                self._interface_manager.get_current_project().delete_specific_protein(tmp_protein.get_molecule_object())
                self._interface_manager.watcher.remove_protein(tmp_protein.get_molecule_object())
                self._interface_manager.remove_protein_from_proteins_model()
                self._interface_manager.refresh_main_view()
        except Exception as e:
            logger.error(f"An error occurred during the protein deletion process: {e}")
            self._interface_manager.status_bar_manager.show_error_message("Protein delete failed!")

    def __slot_save_selected_protein_structure_as_pdb_file(self) -> None:
        """Saves selected protein as pdb file."""
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Export protein' button on the 'Proteins Tab' was clicked.")
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
                    target=protein_async.save_selected_protein_structure_as_pdb_file,
                    args=(
                        tmp_protein,
                        file_path,
                        self._interface_manager.get_current_project().get_database_filepath()
                    ),
                    post_func=self.__await_save_selected_protein_structure_as_pdb_file,
                )
            else:
                self._interface_manager.stop_wait_cursor()
                self._interface_manager.refresh_main_view()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
        else:
            self._active_task.start()
            self._interface_manager.start_wait_cursor()

    def __await_save_selected_protein_structure_as_pdb_file(self, result: tuple) -> None:
        try:
            if result[0] == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
                tmp_dialog = custom_message_box.CustomMessageBoxOk(
                    "Saving the protein as .pdb file failed!",
                    "Save Protein Structure",
                    custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
                )
                tmp_dialog.exec_()
            elif result[0] == exit_codes.EXIT_CODE_ZERO[0]:
                self._interface_manager.status_bar_manager.show_temporary_message("The protein was successfully saved as .pdb file.")
            else:
                tmp_dialog = custom_message_box.CustomMessageBoxOk(
                    "Saving the protein as .pdb file failed with an unexpected error!",
                    "Save Protein Structure",
                    custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
                )
                tmp_dialog.exec_()
            self._interface_manager.refresh_protein_model()
            self._interface_manager.refresh_main_view()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
        finally:
            self._interface_manager.stop_wait_cursor()

    def __slot_clean_protein_update(self) -> None:
        """Cleans the selected protein structure."""
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Clean protein' context menu action was clicked.")
            tmp_dialog = custom_message_box.CustomMessageBoxYesNo(
                "Are you sure you want to clean this protein?\n" "This will remove all organic and solvent components!", "Clean Protein",
                custom_message_box.CustomMessageBoxIcons.WARNING.value
            )
            tmp_dialog.exec_()
            if tmp_dialog.response:
                tmp_main_socket, tmp_general_purpose_socket = self._interface_manager.job_manager.get_general_purpose_socket_pair()
                self._active_task = tasks.Task(
                    target=protein_async.clean_protein_update,
                    args=(
                        self._interface_manager.get_current_active_protein_object(),
                        self._interface_manager.get_current_project().get_database_filepath(),
                        tmp_main_socket,
                        tmp_general_purpose_socket
                    ),
                    post_func=self.__await_clean_protein_update,
                )
                self._active_task.start()
                self._interface_manager.start_wait_cursor()
                self.update_status("Cleaning protein ...")
            else:
                constants.PYSSA_LOGGER.info("No protein has been cleaned.")
                self._interface_manager.stop_wait_cursor()
                self._interface_manager.refresh_main_view()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
            self._interface_manager.stop_wait_cursor()
            self._interface_manager.refresh_main_view()


    def __await_clean_protein_update(self, return_value: tuple) -> None:
        try:
            if return_value[1] is None:
                self._interface_manager.status_bar_manager.show_error_message("Cleaning protein failed!")
            else:
                self.update_status("Cleaning protein finished.")
                if self._interface_manager.pymol_session_manager.is_the_current_protein_in_session(
                    self._interface_manager.get_current_active_protein_object().get_molecule_object()
                ):
                    self._interface_manager.pymol_session_manager.load_protein_session(
                        self._interface_manager.get_current_active_protein_object()
                    )
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
        finally:
            self._interface_manager.stop_wait_cursor()
            self._interface_manager.refresh_main_view()

    def __slot_rename_selected_protein_structure(self) -> None:
        """Opens a new view to rename the selected protein."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Rename protein' context menu action was clicked.")
        self._external_controller = rename_protein_view_controller.RenameProteinViewController(self._interface_manager)
        self._external_controller.user_input.connect(self.post_rename_selected_protein_structure)
        self._external_controller.restore_ui()
        self._interface_manager.get_rename_protein_view().show()

    def post_rename_selected_protein_structure(self, return_value: tuple) -> None:
        """Renames a selected protein structure."""
        if return_value[1] is True:
            self._active_task = tasks.Task(
                target=protein_async.rename_selected_protein_structure,
                args=(
                    self._interface_manager.get_current_protein_tree_index_object(),
                    return_value[0],
                    self._interface_manager.get_current_project().get_database_filepath()
                ),
                post_func=self.__await_post_rename_selected_protein_structure,
            )
            self._active_task.start()
            self.update_status("Renaming protein ...")
            self._interface_manager.start_wait_cursor()
        else:
            pass

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
        self._interface_manager.stop_wait_cursor()
        self.update_status("Renaming protein finished.")

    def __slot_show_protein_chain_sequence(self) -> None:
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Show protein sequence' context menu action was clicked.")
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
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    def __slot_update_protein_scene(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Update protein scene' button on the 'Proteins Tab' was clicked.")
            self._interface_manager.pymol_session_manager.pymol_interface.scene(a_key="auto", an_action="update")
            self._save_protein_pymol_session()
            self._interface_manager.status_bar_manager.show_temporary_message("PyMOL Scene updated.", a_timeout=1500)
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    def _update_scene(self) -> None:
        """Updates the current selected PyMOL scene."""
        self._interface_manager.pymol_session_manager.pymol_interface.scene(
            a_key="auto", an_action="update"
        )
        self._interface_manager.status_bar_manager.show_temporary_message("PyMOL Scene updated.", a_timeout=1500)

    def __slot_save_scene(self) -> None:
        """Saves the current view as a new PyMOL scene."""
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Create pymol scene' button on the 'Proteins or Protein Pairs Tab' was clicked.")
            self._external_controller = add_scene_view_controller.AddSceneViewController(self._interface_manager)
            self._external_controller.user_input.connect(self.post_save_scene)
            self._external_controller.restore_ui()
            self._interface_manager.get_add_scene_view().show()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    def post_save_scene(self, return_value: tuple):
        tmp_scene_name, _ = return_value
        self._interface_manager.pymol_session_manager.pymol_interface.scene(
            a_key=tmp_scene_name, an_action="append"
        )

        if self._interface_manager.current_tab_index == 1:
            self._active_task = tasks.Task(
                target=pymol_session_async.save_protein_pymol_session_to_database,
                args=(
                    self._interface_manager,
                    0
                ),
                post_func=self.__await_save_scene_protein,
            )
            self._active_task.start()
            self._interface_manager.status_bar_manager.show_temporary_message(
                "Adding new scene to protein ...", False)
            self._interface_manager.start_wait_cursor()
            self._interface_manager.add_scene_to_proteins_model(tmp_scene_name)
            self._interface_manager.pymol_session_manager.current_scene_name = tmp_scene_name
            ui_util.set_pymol_scene_name_into_label(self._interface_manager.pymol_session_manager.current_scene_name,
                                                    self._view.ui.lbl_pymol_protein_scene)

        elif self._interface_manager.current_tab_index == 2:
            # The database thread cannot be used here because the session gets loaded again
            # before the new data is in the db
            self._active_task = tasks.Task(
                target=pymol_session_async.save_protein_pair_pymol_session_to_database,
                args=(
                    self._interface_manager,
                    0
                ),
                post_func=self.__await_save_scene_protein_pair,
            )
            self._active_task.start()
            self._interface_manager.status_bar_manager.show_temporary_message(
                "Adding new scene to protein pair ...", False)
            self._interface_manager.start_wait_cursor()
            self._interface_manager.add_scene_to_protein_pairs_model(tmp_scene_name)
            self._interface_manager.pymol_session_manager.current_scene_name = tmp_scene_name
            ui_util.set_pymol_scene_name_into_label(self._interface_manager.pymol_session_manager.current_scene_name,
                                                    self._view.ui.lbl_pymol_protein_pair_scene)
        else:
            logger.warning("The current tab index is not for the proteins nor for the protein pairs tab?!")
            return

    def __await_save_scene_protein(self, return_value: tuple):
        _, exit_flag = return_value
        self._interface_manager.refresh_main_view()
        if exit_flag:
            self._interface_manager.status_bar_manager.show_temporary_message("Adding new scene to protein finished.")
        else:
            self._interface_manager.status_bar_manager.show_error_message("Adding new scene to protein failed!")
        self._interface_manager.stop_wait_cursor()

    def __await_save_scene_protein_pair(self, return_value: tuple):
        _, exit_flag = return_value
        self._interface_manager.refresh_main_view()
        if exit_flag:
            self._interface_manager.status_bar_manager.show_temporary_message("Adding new scene to protein pair finished.")
        else:
            self._interface_manager.status_bar_manager.show_error_message("Adding new scene to protein pair failed!")
        self._interface_manager.stop_wait_cursor()

    def __slot_delete_current_scene(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Delete pymol scene' button on the 'Proteins or Protein Pairs Tab' was clicked.")
            tmp_dialog = custom_message_box.CustomMessageBoxDelete(
                "Are you sure you want to delete this scene?", "Delete PyMOL Scene",
                custom_message_box.CustomMessageBoxIcons.WARNING.value
            )
            tmp_dialog.exec_()
            response: bool = tmp_dialog.response
            if response:
                self._interface_manager.pymol_session_manager.pymol_interface.scene(
                    a_key=self._interface_manager.pymol_session_manager.current_scene_name, an_action="clear"
                )
                if self._interface_manager.current_tab_index == 1:
                    self._save_protein_pymol_session()
                    self._interface_manager.remove_scene_from_proteins_model(self._interface_manager.get_current_protein_tree_index())
                    self._interface_manager.refresh_main_view()
                    self._view.ui.lbl_pymol_protein_scene.setText("PyMOL Scene: No Scene Selected")
                elif self._interface_manager.current_tab_index == 2:
                    # The database thread cannot be used here because the session gets loaded again
                    # before the new data is in the db
                    self._active_task = tasks.Task(
                        target=pymol_session_async.save_protein_pair_pymol_session_to_database,
                        args=(
                            self._interface_manager,
                            0
                        ),
                        post_func=self.__await_delete_current_scene,
                    )
                    self._active_task.start()
                    self._interface_manager.status_bar_manager.show_temporary_message(
                        "Deleting selected scene ...", False
                    )
                    self._interface_manager.start_wait_cursor()
                    self._interface_manager.remove_scene_from_protein_pairs_model(
                        self._interface_manager.get_current_protein_pair_tree_index()
                    )
                    self._view.ui.lbl_pymol_protein_pair_scene.setText("PyMOL Scene: No Scene Selected")
                else:
                    logger.warning("The current tab index is not for the proteins nor for the protein pairs tab?!")
                    return
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    def __await_delete_current_scene(self, return_value: tuple):
        try:
            _, exit_flag = return_value
            self._interface_manager.refresh_main_view()
            if exit_flag:
                self._interface_manager.status_bar_manager.show_temporary_message("Deleted the scene successfully.")
            else:
                self._interface_manager.status_bar_manager.show_error_message("Deleting the scene failed!")
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
        finally:
            self._interface_manager.stop_wait_cursor()

    def _save_protein_pymol_session(self):
        """Saves the session as base64 string and updates the database"""
        tmp_protein = self._interface_manager.get_current_active_protein_object()
        tmp_protein.pymol_session = self._interface_manager.pymol_session_manager.save_current_session_as_base64()
        tmp_database_operation = database_operation.DatabaseOperation(
            enums.SQLQueryType.UPDATE_PYMOL_SESSION_PROTEIN,
            (0, self._interface_manager.get_current_active_protein_object())
        )
        self._database_thread.put_database_operation_into_queue(tmp_database_operation)

    # </editor-fold>

    # <editor-fold desc="Protein Pairs tab methods">
    def open_context_menu_for_protein_pairs(self, position):
        tmp_context_menu = self._protein_pair_tree_context_menu.get_context_menu(
            self._view.ui.protein_pairs_tree_view.selectedIndexes(),
            self._interface_manager.pymol_session_manager.is_the_current_protein_pair_in_session(self._interface_manager.get_current_active_protein_pair_object().name)
        )
        tmp_context_menu.exec_(self._view.ui.protein_pairs_tree_view.viewport().mapToGlobal(position))

    def __slot_open_protein_pair_pymol_session(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Open protein pair pymol session' button on the 'Protein Pairs Tab' was clicked.")
            self._view.tg_protein_pair_white_bg.toggle_button.setCheckState(False)
            tmp_protein_pair: "protein_pair.ProteinPair" = self._interface_manager.get_current_active_protein_pair_object()
            print(tmp_protein_pair.protein_1.get_molecule_object())
            print(tmp_protein_pair.protein_1)
            print(tmp_protein_pair.protein_2.get_molecule_object())
            print(tmp_protein_pair.protein_2)
            # fixme: I am no sure if the code below is needed
            # if not self._interface_manager.pymol_session_manager.is_the_current_session_empty():
            #     tmp_flag = True  # Session is NOT empty and needs reinitialization
            # else:
            #     tmp_flag = False  # Session is empty

            tmp_flag = False
            self._active_task = tasks.Task(
                target=pymol_session_async.load_protein_pair_pymol_session,
                args=(
                    tmp_protein_pair,
                    self._interface_manager.pymol_session_manager,
                    tmp_flag
                ),
                post_func=self.__await_open_protein_pair_pymol_session,
            )
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
        else:
            self._active_task.start()
            self._interface_manager.start_wait_cursor()
            self._interface_manager.status_bar_manager.show_temporary_message(
                f"Loading PyMOL session of {tmp_protein_pair.name} ...", False
            )

    def __await_open_protein_pair_pymol_session(self, return_value: tuple):
        try:
            _, exit_boolean = return_value
            self._view.ui.action_protein_regions.setEnabled(False)
            if exit_boolean:
                self._view.ui.action_protein_regions.setEnabled(True)
                self._view.ui.btn_create_protein_pair_scene.setEnabled(True)
                self._view.ui.btn_update_protein_pair_scene.setEnabled(True)
                self._view.ui.lbl_session_name.setText(f"Session Name: {self._interface_manager.pymol_session_manager.session_name}")
                self._interface_manager.pymol_session_manager.current_scene_name = self._view.ui.protein_pairs_tree_view.currentIndex().child(0, 0).child(1, 0).data(Qt.DisplayRole)
                self._interface_manager.pymol_session_manager.load_current_scene()
                ui_util.set_pymol_scene_name_into_label(self._interface_manager.pymol_session_manager.current_scene_name,
                                                        self._view.ui.lbl_pymol_protein_pair_scene)
                logger.info("Successfully opened protein pair session.")
                self._interface_manager.status_bar_manager.show_temporary_message("Loading the PyMOL session was successful.")
                self._view.ui.lbl_info_3.setText("Please select a chain.")
            else:
                logger.error("The protein name could not be found in the object list in PyMOL!")
                self._view.ui.btn_create_protein_pair_scene.setEnabled(False)
                self._view.ui.btn_update_protein_pair_scene.setEnabled(False)
                self._interface_manager.status_bar_manager.show_error_message("Loading the PyMOL session failed!")
                self._view.ui.lbl_info_3.setText("Please load the PyMOL session of the selected protein.")
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
        finally:
            self._interface_manager.refresh_main_view()
            self._interface_manager.stop_wait_cursor()

    def _save_protein_pair_pymol_session(self):
        tmp_protein_pair = self._interface_manager.get_current_active_protein_pair_object()
        tmp_protein_pair.pymol_session = self._interface_manager.pymol_session_manager.save_current_session_as_base64()
        tmp_database_operation = database_operation.DatabaseOperation(
            enums.SQLQueryType.UPDATE_PYMOL_SESSION_PROTEIN_PAIR,
            (0, tmp_protein_pair.get_id(), tmp_protein_pair)
        )
        self._database_thread.put_database_operation_into_queue(tmp_database_operation)

    def __slot_delete_protein_pair_from_project(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Delete protein pair' button on the 'Protein Pairs Tab' was clicked.")
            tmp_dialog = custom_message_box.CustomMessageBoxDelete(
                "Are you sure you want to delete this protein pair?","Delete Protein Pair",
                custom_message_box.CustomMessageBoxIcons.WARNING.value
            )
            tmp_dialog.exec_()
            response: bool = tmp_dialog.response
            if response:
                tmp_protein_pair: "protein_pair.ProteinPair" = self._interface_manager.get_current_active_protein_pair_object()
                if self._interface_manager.pymol_session_manager.is_the_current_protein_pair_in_session(tmp_protein_pair.name):
                    self._interface_manager.pymol_session_manager.reinitialize_session()
                tmp_database_operation = database_operation.DatabaseOperation(
                    enums.SQLQueryType.DELETE_EXISTING_PROTEIN_PAIR, (0, tmp_protein_pair.get_id())
                )
                self._database_thread.put_database_operation_into_queue(tmp_database_operation)
                self._interface_manager.get_current_project().delete_specific_protein_pair(tmp_protein_pair.name)
                self._interface_manager.watcher.remove_protein_pair(tmp_protein_pair.name)
                self._interface_manager.remove_protein_pair_from_protein_pairs_model()
                self._interface_manager.refresh_main_view()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    def __slot_get_information_about_selected_object_in_protein_pair_branch(self):
        try:
            tmp_type = self._interface_manager.get_current_protein_pair_tree_index_type()

            if tmp_type == "protein_pair":
                logger.log(
                    log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                    f"The protein pair object '{self._view.ui.protein_pairs_tree_view.currentIndex().data(Qt.DisplayRole)}' on the 'Protein Pairs Tab' was clicked."
                )
            elif tmp_type == "protein":
                logger.log(
                    log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                    f"The protein object '{self._view.ui.protein_pairs_tree_view.currentIndex().data(Qt.DisplayRole)}' on the 'Protein Pairs Tab' was clicked."
                )
            elif tmp_type == "scene":
                tmp_scene_name = self._view.ui.protein_pairs_tree_view.currentIndex().data(Qt.DisplayRole)
                tmp_protein_pair_name = self._view.ui.protein_pairs_tree_view.currentIndex().parent().parent().data(Qt.DisplayRole)
                logger.log(
                    log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                    f"The scene '{tmp_scene_name}' of the protein pair '{tmp_protein_pair_name}' on the 'Protein Pairs Tab' was clicked."
                )
                if self._interface_manager.pymol_session_manager.is_the_current_protein_pair_in_session(self._interface_manager.get_current_active_protein_pair_object().name):
                    tmp_scene_name = self._interface_manager.get_current_active_scene_name_of_protein_pair()
                    self._interface_manager.pymol_session_manager.current_scene_name = tmp_scene_name
                    ui_util.set_pymol_scene_name_into_label(self._interface_manager.pymol_session_manager.current_scene_name,
                                                            self._view.ui.lbl_pymol_protein_pair_scene)
                    self._interface_manager.pymol_session_manager.load_scene(tmp_scene_name)

            elif tmp_type == "chain":
                tmp_chain_letter = self._view.ui.protein_pairs_tree_view.currentIndex().data(Qt.DisplayRole)
                tmp_protein_name = self._view.ui.protein_pairs_tree_view.currentIndex().parent().parent().data(Qt.DisplayRole)
                tmp_protein_pair_name = self._view.ui.protein_pairs_tree_view.currentIndex().parent().parent().parent().data(Qt.DisplayRole)
                logger.log(
                    log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                    f"The chain object '{tmp_chain_letter}' of the protein '{tmp_protein_name}' of the protein pair '{tmp_protein_pair_name}' on the 'Protein Pairs Tab' was clicked."
                )
                if self._interface_manager.pymol_session_manager.current_scene_name != "" and self._interface_manager.pymol_session_manager.is_the_current_protein_pair_in_session(self._interface_manager.get_current_active_protein_pair_object().name):
                    self.set_icon_for_current_color_in_protein_pairs_tab()
                    self._interface_manager.show_chain_pymol_parameter_for_protein_pairs(self._interface_manager.pymol_session_manager)
                    self._interface_manager.set_repr_state_in_ui_for_protein_pair_chain(self._interface_manager.pymol_session_manager)
                    self._main_view_state.selected_chain_protein_pairs = self._interface_manager.get_current_protein_pair_tree_index()

            elif tmp_type == "header":
                logger.log(
                    log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                    f"The header '{self._view.ui.protein_pairs_tree_view.currentIndex().data(Qt.DisplayRole)}' on the 'Protein Pairs Tab' was clicked."
                )
            else:
                logger.warning("Unknown object type occurred in Protein Pairs tab.")
                return
            self._interface_manager.manage_ui_of_protein_pairs_tab(tmp_type, self._interface_manager.pymol_session_manager)
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    def __slot_color_protein_pair_by_rmsd(self) -> None:
        """Colors the residues in 5 colors depending on their distance to the reference."""
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Color protein pair by rmsd' context menu action was clicked.")
            self._active_task = tasks.Task(
                target=protein_pair_async.color_protein_pair_by_rmsd_value,
                args=(
                    self._interface_manager.get_current_protein_pair_tree_index_object(),
                    self._interface_manager.pymol_session_manager
                ),
                post_func=self.__await_color_protein_pair_by_rmsd,
            )
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
        else:
            self._active_task.start()
            self._interface_manager.start_wait_cursor()

    def __await_color_protein_pair_by_rmsd(self, result: tuple) -> None:
        self._interface_manager.stop_wait_cursor()
        self._interface_manager.refresh_main_view()

    def __slot_change_chain_color_protein_pairs(self) -> None:
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "The 'Color' attribute of a protein chain on the 'Protein Pairs Tab' changed.")
            tmp_protein_pair = self._interface_manager.get_current_active_protein_pair_object()
            tmp_protein = self._interface_manager.get_current_active_protein_object_of_protein_pair()
            tmp_chain = self._interface_manager.get_current_active_chain_object_of_protein_pair()
            if self._interface_manager._settings_manager.settings.protein_pairs_tab_use_combobox_for_colors == 1:
                tmp_color = self._view.ui.box_protein_pair_color.currentText()
            else:
                tmp_color = self._view.ui.lbl_protein_pair_current_color.text().strip()
            if self._interface_manager.pymol_session_manager.session_object_type == "protein_pair" and self._interface_manager.pymol_session_manager.session_name == tmp_protein_pair.name:
                # Update pymol parameter in PyMOL
                tmp_protein.pymol_selection.set_selection_for_a_single_chain(tmp_chain.chain_letter)
                self._interface_manager.pymol_session_manager.color_protein(tmp_color,
                                                                            tmp_protein.pymol_selection.selection_string)
                self._update_scene()
                self._save_protein_pair_pymol_session()
            else:
                logger.warning("The color of a protein chain could not be changed. This can be due to UI setup reasons.")
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    def __slot_change_chain_color_protein_pairs_atoms(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Color atoms by element' button on the 'Protein Pairs Tab' was clicked.")
            tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
            tmp_selection.set_selection_for_a_single_chain(
                self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
            self._interface_manager.pymol_session_manager.color_protein(
                'atomic', f"{tmp_selection.selection_string} and not elem C"
            )
            self._interface_manager.pymol_session_manager.color_protein(
                'grey70', f"{tmp_selection.selection_string} and elem C"
            )
            self.reset_icon_for_last_color_in_protein_pairs_tab()
            self._view.ui.lbl_protein_pair_current_color.setText("grey70    ")
            self._view.color_grid_protein_pairs.c_grey_70.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
            self._view.color_grid_protein_pairs.c_grey_70.setIconSize(
                self._view.color_grid_protein_pairs.c_grey_70.icon().actualSize(QtCore.QSize(14, 14)))
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    def __slot_change_chain_reset_protein_pairs_atoms(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Reset color atoms by element' button on the 'Protein Pairs Tab' was clicked.")
            tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
            tmp_selection.set_selection_for_a_single_chain(
                self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
            self._interface_manager.pymol_session_manager.color_protein(
                self._view.color_grid_protein_pairs.last_clicked_color, f"{tmp_selection.selection_string}"
            )
            self.set_icon_for_current_color_in_protein_pairs_tab()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    def __slot_protein_pair_change_background_color(self):
        try:
            if self._view.tg_protein_pair_white_bg.toggle_button.isChecked():
                self._interface_manager.pymol_session_manager.pymol_interface.set_background_color("white")
            else:
                self._interface_manager.pymol_session_manager.pymol_interface.set_background_color("black")
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # <editor-fold desc="Color Grid slot methods">
    def set_icon_for_current_color_in_protein_pairs_tab(self):
        try:
            color_index_functions = {
                "red": self.set_color_name_in_label_red_in_protein_pairs_tab,
                "tv_red": self.set_color_name_in_label_tv_red_in_protein_pairs_tab,
                "salmon": self.set_color_name_in_label_salmon_in_protein_pairs_tab,
                "raspberry": self.set_color_name_in_label_raspberry_in_protein_pairs_tab,
                "green": self.set_color_name_in_label_green_in_protein_pairs_tab,
                "tv_green": self.set_color_name_in_label_tv_green_in_protein_pairs_tab,
                "palegreen": self.set_color_name_in_label_palegreen_in_protein_pairs_tab,
                "forest": self.set_color_name_in_label_forest_in_protein_pairs_tab,
                "blue": self.set_color_name_in_label_blue_in_protein_pairs_tab,
                "tv_blue": self.set_color_name_in_label_tv_blue_in_protein_pairs_tab,
                "lightblue": self.set_color_name_in_label_lightblue_in_protein_pairs_tab,
                "skyblue": self.set_color_name_in_label_skyblue_in_protein_pairs_tab,
                "yellow": self.set_color_name_in_label_yellow_in_protein_pairs_tab,
                "tv_yellow": self.set_color_name_in_label_tv_yellow_in_protein_pairs_tab,
                "paleyellow": self.set_color_name_in_label_paleyellow_in_protein_pairs_tab,
                "sand": self.set_color_name_in_label_sand_in_protein_pairs_tab,
                "magenta": self.set_color_name_in_label_magenta_in_protein_pairs_tab,
                "purple": self.set_color_name_in_label_purple_in_protein_pairs_tab,
                "pink": self.set_color_name_in_label_pink_in_protein_pairs_tab,
                "hotpink": self.set_color_name_in_label_hotpink_in_protein_pairs_tab,
                "cyan": self.set_color_name_in_label_cyan_in_protein_pairs_tab,
                "aquamarine": self.set_color_name_in_label_aquamarine_in_protein_pairs_tab,
                "palecyan": self.set_color_name_in_label_palecyan_in_protein_pairs_tab,
                "teal": self.set_color_name_in_label_teal_in_protein_pairs_tab,
                "orange": self.set_color_name_in_label_orange_in_protein_pairs_tab,
                "tv_orange": self.set_color_name_in_label_tv_orange_in_protein_pairs_tab,
                "lightorange": self.set_color_name_in_label_lightorange_in_protein_pairs_tab,
                "olive": self.set_color_name_in_label_olive_in_protein_pairs_tab,
                "white": self.set_color_name_in_label_white_in_protein_pairs_tab,
                "grey70": self.set_color_name_in_label_grey_70_in_protein_pairs_tab,
                "grey30": self.set_color_name_in_label_grey_30_in_protein_pairs_tab,
                "black": self.set_color_name_in_label_black_in_protein_pairs_tab
            }
            tmp_protein = self._interface_manager.get_current_active_protein_object_of_protein_pair()
            tmp_chain = self._interface_manager.get_current_active_chain_object_of_protein_pair()
            tmp_protein.pymol_selection.set_selection_for_a_single_chain(tmp_chain.chain_letter)
            color_index_functions[
                tmp_chain.get_color(tmp_protein.pymol_selection.selection_string,
                                    self._interface_manager.pymol_session_manager)[0]
            ]()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    def reset_icon_for_last_color_in_protein_pairs_tab(self):
        try:
            tmp_color_name = self._view.ui.lbl_protein_pair_current_color.text()
            if tmp_color_name == "":
                return
            color_index_functions = {
                "red": self.reset_icon_for_red_in_protein_pairs_tab,
                "tv_red": self.reset_icon_for_tv_red_in_protein_pairs_tab,
                "salmon": self.reset_icon_for_salmon_in_protein_pairs_tab,
                "raspberry": self.reset_icon_for_raspberry_in_protein_pairs_tab,
                "green": self.reset_icon_for_green_in_protein_pairs_tab,
                "tv_green": self.reset_icon_for_tv_green_in_protein_pairs_tab,
                "palegreen": self.reset_icon_for_palegreen_in_protein_pairs_tab,
                "forest": self.reset_icon_for_forest_in_protein_pairs_tab,
                "blue": self.reset_icon_for_blue_in_protein_pairs_tab,
                "tv_blue": self.reset_icon_for_tv_blue_in_protein_pairs_tab,
                "lightblue": self.reset_icon_for_lightblue_in_protein_pairs_tab,
                "skyblue": self.reset_icon_for_skyblue_in_protein_pairs_tab,
                "yellow": self.reset_icon_for_yellow_in_protein_pairs_tab,
                "tv_yellow": self.reset_icon_for_tv_yellow_in_protein_pairs_tab,
                "paleyellow": self.reset_icon_for_paleyellow_in_protein_pairs_tab,
                "sand": self.reset_icon_for_sand_in_protein_pairs_tab,
                "magenta": self.reset_icon_for_magenta_in_protein_pairs_tab,
                "purple": self.reset_icon_for_purple_in_protein_pairs_tab,
                "pink": self.reset_icon_for_pink_in_protein_pairs_tab,
                "hotpink": self.reset_icon_for_hotpink_in_protein_pairs_tab,
                "cyan": self.reset_icon_for_cyan_in_protein_pairs_tab,
                "aquamarine": self.reset_icon_for_aquamarine_in_protein_pairs_tab,
                "palecyan": self.reset_icon_for_palecyan_in_protein_pairs_tab,
                "teal": self.reset_icon_for_teal_in_protein_pairs_tab,
                "orange": self.reset_icon_for_orange_in_protein_pairs_tab,
                "tv_orange": self.reset_icon_for_tv_orange_in_protein_pairs_tab,
                "lightorange": self.reset_icon_for_lightorange_in_protein_pairs_tab,
                "olive": self.reset_icon_for_olive_in_protein_pairs_tab,
                "white": self.reset_icon_for_white_in_protein_pairs_tab,
                "grey70": self.reset_icon_for_grey_70_in_protein_pairs_tab,
                "grey30": self.reset_icon_for_grey_30_in_protein_pairs_tab,
                "black": self.reset_icon_for_black_in_protein_pairs_tab
            }
            color_index_functions[tmp_color_name.strip()]()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # <editor-fold desc="Set color and icon">
    def set_color_name_in_label_red_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("red    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "red"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_red.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_red.setIconSize(self._view.color_grid_protein_pairs.c_red.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_tv_red_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("tv_red    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "tv_red"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_tv_red.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_tv_red.setIconSize(self._view.color_grid_protein_pairs.c_tv_red.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_salmon_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("salmon    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "salmon"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_salomon.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_salomon.setIconSize(self._view.color_grid_protein_pairs.c_salomon.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_raspberry_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("raspberry    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "raspberry"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_raspberry.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_raspberry.setIconSize(self._view.color_grid_protein_pairs.c_raspberry.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_green_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("green    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "green"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_green.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_green.setIconSize(
            self._view.color_grid_protein_pairs.c_green.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_tv_green_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("tv_green    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "tv_green"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_tv_green.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_tv_green.setIconSize(
            self._view.color_grid_protein_pairs.c_tv_green.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_palegreen_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("palegreen    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "palegreen"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_palegreen.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_palegreen.setIconSize(
            self._view.color_grid_protein_pairs.c_palegreen.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_forest_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("forest    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "forest"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_forest.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_forest.setIconSize(
            self._view.color_grid_protein_pairs.c_forest.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_blue_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("blue    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "blue"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_blue.setIcon(QtGui.QIcon(
            ":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_blue.setIconSize(self._view.color_grid_protein_pairs.c_blue.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_tv_blue_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("tv_blue    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "tv_blue"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_tv_blue.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_tv_blue.setIconSize(
            self._view.color_grid_protein_pairs.c_tv_blue.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_lightblue_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("lightblue    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "lightblue"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_lightblue.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_lightblue.setIconSize(
            self._view.color_grid_protein_pairs.c_lightblue.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_skyblue_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("skyblue    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "skyblue"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_skyblue.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_skyblue.setIconSize(
            self._view.color_grid_protein_pairs.c_skyblue.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_yellow_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("yellow    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "yellow"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_yellow.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_yellow.setIconSize(
            self._view.color_grid_protein_pairs.c_yellow.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_tv_yellow_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("tv_yellow    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "tv_yellow"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_tv_yellow.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_tv_yellow.setIconSize(
            self._view.color_grid_protein_pairs.c_tv_yellow.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_paleyellow_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("paleyellow    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "paleyellow"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_paleyellow.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_paleyellow.setIconSize(
            self._view.color_grid_protein_pairs.c_paleyellow.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_sand_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("sand    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "sand"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_sand.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_sand.setIconSize(
            self._view.color_grid_protein_pairs.c_sand.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_magenta_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("magenta    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "magenta"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_magenta.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_magenta.setIconSize(
            self._view.color_grid_protein_pairs.c_magenta.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_purple_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("purple    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "purple"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_purple.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_purple.setIconSize(
            self._view.color_grid_protein_pairs.c_purple.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_pink_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("pink    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "pink"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_pink.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_pink.setIconSize(
            self._view.color_grid_protein_pairs.c_pink.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_hotpink_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("hotpink    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "hotpink"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_hotpink.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_hotpink.setIconSize(
            self._view.color_grid_protein_pairs.c_hotpink.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_cyan_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("cyan    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "cyan"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_cyan.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_cyan.setIconSize(
            self._view.color_grid_protein_pairs.c_cyan.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_aquamarine_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("aquamarine    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "aquamarine"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_aquamarine.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_aquamarine.setIconSize(
            self._view.color_grid_protein_pairs.c_aquamarine.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_palecyan_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("palecyan    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "palecyan"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_palecyan.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_palecyan.setIconSize(
            self._view.color_grid_protein_pairs.c_palecyan.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_teal_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("teal    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "teal"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_teal.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_teal.setIconSize(
            self._view.color_grid_protein_pairs.c_teal.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_orange_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("orange    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "orange"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_orange.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_orange.setIconSize(
            self._view.color_grid_protein_pairs.c_orange.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_tv_orange_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("tv_orange    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "tv_orange"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_tv_orange.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_tv_orange.setIconSize(
            self._view.color_grid_protein_pairs.c_tv_orange.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_lightorange_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("lightorange    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "lightorange"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_lightorange.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_lightorange.setIconSize(
            self._view.color_grid_protein_pairs.c_lightorange.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_olive_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("olive    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "olive"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_olive.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_olive.setIconSize(
            self._view.color_grid_protein_pairs.c_olive.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_white_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("white    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "white"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_white.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_white.setIconSize(
            self._view.color_grid_protein_pairs.c_white.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_grey_70_in_protein_pairs_tab(self):
        tmp_protein = self._interface_manager.get_current_active_protein_object_of_protein_pair()
        tmp_chain = self._interface_manager.get_current_active_chain_object_of_protein_pair()
        tmp_protein.pymol_selection.set_selection_for_a_single_chain(tmp_chain.chain_letter)
        tmp_is_colored_by_elements = tmp_chain.get_color(tmp_protein.pymol_selection.selection_string,
                                                         self._interface_manager.pymol_session_manager)[1]
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("grey70    ")
        self.__slot_change_chain_color_protein_pairs()
        if tmp_is_colored_by_elements:
            self.__slot_change_chain_color_protein_pairs_atoms()
        else:
            self._view.color_grid_protein_pairs.last_clicked_color = "grey70"
        self._view.color_grid_protein_pairs.c_grey_70.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_grey_70.setIconSize(
            self._view.color_grid_protein_pairs.c_grey_70.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_grey_30_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("grey30    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "grey30"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_grey_30.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_grey_30.setIconSize(
            self._view.color_grid_protein_pairs.c_grey_30.icon().actualSize(QtCore.QSize(14, 14)))

    def set_color_name_in_label_black_in_protein_pairs_tab(self):
        self.reset_icon_for_last_color_in_protein_pairs_tab()
        self._view.ui.lbl_protein_pair_current_color.setText("black    ")
        self._view.color_grid_protein_pairs.last_clicked_color = "black"
        self.__slot_change_chain_color_protein_pairs()
        self._view.color_grid_protein_pairs.c_black.setIcon(QtGui.QIcon(":icons/done_round_edges_w200_g200.svg"))
        self._view.color_grid_protein_pairs.c_black.setIconSize(
            self._view.color_grid_protein_pairs.c_black.icon().actualSize(QtCore.QSize(14, 14)))
    # </editor-fold>

    # <editor-fold desc="Reset Icon">
    def reset_icon_for_red_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("red")
        self._view.color_grid_protein_pairs.c_red.setIcon(QtGui.QIcon())

    def reset_icon_for_tv_red_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("tv_red")
        self._view.color_grid_protein_pairs.c_tv_red.setIcon(QtGui.QIcon())

    def reset_icon_for_salmon_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("salmon")
        self._view.color_grid_protein_pairs.c_salomon.setIcon(QtGui.QIcon())

    def reset_icon_for_raspberry_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("raspberry")
        self._view.color_grid_protein_pairs.c_raspberry.setIcon(QtGui.QIcon())

    def reset_icon_for_green_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("green")
        self._view.color_grid_protein_pairs.c_green.setIcon(QtGui.QIcon())

    def reset_icon_for_tv_green_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("tv_green")
        self._view.color_grid_protein_pairs.c_tv_green.setIcon(QtGui.QIcon())

    def reset_icon_for_palegreen_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("palegreen")
        self._view.color_grid_protein_pairs.c_palegreen.setIcon(QtGui.QIcon())

    def reset_icon_for_forest_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("forest")
        self._view.color_grid_protein_pairs.c_forest.setIcon(QtGui.QIcon())

    def reset_icon_for_blue_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("blue")
        self._view.color_grid_protein_pairs.c_blue.setIcon(QtGui.QIcon())

    def reset_icon_for_tv_blue_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("tv_blue")
        self._view.color_grid_protein_pairs.c_tv_blue.setIcon(QtGui.QIcon())

    def reset_icon_for_lightblue_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("lightblue")
        self._view.color_grid_protein_pairs.c_lightblue.setIcon(QtGui.QIcon())

    def reset_icon_for_skyblue_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("skyblue")
        self._view.color_grid_protein_pairs.c_skyblue.setIcon(QtGui.QIcon())

    def reset_icon_for_yellow_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("yellow")
        self._view.color_grid_protein_pairs.c_yellow.setIcon(QtGui.QIcon())

    def reset_icon_for_tv_yellow_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("tv_yellow")
        self._view.color_grid_protein_pairs.c_tv_yellow.setIcon(QtGui.QIcon())

    def reset_icon_for_paleyellow_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("paleyellow")
        self._view.color_grid_protein_pairs.c_paleyellow.setIcon(QtGui.QIcon())

    def reset_icon_for_sand_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("sand")
        self._view.color_grid_protein_pairs.c_sand.setIcon(QtGui.QIcon())

    def reset_icon_for_magenta_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("magenta")
        self._view.color_grid_protein_pairs.c_magenta.setIcon(QtGui.QIcon())

    def reset_icon_for_purple_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("purple")
        self._view.color_grid_protein_pairs.c_purple.setIcon(QtGui.QIcon())

    def reset_icon_for_pink_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("pink")
        self._view.color_grid_protein_pairs.c_pink.setIcon(QtGui.QIcon())

    def reset_icon_for_hotpink_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("hotpink")
        self._view.color_grid_protein_pairs.c_hotpink.setIcon(QtGui.QIcon())

    def reset_icon_for_cyan_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("cyan")
        self._view.color_grid_protein_pairs.c_cyan.setIcon(QtGui.QIcon())

    def reset_icon_for_aquamarine_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("aquamarine")
        self._view.color_grid_protein_pairs.c_aquamarine.setIcon(QtGui.QIcon())

    def reset_icon_for_palecyan_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("palecyan")
        self._view.color_grid_protein_pairs.c_palecyan.setIcon(QtGui.QIcon())

    def reset_icon_for_teal_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("teal")
        self._view.color_grid_protein_pairs.c_teal.setIcon(QtGui.QIcon())

    def reset_icon_for_orange_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("orange")
        self._view.color_grid_protein_pairs.c_orange.setIcon(QtGui.QIcon())

    def reset_icon_for_tv_orange_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("tv_orange")
        self._view.color_grid_protein_pairs.c_tv_orange.setIcon(QtGui.QIcon())

    def reset_icon_for_lightorange_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("lightorange")
        self._view.color_grid_protein_pairs.c_lightorange.setIcon(QtGui.QIcon())

    def reset_icon_for_olive_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("olive")
        self._view.color_grid_protein_pairs.c_olive.setIcon(QtGui.QIcon())

    def reset_icon_for_white_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("white")
        self._view.color_grid_protein_pairs.c_white.setIcon(QtGui.QIcon())

    def reset_icon_for_grey_70_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("grey70")
        self._view.color_grid_protein_pairs.c_grey_70.setIcon(QtGui.QIcon())

    def reset_icon_for_grey_30_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("grey30")
        self._view.color_grid_protein_pairs.c_grey_30.setIcon(QtGui.QIcon())

    def reset_icon_for_black_in_protein_pairs_tab(self):
        self._view.ui.lbl_protein_pair_current_color.setText("black")
        self._view.color_grid_protein_pairs.c_black.setIcon(QtGui.QIcon())
    # </editor-fold>
    # </editor-fold>

    # <editor-fold desc="Representations">
    # hydrogens
    def __slot_protein_pair_chain_with_hydrogens(self):
        pass


    def __slot_show_protein_pair_chain_with_hydrogens(self):
        try:
            tmp_protein = self._interface_manager.get_current_active_protein_object_of_protein_pair()
            tmp_chain = self._interface_manager.get_current_active_chain_object_of_protein_pair()
            if self._view.tg_protein_pair_sticks.toggle_button.isChecked():
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    "sticks", f"h. and chain {tmp_chain.chain_letter} and {tmp_protein.get_molecule_object()}"
                )
            if self._view.tg_protein_pair_lines.toggle_button.isChecked():
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    "lines", f"h. and chain {tmp_chain.chain_letter} and {tmp_protein.get_molecule_object()}"
                )
            if self._view.tg_protein_pair_spheres.toggle_button.isChecked():
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    "spheres", f"h. and chain {tmp_chain.chain_letter} and {tmp_protein.get_molecule_object()}"
                )
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    def __slot_hide_protein_pair_chain_with_hydrogens(self):
        try:
            tmp_protein = self._interface_manager.get_current_active_protein_object_of_protein_pair()
            tmp_chain = self._interface_manager.get_current_active_chain_object_of_protein_pair()
            if self._view.tg_protein_pair_sticks.toggle_button.isChecked():
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    "sticks", f"h. and chain {tmp_chain.chain_letter} and {tmp_protein.get_molecule_object()}"
                )
            if self._view.tg_protein_pair_lines.toggle_button.isChecked():
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    "lines", f"h. and chain {tmp_chain.chain_letter} and {tmp_protein.get_molecule_object()}"
                )
            if self._view.tg_protein_pair_spheres.toggle_button.isChecked():
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    "spheres", f"h. and chain {tmp_chain.chain_letter} and {tmp_protein.get_molecule_object()}"
                )
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # cartoon
    def __slot_protein_pair_chain_as_cartoon(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Cartoon' toggle on the 'Protein Pairs Tab' was clicked.")
            tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
            tmp_selection.set_selection_for_a_single_chain(
                self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
            if self._view.ui.cb_protein_pair_cartoon.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 0:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.CARTOON.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.CARTOON.value, tmp_selection.selection_string
                )
            elif self._view.tg_protein_pair_cartoon.toggle_button.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 1:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.CARTOON.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.CARTOON.value, tmp_selection.selection_string
                )
            else:
                #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.CARTOON.value)
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    enums.PyMOLRepresentation.CARTOON.value, tmp_selection.selection_string
                )
            self._update_scene()
            self._save_protein_pair_pymol_session()
            self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()
            self._interface_manager.manage_hydrogen_representation_for_protein_pair_chain()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # sticks
    def __slot_protein_pair_chain_as_sticks(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Sticks' toggle on the 'Protein Pairs Tab' was clicked.")
            tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
            tmp_selection.set_selection_for_a_single_chain(
                self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
            if self._view.ui.cb_protein_pair_sticks.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 0:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.STICKS.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.STICKS.value, tmp_selection.selection_string
                )
            elif self._view.tg_protein_pair_sticks.toggle_button.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 1:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.STICKS.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.STICKS.value, tmp_selection.selection_string
                )
            else:
                #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.STICKS.value)
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    enums.PyMOLRepresentation.STICKS.value, tmp_selection.selection_string
                )
            self._update_scene()
            self._save_protein_pair_pymol_session()
            self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()
            self._interface_manager.manage_hydrogen_representation_for_protein_pair_chain()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # ribbon
    def __slot_protein_pair_chain_as_ribbon(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Ribbon' toggle on the 'Protein Pairs Tab' was clicked.")
            tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
            tmp_selection.set_selection_for_a_single_chain(
                self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
            if self._view.ui.cb_protein_pair_ribbon.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 0:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.RIBBON.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.RIBBON.value, tmp_selection.selection_string
                )
            elif self._view.tg_protein_pair_ribbon.toggle_button.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 1:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.RIBBON.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.RIBBON.value, tmp_selection.selection_string
                )
            else:
                #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.RIBBON.value)
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    enums.PyMOLRepresentation.RIBBON.value, tmp_selection.selection_string
                )
            self._update_scene()
            self._save_protein_pair_pymol_session()
            self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()
            self._interface_manager.manage_hydrogen_representation_for_protein_pair_chain()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # lines
    def __slot_protein_pair_chain_as_lines(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Lines' toggle on the 'Protein Pairs Tab' was clicked.")
            tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
            tmp_selection.set_selection_for_a_single_chain(
                self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
            if self._view.ui.cb_protein_pair_lines.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 0:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.LINES.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.LINES.value, tmp_selection.selection_string
                )
            elif self._view.tg_protein_pair_lines.toggle_button.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 1:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.LINES.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.LINES.value, tmp_selection.selection_string
                )
            else:
                #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.LINES.value)
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    enums.PyMOLRepresentation.LINES.value, tmp_selection.selection_string
                )
            self._update_scene()
            self._save_protein_pair_pymol_session()
            self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()
            self._interface_manager.manage_hydrogen_representation_for_protein_pair_chain()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # spheres
    def __slot_protein_pair_chain_as_spheres(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Spheres' toggle on the 'Protein Pairs Tab' was clicked.")
            tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
            tmp_selection.set_selection_for_a_single_chain(
                self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
            if self._view.ui.cb_protein_pair_spheres.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 0:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SPHERES.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.SPHERES.value, tmp_selection.selection_string
                )
            elif self._view.tg_protein_pair_spheres.toggle_button.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 1:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SPHERES.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.SPHERES.value, tmp_selection.selection_string
                )
            else:
                #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.SPHERES.value)
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    enums.PyMOLRepresentation.SPHERES.value, tmp_selection.selection_string
                )
            self._update_scene()
            self._save_protein_pair_pymol_session()
            self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()
            self._interface_manager.manage_hydrogen_representation_for_protein_pair_chain()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # dots
    def __slot_protein_pair_chain_as_dots(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Dots' toggle on the 'Protein Pairs Tab' was clicked.")
            tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
            tmp_selection.set_selection_for_a_single_chain(
                self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
            if self._view.ui.cb_protein_pair_dots.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 0:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.DOTS.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.DOTS.value, tmp_selection.selection_string
                )
            elif self._view.tg_protein_pair_dots.toggle_button.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 1:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.DOTS.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.DOTS.value, tmp_selection.selection_string
                )
            else:
                #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.DOTS.value)
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    enums.PyMOLRepresentation.DOTS.value, tmp_selection.selection_string
                )
            self._update_scene()
            self._save_protein_pair_pymol_session()
            self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()
            self._interface_manager.manage_hydrogen_representation_for_protein_pair_chain()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # mesh
    def __slot_protein_pair_chain_as_mesh(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Mesh' toggle on the 'Protein Pairs Tab' was clicked.")
            tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
            tmp_selection.set_selection_for_a_single_chain(
                self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
            if self._view.ui.cb_protein_pair_mesh.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 0:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.MESH.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.MESH.value, tmp_selection.selection_string
                )
            elif self._view.tg_protein_pair_mesh.toggle_button.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 1:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.MESH.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.MESH.value, tmp_selection.selection_string
                )
            else:
                #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.MESH.value)
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    enums.PyMOLRepresentation.MESH.value, tmp_selection.selection_string
                )
            self._update_scene()
            self._save_protein_pair_pymol_session()
            self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()
            self._interface_manager.manage_hydrogen_representation_for_protein_pair_chain()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # surface
    def __slot_protein_pair_chain_as_surface(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Surface' toggle on the 'Protein Pairs Tab' was clicked.")
            tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
            tmp_selection.set_selection_for_a_single_chain(
                self._interface_manager.get_current_active_chain_object_of_protein_pair().chain_letter)
            if self._view.ui.cb_protein_pair_surface.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 0:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SURFACE.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.SURFACE.value, tmp_selection.selection_string
                )
            elif self._view.tg_protein_pair_surface.toggle_button.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 1:
                #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SURFACE.value)
                self._interface_manager.pymol_session_manager.show_specific_representation(
                    enums.PyMOLRepresentation.SURFACE.value, tmp_selection.selection_string
                )
            else:
                #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.SURFACE.value)
                self._interface_manager.pymol_session_manager.hide_specific_representation(
                    enums.PyMOLRepresentation.SURFACE.value, tmp_selection.selection_string
                )
            self._update_scene()
            self._save_protein_pair_pymol_session()
            self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()
            self._interface_manager.manage_hydrogen_representation_for_protein_pair_chain()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # all
    def __slot_hide_protein_pair_chain_all(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                       "'Hide all' representations button on the 'Protein Pairs Tab' was clicked.")
            self._view.ui.cb_protein_pair_cartoon.setChecked(False)
            self._view.ui.cb_protein_pair_sticks.setChecked(False)
            self._view.ui.cb_protein_pair_ribbon.setChecked(False)
            self._view.ui.cb_protein_pair_lines.setChecked(False)
            self._view.ui.cb_protein_pair_spheres.setChecked(False)
            self._view.ui.cb_protein_pair_dots.setChecked(False)
            self._view.ui.cb_protein_pair_mesh.setChecked(False)
            self._view.ui.cb_protein_pair_surface.setChecked(False)
            self._view.tg_protein_pair_cartoon.toggle_button.setChecked(False)
            self._view.tg_protein_pair_sticks.toggle_button.setChecked(False)
            self._view.tg_protein_pair_ribbon.toggle_button.setChecked(False)
            self._view.tg_protein_pair_lines.toggle_button.setChecked(False)
            self._view.tg_protein_pair_spheres.toggle_button.setChecked(False)
            self._view.tg_protein_pair_dots.toggle_button.setChecked(False)
            self._view.tg_protein_pair_mesh.toggle_button.setChecked(False)
            self._view.tg_protein_pair_surface.toggle_button.setChecked(False)
            self.__slot_protein_pair_chain_as_cartoon()
            self.__slot_protein_pair_chain_as_sticks()
            self.__slot_protein_pair_chain_as_ribbon()
            self.__slot_protein_pair_chain_as_lines()
            self.__slot_protein_pair_chain_as_spheres()
            self.__slot_protein_pair_chain_as_dots()
            self.__slot_protein_pair_chain_as_mesh()
            self.__slot_protein_pair_chain_as_surface()
            self._update_scene()
            self._save_protein_pair_pymol_session()
            self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # </editor-fold>

    def __slot_update_protein_pair_scene(self):
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Update protein scene' button on the 'Protein Pairs Tab' was clicked.")
            self._interface_manager.pymol_session_manager.pymol_interface.scene("auto", "update")
            self._save_protein_pair_pymol_session()
            self._interface_manager.status_bar_manager.show_temporary_message("PyMOL Scene updated.", a_timeout=1500)
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    def _check_for_results(self) -> None:
        if self._view.ui.protein_pairs_tree_view.model().data(self._view.ui.protein_pairs_tree_view.currentIndex(), Qt.DisplayRole).find("_vs_") != -1:
            self._view.ui.action_results_summary.setEnabled(True)
        else:
            self._view.ui.action_results_summary.setEnabled(False)

    def __slot_results_summary(self) -> None:
        try:
            logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Results/Summary' clicked.")
            tmp_protein_pair = self._view.ui.protein_pairs_tree_view.model().data(self._view.ui.protein_pairs_tree_view.currentIndex(),
                                                                                  enums.ModelEnum.OBJECT_ROLE)
            self._external_controller = results_view_controller.ResultsViewController(self._interface_manager,
                                                                                      tmp_protein_pair,
                                                                                      self._interface_manager.pymol_session_manager)
            self._interface_manager.get_results_view().show()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

    # </editor-fold>
