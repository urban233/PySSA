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
"""Module for the main view controller."""
import copy
import logging
import os
import pathlib
import shutil
import subprocess
from typing import Optional, Any

import pygetwindow
from Bio import SeqRecord
from Bio.Seq import Seq

from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from PyQt5 import QtGui

from pyssa.gui.ui import icon_resources  # this import is used for the icons! DO NOT DELETE THIS
from pyssa.gui.ui.custom_context_menus import protein_tree_context_menu, protein_pair_tree_context_menu, \
  sequence_list_context_menu
from pyssa.gui.ui.custom_dialogs import custom_message_box
from pyssa.gui.ui.custom_widgets import job_entry
from pyssa.internal.thread.async_pyssa import util_async, custom_signals, project_async, image_async, \
  pymol_session_async, protein_async, sequence_async, protein_pair_async
from pyssa.controller import results_view_controller, rename_protein_view_controller, use_project_view_controller, \
  pymol_session_manager, add_sequence_view_controller, add_scene_view_controller, add_protein_view_controller, \
  settings_view_controller, predict_protein_view_controller, import_sequence_view_controller, \
  rename_sequence_view_controller
from pyssa.internal.data_structures import chain, job
from pyssa.internal.data_structures.data_classes import residue_color_config
from pyssa.gui.ui.dialogs import dialog_settings_global, dialog_tutorial_videos, dialog_about
from pyssa.internal.data_structures import project, settings, protein, protein_pair
from pyssa.internal.data_structures.data_classes import database_operation
from pyssa.internal.thread import database_thread
from pyssa.io_pyssa import filesystem_io
from pyssa.logging_pyssa import log_handlers, log_levels
from pyssa.util import constants, enums, exit_codes, tools, ui_util, exception
from pyssa.gui.ui.views import main_view
from pyssa.controller import interface_manager, distance_analysis_view_controller, delete_project_view_controller, \
  create_project_view_controller, open_project_view_controller, database_manager
from pyssa.util import globals
from pyssa_pymol import pymol_enums
from tea.thread import tasks, task_manager, task_result_factory, task_result, action
from tea.thread import tasks, task_result_factory, task_result, action

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class MainViewController:
  """Class for main presenter of the pyssa plugin."""

  # <editor-fold desc="Class attributes">
  _view: "main_view.MainView"
  """The main view of the pyssa plugin."""

  _application_settings: "settings.Settings"
  """The application settings object."""

  _interface_manager: "interface_manager.InterfaceManager"
  """A watcher for the project state."""

  _pymol_session_manager: "pymol_session_manager.PymolSessionManager"
  """A manager for the pymol session."""

  _active_task: tasks.LegacyTask
  """The active task of the application."""

  _help_task: tasks.LegacyTask
  """The active task which controls the help window of the application."""

  _database_thread: "database_thread.DatabaseThread"
  """A thread for database related processes."""

  # </editor-fold>

  def __init__(
          self, the_interface_manager: "interface_manager.InterfaceManager"
  ) -> None:
    """Constructor.

    Args:
        the_interface_manager (interface_manager.InterfaceManager): The InterfaceManager object.

    Raises:
        exception.IllegalArgumentError: If `the_interface_manager` is None.
    """
    # <editor-fold desc="Checks">
    if the_interface_manager is None:
      logger.error("the_interface_manager is None.")
      raise exception.IllegalArgumentError("the_interface_manager is None.")

    # </editor-fold>

    self._view: "main_view.MainView" = the_interface_manager.get_main_view()
    self._interface_manager: "interface_manager.InterfaceManager" = (
      the_interface_manager
    )
    self._database_manager = database_manager.DatabaseManager("")
    self._database_manager.set_application_settings(
      self._interface_manager.get_application_settings()
    )
    self._database_thread: "database_thread.DatabaseThread" = (
      database_thread.DatabaseThread("")
    )
    self._external_view = None
    self.active_custom_message_box: "custom_message_box.CustomMessageBoxOk" = (
      None
    )
    # self._main_view_state = main_view_state.MainViewState(
    #     self._view.ui.seqs_list_view,
    #     self.__slot_show_sequence_information,
    #     self._view.ui.proteins_tree_view,
    #     self.__slot_get_information_about_selected_object_in_protein_branch,
    #     self._view.ui.protein_pairs_tree_view,
    #     self.__slot_get_information_about_selected_object_in_protein_pair_branch,
    # )
    self.custom_progress_signal = custom_signals.ProgressSignal()
    self.abort_signal = custom_signals.AbortSignal()
    self._sequence_list_context_menu = (
      sequence_list_context_menu.SequenceListContextMenu()
    )
    self._protein_tree_context_menu = (
      protein_tree_context_menu.ProteinTreeContextMenu()
    )
    self._protein_pair_tree_context_menu = (
      protein_pair_tree_context_menu.ProteinPairTreeContextMenu()
    )

    self.thread_pool = QtCore.QThreadPool()
    self.thread_pool.setMaxThreadCount(os.cpu_count())

    self._setup_statusbar()
    self._init_generic_help_context_menus()
    self._interface_manager.refresh_main_view()
    self._connect_all_ui_elements_with_slot_functions()
    if (
            self._interface_manager.get_application_settings().start_help_at_startup
            == 1
    ):
      self._start_documentation_server()

  def _connect_all_ui_elements_with_slot_functions(self) -> None:
    """Connects all UI elements to their corresponding slot functions in the class."""
    self._view.dialogClosed.connect(self._close_main_window)
    self._interface_manager.get_restart_pymol_view().return_value.connect(
      self._force_close_all
    )
    self.custom_progress_signal.progress.connect(self._update_progress_bar)
    self.abort_signal.abort.connect(self._abort_task)
    self._interface_manager.refresh_after_job_finished_signal.refresh.connect(
      self._update_main_view_ui
    )
    self._view.btn_open_job_overview.clicked.connect(
      self.__slot_open_job_overview_panel
    )
    self._view.btn_open_job_notification.clicked.connect(
      self.__slot_open_notification_panel
    )

    # <editor-fold desc="Menu">
    self._view.ui.action_new_project.triggered.connect(
      self.__slot_create_project
    )
    self._interface_manager.get_create_view().dialogClosed.connect(
      self.__slot_refresh_main_view
    )
    self._view.ui.action_open_project.triggered.connect(
      self.__slot_open_project
    )
    self._interface_manager.get_open_view().dialogClosed.connect(
      self.__slot_refresh_main_view
    )
    self._view.ui.action_use_project.triggered.connect(self.__slot_use_project)
    self._view.ui.action_delete_project.triggered.connect(
      self.__slot_delete_project
    )
    self._interface_manager.get_delete_view().dialogClosed.connect(
      self.__slot_refresh_main_view
    )
    self._view.ui.action_import_project.triggered.connect(
      self.__slot_import_project
    )
    self._view.ui.action_export_project.triggered.connect(
      self.__slot_export_current_project
    )
    self._view.ui.action_close_project.triggered.connect(
      self.__slot_close_project
    )
    self._view.ui.action_exit_application.triggered.connect(
      self.__slot_close_all
    )

    self._view.ui.action_results_summary.triggered.connect(
      self.__slot_results_summary
    )
    self._view.ui.action_preview_image.triggered.connect(
      self.__slot_preview_image
    )
    self._view.ui.action_ray_tracing_image.triggered.connect(
      self.__slot_create_ray_traced_image
    )
    self._view.ui.action_simple_image.triggered.connect(
      self.__slot_create_drawn_image
    )
    self._view.ui.action_protein_regions.triggered.connect(
      self.__slot_hotspots_protein_regions
    )

    self._view.ui.action_edit_settings.triggered.connect(
      self.__slot_open_settings_global
    )
    self._view.ui.action_restore_settings.triggered.connect(
      self.__slot_restore_settings
    )
    self._view.ui.action_show_log_in_explorer.triggered.connect(
      self.__slot_open_logs
    )
    self._view.ui.action_clear_logs.triggered.connect(
      self.__slot_clear_all_log_files
    )
    self._view.ui.action_documentation.triggered.connect(
      self.__slot_open_help_center
    )
    self._view.ui.action_tutorials.triggered.connect(self.__slot_open_tutorial)
    self._view.ui.action_get_demo_projects.triggered.connect(
      self.__slot_get_demo_projects
    )
    self._view.ui.action_restart_pymol.triggered.connect(
      self.__slot_restart_pymol
    )
    self._view.ui.action_about.triggered.connect(self.__slot_open_about)
    self._view.ui.action_predict_monomer.triggered.connect(
      self.__slot_predict_monomer
    )
    self._view.ui.action_predict_multimer.triggered.connect(
      self.__slot_predict_multimer
    )
    self._view.ui.action_abort_prediction.triggered.connect(
      self.__slot_abort_prediction
    )
    self._view.ui.action_distance_analysis.triggered.connect(
      self.__slot_distance_analysis
    )
    self._view.ui.action_arrange_windows.triggered.connect(
      self.__slot_arrange_windows
    )

    self._view.ui.project_tab_widget.currentChanged.connect(self._update_tab)
    # </editor-fold>

    # <editor-fold desc="Sequence Tab">
    self._view.ui.seqs_list_view.customContextMenuRequested.connect(
      self.open_context_menu_for_sequences
    )
    self._view.ui.seqs_list_view.clicked.connect(
      self.__slot_show_sequence_information
    )
    self._view.ui.btn_add_sequence.clicked.connect(self.__slot_add_sequence)
    self._view.ui.btn_import_seq.clicked.connect(self.__slot_import_sequence)
    self._view.ui.btn_save_sequence.clicked.connect(
      self.__slot_save_selected_sequence_as_fasta_file
    )
    self._view.ui.btn_delete_sequence.clicked.connect(
      self.__slot_delete_selected_sequence
    )
    self._view.ui.seqs_table_widget.cellClicked.connect(
      self.__slot_open_text_editor_for_seq
    )
    self._view.line_edit_seq_name.textChanged.connect(
      self._set_new_sequence_name_in_table_item
    )
    # self._view.ui.seqs_table_widget.cellChanged.connect(self._rename_sequence)
    self._view.ui.btn_help.clicked.connect(self.__slot_open_sequences_tab_help)

    # <editor-fold desc="Context menu">
    self._sequence_list_context_menu.connect_rename_sequence_action(
      self.__slot_rename_selected_sequence
    )
    self._sequence_list_context_menu.connect_help_action(
      self.__slot_open_sequences_tab_help
    )
    # </editor-fold>

    # </editor-fold>

    # <editor-fold desc="Proteins Tab">
    self._view.ui.proteins_tree_view.customContextMenuRequested.connect(
      self._open_context_menu_for_proteins
    )
    self._view.ui.proteins_tree_view.clicked.connect(
      self.__slot_get_information_about_selected_object_in_protein_branch
    )
    self._view.ui.btn_protein_tree_view_expand.clicked.connect(
      self.__slot_expand_all_proteins
    )
    self._view.ui.btn_protein_tree_view_collapse.clicked.connect(
      self.__slot_collapse_all_proteins
    )
    self._view.ui.btn_save_protein.clicked.connect(
      self.__slot_save_selected_protein_structure_as_pdb_file
    )
    # import
    self._view.ui.btn_import_protein.clicked.connect(
      self.__slot_import_protein_structure
    )
    self._interface_manager.get_add_protein_view().return_value.connect(
      self._post_import_protein_structure
    )
    self._view.ui.btn_open_protein_session.clicked.connect(
      self.__slot_open_protein_pymol_session
    )
    self._view.ui.btn_create_protein_scene.clicked.connect(
      self.__slot_save_scene
    )
    self._view.ui.btn_delete_protein.clicked.connect(self.__slot_delete_protein)
    self._view.ui.btn_update_protein_scene.clicked.connect(
      self.__slot_update_protein_scene
    )
    self._view.ui.btn_delete_protein_scene.clicked.connect(
      self.__slot_delete_current_scene
    )
    self._view.ui.box_protein_color.currentIndexChanged.connect(
      self.__slot_change_chain_color_proteins
    )
    # self._view.ui.btn_protein_color_atoms.clicked.connect(self.__slot_change_chain_color_proteins_atoms)
    # self._view.ui.btn_protein_reset_atoms.clicked.connect(self.__slot_change_chain_reset_proteins_atoms)
    self._view.tg_protein_white_bg.toggleChanged.connect(
      self.__slot_protein_change_background_color
    )
    # self._view.ui.btn_protein_show_cartoon.clicked.connect(self.__slot_show_protein_chain_as_cartoon)
    # self._view.ui.btn_protein_hide_cartoon.clicked.connect(self.__slot_hide_protein_chain_as_cartoon)
    # self._view.ui.btn_protein_show_sticks.clicked.connect(self.__slot_show_protein_chain_as_sticks)
    # self._view.ui.btn_protein_hide_sticks.clicked.connect(self.__slot_hide_protein_chain_as_sticks)
    # self._view.ui.btn_protein_show_ribbon.clicked.connect(self.__slot_show_protein_chain_as_ribbon)
    # self._view.ui.btn_protein_hide_ribbon.clicked.connect(self.__slot_hide_protein_chain_as_ribbon)
    # self._view.ui.cb_protein_atoms.stateChanged.connect(self.__slot_color_protein_chain_atoms_by_element)
    # self._view.ui.cb_protein_cartoon.stateChanged.connect(self.__slot_protein_chain_as_cartoon)
    # self._view.ui.cb_protein_sticks.stateChanged.connect(self.__slot_protein_chain_as_sticks)
    # self._view.ui.cb_protein_ribbon.stateChanged.connect(self.__slot_protein_chain_as_ribbon)
    # self._view.ui.cb_protein_lines.stateChanged.connect(self.__slot_protein_chain_as_lines)
    # self._view.ui.cb_protein_spheres.stateChanged.connect(self.__slot_protein_chain_as_spheres)
    # self._view.ui.cb_protein_dots.stateChanged.connect(self.__slot_protein_chain_as_dots)
    # self._view.ui.cb_protein_mesh.stateChanged.connect(self.__slot_protein_chain_as_mesh)
    # self._view.ui.cb_protein_surface.stateChanged.connect(self.__slot_protein_chain_as_surface)
    # representation
    self._view.tg_protein_color_atoms.toggleChanged.connect(
      self.__slot_color_protein_atoms_by_element
    )
    # self._view.tg_protein_hydrogen_atoms.toggleChanged.connect(self.__slot_chain_protein_with_hydrogens)  # this could be useful
    self._view.tg_protein_cartoon.toggleChanged.connect(
      self.__slot_protein_chain_as_cartoon
    )
    self._view.tg_protein_sticks.toggleChanged.connect(
      self.__slot_protein_chain_as_sticks
    )
    self._view.tg_protein_ribbon.toggleChanged.connect(
      self.__slot_protein_chain_as_ribbon
    )
    self._view.tg_protein_lines.toggleChanged.connect(
      self.__slot_protein_chain_as_lines
    )
    self._view.tg_protein_spheres.toggleChanged.connect(
      self.__slot_protein_chain_as_spheres
    )
    self._view.tg_protein_dots.toggleChanged.connect(
      self.__slot_protein_chain_as_dots
    )
    self._view.tg_protein_mesh.toggleChanged.connect(
      self.__slot_protein_chain_as_mesh
    )
    self._view.tg_protein_surface.toggleChanged.connect(
      self.__slot_protein_chain_as_surface
    )

    self._view.ui.btn_protein_hide_all_representations.clicked.connect(
      self.__slot_hide_protein_chain_all
    )
    self._view.ui.btn_help_2.clicked.connect(self.__slot_open_proteins_tab_help)

    # <editor-fold desc="Color Grid">
    # It needs to be checked if the lambda function usage is problematic in this context.
    self._view.color_grid_proteins.c_red.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("red")
    )
    self._view.color_grid_proteins.c_tv_red.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("tv_red")
    )
    self._view.color_grid_proteins.c_salomon.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("salmon")
    )
    self._view.color_grid_proteins.c_raspberry.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("raspberry")
    )

    self._view.color_grid_proteins.c_green.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("green")
    )
    self._view.color_grid_proteins.c_tv_green.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("tv_green")
    )
    self._view.color_grid_proteins.c_palegreen.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("palegreen")
    )
    self._view.color_grid_proteins.c_forest.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("forest")
    )

    self._view.color_grid_proteins.c_blue.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("blue")
    )
    self._view.color_grid_proteins.c_tv_blue.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("tv_blue")
    )
    self._view.color_grid_proteins.c_lightblue.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("lightblue")
    )
    self._view.color_grid_proteins.c_skyblue.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("skyblue")
    )

    self._view.color_grid_proteins.c_yellow.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("yellow")
    )
    self._view.color_grid_proteins.c_tv_yellow.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("tv_yellow")
    )
    self._view.color_grid_proteins.c_paleyellow.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("paleyellow")
    )
    self._view.color_grid_proteins.c_sand.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("sand")
    )

    self._view.color_grid_proteins.c_magenta.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("magenta")
    )
    self._view.color_grid_proteins.c_purple.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("purple")
    )
    self._view.color_grid_proteins.c_pink.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("pink")
    )
    self._view.color_grid_proteins.c_hotpink.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("hotpink")
    )

    self._view.color_grid_proteins.c_cyan.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("cyan")
    )
    self._view.color_grid_proteins.c_aquamarine.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("aquamarine")
    )
    self._view.color_grid_proteins.c_palecyan.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("palecyan")
    )
    self._view.color_grid_proteins.c_teal.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("teal")
    )

    self._view.color_grid_proteins.c_orange.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("orange")
    )
    self._view.color_grid_proteins.c_tv_orange.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("tv_orange")
    )
    self._view.color_grid_proteins.c_lightorange.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("lightorange")
    )
    self._view.color_grid_proteins.c_olive.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("olive")
    )

    self._view.color_grid_proteins.c_white.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("white")
    )
    self._view.color_grid_proteins.c_grey_70.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("grey70")
    )
    self._view.color_grid_proteins.c_grey_30.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("grey30")
    )
    self._view.color_grid_proteins.c_black.clicked.connect(
      lambda: self.__slot_change_chain_color_proteins("black")
    )
    # </editor-fold>

    # <editor-fold desc="Protein tree context menu">
    self._protein_tree_context_menu.connect_expand_protein_action(
      self.__slot_expand_protein
    )
    self._protein_tree_context_menu.connect_collapse_protein_action(
      self.__slot_collapse_protein
    )
    self._protein_tree_context_menu.connect_clean_protein_action(
      self.__slot_clean_protein_update
    )
    self._protein_tree_context_menu.connect_rename_protein_action(
      self.__slot_rename_selected_protein_structure
    )
    self._protein_tree_context_menu.connect_show_sequence_action(
      self.__slot_show_protein_chain_sequence
    )
    self._protein_tree_context_menu.connect_help_action(
      self.__slot_open_proteins_tab_help
    )
    # </editor-fold>

    # self._view.ui.btn_protein_sticks_show.clicked.connect(self.__slot_show_protein_regions_resi_sticks)
    # self._view.ui.btn_protein_sticks_hide.clicked.connect(self.__slot_hide_protein_regions_resi_sticks)
    # self._view.ui.btn_protein_disulfide_bonds_show.clicked.connect(
    #     self.__slot_show_protein_regions_disulfide_bonds)
    # self._view.ui.btn_protein_disulfide_bonds_hide.clicked.connect(
    #     self.__slot_hide_protein_regions_disulfide_bonds)
    # self._view.ui.btn_protein_position_zoom.clicked.connect(self.__slot_zoom_protein_regions_resi_position)

    # </editor-fold>

    # <editor-fold desc="Proteins Pair Tab">
    self._view.ui.protein_pairs_tree_view.customContextMenuRequested.connect(
      self.open_context_menu_for_protein_pairs
    )
    self._view.ui.protein_pairs_tree_view.clicked.connect(
      self.__slot_get_information_about_selected_object_in_protein_pair_branch
    )
    self._view.ui.btn_protein_pair_tree_view_expand.clicked.connect(
      self.__slot_expand_all_protein_pairs
    )
    self._view.ui.btn_protein_pair_tree_view_collapse.clicked.connect(
      self.__slot_collapse_all_protein_pairs
    )
    self._view.ui.btn_delete_protein_pair.clicked.connect(
      self.__slot_delete_protein_pair_from_project
    )
    self._view.ui.btn_open_protein_pair_session.clicked.connect(
      self.__slot_open_protein_pair_pymol_session
    )
    self._view.ui.btn_create_protein_pair_scene.clicked.connect(
      self.__slot_save_scene
    )
    self._view.ui.btn_update_protein_pair_scene.clicked.connect(
      self.__slot_update_protein_pair_scene
    )
    self._view.ui.btn_delete_protein_pair_scene.clicked.connect(
      self.__slot_delete_current_scene
    )
    self._view.ui.protein_pairs_tree_view.clicked.connect(
      self.__slot_check_for_results
    )
    self._view.ui.box_protein_pair_color.currentIndexChanged.connect(
      self.__slot_change_chain_color_protein_pairs
    )
    self._view.tg_protein_pair_color_atoms.toggleChanged.connect(
      self.__slot_color_protein_pair_atoms_by_element
    )
    self._view.tg_protein_pair_white_bg.toggleChanged.connect(
      self.__slot_protein_pair_change_background_color
    )
    # self._view.ui.cb_protein_pair_cartoon.stateChanged.connect(self.__slot_protein_pair_chain_as_cartoon)
    # self._view.ui.cb_protein_pair_sticks.stateChanged.connect(self.__slot_protein_pair_chain_as_sticks)
    # self._view.ui.cb_protein_pair_ribbon.stateChanged.connect(self.__slot_protein_pair_chain_as_ribbon)
    # self._view.ui.cb_protein_pair_lines.stateChanged.connect(self.__slot_protein_pair_chain_as_lines)
    # self._view.ui.cb_protein_pair_spheres.stateChanged.connect(self.__slot_protein_pair_chain_as_spheres)
    # self._view.ui.cb_protein_pair_dots.stateChanged.connect(self.__slot_protein_pair_chain_as_dots)
    # self._view.ui.cb_protein_pair_mesh.stateChanged.connect(self.__slot_protein_pair_chain_as_mesh)
    # self._view.ui.cb_protein_pair_surface.stateChanged.connect(self.__slot_protein_pair_chain_as_surface)
    # toggle representation
    # self._view.tg_protein_pair_hydrogen_atoms.toggleChanged.connect(self.__slot_protein_pair_chain_with_hydrogens) # this could be useful
    self._view.tg_protein_pair_cartoon.toggleChanged.connect(
      self.__slot_protein_pair_chain_as_cartoon
    )
    self._view.tg_protein_pair_sticks.toggleChanged.connect(
      self.__slot_protein_pair_chain_as_sticks
    )
    self._view.tg_protein_pair_ribbon.toggleChanged.connect(
      self.__slot_protein_pair_chain_as_ribbon
    )
    self._view.tg_protein_pair_lines.toggleChanged.connect(
      self.__slot_protein_pair_chain_as_lines
    )
    self._view.tg_protein_pair_spheres.toggleChanged.connect(
      self.__slot_protein_pair_chain_as_spheres
    )
    self._view.tg_protein_pair_dots.toggleChanged.connect(
      self.__slot_protein_pair_chain_as_dots
    )
    self._view.tg_protein_pair_mesh.toggleChanged.connect(
      self.__slot_protein_pair_chain_as_mesh
    )
    self._view.tg_protein_pair_surface.toggleChanged.connect(
      self.__slot_protein_pair_chain_as_surface
    )

    self._view.ui.btn_protein_pair_hide_all_representations.clicked.connect(
      self.__slot_hide_protein_pair_chain_all
    )
    self._view.ui.btn_help_3.clicked.connect(self._open_protein_pairs_tab_help)

    # <editor-fold desc="Color Grid">
    self._view.color_grid_protein_pairs.c_red.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("red")
    )
    self._view.color_grid_protein_pairs.c_tv_red.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("tv_red")
    )
    self._view.color_grid_protein_pairs.c_salomon.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("salmon")
    )
    self._view.color_grid_protein_pairs.c_raspberry.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("raspberry")
    )

    self._view.color_grid_protein_pairs.c_green.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("green")
    )
    self._view.color_grid_protein_pairs.c_tv_green.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("tv_green")
    )
    self._view.color_grid_protein_pairs.c_palegreen.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("palegreen")
    )
    self._view.color_grid_protein_pairs.c_forest.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("forest")
    )

    self._view.color_grid_protein_pairs.c_blue.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("blue")
    )
    self._view.color_grid_protein_pairs.c_tv_blue.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("tv_blue")
    )
    self._view.color_grid_protein_pairs.c_lightblue.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("lightblue")
    )
    self._view.color_grid_protein_pairs.c_skyblue.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("skyblue")
    )

    self._view.color_grid_protein_pairs.c_yellow.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("yellow")
    )
    self._view.color_grid_protein_pairs.c_tv_yellow.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("tv_yellow")
    )
    self._view.color_grid_protein_pairs.c_paleyellow.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("paleyellow")
    )
    self._view.color_grid_protein_pairs.c_sand.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("sand")
    )

    self._view.color_grid_protein_pairs.c_magenta.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("magenta")
    )
    self._view.color_grid_protein_pairs.c_purple.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("pink")
    )
    self._view.color_grid_protein_pairs.c_pink.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("pink")
    )
    self._view.color_grid_protein_pairs.c_hotpink.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("hotpink")
    )

    self._view.color_grid_protein_pairs.c_cyan.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("cyan")
    )
    self._view.color_grid_protein_pairs.c_aquamarine.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("aquamarine")
    )
    self._view.color_grid_protein_pairs.c_palecyan.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("palecyan")
    )
    self._view.color_grid_protein_pairs.c_teal.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("teal")
    )

    self._view.color_grid_protein_pairs.c_orange.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("orange")
    )
    self._view.color_grid_protein_pairs.c_tv_orange.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("tv_orange")
    )
    self._view.color_grid_protein_pairs.c_lightorange.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("lightorange")
    )
    self._view.color_grid_protein_pairs.c_olive.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("olive")
    )

    self._view.color_grid_protein_pairs.c_white.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("white")
    )
    self._view.color_grid_protein_pairs.c_grey_70.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("grey_70")
    )
    self._view.color_grid_protein_pairs.c_grey_30.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("grey_30")
    )
    self._view.color_grid_protein_pairs.c_black.clicked.connect(
      lambda: self.__slot_change_chain_color_protein_pairs("black")
    )

    # </editor-fold>

    # <editor-fold desc="Context menu">
    self._protein_pair_tree_context_menu.connect_expand_protein_pair_action(
      self.__slot_expand_protein_pair
    )
    self._protein_pair_tree_context_menu.connect_collapse_protein_pair_action(
      self.__slot_collapse_protein_pair
    )
    self._protein_pair_tree_context_menu.connect_open_results_summary_action(
      self.__slot_results_summary
    )
    self._protein_pair_tree_context_menu.connect_color_based_on_rmsd_action(
      self.__slot_color_protein_pair_by_rmsd
    )
    self._protein_pair_tree_context_menu.connect_help_action(
      self._open_protein_pairs_tab_help
    )
    # </editor-fold>

    # self._view.ui.btn_protein_pair_sticks_show.clicked.connect(self.__slot_show_protein_regions_resi_sticks)
    # self._view.ui.btn_protein_pair_sticks_hide.clicked.connect(self.__slot_hide_protein_regions_resi_sticks)
    # self._view.ui.btn_protein_pair_disulfide_bonds_show.clicked.connect(
    #     self.__slot_show_protein_regions_disulfide_bonds)
    # self._view.ui.btn_protein_pair_disulfide_bonds_hide.clicked.connect(
    #     self.__slot_hide_protein_regions_disulfide_bonds)
    # self._view.ui.btn_protein_pair_position_zoom.clicked.connect(self.__slot_zoom_protein_regions_resi_position)

    # </editor-fold>

  def _close_main_window(self, return_value: tuple[str, list[tuple[bool, tuple]]]) -> None:
    """Cleans after the main window closes.

    Args:
        return_value (tuple[str, list[tuple[bool, tuple]]]): A tuple containing two values: the first value is the main window object, and the second value is the event object.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      return

    # </editor-fold>

    tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
    _, tmp_event = tmp_result
    logger.info("Check if any jobs are running before closing PySSA.")
    if self._interface_manager.job_manager.there_are_jobs_running():
      logger.info("Running jobs found!")
      tmp_dialog = custom_message_box.CustomMessageBoxYesNo(
        "There are still jobs running!\n\nAre you sure you want to close PySSA?\nThis could lead to data loss and a damaged project!",
        "Close PySSA",
        custom_message_box.CustomMessageBoxIcons.WARNING.value,
      )
      tmp_dialog.exec_()
      if not tmp_dialog.response:
        logger.info("PySSA should not be closed.")
        tmp_event.ignore()
        return
      logger.info(
        "There are jobs running but PySSA will closed! (Requested by the user)"
      )
    # Closes the documentation browser if it is still open
    if (
            len(
              pygetwindow.getWindowsWithTitle(
                constants.WINDOW_TITLE_OF_HELP_CENTER
              )
            )
            == 1
    ):
      pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)[
        0
      ].close()
    self._interface_manager.job_manager.stop_auxiliary_pymol()
    self._interface_manager.app_process_manager.close_manager()
    tmp_event.accept()

  def _force_close_all(self) -> None:
    """Forcefully closes all open windows and processes related to PySSA."""
    tmp_number_of_help_windows = len(
      pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)
    )
    # PySSA should be closed
    # if not self._view.ui.lbl_logo.isVisible():
    #     logger.info("A project is currently opened. It will now be saved and the application exists afterwards.")
    #     self.__slot_close_project()
    # Help windows
    if tmp_number_of_help_windows == 1:
      logger.info("The documentation window is open. It will be closed now.")
      pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)[
        0
      ].close()
    elif tmp_number_of_help_windows > 1:
      for tmp_window_index in range(tmp_number_of_help_windows):
        pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)[
          tmp_window_index
        ].close()

    tmp_number_of_pymol_windows = len(
      pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYMOL_PART)
    )
    self._interface_manager.app_process_manager.close_manager()
    # PyMOL windows
    if tmp_number_of_pymol_windows == 1:
      logger.info("PyMOL will be closed now.")
      pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYMOL_PART)[
        0
      ].close()
    elif tmp_number_of_pymol_windows > 1:
      for tmp_window_index in range(tmp_number_of_help_windows + 1):
        pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYMOL_PART)[
          tmp_window_index
        ].close()

    tmp_number_of_pyssa_windows = len(
      pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYSSA)
    )
    tmp_number_of_exact_pyssa_match_windows = 0
    for tmp_window_index in range(tmp_number_of_pyssa_windows - 1):
      if (
              pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYSSA)[
                tmp_window_index
              ].title
              == "PySSA"
      ):
        tmp_number_of_exact_pyssa_match_windows += 1
    # PySSA windows
    if tmp_number_of_exact_pyssa_match_windows == 1:
      logger.info("PySSA will be closed now.")
      pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYSSA)[
        0
      ].close()
    elif tmp_number_of_exact_pyssa_match_windows > 1:
      for tmp_window_index in range(tmp_number_of_help_windows + 1):
        pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYSSA)[
          tmp_window_index
        ].close()

  def __slot_close_all(self) -> None:
    """Close all open slots.

    This method is triggered when the 'Project/Exit Application' menu entry is clicked.
    The method prompts the user with a confirmation dialog box to ensure they want to close the application.
    If the user confirms, the method initiates the process of closing the application by saving any open projects and terminating any running tasks.
    """
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Project/Exit Application' clicked.",
      )
      tmp_dialog = custom_message_box.CustomMessageBoxYesNo(
        "Are you sure you want to close PySSA?",
        "Close PySSA",
        custom_message_box.CustomMessageBoxIcons.WARNING.value,
      )
      tmp_dialog.exec_()
      if tmp_dialog.response:
        # PySSA should be closed
        if not self._view.ui.lbl_logo.isVisible():
          # self._active_task = tasks.LegacyTask(
          #     target=project_async.close_project,
          #     args=(
          #         self._database_thread,
          #         self._interface_manager.pymol_session_manager,
          #     ),
          #     post_func=self.__await_close_project_for_closing_app,
          # )

          self._interface_manager.get_task_manager().append_task_result(
            task_result_factory.TaskResultFactory.run_task_result(
              a_task_result=task_result.TaskResult.from_action(
                an_action=action.Action(
                  a_target=project_async.close_project,
                  args=(
                    self._database_thread,
                    self._interface_manager.pymol_session_manager,
                  ),
                ),
                an_await_function=self.__await_close_project_for_closing_app,
              ),
              a_task_scheduler=self._interface_manager.get_task_scheduler(),
            )
          )

          logger.info(
            "A project is currently opened. It will now be saved and the application exists afterwards."
          )
          # self._active_task.start()
          # self.__slot_close_project()
        else:
          self.__await_close_project_for_closing_app()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def __await_close_project_for_closing_app(self) -> None:
    """Closes the project and then the application."""
    # Help windows
    tmp_number_of_help_windows = len(
      pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)
    )
    if tmp_number_of_help_windows == 1:
      logger.info("The documentation window is open. It will be closed now.")
      pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)[
        0
      ].close()
    elif tmp_number_of_help_windows > 1:
      for tmp_window_index in range(tmp_number_of_help_windows):
        pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)[
          tmp_window_index
        ].close()
    else:
      logger.info("No documentation window is open. Nothing to do.")

    tmp_number_of_pymol_windows = len(
      pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYMOL_PART)
    )
    self._interface_manager.app_process_manager.close_manager()
    # PyMOL windows
    if tmp_number_of_pymol_windows == 1:
      logger.info("PyMOL will be closed now.")
      pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYMOL_PART)[
        0
      ].close()
    elif tmp_number_of_pymol_windows > 1:
      tmp_dialog = custom_message_box.CustomMessageBoxYesNo(
        "There are multiple windows open which contain PyMOL as window title.\nDo you want to close all?",
        "Close PySSA",
        custom_message_box.CustomMessageBoxIcons.WARNING.value,
      )
      tmp_dialog.exec_()
      if tmp_dialog.response:
        for tmp_window_index in range(tmp_number_of_pymol_windows + 1):
          pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYMOL_PART)[
            tmp_window_index
          ].close()

    tmp_number_of_pyssa_windows = len(
      pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYSSA)
    )
    tmp_number_of_exact_pyssa_match_windows = 0
    for tmp_window_index in range(tmp_number_of_pyssa_windows - 1):
      logger.info(
        pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYSSA)[
          tmp_window_index
        ].title
      )
      if (
              pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYSSA)[
                tmp_window_index
              ].title
              == "PySSA"
      ):
        tmp_number_of_exact_pyssa_match_windows += 1
    logger.info(tmp_number_of_exact_pyssa_match_windows)
    # PySSA windows
    if tmp_number_of_exact_pyssa_match_windows == 1:
      logger.info("PySSA will be closed now.")
      pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYSSA)[
        0
      ].close()
    elif tmp_number_of_exact_pyssa_match_windows > 1:
      tmp_dialog = custom_message_box.CustomMessageBoxYesNo(
        "There are multiple windows open which contain PySSA as window title.\nDo you want to close all?",
        "Close PySSA",
        custom_message_box.CustomMessageBoxIcons.WARNING.value,
      )
      tmp_dialog.exec_()
      if tmp_dialog.response:
        for tmp_window_index in range(tmp_number_of_pyssa_windows - 1):
          pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYSSA)[
            tmp_window_index
          ].close()
    else:
      logger.error("No PySSA window found.")

  def _abort_task(self, return_value: tuple[str, list[tuple[bool, tuple]]]) -> None:
    """Aborts a task.

    Args:
        return_value (tuple[str, list[tuple[bool, tuple]]]): A tuple containing two values. The first value indicates whether the task should be aborted (True) or not (False). The second value is a string representing the task.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      return

    # </editor-fold>

    tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
    if tmp_result[0] is True and tmp_result[1] == "ColabFold Prediction":
      self.__slot_abort_prediction()

  def _update_main_view_ui(
          self, refresh_after_job_finished_signal_values: tuple
  ) -> None:
    """Updates the main view UI based on the completed job.

    Args:
        refresh_after_job_finished_signal_values (tuple): A tuple containing the following information:
            - tmp_job_is_for_current_project_flag: A boolean indicating whether the job is for the current project or not.
            - tmp_job_base_information: An object containing information about the job.
            - tmp_job_notification_widget: The notification widget for the job.
    """
    (
      tmp_job_is_for_current_project_flag,
      tmp_job_base_information,
      tmp_job_notification_widget,
    ) = refresh_after_job_finished_signal_values
    self._interface_manager.remove_job_notification_widget(
      tmp_job_notification_widget
    )
    if tmp_job_base_information.job_progress == enums.JobProgress.FAILED:
      return
    if tmp_job_is_for_current_project_flag:
      if tmp_job_base_information.job_type == enums.JobType.PREDICTION:
        # refresh protein model
        # self._active_task = tasks.LegacyTask(
        #     target=util_async.add_proteins_to_project_and_model,
        #     args=(
        #         self._interface_manager,
        #         tmp_job_base_information.protein_names,
        #     ),
        #     post_func=self._post_update_project_and_model,
        # )

        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=util_async.add_proteins_to_project_and_model,
                args=(
                  self._interface_manager,
                  tmp_job_base_information.protein_names,
                ),
              ),
              an_await_function=self._post_update_project_and_model,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task.start()
      elif tmp_job_base_information.job_type == enums.JobType.DISTANCE_ANALYSIS:
        # refresh protein pair model
        # self._active_task = tasks.LegacyTask(
        #     target=util_async.add_protein_pairs_to_project_and_model,
        #     args=(
        #         self._interface_manager,
        #         tmp_job_base_information.protein_pair_names,
        #     ),
        #     post_func=self._post_update_project_and_model,
        # )

        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=util_async.add_protein_pairs_to_project_and_model,
                args=(
                  self._interface_manager,
                  tmp_job_base_information.protein_pair_names,
                ),
              ),
              an_await_function=self._post_update_project_and_model,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task.start()
      elif (
              tmp_job_base_information.job_type
              == enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS
      ):
        # refresh protein and protein pair model
        # self._active_task = tasks.LegacyTask(
        #     target=util_async.add_proteins_and_protein_pairs_to_project_and_model,
        #     args=(
        #         self._interface_manager,
        #         tmp_job_base_information.protein_names,
        #         tmp_job_base_information.protein_pair_names,
        #     ),
        #     post_func=self._post_update_project_and_model,
        # )

        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=util_async.add_proteins_and_protein_pairs_to_project_and_model,
                args=(
                  self._interface_manager,
                  tmp_job_base_information.protein_names,
                  tmp_job_base_information.protein_pair_names,
                ),
              ),
              an_await_function=self._post_update_project_and_model,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task.start()
      elif tmp_job_base_information.job_type == enums.JobType.RAY_TRACING:
        if not os.path.exists(tmp_job_base_information.get_image_filepath()):
          self._interface_manager.status_bar_manager.show_error_message(
            "Could not find image file!"
          )
        os.startfile(tmp_job_base_information.get_image_filepath())
        return
    else:
      # open other project
      tmp_a_project_is_open = self._view.ui.lbl_project_name.isVisible()
      if tmp_a_project_is_open:
        self._disconnect_sequence_selection_model()
        # self._active_task = tasks.LegacyTask(
        #     target=util_async.close_project_automatically,
        #     args=(
        #         tmp_a_project_is_open,
        #         self._database_thread,
        #         self._interface_manager.pymol_session_manager,
        #         tmp_job_base_information.project_name,
        #     ),
        #     post_func=self._post_close_project_automatically,
        # )

        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=util_async.close_project_automatically,
                args=(
                  tmp_a_project_is_open,
                  self._database_thread,
                  self._interface_manager.pymol_session_manager,
                  tmp_job_base_information.project_name,
                ),
              ),
              an_await_function=self._post_update_project_and_model,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
      # self._active_task.start()
      self.update_status("Saving current project ...")
      self._interface_manager.restore_default_main_view()
      self._interface_manager.close_job_notification_panel()
    self._interface_manager.block_gui(with_wait_cursor=True)

  def _post_update_project_and_model(self, return_value: tuple[str, list[tuple[bool, tuple]]]) -> None:
    """Refreshes the main view and reverts the cursor.

    Args:
        return_value (tuple[str, list[tuple[bool, tuple]]]): A tuple containing the return value of the update project and model method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      return

    # </editor-fold>

    tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
    self._interface_manager.refresh_main_view()
    self._interface_manager.stop_wait_cursor()

  def _post_close_project_automatically(self, return_value: tuple) -> None:
    """Tries to close the project automatically without user intervention.

    Args:
        return_value (tuple[str, list[tuple[bool, tuple]]]): The return value from the previous method call.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      return

    # </editor-fold>

    try:
      self._interface_manager.documentation_window = return_value[2]
      tmp_project_name = return_value[0]
      self._interface_manager.status_bar_manager.show_temporary_message(
        enums.StatusMessages.OPENING_PROJECT.value,
        False,
      )
      tmp_project_database_filepath = str(
        pathlib.Path(
          f"{self._interface_manager.get_application_settings().workspace_path}/{tmp_project_name}.db",
        ),
      )
      self._database_thread = database_thread.DatabaseThread(
        tmp_project_database_filepath
      )
      self._database_manager.set_database_filepath(
        tmp_project_database_filepath
      )
      # self._active_task = tasks.LegacyTask(
      #     target=project_async.open_project,
      #     args=(
      #         tmp_project_name,
      #         tmp_project_database_filepath,
      #         self._interface_manager,
      #         self._interface_manager.pymol_session_manager,
      #         self.custom_progress_signal,
      #         self._interface_manager.watcher,
      #     ),
      #     post_func=self.__await_open_project,
      # )
      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=project_async.open_project,
              args=(
                tmp_project_name,
                tmp_project_database_filepath,
                self._interface_manager,
                self._interface_manager.pymol_session_manager,
                self.custom_progress_signal,
                self._interface_manager.watcher,
              ),
            ),
            an_await_function=self.__await_open_help,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )

    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def _update_progress_bar(self, return_value: tuple) -> None:
    """Updates the progress bar.

    Args:
        return_value (tuple[str, list[tuple[bool, tuple]]]): The return value to be passed to the progress bar.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      return

    # </editor-fold>

    self._interface_manager.status_bar_manager.update_progress_bar(return_value)

  def __slot_restart_pymol(self) -> None:
    """Restarts PyMOL.

    This method restarts PyMOL by closing the PyMOL window and displaying a temporary message on the status bar indicating that PyMOL is restarting.
    """
    logger.warning("Recived user request to restart PyMOL.")
    pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_PYMOL_PART)[
      0
    ].close()
    self._interface_manager.status_bar_manager.show_temporary_message(
      "Restarting PyMOL ...",
      a_with_timeout_flag=False,
    )

  # <editor-fold desc="Util methods">
  def update_status(self, message: str) -> None:
    """Updates the status bar of the main view with a custom message.

    Args:
        message: The message to display in the statusbar.

    Raises:
        exception.IllegalArgumentError: If `message` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if message is None or message == "":
      logger.error("message is either None or an empty string.")
      raise exception.IllegalArgumentError(
        "message is either None or an empty string."
      )

    # </editor-fold>

    self._interface_manager.status_bar_manager.show_temporary_message(message)

  def _update_tab(self) -> None:
    """Update the current tab based on the current state of the interface manager."""
    self._interface_manager.current_tab_index = (
      self._view.ui.project_tab_widget.currentIndex()
    )
    if (
            self._interface_manager.pymol_session_manager.session_object_type
            == "protein"
            and self._interface_manager.current_tab_index == 2
    ):
      self._interface_manager.hide_protein_pair_pymol_scene_configuration()
      self._view.ui.lbl_info_3.setText(
        "Please load the PyMOL session of the \nselected protein pair."
      )
    elif (
            self._interface_manager.pymol_session_manager.session_object_type
            == "protein_pair"
            and self._interface_manager.current_tab_index == 1
    ):
      self._interface_manager.hide_protein_pymol_scene_configuration()
      self._view.ui.lbl_info.setText(
        "Please load the PyMOL session of the selected protein."
      )
    elif (
            self._interface_manager.pymol_session_manager.is_the_current_session_empty()
    ):
      self._interface_manager.hide_protein_pymol_scene_configuration()
      self._interface_manager.hide_protein_pair_pymol_scene_configuration()
      self._view.ui.lbl_info.setText(
        "Please load the PyMOL session of the selected protein."
      )
      self._view.ui.lbl_info_3.setText(
        "Please load the PyMOL session of the \nselected protein pair."
      )

  def _setup_statusbar(self) -> None:
    """Sets up the status bar and fills it with the current workspace."""
    self._interface_manager.get_main_view().setStatusBar(
      self._interface_manager.get_main_view().status_bar
    )

  def __slot_refresh_main_view(self) -> None:
    """Refreshes the main view."""
    self._interface_manager.refresh_main_view()

  # <editor-fold desc="Help related methods">
  def _start_documentation_server(self) -> None:
    """Starts the documentation server."""
    # self._help_task = tasks.LegacyTask(
    #     target=util_async.start_documentation_server,
    #     args=(0, 0),
    #     post_func=self.__await_start_documentation_server,
    # )

    self._interface_manager.get_task_manager().append_task_result(
      task_result_factory.TaskResultFactory.run_task_result(
        a_task_result=task_result.TaskResult.from_action(
          an_action=action.Action(
            a_target=util_async.start_documentation_server,
            args=(
              (0, 0)
            ),
          ),
          an_await_function=self.__await_start_documentation_server,
        ),
        a_task_scheduler=self._interface_manager.get_task_scheduler(),
      )
    )

    # self._help_task.start()

  def __await_start_documentation_server(self, return_value: tuple[str, list[tuple[bool, tuple]]]) -> None:
    """Checks if the documentation server started correctly.

    Args:
        return_value (tuple): A tuple containing information about the start of the documentation server.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "Opening help center failed!"
      )
      return
    if return_value[0] == "":
      self._interface_manager.status_bar_manager.show_error_message(
        "Opening help center failed!"
      )
      return

    # </editor-fold>

    tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
    self._interface_manager.documentation_window = tmp_result[2]
    # self._interface_manager.documentation_window = return_value[1]
    self._interface_manager.status_bar_manager.show_temporary_message(
      "Opening help center finished."
    )

  def open_help(self, a_page_name: str) -> None:
    """Opens the pyssa documentation window if it's not already open.

    Args:
        a_page_name (str): a name of a documentation page to display

    Raises:
        exception.IllegalArgumentError: If `a_page_name` is None.
    """
    # <editor-fold desc="Checks">
    if a_page_name is None:
      logger.error("a_page_name is None.")
      raise exception.IllegalArgumentError("a_page_name is None.")

    # </editor-fold>
    try:
      self._interface_manager.status_bar_manager.show_temporary_message(
        "Opening help center ...", False
      )

      if (
              len(
                pygetwindow.getWindowsWithTitle(
                  constants.WINDOW_TITLE_OF_HELP_CENTER
                )
              )
              != 1
      ):
        self._interface_manager.documentation_window = None
      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=util_async.open_documentation_on_certain_page,
              args=(
                a_page_name,
                self._interface_manager.documentation_window,
              ),
            ),
            an_await_function=self.__await_open_help,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )

      # self._help_task = tasks.LegacyTask(
      #     target=util_async.open_documentation_on_certain_page,
      #     args=(a_page_name, self._interface_manager.documentation_window),
      #     post_func=self.__await_open_help,
      # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui(with_wait_cursor=True)
      # self._help_task.start()

  def __await_open_help(self, return_value: tuple[str, list[tuple[bool, tuple]]]) -> None:
    """Opens the help center and performs necessary actions based on the return value.

    Args:
        return_value (tuple[str, list[tuple[bool, tuple]]]): The return value from opening the help center.
    """
    try:
      logger.debug("Method __await_open_help started")
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      self._interface_manager.documentation_window = tmp_result[2]
      if not os.path.exists(constants.HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH):
        tmp_dialog = custom_message_box.CustomMessageBoxOk(
          "The script for bringing the documentation window in front could not be found!",
          "Documentation",
          custom_message_box.CustomMessageBoxIcons.ERROR.value,
        )
        tmp_dialog.exec_()
      else:
        self._interface_manager.documentation_window.restore()
        subprocess.run([constants.HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH])
        self._interface_manager.status_bar_manager.show_temporary_message(
          "Opening help center finished."
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  def _init_generic_help_context_menus(self) -> None:
    """Initializes the generic help context menus for the application."""
    # <editor-fold desc="General context menu setup">
    context_menu = QtWidgets.QMenu()
    self.help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    self.help_context_action.triggered.connect(
      self.__slot_open_sequences_tab_help
    )

    # </editor-fold>

    # Set the context menu for the buttons
    self._view.ui.seqs_table_widget.setContextMenuPolicy(3)
    self._view.ui.seqs_table_widget.customContextMenuRequested.connect(
      self._show_context_menu_for_seq_table
    )
    self._view.ui.btn_import_seq.setContextMenuPolicy(
      3
    )  # 3 corresponds to Qt.CustomContextMenu
    self._view.ui.btn_import_seq.customContextMenuRequested.connect(
      self._show_context_menu_for_seq_import
    )
    self._view.ui.btn_add_sequence.setContextMenuPolicy(
      3
    )  # 3 corresponds to Qt.CustomContextMenu
    self._view.ui.btn_add_sequence.customContextMenuRequested.connect(
      self._show_context_menu_for_seq_add
    )
    self._view.ui.btn_save_sequence.setContextMenuPolicy(
      3
    )  # 3 corresponds to Qt.CustomContextMenu
    self._view.ui.btn_save_sequence.customContextMenuRequested.connect(
      self._show_context_menu_for_seq_save
    )
    self._view.ui.btn_delete_sequence.setContextMenuPolicy(
      3
    )  # 3 corresponds to Qt.CustomContextMenu
    self._view.ui.btn_delete_sequence.customContextMenuRequested.connect(
      self._show_context_menu_for_seq_delete
    )
    self._view.ui.btn_import_protein.setContextMenuPolicy(3)
    self._view.ui.btn_import_protein.customContextMenuRequested.connect(
      self._show_context_menu_for_protein_import,
    )
    self._view.ui.btn_save_protein.setContextMenuPolicy(3)
    self._view.ui.btn_save_protein.customContextMenuRequested.connect(
      self._show_context_menu_for_protein_save,
    )
    self._view.ui.btn_delete_protein.setContextMenuPolicy(3)
    self._view.ui.btn_delete_protein.customContextMenuRequested.connect(
      self._show_context_menu_for_protein_delete,
    )
    self._view.ui.btn_open_protein_session.setContextMenuPolicy(3)
    self._view.ui.btn_open_protein_session.customContextMenuRequested.connect(
      self._show_context_menu_for_protein_load_session,
    )
    self._view.ui.btn_create_protein_scene.setContextMenuPolicy(3)
    self._view.ui.btn_create_protein_scene.customContextMenuRequested.connect(
      self._show_context_menu_for_protein_add_scene,
    )
    self._view.ui.btn_update_protein_scene.setContextMenuPolicy(3)
    self._view.ui.btn_update_protein_scene.customContextMenuRequested.connect(
      self._show_context_menu_for_protein_update_scene,
    )
    self._view.ui.btn_delete_protein_scene.setContextMenuPolicy(3)
    self._view.ui.btn_delete_protein_scene.customContextMenuRequested.connect(
      self._show_context_menu_for_protein_delete_scene,
    )
    self._view.ui.frame_protein_pymol_scene.setContextMenuPolicy(3)
    self._view.ui.frame_protein_pymol_scene.customContextMenuRequested.connect(
      self._show_context_menu_for_protein_pymol_scene_config,
    )
    self._view.ui.btn_delete_protein_pair.setContextMenuPolicy(3)
    self._view.ui.btn_delete_protein_pair.customContextMenuRequested.connect(
      self._show_context_menu_for_protein_pair_delete,
    )
    self._view.ui.btn_open_protein_pair_session.setContextMenuPolicy(3)
    self._view.ui.btn_open_protein_pair_session.customContextMenuRequested.connect(
      self._show_context_menu_for_protein_pair_load_session,
    )
    self._view.ui.btn_create_protein_pair_scene.setContextMenuPolicy(3)
    self._view.ui.btn_create_protein_pair_scene.customContextMenuRequested.connect(
      self._show_context_menu_for_protein_pair_add_scene,
    )
    self._view.ui.btn_update_protein_pair_scene.setContextMenuPolicy(3)
    self._view.ui.btn_update_protein_pair_scene.customContextMenuRequested.connect(
      self._show_context_menu_for_protein_pair_update_scene,
    )
    self._view.ui.btn_delete_protein_pair_scene.setContextMenuPolicy(3)
    self._view.ui.btn_delete_protein_pair_scene.customContextMenuRequested.connect(
      self._show_context_menu_for_protein_pair_delete_scene,
    )
    self._view.ui.frame_protein_pair_pymol_scene.setContextMenuPolicy(3)
    self._view.ui.frame_protein_pair_pymol_scene.customContextMenuRequested.connect(
      self._show_context_menu_for_protein_pair_pymol_scene_config,
    )

  # <editor-fold desc="Help pages">
  def __slot_open_help_center(self) -> None:
    """Opens the help dialog on the homepage."""
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
      "Menu entry 'Help/Documentation' clicked.",
    )
    self.open_help("help/")

  def __slot_open_sequences_tab_help(self) -> None:
    """Opens the help dialog on the sequence tab page."""
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
      "'Help' button on the 'Sequence Tab' was clicked.",
    )
    self.open_help("help/sequences/sequences_tab/")

  def _open_additional_information_table_help(self) -> None:
    """Opens the help dialog on the additional sequence page."""
    self.open_help("help/sequences/additional_sequence_information/")

  def _open_sequence_import_help(self) -> None:
    """Opens the help dialog on the import sequence page."""
    self.open_help("help/sequences/sequence_import/")

  def _open_sequence_add_help(self) -> None:
    """Opens the help dialog on the add sequence page."""
    self.open_help("help/sequences/sequence_add/")

  def _open_sequence_save_help(self) -> None:
    """Opens the help dialog on the save sequence page."""
    self.open_help("help/sequences/sequence_save/")

  def _open_sequence_delete_help(self) -> None:
    """Opens the help dialog on the delete sequence page."""
    self.open_help("help/sequences/sequence_delete/")

  def __slot_open_proteins_tab_help(self) -> None:
    """Opens the help dialog on the proteins tab page."""
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
      "'Help' button on the 'Proteins Tab' was clicked.",
    )
    self.open_help("help/proteins/proteins_tab/")

  def _open_protein_import_help(self) -> None:
    """Opens the help dialog on the import protein page."""
    self.open_help("help/proteins/protein_import/")

  def _open_protein_save_help(self) -> None:
    """Opens the help dialog on the save protein page."""
    self.open_help("help/proteins/protein_save/")

  def _open_protein_delete_help(self) -> None:
    """Opens the help dialog on the delete protein page."""
    self.open_help("help/proteins/protein_delete/")

  def _open_protein_pymol_scene_config_help(self) -> None:
    """Opens the help dialog on the protein pymol scene config page."""
    self.open_help("help/proteins/protein_pymol_scene_configuration/")

  def _open_protein_load_session_help(self) -> None:
    """Opens the help dialog on the open protein pymol session page."""
    self.open_help("help/proteins/protein_load_session/")

  def _open_protein_add_scene_help(self) -> None:
    """Opens the help dialog on the protein add scene page."""
    self.open_help("help/proteins/protein_add_scene/")

  def _open_protein_update_scene_help(self) -> None:
    """Opens the help dialog on the protein update scene page."""
    self.open_help("help/proteins/protein_update_scene/")

  def _open_protein_delete_scene_help(self) -> None:
    """Opens the help dialog on the protein delete scene page."""
    self.open_help("help/proteins/protein_delete_scene/")

  def _open_protein_pairs_tab_help(self) -> None:
    """Opens the help dialog on the protein pairs tab page."""
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
      "'Help' button on the 'Protein Pairs Tab' was clicked.",
    )
    self.open_help("help/protein_pairs/protein_pairs_tab/")

  def _open_protein_pair_delete_help(self) -> None:
    """Opens the help dialog on the delete protein pair page."""
    self.open_help("help/protein_pairs/protein_pair_delete/")

  def _open_protein_pair_pymol_scene_config_help(self) -> None:
    """Opens the help dialog on the protein pair pymol scene config page."""
    self.open_help("help/protein_pairs/protein_pair_pymol_scene_configuration/")

  def _open_protein_pair_load_session_help(self) -> None:
    """Opens the help dialog on the open protein pair pymol session page."""
    self.open_help("help/protein_pairs/protein_pair_load_session/")

  def _open_protein_pair_add_scene_help(self) -> None:
    """Opens the help dialog on the protein pair add scene page."""
    self.open_help("help/protein_pairs/protein_pair_add_scene/")

  def _open_protein_pair_update_scene_help(self) -> None:
    """Opens the help dialog on the protein pair update scene page."""
    self.open_help("help/protein_pairs/protein_pair_update_scene/")

  def _open_protein_pair_delete_scene_help(self) -> None:
    """Opens the help dialog on the protein pair delete scene page."""
    self.open_help("help/protein_pairs/protein_pair_delete_scene/")

  # </editor-fold>

  # <editor-fold desc="Context menu connections">
  def _show_context_menu_for_seq_table(self, a_point: QtCore.QPoint) -> None:
    """Shows a context menu on the specific point.

    Args:
        a_point: a QPoint object representing the coordinates of the point where the context menu is to be shown
    """
    context_menu = QtWidgets.QMenu()
    help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    help_context_action.triggered.connect(
      self._open_additional_information_table_help
    )
    context_menu.exec_(self._view.ui.seqs_table_widget.mapToGlobal(a_point))

  def _show_context_menu_for_seq_import(self, a_point: QtCore.QPoint) -> None:
    """Shows a context menu on the specific point.

    Args:
        a_point: a QPoint object representing the coordinates of the point where the context menu is to be shown
    """
    context_menu = QtWidgets.QMenu()
    help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    help_context_action.triggered.connect(self._open_sequence_import_help)
    context_menu.exec_(self._view.ui.btn_import_seq.mapToGlobal(a_point))

  def _show_context_menu_for_seq_add(self, a_point: QtCore.QPoint) -> None:
    """Shows a context menu on the specific point.

    Args:
        a_point: a QPoint object representing the coordinates of the point where the context menu is to be shown
    """
    context_menu = QtWidgets.QMenu()
    help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    help_context_action.triggered.connect(self._open_sequence_add_help)
    context_menu.exec_(self._view.ui.btn_add_sequence.mapToGlobal(a_point))

  def _show_context_menu_for_seq_save(self, a_point: QtCore.QPoint) -> None:
    """Shows a context menu on the specific point.

    Args:
        a_point: a QPoint object representing the coordinates of the point where the context menu is to be shown
    """
    context_menu = QtWidgets.QMenu()
    help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    help_context_action.triggered.connect(self._open_sequence_save_help)
    context_menu.exec_(self._view.ui.btn_save_sequence.mapToGlobal(a_point))

  def _show_context_menu_for_seq_delete(self, a_point: QtCore.QPoint) -> None:
    """Shows a context menu on the specific point.

    Args:
        a_point: a QPoint object representing the coordinates of the point where the context menu is to be shown
    """
    context_menu = QtWidgets.QMenu()
    help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    help_context_action.triggered.connect(self._open_sequence_delete_help)
    context_menu.exec_(self._view.ui.btn_delete_sequence.mapToGlobal(a_point))

  def _show_context_menu_for_protein_import(
          self, a_point: QtCore.QPoint
  ) -> None:
    """Shows a context menu on the specific point.

    Args:
        a_point: a QPoint object representing the coordinates of the point where the context menu is to be shown
    """
    context_menu = QtWidgets.QMenu()
    help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    help_context_action.triggered.connect(self._open_protein_import_help)
    context_menu.exec_(self._view.ui.btn_import_protein.mapToGlobal(a_point))

  def _show_context_menu_for_protein_save(self, a_point: QtCore.QPoint) -> None:
    """Shows a context menu on the specific point.

    Args:
        a_point: a QPoint object representing the coordinates of the point where the context menu is to be shown
    """
    context_menu = QtWidgets.QMenu()
    help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    help_context_action.triggered.connect(self._open_protein_save_help)
    context_menu.exec_(self._view.ui.btn_save_protein.mapToGlobal(a_point))

  def _show_context_menu_for_protein_delete(
          self, a_point: QtCore.QPoint
  ) -> None:
    """Shows a context menu on the specific point.

    Args:
        a_point: a QPoint object representing the coordinates of the point where the context menu is to be shown
    """
    context_menu = QtWidgets.QMenu()
    help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    help_context_action.triggered.connect(self._open_protein_delete_help)
    context_menu.exec_(self._view.ui.btn_delete_protein.mapToGlobal(a_point))

  def _show_context_menu_for_protein_load_session(
          self, a_point: QtCore.QPoint
  ) -> None:
    """Shows a context menu on the specific point.

    Args:
        a_point: a QPoint object representing the coordinates of the point where the context menu is to be shown
    """
    context_menu = QtWidgets.QMenu()
    help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    help_context_action.triggered.connect(self._open_protein_load_session_help)
    context_menu.exec_(
      self._view.ui.btn_open_protein_session.mapToGlobal(a_point)
    )

  def _show_context_menu_for_protein_add_scene(
          self, a_point: QtCore.QPoint
  ) -> None:
    """Shows a context menu on the specific point.

    Args:
        a_point: a QPoint object representing the coordinates of the point where the context menu is to be shown
    """
    context_menu = QtWidgets.QMenu()
    help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    help_context_action.triggered.connect(self._open_protein_add_scene_help)
    context_menu.exec_(
      self._view.ui.btn_create_protein_scene.mapToGlobal(a_point)
    )

  def _show_context_menu_for_protein_update_scene(
          self, a_point: QtCore.QPoint
  ) -> None:
    """Shows a context menu on the specific point.

    Args:
        a_point: a QPoint object representing the coordinates of the point where the context menu is to be shown
    """
    context_menu = QtWidgets.QMenu()
    help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    help_context_action.triggered.connect(self._open_protein_update_scene_help)
    context_menu.exec_(
      self._view.ui.btn_update_protein_scene.mapToGlobal(a_point)
    )

  def _show_context_menu_for_protein_delete_scene(
          self, a_point: QtCore.QPoint
  ) -> None:
    """Shows a context menu on the specific point.

    Args:
        a_point: a QPoint object representing the coordinates of the point where the context menu is to be shown
    """
    context_menu = QtWidgets.QMenu()
    help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    help_context_action.triggered.connect(self._open_protein_delete_scene_help)
    context_menu.exec_(
      self._view.ui.btn_delete_protein_scene.mapToGlobal(a_point)
    )

  def _show_context_menu_for_protein_pymol_scene_config(
          self, a_point: QtCore.QPoint
  ) -> None:
    """Shows a context menu on the specific point.

    Args:
        a_point: a QPoint object representing the coordinates of the point where the context menu is to be shown
    """
    context_menu = QtWidgets.QMenu()
    help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    help_context_action.triggered.connect(
      self._open_protein_pymol_scene_config_help
    )
    context_menu.exec_(
      self._view.ui.frame_protein_pymol_scene.mapToGlobal(a_point)
    )

  def _show_context_menu_for_protein_pair_delete(
          self, a_point: QtCore.QPoint
  ) -> None:
    """Shows a context menu on the specific point.

    Args:
        a_point: a QPoint object representing the coordinates of the point where the context menu is to be shown
    """
    context_menu = QtWidgets.QMenu()
    help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    help_context_action.triggered.connect(self._open_protein_pair_delete_help)
    context_menu.exec_(
      self._view.ui.btn_delete_protein_pair.mapToGlobal(a_point)
    )

  def _show_context_menu_for_protein_pair_load_session(
          self, a_point: QtCore.QPoint
  ) -> None:
    """Shows a context menu on the specific point.

    Args:
        a_point: a QPoint object representing the coordinates of the point where the context menu is to be shown
    """
    context_menu = QtWidgets.QMenu()
    help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    help_context_action.triggered.connect(
      self._open_protein_pair_load_session_help
    )
    context_menu.exec_(
      self._view.ui.btn_open_protein_pair_session.mapToGlobal(a_point)
    )

  def _show_context_menu_for_protein_pair_add_scene(
          self, a_point: QtCore.QPoint
  ) -> None:
    """Shows a context menu on the specific point.

    Args:
        a_point: a QPoint object representing the coordinates of the point where the context menu is to be shown
    """
    context_menu = QtWidgets.QMenu()
    help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    help_context_action.triggered.connect(
      self._open_protein_pair_add_scene_help
    )
    context_menu.exec_(
      self._view.ui.btn_create_protein_pair_scene.mapToGlobal(a_point)
    )

  def _show_context_menu_for_protein_pair_update_scene(
          self, a_point: QtCore.QPoint
  ) -> None:
    """Shows a context menu on the specific point.

    Args:
        a_point: a QPoint object representing the coordinates of the point where the context menu is to be shown
    """
    context_menu = QtWidgets.QMenu()
    help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    help_context_action.triggered.connect(
      self._open_protein_pair_update_scene_help
    )
    context_menu.exec_(
      self._view.ui.btn_update_protein_pair_scene.mapToGlobal(a_point)
    )

  def _show_context_menu_for_protein_pair_delete_scene(
          self, a_point: QtCore.QPoint
  ) -> None:
    """Shows a context menu on the specific point.

    Args:
        a_point: a QPoint object representing the coordinates of the point where the context menu is to be shown
    """
    context_menu = QtWidgets.QMenu()
    help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    help_context_action.triggered.connect(
      self._open_protein_pair_delete_scene_help
    )
    context_menu.exec_(
      self._view.ui.btn_delete_protein_pair_scene.mapToGlobal(a_point)
    )

  def _show_context_menu_for_protein_pair_pymol_scene_config(
          self, a_point: QtCore.QPoint
  ) -> None:
    """Shows a context menu on the specific point.

    Args:
        a_point: a QPoint object representing the coordinates of the point where the context menu is to be shown
    """
    context_menu = QtWidgets.QMenu()
    help_context_action = context_menu.addAction(self._view.tr("Get Help"))
    help_context_action.triggered.connect(
      self._open_protein_pair_pymol_scene_config_help
    )
    context_menu.exec_(
      self._view.ui.frame_protein_pair_pymol_scene.mapToGlobal(a_point)
    )

  # </editor-fold>

  # </editor-fold>

  # </editor-fold>

  # <editor-fold desc="Job panels">
  def __slot_open_job_overview_panel(self) -> None:
    """Hides or shows the job overview panel based on its current state."""
    if self._view.ui.frame_job_overview.isVisible():
      self._view.ui.frame_job_overview.hide()
    elif not self._view.ui.frame_job_overview.isVisible():
      self._view.ui.frame_job_notification.hide()
      self._view.ui.frame_job_overview.show()

  def __slot_open_notification_panel(self) -> None:
    """Hides or shows the job notification panel based on its current state."""
    if self._view.ui.frame_job_notification.isVisible():
      self._view.ui.frame_job_notification.hide()
      self._view.btn_open_job_notification.setIcon(self._view.icon_notify)
    elif not self._view.ui.frame_job_notification.isVisible():
      self._view.ui.frame_job_overview.hide()
      self._view.ui.frame_job_notification.show()
      self._view.btn_open_job_notification.setIcon(self._view.icon_notify)

  # </editor-fold>

  # <editor-fold desc="Project menu">
  def __slot_close_project(self) -> None:
    """Closes the current project."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Project/Close' clicked.",
      )
      # self._active_task = tasks.LegacyTask(
      #     target=project_async.close_project,
      #     args=(
      #         self._database_thread,
      #         self._interface_manager.pymol_session_manager,
      #     ),
      #     post_func=self.__await_close_project,
      # )

      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=project_async.close_project,
              args=(
                self._database_thread,
                self._interface_manager.pymol_session_manager,
              ),
            ),
            an_await_function=self.__await_close_project,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )

      self._interface_manager.restore_default_main_view()
      self._interface_manager.close_job_notification_panel()
      self._interface_manager.close_job_overview_panel()
      self._disconnect_sequence_selection_model()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui(with_wait_cursor=True)
      self.update_status("Saving current project ...")
      # self._active_task.start()

  def __await_close_project(self, return_value: tuple[str, list[tuple[bool, tuple]]]) -> None:
    """Await the async closing process.

    Args:
        return_value (tuple[str, list[tuple[bool, tuple]]]): A tuple of the return value of the async closing process.
    """
    # <editor-fold desc="Checks">
    if return_value[0] == "":
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.status_bar_manager.show_error_message(
        "Closing the project failed!"
      )
      return
    # </editor-fold>

    try:
      self._interface_manager.set_new_project(project.Project())
      self._view.ui.project_tab_widget.setCurrentIndex(0)
      self.update_status("Closing project finished.")
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()

  def __slot_create_project(self) -> None:
    """Initializes the CreateProjectViewController and shows the create project dialog."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Project/Create' clicked.",
      )
      self._external_controller = (
        create_project_view_controller.CreateProjectViewController(
          self._interface_manager
        )
      )
      self._external_controller.user_input.connect(self._post_create_project)
      self._interface_manager.get_create_view().show()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def _post_create_project(self, user_input: tuple) -> None:
    """Creates a new project based on the user's input.

    Args:
        user_input: A tuple containing the user input for creating a new project. The tuple must have two elements:
                    - tmp_project_name: A string representing the name of the new project.
                    - tmp_protein_name: A string representing the name of the protein associated with the new project.

    Notes:
        This method is called after the user provides input to create a new project. It performs the following steps:
    """
    # <editor-fold desc="Checks">
    if user_input is None:
      logger.error("user_input is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      return

    # </editor-fold>

    try:
      self._interface_manager.block_gui(with_wait_cursor=True)
      tmp_project_name, tmp_protein_name = user_input
      tmp_project_database_filepath = str(
        pathlib.Path(
          f"{self._interface_manager.get_application_settings().workspace_path}/{tmp_project_name}.db"
        )
      )
      with database_manager.DatabaseManager(
              tmp_project_database_filepath
      ) as db_manager:
        db_manager.build_new_database()
      self._database_thread = database_thread.DatabaseThread(
        tmp_project_database_filepath
      )
      # self._database_thread.start()
      # self._active_task = tasks.LegacyTask(
      #     target=project_async.create_new_project,
      #     args=(
      #         tmp_project_name,
      #         self._interface_manager.get_application_settings().get_workspace_path(),
      #         self._interface_manager.watcher,
      #         self._interface_manager,
      #     ),
      #     post_func=self.__await_create_project,
      # )

      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=project_async.create_new_project,
              args=(
                tmp_project_name,
                self._interface_manager.get_application_settings().get_workspace_path(),
                self._interface_manager.watcher,
                self._interface_manager,
              ),
            ),
            an_await_function=self.__await_create_project,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def __await_create_project(self, return_value: tuple[str, list[tuple[bool, tuple]]]) -> None:
    """Awaits the async method that creates the new project.

    Args:
        return_value (tuple[str, list[tuple[bool, tuple]]]): A tuple containing the return value of the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      return
    if return_value[0] == "":
      self._interface_manager.status_bar_manager.show_error_message(
        "Creating the project failed!"
      )
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      return
    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      (
        _,
        tmp_project,
        self._interface_manager.watcher,
        self._interface_manager,
      ) = tmp_result
      self._interface_manager.set_new_project(tmp_project)
      self._interface_manager.refresh_workspace_model()
      self._interface_manager.pymol_session_manager.reinitialize_session()
      self._connect_sequence_selection_model()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  def __slot_open_project(self) -> None:
    """Opens the dialog to open a project."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Project/Open' clicked.",
      )
      self._external_controller = (
        open_project_view_controller.OpenProjectViewController(
          self._interface_manager
        )
      )
      self._external_controller.return_value.connect(self._post_open_project)
      self._interface_manager.get_open_view().show()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def _post_open_project(self, return_value: tuple) -> None:
    """Post method that gets executed after the open project dialog closes.

    Args:
        return_value (tuple[str, list[tuple[bool, tuple]]]): The return value from the method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      return
    if return_value[1] is False:
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      self._interface_manager.status_bar_manager.show_temporary_message(
        enums.StatusMessages.OPENING_PROJECT.value,
        False,
      )
      self._interface_manager.block_gui(with_wait_cursor=True)
      tmp_project_name = return_value
      tmp_project_database_filepath = str(
        pathlib.Path(
          f"{self._interface_manager.get_application_settings().workspace_path}/{tmp_project_name}.db",
        ),
      )
      self._database_thread = database_thread.DatabaseThread(
        tmp_project_database_filepath
      )
      # self._database_thread.start()
      self._database_manager.set_database_filepath(
        tmp_project_database_filepath
      )
      # self._active_task = tasks.LegacyTask(
      #     target=project_async.open_project,
      #     args=(
      #         tmp_project_name,
      #         tmp_project_database_filepath,
      #         self._interface_manager,
      #         self._interface_manager.pymol_session_manager,
      #         self.custom_progress_signal,
      #         self._interface_manager.watcher,
      #     ),
      #     post_func=self.__await_open_project,
      # )

      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=project_async.open_project,
              args=(
                tmp_project_name,
                tmp_project_database_filepath,
                self._interface_manager,
                self._interface_manager.pymol_session_manager,
                self.custom_progress_signal,
                self._interface_manager.watcher,
              ),
            ),
            an_await_function=self.__await_open_project,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def __await_open_project(self, return_value: tuple[str, list[tuple[bool, tuple]]]) -> None:
    """Finishes the project opening process.

    Args:
        return_value (tuple[str, list[tuple[bool, tuple]]]): The return value from the async method.
    """
    self._interface_manager.status_bar_manager.hide_progress_bar()
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      return
    if return_value[0] == "":
      self._interface_manager.status_bar_manager.show_error_message(
        enums.StatusMessages.OPENING_PROJECT_FAILED.value,
      )
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      _, tmp_project, tmp_interface_manager, tmp_watcher = tmp_result
      self._interface_manager = tmp_interface_manager
      self._interface_manager.watcher = tmp_watcher
      self._interface_manager.refresh_main_view()
      self._interface_manager.hide_progress_bar()
      self._interface_manager.status_bar_manager.show_temporary_message(
        enums.StatusMessages.OPENING_PROJECT_FINISHED.value,
      )
      self._connect_sequence_selection_model()
      # Expand all available proteins
      self.__slot_expand_all_proteins()
      # Expand all available protein pairs
      self.__slot_expand_all_protein_pairs()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()

  def __slot_use_project(self) -> None:
    """Opens the use project dialog."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Project/Use' clicked.",
      )
      self._external_controller = (
        use_project_view_controller.UseProjectViewController(
          self._interface_manager
        )
      )
      self._external_controller.user_input.connect(self._post_use_project)
      self._interface_manager.get_use_project_view().show()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def _post_use_project(self, user_input: tuple) -> None:
    """Starts the use project process.

    Args:
        user_input (tuple): The user inputs from the use project dialog.
    """
    # <editor-fold desc="Checks">
    if user_input is None:
      logger.error("user_input is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      return

    # </editor-fold>

    try:
      tmp_project_database_filepath = str(
        pathlib.Path(
          f"{self._interface_manager.get_application_settings().get_workspace_path()}/{user_input[0]}.db"
        )
      )
      with database_manager.DatabaseManager(
              tmp_project_database_filepath
      ) as db_manager:
        db_manager.build_new_database()

      # self._active_task = tasks.LegacyTask(
      #   target=project_async.create_use_project,
      #   args=(
      #     user_input[0],
      #     self._interface_manager.get_application_settings().get_workspace_path(),
      #     user_input[1],
      #     self._interface_manager.watcher,
      #     self._interface_manager,
      #   ),
      #   post_func=self.__await_use_project,
      # )

      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=project_async.create_use_project,
              args=(
                user_input[0],
                self._interface_manager.get_application_settings().get_workspace_path(),
                user_input[1],
                self._interface_manager.watcher,
                self._interface_manager,
              ),
            ),
            an_await_function=self.__await_use_project,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui(with_wait_cursor=True)
      # self._active_task.start()

  def __await_use_project(self, return_value: tuple[str, list[tuple[bool, tuple]]]) -> None:
    """Finishes the use project process.

    Args:
        return_value (tuple[str, list[tuple[bool, tuple]]]): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      return
    if return_value[0] == "":
      self._interface_manager.status_bar_manager.show_error_message(
        "Using the project failed!"
      )
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      (
        _,
        tmp_project,
        self._interface_manager.watcher,
        self._interface_manager,
      ) = tmp_result
      self._interface_manager.set_new_project(tmp_project)
      self._interface_manager.add_project_to_workspace_model(
        tmp_project.get_project_name()
      )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._connect_sequence_selection_model()
      self._interface_manager.status_bar_manager.show_temporary_message(
        "Use process finished."
      )
    finally:
      self._interface_manager.pymol_session_manager.reinitialize_session()
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()

  def __slot_delete_project(self) -> None:
    """Opens the delete project dialog."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Project/Delete' clicked.",
      )
      self._external_controller = (
        delete_project_view_controller.DeleteProjectViewController(
          self._interface_manager
        )
      )
      self._interface_manager.get_delete_view().show()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def _post_delete_project(self) -> None:
    """Refreshes the main view after the delete project dialog closed."""
    self._interface_manager.refresh_main_view()

  def __slot_import_project(self) -> None:
    """Imports a project into the current workspace."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Project/Import' clicked.",
      )
      file_dialog = QtWidgets.QFileDialog()
      desktop_path = QtCore.QStandardPaths.standardLocations(
        QtCore.QStandardPaths.DesktopLocation
      )[0]
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
          text=tmp_import_filepath.name.replace(".db", ""),
        )[0]
        # self._active_task = tasks.LegacyTask(
        #   target=project_async.import_project,
        #   args=(
        #     tmp_new_project_name,
        #     tmp_import_filepath,
        #     self._interface_manager,
        #   ),
        #   post_func=self._await__slot_import_project,
        # )

        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=project_async.import_project,
                args=(
                  tmp_new_project_name,
                  tmp_import_filepath,
                  self._interface_manager,
                ),
              ),
              an_await_function=self._await__slot_import_project,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui(with_wait_cursor=True)
      self._interface_manager.status_bar_manager.show_temporary_message(
        "Importing project ...", a_with_timeout_flag=False
      )
      # self._active_task.start()

  def _await__slot_import_project(self, return_value: tuple[str, list[tuple[bool, tuple]]]) -> None:
    """Finishes the import project process.

    Args:
        return_value (tuple[str, list[tuple[bool, tuple]]]): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      self._database_thread = database_thread.DatabaseThread(tmp_result[1])
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      self._interface_manager.status_bar_manager.show_temporary_message(
        "Importing project finished."
      )

  def __slot_export_current_project(self) -> None:
    """Exports the current project to an importable format."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Project/Export' clicked.",
      )
      file_dialog = QtWidgets.QFileDialog()
      desktop_path = QtCore.QStandardPaths.standardLocations(
        QtCore.QStandardPaths.DesktopLocation
      )[0]
      file_dialog.setDirectory(desktop_path)
      file_path, _ = file_dialog.getSaveFileName(
        self._view,
        "Export current project",
        "",
        "Project Database File (*.db)",
      )
      if file_path:
        shutil.copyfile(
          self._interface_manager.get_current_project().get_database_filepath(),
          file_path,
        )
        # tmp_dialog = custom_message_box.CustomMessageBoxOk(
        #     "The project was successfully exported.", "Export Project",
        #     custom_message_box.CustomMessageBoxIcons.INFORMATION.value
        # )
        # tmp_dialog.exec_()
        self._interface_manager.status_bar_manager.show_temporary_message(
          "The project was successfully exported."
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  # </editor-fold>

  # <editor-fold desc="Analysis menu">
  def __slot_distance_analysis(self) -> None:
    """Opens the distance analysis dialog."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Analysis/Distance' clicked.",
      )
      self._external_controller = (
        distance_analysis_view_controller.DistanceAnalysisViewController(
          self._interface_manager,
          self._interface_manager.watcher,
        )
      )
      self._external_controller.job_input.connect(self._post_distance_analysis)
      self._interface_manager.get_distance_analysis_view().show()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def _post_distance_analysis(self, job_input: tuple) -> None:
    """Sets up the worker for the analysis task.

    Args:
        job_input (tuple): A tuple with the data needed for the distance analysis.
    """
    # <editor-fold desc="Checks">
    if job_input is None:
      logger.error("job_input is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      return

    # </editor-fold>

    try:
      _, tmp_raw_analysis_run_names, tmp_checkbox_state = job_input
      # --- New job approach
      self._interface_manager.watcher.add_protein_pairs_from_new_job(
        tmp_raw_analysis_run_names
      )
      tmp_distance_analysis_job, tmp_distance_analysis_entry_widget = (
        self._interface_manager.job_manager.create_distance_analysis_job(
          self._interface_manager.get_current_project(),
          self._interface_manager.project_lock,
          self._interface_manager,
          tmp_raw_analysis_run_names,
          self._interface_manager.get_settings_manager().settings.cutoff,
          self._interface_manager.get_settings_manager().settings.cycles,
        )
      )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.job_manager.put_job_into_queue(
        tmp_distance_analysis_job
      )
      self._interface_manager.add_job_entry_to_job_overview_layout(
        tmp_distance_analysis_entry_widget
      )

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

  # def _add_new_protein_pairs_to_protein_pair_model(self):
  #     """Adds the new protein pairs to the interface manager's protein pair model."""
  #     tmp_protein_pairs_to_add = self._main_view_state.get_not_matching_protein_pairs(
  #         self._interface_manager.get_current_project().protein_pairs,
  #     )
  #     for tmp_protein_pair in tmp_protein_pairs_to_add:
  #         self._interface_manager.add_protein_pair_to_protein_pairs_model(tmp_protein_pair)

  # def __await_unfreeze_pymol_session_after_analysis(self):
  #     self._interface_manager.main_tasks_manager.distance_analysis_task = None
  #     self._interface_manager.refresh_main_view()
  #     # self.active_custom_message_box = custom_message_box.CustomMessageBoxOk(
  #     #     "All structure analysis' are done. \nGo to the Protein Pairs tab to view the new results.",
  #     #     "Distance Analysis",
  #     #     custom_message_box.CustomMessageBoxIcons.INFORMATION.value
  #     # )
  #     # self.active_custom_message_box.exec_()
  #     constants.PYSSA_LOGGER.info("All structure analysis' are done.")
  #     self._interface_manager.status_bar_manager.show_temporary_message("All structure analysis' are done.")

  # </editor-fold>

  # <editor-fold desc="Prediction menu">
  def _connect_sequence_selection_model(self) -> None:
    """Connects the selection model of the sequence list view to the _check_options_for_sequence_selection method."""
    try:
      self._view.ui.seqs_list_view.selectionModel().selectionChanged.connect(
        self._check_options_for_sequence_selection
      )
    except Exception as e:
      logger.warning(
        f"The self._view.ui.seqs_list_view.selectionModel() is None. Error message: {e}"
      )

  def _disconnect_sequence_selection_model(self) -> None:
    """Disconnects the selection model of the sequence list view to the _check_options_for_sequence_selection method."""
    try:
      self._view.ui.seqs_list_view.selectionModel().selectionChanged.disconnect(
        self._check_options_for_sequence_selection
      )
    except TypeError:
      logger.warning(
        "Catching error because no sequences exist in the project. Therefore the selectionChanged signal does not need to be disconnected."
      )

  def _check_options_for_sequence_selection(self) -> None:
    """Checks prediction options for the sequence selection."""
    if len(self._view.ui.seqs_list_view.selectedIndexes()) > 0:
      tmp_enable_monomer_flag = False
      tmp_enable_multimer_flag = False
      for tmp_model_index in self._view.ui.seqs_list_view.selectedIndexes():
        tmp_type = tmp_model_index.data(enums.ModelEnum.TYPE_ROLE)
        tmp_sequence_name: SeqRecord.SeqRecord = tmp_model_index.data(
          enums.ModelEnum.OBJECT_ROLE
        ).name
        if (
                tmp_type == enums.ModelTypeEnum.MONOMER_SEQ
                and tmp_enable_monomer_flag is False
                and not self._interface_manager.get_current_project().is_sequence_as_protein_in_project(
          tmp_sequence_name
        )
        ):
          tmp_enable_monomer_flag = True
        elif (
                tmp_model_index.data(enums.ModelEnum.TYPE_ROLE)
                == enums.ModelTypeEnum.MULTIMER_SEQ
                and tmp_enable_multimer_flag is False
                and not self._interface_manager.get_current_project().is_sequence_as_protein_in_project(
          tmp_sequence_name
        )
        ):
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
  def __slot_predict_monomer(self) -> None:
    """Opens the predict protein dialog."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Prediction/Monomer' clicked.",
      )
      tmp_indexes = []
      if len(self._view.ui.seqs_list_view.selectedIndexes()) == 0:
        tmp_model = (
          self._interface_manager.get_main_view().ui.seqs_list_view.model()
        )
        for tmp_row_no in range(tmp_model.rowCount()):
          tmp_index = tmp_model.index(tmp_row_no, 0)
          if (
                  tmp_index.data(enums.ModelEnum.TYPE_ROLE)
                  == enums.ModelTypeEnum.MONOMER_SEQ
          ):
            tmp_indexes.append(tmp_index)
      else:
        tmp_indexes = self._view.ui.seqs_list_view.selectedIndexes()
      self._external_controller = (
        predict_protein_view_controller.PredictProteinViewController(
          self._interface_manager,
          self._interface_manager.watcher,
          tmp_indexes,
          "monomer",
        )
      )
      if self._external_controller.has_internet_connection is False:
        return
      self._external_controller.job_input.connect(self._post_predict_protein)
      self._interface_manager.get_predict_protein_view().show()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def __slot_predict_multimer(self) -> None:
    """Opens the predict protein dialog."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Prediction/Multimer' clicked.",
      )
      tmp_indexes = []
      if len(self._view.ui.seqs_list_view.selectedIndexes()) == 0:
        tmp_model = (
          self._interface_manager.get_main_view().ui.seqs_list_view.model()
        )
        for tmp_row_no in range(tmp_model.rowCount()):
          tmp_index = tmp_model.index(tmp_row_no, 0)
          if (
                  tmp_index.data(enums.ModelEnum.TYPE_ROLE)
                  == enums.ModelTypeEnum.MULTIMER_SEQ
          ):
            tmp_indexes.append(tmp_index)
      else:
        tmp_indexes = self._view.ui.seqs_list_view.selectedIndexes()
      self._external_controller = (
        predict_protein_view_controller.PredictProteinViewController(
          self._interface_manager,
          self._interface_manager.watcher,
          tmp_indexes,
          "multimer",
        )
      )
      if self._external_controller.has_internet_connection is False:
        return
      self._external_controller.job_input.connect(self._post_predict_protein)
      self._interface_manager.get_predict_protein_view().show()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def _setup_prediction_job(
          self, result: tuple
  ) -> Optional[tuple["job.PredictionJob", "job_entry.JobEntryWidget"]]:
    """Sets up the prediction job with the data from the predict protein dialog.

    Args:
        result (tuple): The user inputs form the predict protein dialog.

    Returns:
        A tuple with the job and job entry widget or None if the argument is None.
    """
    # <editor-fold desc="Checks">
    if result is None:
      logger.error("result is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      return

    # </editor-fold>

    _, tmp_prediction_protein_infos, tmp_prediction_configuration, _ = result
    self._interface_manager.watcher.add_proteins_from_new_job(
      tmp_prediction_protein_infos
    )
    return self._interface_manager.job_manager.create_prediction_job(
      self._interface_manager.get_current_project(),
      tmp_prediction_protein_infos,
      tmp_prediction_configuration,
      self._interface_manager.project_lock,
      self._interface_manager,
    )

  def _setup_prediction_and_analysis_job(
          self, result: tuple
  ) -> Optional[
    tuple["job.PredictionAndDistanceAnalysisJob", "job_entry.JobEntryWidget"]
  ]:
    """Sets up the prediction and analysis job with the data from the predict protein dialog.

    Args:
        result (tuple): The user inputs form the predict protein dialog.

    Returns:
        A tuple with the job and job entry widget or None if the argument is None.
    """
    # <editor-fold desc="Checks">
    if result is None:
      logger.error("result is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      return

    # </editor-fold>

    _, tmp_prediction_protein_infos, tmp_prediction_configuration, _ = result
    self._interface_manager.watcher.add_proteins_from_new_job(
      tmp_prediction_protein_infos
    )
    tmp_prediction_job, _ = (
      self._interface_manager.job_manager.create_prediction_job(
        self._interface_manager.get_current_project(),
        tmp_prediction_protein_infos,
        tmp_prediction_configuration,
        self._interface_manager.project_lock,
        self._interface_manager,
      )
    )
    tmp_raw_analysis_run_names: list = []
    for row_no in range(
            self._interface_manager.get_predict_protein_view().ui.list_analysis_overview.count()
    ):
      tmp_raw_analysis_run_names.append(
        self._interface_manager.get_predict_protein_view()
        .ui.list_analysis_overview.item(row_no)
        .text()
      )

    self._interface_manager.watcher.add_protein_pairs_from_new_job(
      tmp_raw_analysis_run_names
    )
    tmp_distance_analysis_job, _ = (
      self._interface_manager.job_manager.create_distance_analysis_job(
        self._interface_manager.get_current_project(),
        self._interface_manager.project_lock,
        self._interface_manager,
        tmp_raw_analysis_run_names,
        self._interface_manager.get_settings_manager().settings.cutoff,
        self._interface_manager.get_settings_manager().settings.cycles,
      )
    )
    tmp_job = self._interface_manager.job_manager.create_prediction_and_distance_analysis_job(
      tmp_prediction_job,
      tmp_distance_analysis_job,
      self._interface_manager,
    )
    return tmp_job

  def _post_predict_protein(self, result: tuple) -> None:
    """Sets up and starts the job with the data from the predict protein dialog.

    Args:
        result (tuple): A tuple with the user input needed for the job.
    """
    # <editor-fold desc="Checks">
    if result is None:
      logger.error("result is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      # <editor-fold desc="Check if WSL2 and ColabFold are installed">
      if globals.g_os == "win32":
        constants.PYSSA_LOGGER.info("Checking if WSL2 is installed ...")
        if not dialog_settings_global.is_wsl2_installed():
          constants.PYSSA_LOGGER.warning("WSL2 is NOT installed.")
          self._interface_manager.get_application_settings().wsl_install = 0
          tmp_dialog = custom_message_box.CustomMessageBoxOk(
            "Prediction failed because the WSL2 environment is not installed!",
            "Structure Prediction",
            custom_message_box.CustomMessageBoxIcons.DANGEROUS.value,
          )
          tmp_dialog.exec_()
          return
        constants.PYSSA_LOGGER.info(
          "Checking if Local Colabfold is installed ..."
        )
        if not dialog_settings_global.is_local_colabfold_installed():
          constants.PYSSA_LOGGER.warning("Local Colabfold is NOT installed.")
          self._interface_manager.get_application_settings().local_colabfold = 0
          tmp_dialog = custom_message_box.CustomMessageBoxOk(
            "Prediction failed because the ColabFold is not installed!",
            "Structure Prediction",
            custom_message_box.CustomMessageBoxIcons.DANGEROUS.value,
          )
          tmp_dialog.exec_()
          return

      # </editor-fold>

      constants.PYSSA_LOGGER.info("Begin prediction process.")
      if result[3] is True:
        constants.PYSSA_LOGGER.info(
          "Running prediction with subsequent analysis."
        )
        tmp_job, tmp_job_widget = self._setup_prediction_and_analysis_job(
          result
        )
      else:
        constants.PYSSA_LOGGER.info(
          "Running prediction without subsequent analysis."
        )
        tmp_job, tmp_job_widget = self._setup_prediction_job(result)
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.job_manager.put_job_into_queue(tmp_job)
      self._interface_manager.add_job_entry_to_job_overview_layout(
        tmp_job_widget
      )
    finally:
      self._interface_manager.refresh_main_view()

  def __slot_abort_prediction(self) -> None:
    """Aborts the running prediction."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Prediction/Abort' clicked.",
      )
      constants.PYSSA_LOGGER.info(
        "Structure prediction process was aborted manually."
      )
      subprocess.run(["wsl", "--shutdown"])
      constants.PYSSA_LOGGER.info("Shutdown of wsl environment.")
      filesystem_io.FilesystemCleaner.clean_prediction_scratch_folder()
      constants.PYSSA_LOGGER.info("Cleaned scratch directory.")
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      tmp_dialog = custom_message_box.CustomMessageBoxOk(
        "The structure prediction was aborted.",
        "Abort Structure Prediction",
        custom_message_box.CustomMessageBoxIcons.INFORMATION.value,
      )
      tmp_dialog.exec_()
    finally:
      self._interface_manager.status_bar_manager.hide_progress_bar()
      self._interface_manager.refresh_main_view()

  # </editor-fold>

  # </editor-fold>

  # <editor-fold desc="Hotspots">
  def __slot_hotspots_protein_regions(self) -> None:
    """Shows the selected hotspots protein region."""
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
      "'Highlight Region' menu entry was clicked.",
    )
    self._highlight_protein_region()

  def _highlight_protein_region(self) -> None:
    """Highlights the selected protein region in the PyMOL viewer window."""
    try:
      if self._interface_manager.pymol_session_manager.check_if_sele_is_empty():
        return

      if self._interface_manager.current_tab_index == 1:
        self._cut_representation_to_selected_protein_region(
          self._interface_manager.get_current_protein_representation_states(),
        )
      elif self._interface_manager.current_tab_index == 2:
        self._cut_representation_to_selected_protein_region(
          self._interface_manager.get_current_protein_pair_representation_states(),
        )
      self._interface_manager.pymol_session_manager.user_pymol_connector.set_custom_setting(
        "valence",
        0,
      )
      self._interface_manager.pymol_session_manager.zoom_to_residue_in_protein_position(
        "sele"
      )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def _cut_representation_to_selected_protein_region(
          self,
          all_representation_toggle_states: list[
            tuple[enums.PyMOLRepresentation, bool]
          ],
  ) -> None:
    """Removes the representation to the selected protein region based on the given toggle states.

    Args:
        all_representation_toggle_states (list[tuple[enums.PyMOLRepresentation, bool]]): A list of tuples containing the representation and toggle states for all representations.

    Raises:
        exception.IllegalArgumentError: If `all_representation_toggle_states` is None.
    """
    if all_representation_toggle_states is None:
      logger.error("all_representation_toggle_states is None.")
      raise exception.IllegalArgumentError(
        "all_representation_toggle_states is None."
      )

    self._interface_manager.pymol_session_manager.user_pymol_connector.select(
      "sele",
      "sele and not hydrogens",
    )
    for tmp_toggle_state in all_representation_toggle_states:
      tmp_representation, tmp_toggle_check_state = tmp_toggle_state
      if tmp_toggle_check_state:
        self._interface_manager.pymol_session_manager.hide_specific_representation(
          tmp_representation.value,
          "not sele",
        )
        self._interface_manager.pymol_session_manager.show_specific_representation(
          tmp_representation.value,
          "sele and not hydrogens",
        )

  # </editor-fold>

  # <editor-fold desc="Settings menu methods">
  def __slot_open_settings_global(self) -> None:
    """Opens the dialog for the global settings."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Settings/Edit' clicked.",
      )
      self._external_controller = (
        settings_view_controller.SettingsViewController(
          self._interface_manager
        )
      )
      self._external_controller.user_input.connect(
        self.post_open_settings_global
      )
      self._external_controller.restore_ui()
      self._interface_manager.get_settings_view().show()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    # dialog = dialog_settings_global.DialogSettingsGlobal(self._interface_manager)
    # dialog.exec_()
    # self._interface_manager.update_settings()
    # self._workspace_label = QtWidgets.QLabel(f"Current Workspace: {self._workspace_path}")

  def post_open_settings_global(self, return_value: tuple) -> None:
    """Refreshes the workspae model and the main view after the settings dialog closed."""
    try:
      self._interface_manager.refresh_workspace_model()
      self._interface_manager.refresh_main_view()
      try:
        tmp_type = self._interface_manager.get_current_protein_tree_index_type()
      except AttributeError:
        return
      if tmp_type == "chain":
        self._interface_manager.set_current_chain_color_for_ui_for_proteins(
          self._interface_manager.pymol_session_manager
        )
        self._interface_manager.set_repr_state_in_ui_for_protein_chain(
          self._interface_manager.pymol_session_manager
        )
        self._interface_manager.show_protein_pymol_scene_configuration()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def __slot_restore_settings(self) -> None:
    """Restores the settings.xml file to the default values."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Settings/Restore' clicked.",
      )
      tmp_dialog = custom_message_box.CustomMessageBoxYesNo(
        "Are you sure you want to restore all settings?",
        "Restore Settings",
        custom_message_box.CustomMessageBoxIcons.INFORMATION.value,
      )
      tmp_dialog.exec_()
      if tmp_dialog.response:
        tools.restore_default_settings(
          self._interface_manager.get_application_settings()
        )
        self._view.status_bar.showMessage(
          "Settings were successfully restored."
        )
        logging.info("Settings were successfully restored.")
      else:
        self._view.status_bar.showMessage("Settings were not modified.")
        logging.info("Settings were not modified.")
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  # </editor-fold>

  # <editor-fold desc="Help menu methods">
  def __slot_arrange_windows(self) -> None:
    """Arranges the PySSA and PyMOL window in a vertical split."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Help/Arrange Windows' clicked.",
      )
      if not os.path.exists(constants.ARRANGE_WINDOWS_EXE_FILEPATH):
        tmp_dialog = custom_message_box.CustomMessageBoxOk(
          "The script for arranging the windows could not be found!",
          "Arrange Windows",
          custom_message_box.CustomMessageBoxIcons.ERROR.value,
        )
        tmp_dialog.exec_()
      else:
        logger.debug("Started script to arrange window ...")
        subprocess.Popen(
          [constants.ARRANGE_WINDOWS_EXE_FILEPATH],
          creationflags=subprocess.CREATE_NO_WINDOW,
        )
        logger.debug("Script to arrange windows finished.")
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def __slot_open_logs(self) -> None:
    """Opens a file explorer with all log files and can open a log file in the default application."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Help/Show Logs in Explorer' clicked.",
      )
      file_dialog = QtWidgets.QFileDialog()
      log_path = str(constants.LOG_PATH)
      file_dialog.setDirectory(log_path)
      file_path, _ = file_dialog.getOpenFileName(
        self._view, "Select a log file to open", "", "LOG File (*.log)"
      )
      if file_path:
        os.startfile(file_path)
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def __slot_clear_all_log_files(self) -> None:
    """Clears all log files generated under .pyssa/logs."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Help/Clear All Logs' clicked.",
      )
      tmp_dialog = custom_message_box.CustomMessageBoxYesNo(
        "Are you sure you want to delete all log files?",
        "Clear Log Files",
        custom_message_box.CustomMessageBoxIcons.WARNING.value,
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
          self._interface_manager.status_bar_manager.show_temporary_message(
            "All log files could be deleted."
          )
          constants.PYSSA_LOGGER.info("All log files were deleted.")
        else:
          tmp_dialog = custom_message_box.CustomMessageBoxOk(
            "Not all log files could be deleted.",
            "Clear Log Files",
            custom_message_box.CustomMessageBoxIcons.WARNING.value,
          )
          tmp_dialog.exec_()
          constants.PYSSA_LOGGER.warning("Not all log files were deleted!")
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def __slot_open_tutorial(self) -> None:
    """Opens the official tutorial pdf file."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Help/Tutorials' clicked.",
      )
      tmp_dialog = dialog_tutorial_videos.TutorialVideosDialog()
      tmp_dialog.exec_()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def open_documentation(self) -> None:
    """Opens the official plugin documentation as PDF."""
    try:
      os.startfile(constants.DOCS_PATH)
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def __slot_open_about(self) -> None:
    """Opens the About dialog."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Help/About' clicked.",
      )
      dialog = dialog_about.DialogAbout()
      dialog.exec_()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def __slot_get_demo_projects(self) -> None:
    """Starts download of demo projects."""
    try:
      if not tools.check_internet_connectivity():
        tmp_dialog = custom_message_box.CustomMessageBoxOk(
          "You do not have a working internet connection\nbut that is necessary for this operation!",
          "Internet Connection",
          custom_message_box.CustomMessageBoxIcons.ERROR.value,
        )
        tmp_dialog.exec_()
        return
      # self._active_task = tasks.LegacyTask(
      #   target=util_async.download_demo_projects,
      #   args=(
      #     self._interface_manager.get_application_settings().workspace_path,
      #     0,
      #   ),
      #   post_func=self.__await_download_demo_projects,
      # )

      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=util_async.download_demo_projects,
              args=(
                self._interface_manager.get_application_settings().workspace_path,
                0,
              ),
            ),
            an_await_function=self.__await_download_demo_projects,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.status_bar_manager.show_temporary_message(
        "Getting demo projects ...", a_with_timeout_flag=False
      )
      self._interface_manager.block_gui()

    # try:
    #     logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "Menu entry 'Help/Get Demo Projects' clicked.")
    #     self._interface_manager.status_bar_manager.show_temporary_message(
    #         "Getting demo projects ...", False)
    #     import zipfile
    #     download_dest = pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects.zip")
    #     if os.path.exists(download_dest):
    #         os.remove(download_dest)
    #     if not os.path.exists(download_dest):
    #         # download demo projects
    #         url = f'https://w-hs.sciebo.de/s/ZHJa6XB9SKWtqGi/download'
    #         tmp_error_flag = False
    #         try:
    #             response = requests.get(url)
    #             response.raise_for_status()  # Check for errors
    #             zipfile = zipfile.ZipFile(BytesIO(response.content))
    #             zipfile.extractall(pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects"))
    #         except requests.exceptions.HTTPError as errh:
    #             constants.PYSSA_LOGGER.error(f"HTTP Error: {errh}")
    #             tmp_error_flag = True
    #         except requests.exceptions.ConnectionError as errc:
    #             constants.PYSSA_LOGGER.error(f"Error Connecting: {errc}")
    #             tmp_error_flag = True
    #         except requests.exceptions.Timeout as errt:
    #             constants.PYSSA_LOGGER.error(f"Timeout Error: {errt}")
    #             tmp_error_flag = True
    #         except requests.exceptions.RequestException as err:
    #             constants.PYSSA_LOGGER.error(f"Error: {err}")
    #             tmp_error_flag = True
    #         else:
    #             constants.PYSSA_LOGGER.info(f"Demo projects downloaded and extracted successfully.")
    #
    #         if tmp_error_flag:
    #             tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #                 "The download of the demo projects failed. Please try again later.",
    #                 "Get Demo Projects",
    #                 custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #             )
    #             tmp_dialog.exec_()
    #             self._interface_manager.status_bar_manager.show_error_message("The download of the demo projects failed.")
    #             return
    #     else:
    #         constants.PYSSA_LOGGER.info("Demo projects are getting extracted ...")
    #         try:
    #             with zipfile.ZipFile(pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects.zip"), "r") as zip_ref:
    #                 zip_ref.extractall(pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects"))
    #             constants.PYSSA_LOGGER.info(
    #                 "Demo projects are downloaded and extracted.\n Import of demo projects started ...",
    #             )
    #         except Exception as e:
    #             constants.PYSSA_LOGGER.error(f"Extraction process of demo projects finished with the error: {e}.")
    #             tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #                 "Extraction process of demo projects finished with an error. Check the logs to get more information.",
    #                 "Get Demo Projects",
    #                 custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #             )
    #             tmp_dialog.exec_()
    #             self._interface_manager.status_bar_manager.show_error_message("Extraction process of demo projects finished with an error.")
    #             return
    #     try:
    #         path_of_demo_projects = pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects")
    #         for tmp_filename in os.listdir(path_of_demo_projects):
    #             # Copy db file into new workspace
    #             tmp_project_database_filepath = str(
    #                 pathlib.Path(
    #                     f"{self._interface_manager.get_application_settings().workspace_path}/{tmp_filename}"
    #                 )
    #             )
    #             tmp_src_filepath = str(pathlib.Path(f"{path_of_demo_projects}/{tmp_filename}"))
    #             shutil.copyfile(tmp_src_filepath, tmp_project_database_filepath)
    #         constants.PYSSA_LOGGER.info("Import process of demo projects finished.")
    #     except Exception as e:
    #         constants.PYSSA_LOGGER.error(f"Import process of demo projects finished with the error: {e}.")
    #         tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #             "Import process of demo projects finished with an error. Check the logs to get more information.",
    #             "Get Demo Projects",
    #             custom_message_box.CustomMessageBoxIcons.DANGEROUS.value
    #         )
    #         tmp_dialog.exec_()
    #         self._interface_manager.status_bar_manager.show_error_message("Import process of demo projects finished with an error.")
    #     else:
    #         self._interface_manager.refresh_workspace_model()
    #         self._interface_manager.refresh_main_view()
    #         self._interface_manager.status_bar_manager.show_temporary_message("Getting demo projects finished successfully.")
    # except Exception as e:
    #     logger.error(f"An error occurred: {e}")
    #     self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

  def __await_download_demo_projects(self, return_value: tuple[str, list[tuple[bool, tuple]]]) -> None:
    """Integrates demo projects in the user's workspace.

    Args:
        return_value (tuple[str, list[tuple[bool, tuple]]]): The success flag whether the operation was successful.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      if tmp_result[0] is True:
        self._interface_manager.refresh_workspace_model()
        self._interface_manager.status_bar_manager.show_temporary_message(
          "Getting demo projects finished successfully."
        )
      else:
        self._interface_manager.status_bar_manager.show_error_message(
          "The download of the demo projects failed."
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  # </editor-fold>

  # <editor-fold desc="Image menu methods">
  def __slot_preview_image(self) -> None:
    """Starts to create an image preview."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Image/Preview' clicked.",
      )
      # self._active_task = tasks.LegacyTask(
      #   target=image_async.preview_image,
      #   args=(self._interface_manager.pymol_session_manager, 0),
      #   post_func=self.__await_preview_image,
      # )

      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=image_async.preview_image,
              args=(
                self._interface_manager.pymol_session_manager, 0,
              ),
            ),
            an_await_function=self.__await_preview_image,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui(with_wait_cursor=True)
      self.update_status("Creating preview of image ...")

  def __await_preview_image(self, return_value: tuple) -> None:
    """Refreshes the main view and reverts the cursor after image preview."""
    self._interface_manager.stop_wait_cursor()
    self._interface_manager.refresh_main_view()
    self.update_status("Preview finished.")

  def __slot_create_ray_traced_image(self) -> None:
    """Creates a ray-tracing job for a ray-traced image."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Image/Ray-Traced' clicked.",
      )
      save_dialog = QtWidgets.QFileDialog()
      full_file_name = save_dialog.getSaveFileName(
        caption="Save Image", filter="Image (*.png)"
      )
      if full_file_name == ("", ""):
        logger.info("No file has been selected.")
        self.update_status("No file has been selected.")
        return

      # --- New job approach
      tmp_session_filepath = (
        self._interface_manager.pymol_session_manager.save_current_pymol_session_as_pse_cache_file()
      )
      tmp_ray_tracing_job, tmp_ray_tracing_entry_widget = (
        self._interface_manager.job_manager.create_ray_tracing_job(
          full_file_name[0],
          tmp_session_filepath,
          self._interface_manager.get_application_settings().image_ray_trace_mode,
          self._interface_manager.get_application_settings().image_ray_texture,
          self._interface_manager.get_application_settings().image_renderer,
          self._interface_manager,
          self._interface_manager.get_current_project().get_project_name(),
        )
      )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.job_manager.put_job_into_queue(
        tmp_ray_tracing_job
      )
      self._interface_manager.add_job_entry_to_job_overview_layout(
        tmp_ray_tracing_entry_widget
      )

  def __slot_create_drawn_image(self) -> None:
    """Starts to create a drawn image."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Image/Simple' clicked.",
      )
      save_dialog = QtWidgets.QFileDialog()
      full_file_name = save_dialog.getSaveFileName(
        caption="Save Image", filter="Image (*.png)"
      )
      if full_file_name == ("", ""):
        logger.info("No file has been selected.")
        self.update_status("No file has been selected.")
        return

      self._interface_manager.pymol_session_manager.async_cmds(
        self._interface_manager.get_task_manager(),
        self._interface_manager.get_task_scheduler(),
        (pymol_enums.CommandEnum.DRAW, pymol_enums.CommandEnum.PNG),
        ((2400, 2400, 2), (full_file_name[0], 300)),
        self.__await_create_drawn_image
      )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self.update_status("Creating simple image ...")
      self._interface_manager.block_gui(with_wait_cursor=True)

  def __await_create_drawn_image(self, return_value: tuple) -> None:
    """Refreshes the main view and reverts the cursor after image creation."""
    self._interface_manager.stop_wait_cursor()
    self._interface_manager.refresh_main_view()
    self.update_status("Image creation finished.")

  # </editor-fold>

  # <editor-fold desc="Sequences tab methods">
  def __slot_open_text_editor_for_seq(self) -> None:
    """Opens a dialog that displays the sequence."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "A sequence in the 'Addition Information' table was clicked.",
      )
      if (
              self._view.ui.seqs_table_widget.currentColumn() == 1
              and self._view.ui.seqs_table_widget.currentRow() == 0
      ):
        self.__slot_rename_selected_sequence()
      elif (
              self._view.ui.seqs_table_widget.currentColumn() == 1
              and self._view.ui.seqs_table_widget.currentRow() == 1
      ):
        self.tmp_txt_browser = QtWidgets.QTextBrowser()
        try:
          tmp_seq = (
            self._view.ui.seqs_table_widget.currentItem()
            .data(enums.ModelEnum.OBJECT_ROLE)
            .seq
          )
          tmp_seqs = tmp_seq.split(",")
          tmp_seq = ",\n\n".join(tmp_seqs)
          self.tmp_txt_browser.setText(tmp_seq)
        except AttributeError:
          return
        else:
          self.tmp_txt_browser.setWindowTitle("View Sequence")
          self.tmp_txt_browser.setWindowIcon(
            QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH)
          )
          self.tmp_txt_browser.resize(500, 150)
          self.tmp_txt_browser.show()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def _set_new_sequence_name_in_table_item(self) -> None:
    """Sets the new sequence name in the table widget item."""
    try:
      tmp_new_seq = self._view.line_edit_seq_name.text()
      self._view.ui.seqs_table_widget.item(0, 1).setText(tmp_new_seq)
    except AttributeError:
      return

  def _rename_sequence(self) -> None:
    """Renames a selected sequence in the table widget."""
    tmp_old_name = (
      self._view.ui.seqs_list_view.currentIndex()
      .data(enums.ModelEnum.OBJECT_ROLE)
      .name
    )
    try:
      # this is needed because the signal is fired even if the current item is None
      tmp_new_name = self._view.ui.seqs_table_widget.item(0, 1).text()
    except AttributeError:
      return
    else:
      tmp_seq = (
        self._view.ui.seqs_list_view.currentIndex()
        .data(enums.ModelEnum.OBJECT_ROLE)
        .seq
      )
      self._view.ui.seqs_list_view.currentIndex().data(
        enums.ModelEnum.OBJECT_ROLE
      ).name = tmp_new_name
      self._view.ui.seqs_list_view.model().setData(
        self._view.ui.seqs_list_view.currentIndex(),
        tmp_new_name,
        Qt.DisplayRole,
      )
      tmp_database_operation = database_operation.DatabaseOperation(
        enums.SQLQueryType.UPDATE_SEQUENCE_NAME,
        (0, tmp_new_name, tmp_old_name, tmp_seq),
      )
      self._database_thread.put_database_operation_into_queue(
        tmp_database_operation
      )

  def __slot_show_sequence_information(self) -> None:
    """Shows the table widget with additional information about the selected sequence."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        f"The sequence '{self._view.ui.seqs_list_view.currentIndex().data(Qt.DisplayRole)}' on the 'Sequence Tab' was clicked.",
      )
      self._interface_manager.show_sequence_parameters(
        self._view.ui.seqs_list_view.currentIndex(),
      )
      self._view.ui.btn_save_sequence.setEnabled(True)
      self._view.ui.btn_delete_sequence.setEnabled(True)
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def __slot_import_sequence(self) -> None:
    """Opens the import seqeunce dialog."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Import sequence' button on the 'Sequence Tab' was clicked.",
      )
      self._external_controller = (
        import_sequence_view_controller.ImportSequenceViewController(
          self._interface_manager
        )
      )
      self._external_controller.user_input.connect(self._post_import_sequence)
      self._external_controller.restore_ui()
      self._interface_manager.get_import_sequence_view().show()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def _post_import_sequence(self, return_value: tuple) -> None:
    """Inserts the new sequence in the database and refreshes the sequence model.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      for tmp_seq_record in return_value[1]:
        logger.info(
          f"Adding new sequence {tmp_seq_record.name} with {tmp_seq_record.seq} to the current project."
        )
        self._interface_manager.get_current_project().sequences.append(
          tmp_seq_record
        )
        tmp_database_operation = database_operation.DatabaseOperation(
          enums.SQLQueryType.INSERT_NEW_SEQUENCE,
          (0, tmp_seq_record),
        )
        self._database_thread.put_database_operation_into_queue(
          tmp_database_operation
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.refresh_sequence_model()
      self._interface_manager.show_menu_options_with_seq()
    finally:
      self._interface_manager.refresh_main_view()

  def __slot_add_sequence(self) -> None:
    """Opens the add sequence dialog."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Add sequence' button on the 'Sequence Tab' was clicked.",
      )
      self._external_controller = (
        add_sequence_view_controller.AddSequenceViewController(
          self._interface_manager
        )
      )
      self._external_controller.return_value.connect(self._post_add_sequence)
      self._external_controller.restore_default_view()
      self._interface_manager.get_add_sequence_view().show()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def _post_add_sequence(self, return_value: tuple) -> None:
    """Inserts the new sequence in the database and refreshes the sequence model.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      logger.info(
        f"Adding new sequence {return_value[0]} with {return_value[1]} to the current project."
      )
      tmp_seq_name = return_value[0]
      tmp_sequence = return_value[1]
      tmp_seq_record = SeqRecord.SeqRecord(tmp_sequence, name=tmp_seq_name)
      self._interface_manager.get_current_project().sequences.append(
        tmp_seq_record
      )
      tmp_database_operation = database_operation.DatabaseOperation(
        enums.SQLQueryType.INSERT_NEW_SEQUENCE,
        (0, tmp_seq_record),
      )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._database_thread.put_database_operation_into_queue(
        tmp_database_operation
      )
      self._interface_manager.refresh_sequence_model()
    finally:
      self._interface_manager.refresh_main_view()

  def __slot_save_selected_sequence_as_fasta_file(self) -> None:
    """Opens a QFileDialog to choose a filepath and saves the sequence as fasta file to that filepath."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Export sequence' button on the 'Sequence Tab' was clicked.",
      )

      file_dialog = QtWidgets.QFileDialog()  # noqa: ANN001
      desktop_path = QtCore.QStandardPaths.standardLocations(
        QtCore.QStandardPaths.DesktopLocation
      )[0]
      file_dialog.setDirectory(desktop_path)
      file_path, _ = file_dialog.getSaveFileName(
        self._view,
        "Save Protein Sequence",
        "",
        "FASTA File (*.fasta)",
      )
      if file_path:
        tmp_seq_record = (
          self._interface_manager.get_current_sequence_list_index_object()
        )
        # pre-process seq record object for the SeqIO module
        if tmp_seq_record.id == "<unknown id>":
          tmp_seq_record.id = tmp_seq_record.name
        tmp_seq_record.seq = Seq(tmp_seq_record.seq)
        # defines the task to save the sequence as .fasta file
        # self._active_task = tasks.LegacyTask(
        #   target=sequence_async.save_selected_protein_sequence_as_fasta_file,
        #   args=(
        #     tmp_seq_record,
        #     file_path,
        #   ),
        #   post_func=self.__await_save_selected_sequence_as_fasta_file,
        # )

        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=sequence_async.save_selected_protein_sequence_as_fasta_file,
                args=(
                  tmp_seq_record,
                  file_path,
                ),
              ),
              an_await_function=self.__await_save_selected_sequence_as_fasta_file,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui(with_wait_cursor=True)
    finally:
      self._interface_manager.refresh_main_view()

  def __await_save_selected_sequence_as_fasta_file(self, return_value: tuple[str, list[tuple[bool, tuple]]]) -> None:
    """Finishes the saving sequence to fasta file process.

    Args:
        return_value (tuple[str, list[tuple[bool, tuple]]]): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("result is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      self._interface_manager.stop_wait_cursor()
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      if tmp_result[0] == exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0]:
        tmp_dialog = custom_message_box.CustomMessageBoxOk(
          "Saving the sequence as .fasta file failed!",
          "Save Protein Sequence",
          custom_message_box.CustomMessageBoxIcons.DANGEROUS.value,
        )
        tmp_dialog.exec_()
      elif tmp_result[0] == exit_codes.EXIT_CODE_ZERO[0]:
        # tmp_dialog = custom_message_box.CustomMessageBoxOk(
        #     "The sequence was successfully saved as .fasta file.",
        #     "Save Protein Sequence",
        #     custom_message_box.CustomMessageBoxIcons.INFORMATION.value
        # )
        # tmp_dialog.exec_()
        self._interface_manager.status_bar_manager.show_temporary_message(
          "The sequence was successfully saved as .fasta file."
        )
      else:  # noqa: ANN001
        tmp_dialog = custom_message_box.CustomMessageBoxOk(
          "Saving the sequence as .fasta file failed with an unexpected error!",
          "Save Protein Sequence",
          custom_message_box.CustomMessageBoxIcons.DANGEROUS.value,
        )
        tmp_dialog.exec_()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.refresh_sequence_model()
    finally:
      self._interface_manager.refresh_main_view()

  def __slot_delete_selected_sequence(self) -> None:
    """Deletes a selected sequence from the project."""
    # popup message which warns the user that the selected sequence gets deleted
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
      "'Delete sequence' button on the 'Sequence Tab' was clicked.",
    )
    tmp_dialog = custom_message_box.CustomMessageBoxDelete(
      "Are you sure you want to delete this sequence?",
      "Delete Sequence",
      custom_message_box.CustomMessageBoxIcons.WARNING.value,
    )
    tmp_dialog.exec_()
    if tmp_dialog.response:
      try:
        tmp_seq_record: "SeqRecord.SeqRecord" = (
          self._interface_manager.get_current_sequence_list_index_object()
        )
        tmp_database_operation = database_operation.DatabaseOperation(
          enums.SQLQueryType.DELETE_EXISTING_SEQUENCE, (0, tmp_seq_record)
        )
        self._database_thread.put_database_operation_into_queue(
          tmp_database_operation
        )
        self._interface_manager.get_current_project().delete_specific_sequence(
          tmp_seq_record.name
        )
        self._interface_manager.status_bar_manager.show_temporary_message(
          "The sequence was successfully deleted."
        )
      except Exception as e:
        logger.error(f"An error occurred: {e}")
        self._interface_manager.status_bar_manager.show_error_message(
          "An unknown error occurred!"
        )
      else:
        self._interface_manager.refresh_sequence_model()
        self._interface_manager.refresh_main_view()
        # extra ui changes
        self._view.ui.seqs_table_widget.setRowCount(0)
        self._view.build_sequence_table()

  def __slot_rename_selected_sequence(self) -> None:
    """Opens a new view to rename the selected sequence."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Rename sequence' context menu action was clicked.",
      )
      self._external_controller = (
        rename_sequence_view_controller.RenameSequenceViewController(
          self._interface_manager
        )
      )
      self._external_controller.user_input.connect(
        self.post_rename_selected_sequence_structure
      )
      self._external_controller.restore_ui()
      self._interface_manager.get_rename_sequence_view().show()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def post_rename_selected_sequence_structure(
          self, return_value: tuple
  ) -> None:
    """Finishes the rename sequence process.

    Args:
        return_value (tuple): The result data from the async method.
    """
    try:
      tmp_new_name = return_value[0]
      tmp_old_name = (
        self._view.ui.seqs_list_view.currentIndex()
        .data(enums.ModelEnum.OBJECT_ROLE)
        .name
      )
      tmp_seq = (
        self._view.ui.seqs_list_view.currentIndex()
        .data(enums.ModelEnum.OBJECT_ROLE)
        .seq
      )
      self._view.ui.seqs_list_view.currentIndex().data(
        enums.ModelEnum.OBJECT_ROLE
      ).name = tmp_new_name
      self._view.ui.seqs_list_view.model().setData(
        self._view.ui.seqs_list_view.currentIndex(),
        tmp_new_name,
        Qt.DisplayRole,
      )
      tmp_database_operation = database_operation.DatabaseOperation(
        enums.SQLQueryType.UPDATE_SEQUENCE_NAME,
        (0, tmp_new_name, tmp_old_name, tmp_seq),
      )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._database_thread.put_database_operation_into_queue(
        tmp_database_operation
      )
      self._view.ui.seqs_table_widget.item(0, 1).setText(tmp_new_name)
    finally:
      self._interface_manager.refresh_main_view()

  def open_context_menu_for_sequences(self, position) -> None:
    """Opens the context menu for the sequences tab."""
    tmp_context_menu = self._sequence_list_context_menu.get_context_menu(
      self._view.ui.seqs_list_view.selectedIndexes(),
    )
    tmp_context_menu.exec_(
      self._view.ui.seqs_list_view.viewport().mapToGlobal(position)
    )
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
  def __slot_expand_protein(self) -> None:
    """Expands a protein branch."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "A protein of the tree view was expanded.",
      )
      tmp_type = self._interface_manager.get_current_protein_tree_index().data(
        enums.ModelEnum.TYPE_ROLE
      )
      if tmp_type == "protein":
        # protein
        self._view.ui.proteins_tree_view.setExpanded(
          self._interface_manager.get_current_protein_tree_index(),
          True,
        )
        # scenes
        self._view.ui.proteins_tree_view.setExpanded(
          self._interface_manager.get_current_protein_tree_index().child(
            0, 0
          ),
          True,
        )
        # chains
        self._view.ui.proteins_tree_view.setExpanded(
          self._interface_manager.get_current_protein_tree_index().child(
            1, 0
          ),
          True,
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def __slot_collapse_protein(self) -> None:
    """Collapses a protein branch."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "A protein of the tree view was collapsed.",
      )
      tmp_type = self._interface_manager.get_current_protein_tree_index().data(
        enums.ModelEnum.TYPE_ROLE
      )
      if tmp_type == "protein":
        # protein
        self._view.ui.proteins_tree_view.collapse(
          self._interface_manager.get_current_protein_tree_index(),
        )
        # scenes
        self._view.ui.proteins_tree_view.collapse(
          self._interface_manager.get_current_protein_tree_index().child(
            0, 0
          ),
        )
        # chains
        self._view.ui.proteins_tree_view.collapse(
          self._interface_manager.get_current_protein_tree_index().child(
            1, 0
          ),
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def __slot_expand_all_proteins(self) -> None:
    """Expands all protein tree branches."""
    self._view.ui.proteins_tree_view.expandAll()

  def __slot_collapse_all_proteins(self) -> None:
    """Collapses all protein tree branches."""
    self._view.ui.proteins_tree_view.collapseAll()

  def _open_context_menu_for_proteins(self, position) -> None:
    """Opens the context menu for the proteins tab."""
    try:
      tmp_protein = self._interface_manager.get_current_active_protein_object()
    except ValueError:
      tmp_is_protein_in_any_pair_flag = True
      tmp_is_protein_in_session_flag = False
    else:
      tmp_is_protein_in_any_pair_flag = self._interface_manager.get_current_project().check_if_protein_is_in_any_protein_pair(
        tmp_protein.get_molecule_object(),
      )
      tmp_is_protein_in_session_flag = self._interface_manager.pymol_session_manager.is_the_current_protein_in_session(
        self._interface_manager.get_current_active_protein_object().get_molecule_object()
      )
    tmp_is_protein_expanded_flag: bool = False
    try:
      if (
              self._interface_manager.get_current_protein_tree_index().data(
                enums.ModelEnum.TYPE_ROLE
              )
              == "protein"
      ):
        if self._view.ui.proteins_tree_view.isExpanded(
                self._interface_manager.get_current_protein_tree_index()
        ):
          tmp_is_protein_expanded_flag: bool = True
    except Exception as e:
      logger.error(e)
    else:
      tmp_context_menu = self._protein_tree_context_menu.get_context_menu(
        self._view.ui.proteins_tree_view.selectedIndexes(),
        self._interface_manager.get_current_protein_tree_index_type(),
        tmp_is_protein_in_any_pair_flag,
        tmp_is_protein_in_session_flag,
        tmp_is_protein_expanded_flag,
      )
      tmp_context_menu.exec_(
        self._view.ui.proteins_tree_view.viewport().mapToGlobal(position)
      )
      self.__slot_get_information_about_selected_object_in_protein_branch()  # fixme: This should be done in a better way than this!

  # <editor-fold desc="PyMOL session">
  def __slot_open_protein_pymol_session(self) -> None:
    """Starts an async method that opens a protein pymol session."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Open protein pymol session' button on the 'Proteins Tab' was clicked.",
      )
      self._view.tg_protein_white_bg.toggle_button.setCheckState(False)
      tmp_protein: "protein.Protein" = (
        self._interface_manager.get_current_active_protein_object()
      )
      tmp_flag = False

      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=pymol_session_async.load_protein_pymol_session,
              args=(
                tmp_protein,
                self._interface_manager.pymol_session_manager,
                tmp_flag,
              ),
            ),
            an_await_function=self.__await_open_protein_pymol_session,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
      self._interface_manager.status_bar_manager.show_temporary_message(
        f"Loading PyMOL session of {tmp_protein.get_molecule_object()} ...",
        False,
      )
    except Exception as e:
      logger.error(f"The error {e} occurred during the pymol session loading!")
      self._interface_manager.status_bar_manager.show_error_message(
        "Loading the PyMOL session failed!"
      )
      self._interface_manager.stop_wait_cursor()
    else:
      self._interface_manager.block_gui(with_wait_cursor=True)

  def __await_open_protein_pymol_session(self, return_value: tuple) -> None:
    """Finishes the open protein pymol session process.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      return

    # </editor-fold>

    try:
      logger.debug("Returning from async function.")
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      tmp_pymol_session_manager, exit_boolean = tmp_result
      self._interface_manager.pymol_session_manager = tmp_pymol_session_manager
      self._view.ui.action_protein_regions.setEnabled(False)
      if exit_boolean:
        self._view.cb_chain_color.setEnabled(True)
        self._view.cb_chain_representation.setEnabled(True)
        self._view.ui.action_protein_regions.setEnabled(True)
        self._view.ui.btn_create_protein_scene.setEnabled(True)
        self._view.ui.btn_update_protein_scene.setEnabled(True)
        self._view.ui.lbl_session_name.setText(
          f"Session Name: {self._interface_manager.pymol_session_manager.session_name}"
        )
        self._view.ui.lbl_pymol_protein_scene.setText("PyMOL Scene: base")
        self._view.ui.lbl_info.setText("Please select a chain.")
        logger.info("Successfully opened protein session.")
        self._interface_manager.status_bar_manager.show_temporary_message(
          "Loading the PyMOL session was successful."
        )
      else:
        logger.error(
          "The protein name could not be found in the object list in PyMOL!"
        )
        self._view.cb_chain_color.setEnabled(False)
        self._view.cb_chain_representation.setEnabled(False)
        self._view.ui.btn_create_protein_scene.setEnabled(False)
        self._view.ui.btn_update_protein_scene.setEnabled(False)
        self._interface_manager.status_bar_manager.show_error_message(
          "Loading the PyMOL session failed!"
        )
        self._view.ui.lbl_info.setText(
          "Please load the PyMOL session of the selected protein."
        )
    except Exception as e:
      logger.error(f"The error {e} occurred during the pymol session loading!")
      self._interface_manager.status_bar_manager.show_error_message(
        "Loading the PyMOL session failed!"
      )
    finally:
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()

  def _save_protein_pymol_session(self) -> None:
    """Saves the session as base64 string and updates the database."""
    try:
      tmp_protein = self._interface_manager.get_current_active_protein_object()
      tmp_protein.pymol_session = (
        self._interface_manager.pymol_session_manager.save_current_session_as_base64()
      )
      tmp_database_operation = database_operation.DatabaseOperation(
        enums.SQLQueryType.UPDATE_PYMOL_SESSION_PROTEIN,
        (0, self._interface_manager.get_current_active_protein_object()),
      )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._database_thread.put_database_operation_into_queue(
        tmp_database_operation
      )

  # </editor-fold>

  def __slot_get_information_about_selected_object_in_protein_branch(
          self,
  ) -> None:
    """Modifies GUI that all available options are seen for the pymol scene configuration panel."""
    try:
      tmp_type = self._interface_manager.get_current_protein_tree_index_type()
      if tmp_type == "protein":
        logger.log(
          log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
          f"The protein object '{self._view.ui.proteins_tree_view.currentIndex().data(Qt.DisplayRole)}' on the 'Proteins Tab' was clicked.",
        )
      elif tmp_type == "scene":
        self._setup_protein_scene()
      elif tmp_type == "chain":
        self._setup_protein_pymol_scene_config()
      elif tmp_type == "header":
        logger.log(
          log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
          f"The header '{self._view.ui.proteins_tree_view.currentIndex().data(Qt.DisplayRole)}' on the 'Proteins Tab' was clicked.",
        )
      else:
        logger.warning("Unknown object type occurred in Protein tab.")
        return
      self._interface_manager.manage_ui_of_protein_tab(
        tmp_type,
        self._interface_manager.get_current_project().check_if_protein_is_in_any_protein_pair(
          self._interface_manager.get_current_active_protein_object().get_molecule_object(),
        ),
        self._interface_manager.pymol_session_manager,
      )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  # <editor-fold desc="Setup protein scene of proteins tree">
  def _setup_protein_scene(self) -> None:
    """Loads the default pymol scene."""
    try:
      tmp_scene_name = self._view.ui.proteins_tree_view.currentIndex().data(
        Qt.DisplayRole
      )
      tmp_protein_name = (
        self._view.ui.proteins_tree_view.currentIndex()
        .parent()
        .parent()
        .data(Qt.DisplayRole)
      )
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        f"The scene '{tmp_scene_name}' of the protein '{tmp_protein_name}' on the 'Proteins Tab' was clicked.",
      )
      if not self._interface_manager.pymol_session_manager.is_the_current_protein_in_session(
              tmp_protein_name
      ):
        # self._interface_manager.get_current_active_protein_object().get_molecule_object()
        return
      tmp_scene_name = self._interface_manager.get_current_active_scene_name()
      self._interface_manager.pymol_session_manager.current_scene_name = (
        tmp_scene_name
      )

      self._interface_manager.pymol_session_manager.async_cmd(
        self._interface_manager.get_task_manager(),
        self._interface_manager.get_task_scheduler(),
        pymol_enums.CommandEnum.LOAD_SCENE,
        (tmp_scene_name,),
        self.__await_load_scene_protein
      )
      # self._task_result = tasks.TaskResult.run_action(  # fixme: raises currently an error if the unlock function is called
      #     tasks.Action(
      #         a_target=pymol_session_async.load_scene,
      #         args=(
      #             self._interface_manager.pymol_session_manager,
      #             tmp_scene_name,
      #         ),
      #     ),
      #     self.thread_pool,
      #     an_await_function=self.__await_load_scene_protein,
      # )

      # self._active_task = tasks.LegacyTask(
      #     target=pymol_session_async.load_scene,
      #     args=(self._interface_manager.pymol_session_manager, tmp_scene_name),
      #     post_func=self.__await_load_scene_protein
      # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.status_bar_manager.show_temporary_message(
        "Loading PyMOL scene ...", a_with_timeout_flag=False
      )
      self._interface_manager.block_gui()

  def __await_load_scene_protein(self, return_value: tuple[bool]) -> None:
    """Finishes the load default pymol scene process.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      ui_util.set_pymol_scene_name_into_label(
        self._interface_manager.pymol_session_manager.current_scene_name,
        self._view.ui.lbl_pymol_protein_scene,
      )
      if return_value[0]:
        self._interface_manager.status_bar_manager.show_temporary_message(
          "Loading PyMOL scene was successful."
        )
      else:
        self._interface_manager.status_bar_manager.show_temporary_message(
          "Loading PyMOL scene failed! Please try again."
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  # </editor-fold>

  # <editor-fold desc="Setup PyMOL scene configuration">
  def _setup_protein_pymol_scene_config(self) -> None:
    """Sets up the color grid and representation section of the pymol scene configuration panel."""
    logger.debug("Running method '_setup_protein_pymol_scene_config'.")
    try:
      tmp_chain_letter = self._view.ui.proteins_tree_view.currentIndex().data(
        Qt.DisplayRole
      )
      tmp_protein_name = (
        self._view.ui.proteins_tree_view.currentIndex()
        .parent()
        .parent()
        .data(Qt.DisplayRole)
      )
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        f"The chain object '{tmp_chain_letter}' of the protein '{tmp_protein_name}' on the 'Proteins Tab' was clicked.",
      )
      if (
              self._interface_manager.pymol_session_manager.current_scene_name == ""
              or not self._interface_manager.pymol_session_manager.is_the_current_protein_in_session(
        tmp_protein_name
      )
      ):
        return
      # Set icon for color grid
      self._view.color_grid_proteins.reset_icon_for_selected_color()
      tmp_protein = self._interface_manager.get_current_active_protein_object()
      tmp_chain = self._interface_manager.get_current_active_chain_object()

      self._interface_manager.pymol_session_manager.async_cmds(
        self._interface_manager.get_task_manager(),
        self._interface_manager.get_task_scheduler(),
        (
          pymol_enums.CommandEnum.GET_RESIDUE_COLOR_CONFIG,
          pymol_enums.CommandEnum.GET_CHAIN_REPR_STATE
        ),
        (
          (tmp_protein.get_molecule_object(), tmp_chain.chain_letter),
          (tmp_protein.get_molecule_object(), tmp_chain.chain_letter)
        ),
        self.__await_setup_protein_pymol_scene_config
      )

      # self._task_result = tasks.TaskResult.run_actions(
      #     the_actions=(
      #         tasks.Action(
      #             a_target=pymol_session_async.get_residue_color_config_of_a_given_protein_chain,
      #             args=(
      #                 tmp_protein.get_molecule_object(),
      #                 tmp_chain.chain_letter,
      #                 self._interface_manager.pymol_session_manager,
      #             ),
      #         ),
      #         tasks.Action(
      #             a_target=pymol_session_async.get_representation_config_of_a_given_protein_chain,
      #             args=(
      #                 tmp_protein.get_molecule_object(),
      #                 tmp_chain.chain_letter,
      #                 self._interface_manager.pymol_session_manager,
      #             ),
      #         ),
      #     ),
      #     the_threadpool=self.thread_pool,
      #     an_await_function=self.__await_setup_protein_pymol_scene_config,
      # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __await_setup_protein_pymol_scene_config(
          self, a_t_result: tuple[str, list[tuple[bool, Any]]]
  ) -> None:
    """Finishes the setup of the color grid and representation section process.

    Args:
        a_t_result (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if a_t_result is None:
      logger.error("a_t_result is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    logger.debug("Running method '__await_setup_protein_pymol_scene_config'.")
    logger.debug(a_t_result)

    try:
      # <editor-fold desc="Color config">
      tmp_success_flag, tmp_raw_residue_color_config = a_t_result[1][0]
      tmp_residue_color_config = residue_color_config.ResidueColorConfig(
        tmp_raw_residue_color_config["data"][0], tmp_raw_residue_color_config["data"][1],
        tmp_raw_residue_color_config["data"][2]
      )

      logger.debug(f"The return_value is: {a_t_result[0]}.")
      if not tmp_success_flag:
        logger.error("Retrieving color information failed!")
        self._interface_manager.status_bar_manager.show_error_message(
          "Retrieving color information failed!"
        )
        self._interface_manager.stop_wait_cursor()
        self._interface_manager.refresh_main_view()
        return

      logger.debug("Check if chain is colored by element.")
      if tmp_residue_color_config.atoms_are_colored_by_elements():
        ui_util.set_checked_async(
          self._view.tg_protein_color_atoms.toggle_button, True
        )
        self._view.ui.lbl_protein_current_color.setText("By Element    ")
      else:
        ui_util.set_checked_async(
          self._view.tg_protein_color_atoms.toggle_button, False
        )
        self._interface_manager.set_current_active_chain_color_of_protein(
          tmp_residue_color_config.carbon_color
        )
        self._view.color_grid_proteins.set_icon_for_selected_color(
          tmp_residue_color_config.carbon_color
        )
        self._view.ui.lbl_protein_current_color.setText(
          f"{tmp_residue_color_config.carbon_color}    "
        )
      # Set representation toggle states for selected chain
      tmp_protein = self._interface_manager.get_current_active_protein_object()
      tmp_chain = self._interface_manager.get_current_active_chain_object()
      tmp_protein.pymol_selection.selection_string = (
        f"first chain {tmp_chain.chain_letter}"
      )

      # </editor-fold>

      # <editor-fold desc="Representation config">
      logger.debug(
        "Returned from method 'pymol_session_async.get_representation_config_of_a_given_protein_chain'."
      )
      tmp_success_flag, tmp_raw_representation_config = a_t_result[1][1]
      tmp_representation_config = tmp_raw_representation_config["data"]
      logger.debug(f"The return_value is: {a_t_result[1]}.")
      if tmp_success_flag:
        tmp_chain = self._interface_manager.get_current_active_chain_object()
        if tmp_chain.chain_type == "protein_chain":
          self._view.ui.frame_protein_repr.setEnabled(True)
          self._interface_manager.manage_toggle_state_of_protein_repr(
            tmp_representation_config
          )
        else:
          self._view.ui.frame_protein_repr.setEnabled(False)
          ui_util.set_checked_async(
            self._view.tg_protein_cartoon.toggle_button, False
          )
          ui_util.set_checked_async(
            self._view.tg_protein_ribbon.toggle_button, False
          )
          ui_util.set_checked_async(
            self._view.tg_protein_sticks.toggle_button, False
          )
          ui_util.set_checked_async(
            self._view.tg_protein_lines.toggle_button, False
          )
          ui_util.set_checked_async(
            self._view.tg_protein_spheres.toggle_button, False
          )
          ui_util.set_checked_async(
            self._view.tg_protein_dots.toggle_button, False
          )
          ui_util.set_checked_async(
            self._view.tg_protein_mesh.toggle_button, False
          )
          ui_util.set_checked_async(
            self._view.tg_protein_surface.toggle_button, False
          )

      # </editor-fold>

    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  # def legacy_setup_protein_pymol_scene_config(self):
  #     try:
  #         tmp_chain_letter = self._view.ui.proteins_tree_view.currentIndex().data(Qt.DisplayRole)
  #         tmp_protein_name = self._view.ui.proteins_tree_view.currentIndex().parent().parent().data(Qt.DisplayRole)
  #         logger.log(
  #             log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
  #             f"The chain object '{tmp_chain_letter}' of the protein '{tmp_protein_name}' on the 'Proteins Tab' was clicked."
  #         )
  #         if self._interface_manager.pymol_session_manager.current_scene_name == "" or not self._interface_manager.pymol_session_manager.is_the_current_protein_in_session(tmp_protein_name):
  #             return
  #         # Set icon for color grid
  #         self._view.color_grid_proteins.reset_icon_for_selected_color()
  #         tmp_protein = self._interface_manager.get_current_active_protein_object()
  #         tmp_chain = self._interface_manager.get_current_active_chain_object()
  #         self._active_task = tasks.LegacyTask(
  #             target=pymol_session_async.get_residue_color_config_of_a_given_protein_chain,
  #             args=(
  #                 tmp_protein.get_molecule_object(),
  #                 tmp_chain.chain_letter,
  #                 self._interface_manager.pymol_session_manager
  #             ),
  #             post_func=self.__await_get_residue_color_config_of_a_given_protein_chain_of_a_protein
  #         )
  #     except Exception as e:
  #         logger.error(f"An error occurred: {e}")
  #         self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
  #     else:
  #         self._interface_manager.block_gui()
  #         self._active_task.start()
  #
  # def __await_get_residue_color_config_of_a_given_protein_chain_of_a_protein(
  #         self,
  #         return_value: tuple[bool, Optional["residue_color_config.ResidueColorConfig"]]
  # ) -> None:
  #     try:
  #         logger.debug("Returned from method 'pymol_session_async.get_residue_color_config_of_a_given_protein_chain'.")
  #         tmp_success_flag, tmp_residue_color_config = return_value
  #         logger.debug(f"The return_value is: {return_value}.")
  #         if not tmp_success_flag:
  #             logger.error("Retrieving color information failed!")
  #             self._interface_manager.status_bar_manager.show_error_message("Retrieving color information failed!")
  #             self._interface_manager.stop_wait_cursor()
  #             self._interface_manager.refresh_main_view()
  #             return
  #
  #         logger.debug("Check if chain is colored by element.")
  #         if tmp_residue_color_config.atoms_are_colored_by_elements():
  #             ui_util.set_checked_async(self._view.tg_protein_color_atoms.toggle_button, True)
  #             self._view.ui.lbl_protein_current_color.setText("By Element    ")
  #         else:
  #             ui_util.set_checked_async(self._view.tg_protein_color_atoms.toggle_button, False)
  #             self._interface_manager.set_current_active_chain_color_of_protein(tmp_residue_color_config.carbon_color)
  #             self._view.color_grid_proteins.set_icon_for_selected_color(tmp_residue_color_config.carbon_color)
  #             self._view.ui.lbl_protein_current_color.setText(f"{tmp_residue_color_config.carbon_color}    ")
  #         # Set representation toggle states for selected chain
  #         tmp_protein = self._interface_manager.get_current_active_protein_object()
  #         tmp_chain = self._interface_manager.get_current_active_chain_object()
  #         tmp_protein.pymol_selection.selection_string = f"first chain {tmp_chain.chain_letter}"
  #         self._active_task = tasks.LegacyTask(
  #             target=pymol_session_async.get_representation_config_of_a_given_protein_chain,
  #             args=(
  #                 tmp_protein.get_molecule_object(),
  #                 tmp_chain.chain_letter,
  #                 self._interface_manager.pymol_session_manager
  #             ),
  #             post_func=self.__await_get_representation_config_of_a_given_protein_chain_of_a_protein
  #         )
  #     except Exception as e:
  #         logger.error(f"An error occurred: {e}")
  #         self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
  #         self._interface_manager.stop_wait_cursor()
  #         self._interface_manager.refresh_main_view()
  #     else:
  #         self._active_task.start()
  #
  # def __await_get_representation_config_of_a_given_protein_chain_of_a_protein(
  #         self, return_value: tuple[bool, Optional[dict]]
  # ) -> None:
  #     try:
  #         logger.debug(
  #             "Returned from method 'pymol_session_async.get_representation_config_of_a_given_protein_chain'.")
  #         tmp_success_flag, tmp_representation_config = return_value
  #         logger.debug(f"The return_value is: {return_value}.")
  #         if tmp_success_flag:
  #             tmp_chain = self._interface_manager.get_current_active_chain_object()
  #             if tmp_chain.chain_type == "protein_chain":
  #                 self._view.ui.frame_protein_repr.setEnabled(True)
  #                 self._interface_manager.manage_toggle_state_of_protein_repr(tmp_representation_config)
  #             else:
  #                 self._view.ui.frame_protein_repr.setEnabled(False)
  #                 ui_util.set_checked_async(self._view.tg_protein_cartoon.toggle_button, False)
  #                 ui_util.set_checked_async(self._view.tg_protein_ribbon.toggle_button, False)
  #                 ui_util.set_checked_async(self._view.tg_protein_sticks.toggle_button, False)
  #                 ui_util.set_checked_async(self._view.tg_protein_lines.toggle_button, False)
  #                 ui_util.set_checked_async(self._view.tg_protein_spheres.toggle_button, False)
  #                 ui_util.set_checked_async(self._view.tg_protein_dots.toggle_button, False)
  #                 ui_util.set_checked_async(self._view.tg_protein_mesh.toggle_button, False)
  #                 ui_util.set_checked_async(self._view.tg_protein_surface.toggle_button, False)
  #                 # self._view.tg_protein_cartoon.toggle_button.setChecked(False)
  #                 # self._view.tg_protein_ribbon.toggle_button.setChecked(False)
  #                 # self._view.tg_protein_sticks.toggle_button.setChecked(False)
  #                 # self._view.tg_protein_lines.toggle_button.setChecked(False)
  #                 # self._view.tg_protein_spheres.toggle_button.setChecked(False)
  #                 # self._view.tg_protein_dots.toggle_button.setChecked(False)
  #                 # self._view.tg_protein_mesh.toggle_button.setChecked(False)
  #                 # self._view.tg_protein_surface.toggle_button.setChecked(False)
  #     except Exception as e:
  #         logger.error(f"An error occurred: {e}")
  #         self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
  #     finally:
  #         self._interface_manager.refresh_main_view()
  #         self._interface_manager.stop_wait_cursor()
  #         QtWidgets.QApplication.restoreOverrideCursor()  # Due to unknown reasons this line is needed to reset the cursor

  # </editor-fold>

  # <editor-fold desc="Change chain color of protein">
  def __slot_change_chain_color_proteins(self, a_color: str) -> None:
    """Changes the 'Color' attribute of a protein chain on the 'Proteins Tab' and colors the protein in User PyMOL.

    Args:
        a_color (str): The new color for the protein chain.

    Notes:
        This slot method is used if a button of the color grid is clicked.

    Returns:
        None
    """
    # <editor-fold desc="Checks">
    if a_color is None or a_color == "":
      logger.error("a_color is either None or an empty string!")
      self._interface_manager.status_bar_manager.show_error_message(
        "a_color is either None or an empty string!"
      )
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      return
    if a_color not in constants.PYMOL_COLORS_WITH_INDICES.values():
      logger.error("a_color is not part of the PYMOL_COLORS_WITH_INDICES dict!")
      self._interface_manager.status_bar_manager.show_error_message(
        "a_color is not part of the PYMOL_COLORS_WITH_INDICES dict!"
      )
      return

    # </editor-fold>

    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "The 'Color' attribute of a protein chain on the 'Proteins Tab' changed.",
      )
      tmp_protein = self._interface_manager.get_current_active_protein_object()
      if not self._interface_manager.pymol_session_manager.is_the_current_protein_in_session(
              tmp_protein.get_molecule_object()
      ):
        logger.warning(
          f"Protein {tmp_protein.get_molecule_object()} is not in the current session."
        )
        self._interface_manager.status_bar_manager.show_error_message(
          f"Protein {tmp_protein.get_molecule_object()} is not in the current session.",
        )
        return

      tmp_chain = self._interface_manager.get_current_active_chain_object()
      tmp_protein.pymol_selection.set_selection_for_a_single_chain(
        tmp_chain.chain_letter
      )
      # Color protein in User PyMOL
      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=pymol_session_async.color_pymol_selection,
              args=(
                a_color,
                tmp_protein.pymol_selection.selection_string,
                self._interface_manager.pymol_session_manager,
              ),
            ),
            an_await_function=self.__await_color_pymol_selection_for_protein,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
      # self._active_task = tasks.LegacyTask(
      #     target=pymol_session_async.color_pymol_selection,
      #     args=(
      #         a_color,
      #         tmp_protein.pymol_selection.selection_string,
      #         self._interface_manager.pymol_session_manager,
      #     ),
      #     post_func=self.__await_color_pymol_selection_for_protein,
      # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()  # fixme: This blocks the entire UI which leads to quick flashes. This might be a problem.

  def __await_color_pymol_selection_for_protein(
          self, return_value: tuple[bool, str]
  ) -> None:
    """Updates the color of the protein chain based on the return value provided in the object and data model.

    Args:
        return_value (tuple): A tuple containing two values, a boolean flag indicating the success of the operation and a string representing the color.
    """
    # <editor-fold desc="Checks">
    if return_value is None or len(return_value) == 0:
      logger.error("return_value is either None or has a length of 0.")
      self._interface_manager.status_bar_manager.show_error_message(
        "return_value is either None or has a length of 0."
      )
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      tmp_success_flag, tmp_chain_color = tmp_result
      if tmp_success_flag:
        tmp_chain = self._interface_manager.get_current_active_chain_object()

        # <editor-fold desc="Setup color grid icon">
        self._view.color_grid_proteins.reset_icon_for_selected_color()
        if tmp_chain_color != "By Element":
          tmp_chain.pymol_parameters["chain_color"] = tmp_chain_color
          self._interface_manager.set_current_active_chain_color_of_protein(
            tmp_chain_color
          )
          self._view.color_grid_proteins.set_icon_for_selected_color(
            tmp_chain_color
          )
          ui_util.set_checked_async(
            self._view.tg_protein_color_atoms.toggle_button, False
          )
          # self._view.tg_protein_color_atoms.toggle_button.setChecked(False)

        # </editor-fold>

        self._view.ui.lbl_protein_current_color.setText(
          f"{tmp_chain_color}    "
        )
        self._update_protein_scene_legacy()
        self._save_protein_pymol_session()
      else:
        logger.error("The operation failed!")
        self._interface_manager.status_bar_manager.show_error_message(
          "The operation failed!"
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()

  # </editor-fold>

  # <editor-fold desc="Color protein chain atoms by element">
  def __slot_color_protein_atoms_by_element(self) -> None:
    """Color protein atoms by element.

    Notes:
        This method is used to color the protein atoms based on their element.
        If the toggle button is checked, the method will color the atoms based on their element:
        grey70 for carbon atoms and atomic color for non-carbon atoms.
        It will also update the protein chain color in the database.
        If the toggle button is not checked, the method will reset the color to the chain color
        or green if the chain color is set to "By Element".

    """
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
      "Toggle 'By Elements' button on the 'Proteins Tab' was clicked.",
    )
    try:
      # Create selection for atom coloring
      tmp_selection = (
        self._interface_manager.get_current_active_protein_object().pymol_selection
      )
      tmp_selection.set_selection_for_a_single_chain(
        self._interface_manager.get_current_active_chain_object().chain_letter
      )

      if self._view.tg_protein_color_atoms.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.color_pymol_selection_atoms_by_element,
                args=(
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_color_pymol_selection_atoms_by_element_for_protein,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.color_pymol_selection_atoms_by_element,
        #     args=(
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_color_pymol_selection_atoms_by_element_for_protein,
        # )
      else:
        if (
                self._interface_manager.get_current_active_chain_color_of_protein()
                is None
        ):
          tmp_current_active_chain_color = "By Element"
        else:
          tmp_current_active_chain_color = (
            self._interface_manager.get_current_active_chain_color_of_protein()
          )

        tmp_protein = (
          self._interface_manager.get_current_active_protein_object()
        )
        tmp_chain = self._interface_manager.get_current_active_chain_object()
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.reset_color_pymol_selection_atoms_by_element,
                args=(
                  tmp_protein.get_molecule_object(),
                  tmp_chain.chain_letter,
                  tmp_current_active_chain_color,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_reset_color_pymol_selection_atoms_by_element_for_protein,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.reset_color_pymol_selection_atoms_by_element,
        #     args=(
        #         tmp_protein.get_molecule_object(),
        #         tmp_chain.chain_letter,
        #         tmp_current_active_chain_color,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_reset_color_pymol_selection_atoms_by_element_for_protein,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __await_color_pymol_selection_atoms_by_element_for_protein(
          self, return_value: tuple
  ) -> None:
    """Await method for the coloring of atom by their element.

    Args:
        return_value (tuple): A tuple containing the result of the operation. The first element of the tuple determines whether the operation was successful or not.
    """
    # <editor-fold desc="Checks">
    if return_value is None or len(return_value) == 0:
      logger.error("return_value is either None or has a length of 0.")
      self._interface_manager.status_bar_manager.show_error_message(
        "return_value is either None or has a length of 0."
      )
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      if tmp_result[0]:
        self._view.color_grid_proteins.reset_icon_for_selected_color()
        self._view.ui.lbl_protein_current_color.setText("By Element    ")
        self._update_protein_scene_legacy()
        self._save_protein_pymol_session()
      else:
        logger.error("The operation failed!")
        self._interface_manager.status_bar_manager.show_error_message(
          "The operation failed!"
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()

  def __await_reset_color_pymol_selection_atoms_by_element_for_protein(
          self, return_value: tuple
  ) -> None:
    """Await method for the reset coloring of atom by their element.

    Args:
        return_value (tuple): A tuple containing the result of the operation and the chain color.
    """
    # <editor-fold desc="Checks">
    if return_value is None or len(return_value) == 0:
      logger.error("return_value is either None or has a length of 0.")
      self._interface_manager.status_bar_manager.show_error_message(
        "return_value is either None or has a length of 0."
      )
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      return

    # </editor-fold>

    tmp_success_flag, tmp_chain_color = return_value
    try:
      if tmp_success_flag:
        self._view.color_grid_proteins.reset_icon_for_selected_color()

        # <editor-fold desc="Update chain color">
        # Updates chain color in chain object and data model
        tmp_chain = self._interface_manager.get_current_active_chain_object()
        tmp_chain.pymol_parameters["chain_color"] = tmp_chain_color
        self._interface_manager.set_current_active_chain_color_of_protein(
          tmp_chain_color
        )

        # </editor-fold>

        self._view.color_grid_proteins.set_icon_for_selected_color(
          tmp_chain_color
        )
        self._view.ui.lbl_protein_current_color.setText(
          f"{tmp_chain_color}    "
        )

        self._update_protein_scene_legacy()
        self._save_protein_pymol_session()
      else:
        logger.error("The operation failed!")
        self._interface_manager.status_bar_manager.show_error_message(
          "The operation failed!"
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()

  # </editor-fold>

  # <editor-fold desc="Set background color">
  def __slot_protein_change_background_color(self) -> None:
    """Sets the background color for the protein pymol session."""
    try:
      if self._view.tg_protein_white_bg.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.set_background_color,
                args=("white", self._interface_manager.pymol_session_manager),
              ),
              an_await_function=self.__await_set_background_color_for_protein_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.set_background_color,
        #     args=("white", self._interface_manager.pymol_session_manager),
        #     post_func=self.__await_set_background_color_for_protein_session,
        # )
      else:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.set_background_color,
                args=("black", self._interface_manager.pymol_session_manager),
              ),
              an_await_function=self.__await_set_background_color_for_protein_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.set_background_color,
        #     args=("black", self._interface_manager.pymol_session_manager),
        #     post_func=self.__await_set_background_color_for_protein_session,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __await_set_background_color_for_protein_session(
          self, return_value: tuple[bool]
  ) -> None:
    """Await method for setting the background color.

    Args:
        return_value (tuple): A tuple containing the result of the operation.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
    try:
      if tmp_result[0]:
        self._interface_manager.status_bar_manager.show_temporary_message(
          "Background color updated."
        )
      else:
        self._interface_manager.status_bar_manager.show_error_message(
          "Updating background color failed!"
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  # </editor-fold>

  # <editor-fold desc="Show/Hide representations">
  def __await_set_representation_for_protein_session(
          self, return_value: tuple[bool]
  ) -> None:
    """Saves the pymol session after changing the representation.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      if tmp_result[0]:
        self._update_protein_scene_legacy()
        self._save_protein_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_chain()
        self._interface_manager.manage_hydrogen_representation_for_protein_chain()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  # def __slot_set_representation_for_protein_chain(self, a_representation: "enums.PyMOLRepresentation"):
  #     try:
  #         logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
  #                    "'Cartoon' toggle on the 'Proteins Tab' was clicked.")
  #         tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
  #         tmp_selection.set_selection_for_a_single_chain(
  #             self._interface_manager.get_current_active_chain_object().chain_letter)
  #
  #         if self._view.tg_protein_cartoon.toggle_button.isChecked():
  #             self._active_task = tasks.Task(
  #                 target=pymol_session_async.show_specific_representation,
  #                 args=(
  #                     a_representation,
  #                     tmp_selection.selection_string,
  #                     self._interface_manager.pymol_session_manager
  #                 ),
  #                 post_func=self.__await_set_representation_for_protein_session
  #             )
  #         else:
  #             self._active_task = tasks.Task(
  #                 target=pymol_session_async.hide_specific_representation,
  #                 args=(
  #                     a_representation,
  #                     tmp_selection.selection_string,
  #                     self._interface_manager.pymol_session_manager
  #                 ),
  #                 post_func=self.__await_set_representation_for_protein_session
  #             )
  #     except Exception as e:
  #         logger.error(f"An error occurred: {e}")
  #         self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
  #     else:
  #         self._interface_manager.start_wait_cursor()
  #         self._active_task.start()

  def __slot_protein_chain_as_cartoon(self) -> None:
    """Changes the representation based on the toggle state."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Cartoon' toggle on the 'Proteins Tab' was clicked.",
      )
      tmp_selection = (
        self._interface_manager.get_current_active_protein_object().pymol_selection
      )
      tmp_selection.set_selection_for_a_single_chain(
        self._interface_manager.get_current_active_chain_object().chain_letter
      )

      if self._view.tg_protein_cartoon.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.show_specific_representation,
                args=(
                  enums.PyMOLRepresentation.CARTOON,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.show_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.CARTOON,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_session,
        # )
      else:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.hide_specific_representation,
                args=(
                  enums.PyMOLRepresentation.CARTOON,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.hide_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.CARTOON,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_session,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __slot_protein_chain_as_sticks(self) -> None:
    """Changes the representation based on the toggle state."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Cartoon' toggle on the 'Proteins Tab' was clicked.",
      )
      tmp_selection = (
        self._interface_manager.get_current_active_protein_object().pymol_selection
      )
      tmp_selection.set_selection_for_a_single_chain(
        self._interface_manager.get_current_active_chain_object().chain_letter
      )

      if self._view.tg_protein_sticks.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.show_specific_representation,
                args=(
                  enums.PyMOLRepresentation.STICKS,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.show_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.STICKS,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_session,
        # )
      else:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.hide_specific_representation,
                args=(
                  enums.PyMOLRepresentation.STICKS,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.hide_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.STICKS,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_session,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __slot_protein_chain_as_ribbon(self) -> None:
    """Changes the representation based on the toggle state."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Cartoon' toggle on the 'Proteins Tab' was clicked.",
      )
      tmp_selection = (
        self._interface_manager.get_current_active_protein_object().pymol_selection
      )
      tmp_selection.set_selection_for_a_single_chain(
        self._interface_manager.get_current_active_chain_object().chain_letter
      )

      if self._view.tg_protein_ribbon.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.show_specific_representation,
                args=(
                  enums.PyMOLRepresentation.RIBBON,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.show_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.RIBBON,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_session,
        # )
      else:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.hide_specific_representation,
                args=(
                  enums.PyMOLRepresentation.RIBBON,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.hide_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.RIBBON,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_session,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __slot_protein_chain_as_lines(self) -> None:
    """Changes the representation based on the toggle state."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Cartoon' toggle on the 'Proteins Tab' was clicked.",
      )
      tmp_selection = (
        self._interface_manager.get_current_active_protein_object().pymol_selection
      )
      tmp_selection.set_selection_for_a_single_chain(
        self._interface_manager.get_current_active_chain_object().chain_letter
      )

      if self._view.tg_protein_lines.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.show_specific_representation,
                args=(
                  enums.PyMOLRepresentation.LINES,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.show_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.LINES,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_session,
        # )
      else:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.hide_specific_representation,
                args=(
                  enums.PyMOLRepresentation.LINES,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.hide_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.LINES,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_session,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __slot_protein_chain_as_spheres(self) -> None:
    """Changes the representation based on the toggle state."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Cartoon' toggle on the 'Proteins Tab' was clicked.",
      )
      tmp_selection = (
        self._interface_manager.get_current_active_protein_object().pymol_selection
      )
      tmp_selection.set_selection_for_a_single_chain(
        self._interface_manager.get_current_active_chain_object().chain_letter
      )

      if self._view.tg_protein_spheres.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.show_specific_representation,
                args=(
                  enums.PyMOLRepresentation.SPHERES,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.show_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.SPHERES,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_session,
        # )
      else:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.hide_specific_representation,
                args=(
                  enums.PyMOLRepresentation.SPHERES,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.hide_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.SPHERES,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_session,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __slot_protein_chain_as_dots(self) -> None:
    """Changes the representation based on the toggle state."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Cartoon' toggle on the 'Proteins Tab' was clicked.",
      )
      tmp_selection = (
        self._interface_manager.get_current_active_protein_object().pymol_selection
      )
      tmp_selection.set_selection_for_a_single_chain(
        self._interface_manager.get_current_active_chain_object().chain_letter
      )

      if self._view.tg_protein_dots.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.show_specific_representation,
                args=(
                  enums.PyMOLRepresentation.DOTS,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.show_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.DOTS,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_session,
        # )
      else:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.hide_specific_representation,
                args=(
                  enums.PyMOLRepresentation.DOTS,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.hide_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.DOTS,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_session,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __slot_protein_chain_as_mesh(self) -> None:
    """Changes the representation based on the toggle state."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Cartoon' toggle on the 'Proteins Tab' was clicked.",
      )
      tmp_selection = (
        self._interface_manager.get_current_active_protein_object().pymol_selection
      )
      tmp_selection.set_selection_for_a_single_chain(
        self._interface_manager.get_current_active_chain_object().chain_letter
      )

      if self._view.tg_protein_mesh.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.show_specific_representation,
                args=(
                  enums.PyMOLRepresentation.MESH,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.show_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.MESH,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_session,
        # )
      else:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.hide_specific_representation,
                args=(
                  enums.PyMOLRepresentation.MESH,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.hide_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.MESH,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_session,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __slot_protein_chain_as_surface(self) -> None:
    """Changes the representation based on the toggle state."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Cartoon' toggle on the 'Proteins Tab' was clicked.",
      )
      tmp_selection = (
        self._interface_manager.get_current_active_protein_object().pymol_selection
      )
      tmp_selection.set_selection_for_a_single_chain(
        self._interface_manager.get_current_active_chain_object().chain_letter
      )

      if self._view.tg_protein_surface.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.show_specific_representation,
                args=(
                  enums.PyMOLRepresentation.SURFACE,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.show_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.SURFACE,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_session,
        # )
      else:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.hide_specific_representation,
                args=(
                  enums.PyMOLRepresentation.SURFACE,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.hide_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.SURFACE,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_session,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  # def __slot_protein_chain_as_sticks(self):
  #     try:
  #         logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
  #                    "'Sticks' toggle on the 'Proteins Tab' was clicked.")
  #         tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
  #         tmp_selection.set_selection_for_a_single_chain(
  #             self._interface_manager.get_current_active_chain_object().chain_letter)
  #         if self._view.ui.cb_protein_sticks.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.STICKS.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.STICKS.value, tmp_selection.selection_string
  #             )
  #         elif self._view.tg_protein_sticks.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.STICKS.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.STICKS.value, tmp_selection.selection_string
  #             )
  #         else:
  #             self._interface_manager.pymol_session_manager.hide_specific_representation(
  #                 enums.PyMOLRepresentation.STICKS.value, tmp_selection.selection_string
  #             )
  #         #self.__slot_chain_protein_with_hydrogens()
  #         self._update_protein_scene_legacy()
  #         self._save_protein_pymol_session()
  #         self._interface_manager.manage_coloring_by_element_option_for_protein_chain()
  #         self._interface_manager.manage_hydrogen_representation_for_protein_chain()
  #     except Exception as e:
  #         logger.error(f"An error occurred: {e}")
  #         self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
  #
  # def __slot_protein_chain_as_ribbon(self):
  #     try:
  #         logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
  #                    "'Ribbon' toggle on the 'Proteins Tab' was clicked.")
  #         tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
  #         tmp_selection.set_selection_for_a_single_chain(
  #             self._interface_manager.get_current_active_chain_object().chain_letter)
  #         if self._view.ui.cb_protein_ribbon.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.RIBBON.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.RIBBON.value, tmp_selection.selection_string
  #             )
  #         elif self._view.tg_protein_ribbon.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.RIBBON.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.RIBBON.value, tmp_selection.selection_string
  #             )
  #         else:
  #             #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.RIBBON.value)
  #             self._interface_manager.pymol_session_manager.hide_specific_representation(
  #                 enums.PyMOLRepresentation.RIBBON.value, tmp_selection.selection_string
  #             )
  #         self._update_protein_scene_legacy()
  #         self._save_protein_pymol_session()
  #         self._interface_manager.manage_coloring_by_element_option_for_protein_chain()
  #         self._interface_manager.manage_hydrogen_representation_for_protein_chain()
  #     except Exception as e:
  #         logger.error(f"An error occurred: {e}")
  #         self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
  #
  # def __slot_protein_chain_as_lines(self):
  #     try:
  #         logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
  #                    "'Lines' toggle on the 'Proteins Tab' was clicked.")
  #         tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
  #         tmp_selection.set_selection_for_a_single_chain(
  #             self._interface_manager.get_current_active_chain_object().chain_letter)
  #         if self._view.ui.cb_protein_lines.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.LINES.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.LINES.value, tmp_selection.selection_string
  #             )
  #         elif self._view.tg_protein_lines.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.LINES.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.LINES.value, tmp_selection.selection_string
  #             )
  #         else:
  #             #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.LINES.value)
  #             self._interface_manager.pymol_session_manager.hide_specific_representation(
  #                 enums.PyMOLRepresentation.LINES.value, tmp_selection.selection_string
  #             )
  #         #self.__slot_chain_protein_with_hydrogens()
  #         self._update_protein_scene_legacy()
  #         self._save_protein_pymol_session()
  #         self._interface_manager.manage_coloring_by_element_option_for_protein_chain()
  #         self._interface_manager.manage_hydrogen_representation_for_protein_chain()
  #     except Exception as e:
  #         logger.error(f"An error occurred: {e}")
  #         self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
  #
  # def __slot_protein_chain_as_spheres(self):
  #     try:
  #         logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
  #                    "'Spheres' toggle on the 'Proteins Tab' was clicked.")
  #         tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
  #         tmp_selection.set_selection_for_a_single_chain(
  #             self._interface_manager.get_current_active_chain_object().chain_letter)
  #         if self._view.ui.cb_protein_spheres.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SPHERES.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.SPHERES.value, tmp_selection.selection_string
  #             )
  #         elif self._view.tg_protein_spheres.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SPHERES.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.SPHERES.value, tmp_selection.selection_string
  #             )
  #         else:
  #             #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.SPHERES.value)
  #             self._interface_manager.pymol_session_manager.hide_specific_representation(
  #                 enums.PyMOLRepresentation.SPHERES.value, tmp_selection.selection_string
  #             )
  #         #self.__slot_chain_protein_with_hydrogens()
  #         self._update_protein_scene_legacy()
  #         self._save_protein_pymol_session()
  #         self._interface_manager.manage_coloring_by_element_option_for_protein_chain()
  #         self._interface_manager.manage_hydrogen_representation_for_protein_chain()
  #     except Exception as e:
  #         logger.error(f"An error occurred: {e}")
  #         self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
  #
  # def __slot_protein_chain_as_dots(self):
  #     try:
  #         logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
  #                    "'Dots' toggle on the 'Proteins Tab' was clicked.")
  #         tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
  #         tmp_selection.set_selection_for_a_single_chain(
  #             self._interface_manager.get_current_active_chain_object().chain_letter)
  #         if self._view.ui.cb_protein_dots.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.DOTS.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.DOTS.value, tmp_selection.selection_string
  #             )
  #         elif self._view.tg_protein_dots.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.DOTS.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.DOTS.value, tmp_selection.selection_string
  #             )
  #         else:
  #             #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.DOTS.value)
  #             self._interface_manager.pymol_session_manager.hide_specific_representation(
  #                 enums.PyMOLRepresentation.DOTS.value, tmp_selection.selection_string
  #             )
  #         self._update_protein_scene_legacy()
  #         self._save_protein_pymol_session()
  #         self._interface_manager.manage_coloring_by_element_option_for_protein_chain()
  #         self._interface_manager.manage_hydrogen_representation_for_protein_chain()
  #     except Exception as e:
  #         logger.error(f"An error occurred: {e}")
  #         self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
  #
  # def __slot_protein_chain_as_mesh(self):
  #     try:
  #         logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
  #                    "'Mesh' toggle on the 'Proteins Tab' was clicked.")
  #         tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
  #         tmp_selection.set_selection_for_a_single_chain(
  #             self._interface_manager.get_current_active_chain_object().chain_letter)
  #         if self._view.ui.cb_protein_mesh.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.MESH.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.MESH.value, tmp_selection.selection_string
  #             )
  #         elif self._view.tg_protein_mesh.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.MESH.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.MESH.value, tmp_selection.selection_string
  #             )
  #         else:
  #             #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.MESH.value)
  #             self._interface_manager.pymol_session_manager.hide_specific_representation(
  #                 enums.PyMOLRepresentation.MESH.value, tmp_selection.selection_string
  #             )
  #         self._update_protein_scene_legacy()
  #         self._save_protein_pymol_session()
  #         self._interface_manager.manage_coloring_by_element_option_for_protein_chain()
  #         self._interface_manager.manage_hydrogen_representation_for_protein_chain()
  #     except Exception as e:
  #         logger.error(f"An error occurred: {e}")
  #         self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
  #
  # def __slot_protein_chain_as_surface(self):
  #     try:
  #         logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
  #                    "'Surface' toggle on the 'Proteins Tab' was clicked.")
  #         tmp_selection = self._interface_manager.get_current_active_protein_object().pymol_selection
  #         tmp_selection.set_selection_for_a_single_chain(
  #             self._interface_manager.get_current_active_chain_object().chain_letter)
  #         if self._view.ui.cb_protein_surface.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 0:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SURFACE.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.SURFACE.value, tmp_selection.selection_string
  #             )
  #         elif self._view.tg_protein_surface.toggle_button.isChecked() and self._interface_manager.get_protein_repr_toggle_flag() == 1:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SURFACE.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.SURFACE.value, tmp_selection.selection_string
  #             )
  #         else:
  #             #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.SURFACE.value)
  #             self._interface_manager.pymol_session_manager.hide_specific_representation(
  #                 enums.PyMOLRepresentation.SURFACE.value, tmp_selection.selection_string
  #             )
  #         self._update_protein_scene_legacy()
  #         self._save_protein_pymol_session()
  #         self._interface_manager.manage_coloring_by_element_option_for_protein_chain()
  #         self._interface_manager.manage_hydrogen_representation_for_protein_chain()
  #     except Exception as e:
  #         logger.error(f"An error occurred: {e}")
  #         self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

  def __slot_hide_protein_chain_all(self) -> None:
    """Hides all representations for the selected chain."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Hide all' representations button on the 'Proteins Tab' was clicked.",
      )
      tmp_selection = (
        self._interface_manager.get_current_active_protein_object().pymol_selection
      )
      tmp_selection.set_selection_for_a_single_chain(
        self._interface_manager.get_current_active_chain_object().chain_letter
      )
      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=pymol_session_async.hide_all_representations,
              args=(
                tmp_selection.selection_string,
                self._interface_manager.pymol_session_manager,
              ),
            ),
            an_await_function=self.__await_hide_all_representations_of_protein_chain_of_a_protein,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
      # self._active_task = tasks.LegacyTask(
      #     target=pymol_session_async.hide_all_representations,
      #     args=(
      #         tmp_selection.selection_string,
      #         self._interface_manager.pymol_session_manager,
      #     ),
      #     post_func=self.__await_hide_all_representations_of_protein_chain_of_a_protein,
      # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __await_hide_all_representations_of_protein_chain_of_a_protein(
          self, return_value: tuple[bool]
  ) -> None:
    """Finishes hide all representations for the selected chain process.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      QtWidgets.QApplication.restoreOverrideCursor()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      if tmp_result[0]:
        ui_util.set_checked_async(
          self._view.tg_protein_cartoon.toggle_button, False
        )
        ui_util.set_checked_async(
          self._view.tg_protein_ribbon.toggle_button, False
        )
        ui_util.set_checked_async(
          self._view.tg_protein_sticks.toggle_button, False
        )
        ui_util.set_checked_async(
          self._view.tg_protein_lines.toggle_button, False
        )
        ui_util.set_checked_async(
          self._view.tg_protein_spheres.toggle_button, False
        )
        ui_util.set_checked_async(
          self._view.tg_protein_dots.toggle_button, False
        )
        ui_util.set_checked_async(
          self._view.tg_protein_mesh.toggle_button, False
        )
        ui_util.set_checked_async(
          self._view.tg_protein_surface.toggle_button, False
        )
        # self._view.tg_protein_cartoon.toggle_button.setChecked(False)
        # self._view.tg_protein_sticks.toggle_button.setChecked(False)
        # self._view.tg_protein_ribbon.toggle_button.setChecked(False)
        # self._view.tg_protein_lines.toggle_button.setChecked(False)
        # self._view.tg_protein_spheres.toggle_button.setChecked(False)
        # self._view.tg_protein_dots.toggle_button.setChecked(False)
        # self._view.tg_protein_mesh.toggle_button.setChecked(False)
        # self._view.tg_protein_surface.toggle_button.setChecked(False)
        self._update_protein_scene_legacy()
        self._save_protein_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_chain()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      QtWidgets.QApplication.restoreOverrideCursor()

  # </editor-fold>

  # <editor-fold desc="Import protein structure">
  def __slot_import_protein_structure(self) -> None:
    """Shows the import protein dialog."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Import protein' button on the 'Proteins Tab' was clicked.",
      )
      self._external_controller = (
        add_protein_view_controller.AddProteinViewController(
          self._interface_manager
        )
      )
      self._external_controller.user_input.connect(
        self._post_import_protein_structure
      )
      self._external_controller.restore_ui()
      self._interface_manager.get_add_protein_view().show()
    except Exception as e:
      logger.error(f"An error occurred during the protein import: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def _post_import_protein_structure(self, return_value: tuple) -> None:
    """Starts the import process with the user's inputs.

    Args:
        return_value (tuple): The result data from the dialog signal.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      return

    # </editor-fold>

    try:
      tmp_protein_name, tmp_name_len = return_value
      if tmp_name_len == 4:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=project_async.add_protein_from_pdb_to_project,
                args=(
                  tmp_protein_name,
                  self._interface_manager,
                ),
              ),
              an_await_function=self.__await_post_import_protein_structure,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=project_async.add_protein_from_pdb_to_project,
        #     args=(
        #         tmp_protein_name,
        #         self._interface_manager,
        #     ),
        #     post_func=self.__await_post_import_protein_structure,
        # )
      elif tmp_name_len > 0:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=project_async.add_protein_from_local_filesystem_to_project,
                args=(
                  tmp_protein_name,
                  self._interface_manager,
                ),
              ),
              an_await_function=self.__await_post_import_protein_structure,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=project_async.add_protein_from_local_filesystem_to_project,
        #     args=(
        #         tmp_protein_name,
        #         self._interface_manager,
        #     ),
        #     post_func=self.__await_post_import_protein_structure,
        # )
      else:
        logger.warning("No protein object was created.")
        return
      self._interface_manager.status_bar_manager.show_temporary_message(
        "Importing protein structure ...",
        False,
      )
    except Exception as e:
      logger.error(f"An error occurred during the protein import: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "Protein import failed!"
      )
    else:
      self._interface_manager.block_gui(with_wait_cursor=True)

  def __await_post_import_protein_structure(self, return_value: tuple) -> None:
    """Finishes the import protein process.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      return
    if return_value[0] == "":
      self._interface_manager.status_bar_manager.show_error_message(
        "Protein import failed!"
      )
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      return
    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      tmp_protein: "protein.Protein" = tmp_result[1]
      self._interface_manager.get_current_project().add_existing_protein(
        tmp_protein
      )
      self._interface_manager.watcher.add_protein(
        tmp_protein.get_molecule_object()
      )
      self._database_thread.put_database_operation_into_queue(
        database_operation.DatabaseOperation(
          enums.SQLQueryType.INSERT_NEW_PROTEIN, (0, tmp_protein)
        )
      )
      self._interface_manager.refresh_main_view()
      self._interface_manager.status_bar_manager.show_temporary_message(
        "Importing protein structure finished."
      )
    except Exception as e:
      logger.error(f"An error occurred during the protein import: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "Protein import failed!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()

  # </editor-fold>

  # <editor-fold desc="Delete protein">
  def __slot_delete_protein(self) -> None:
    """Deletes the selected protein."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Delete protein' button on the 'Proteins Tab' was clicked.",
      )
      tmp_dialog = custom_message_box.CustomMessageBoxDelete(
        "Are you sure you want to delete this protein?",
        "Delete Protein",
        custom_message_box.CustomMessageBoxIcons.WARNING.value,
      )
      tmp_dialog.exec_()
      response: bool = tmp_dialog.response
      if not response:
        return

      tmp_protein: "protein.Protein" = (
        self._interface_manager.get_current_active_protein_object()
      )
      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=pymol_session_async.reinitialize_session,
              args=(
                self._interface_manager.pymol_session_manager,
                tmp_protein.get_molecule_object(),
              ),
            ),
            an_await_function=self.__await_reinitialize_session_before_delete_protein,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
      # self._active_task = tasks.LegacyTask(
      #     target=pymol_session_async.reinitialize_session,
      #     args=(
      #         self._interface_manager.pymol_session_manager,
      #         tmp_protein.get_molecule_object(),
      #     ),
      #     post_func=self.__await_reinitialize_session_before_delete_protein,
      # )
    except Exception as e:
      logger.error(
        f"An error occurred during the protein deletion process: {e}"
      )
      self._interface_manager.status_bar_manager.show_error_message(
        "Protein delete failed!"
      )
    else:
      self._interface_manager.block_gui()

  def __await_reinitialize_session_before_delete_protein(
          self, return_value: tuple[bool]
  ) -> None:
    """Removes protein from database and refreshes the main view.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      if tmp_result[0]:
        tmp_protein: "protein.Protein" = (
          self._interface_manager.get_current_active_protein_object()
        )
        tmp_database_operation = database_operation.DatabaseOperation(
          enums.SQLQueryType.DELETE_EXISTING_PROTEIN,
          (0, tmp_protein.get_id()),
        )
        self._database_thread.put_database_operation_into_queue(
          tmp_database_operation
        )
        self._interface_manager.get_current_project().delete_specific_protein(
          tmp_protein.get_molecule_object()
        )
        self._interface_manager.watcher.remove_protein(
          tmp_protein.get_molecule_object()
        )
        self._interface_manager.remove_protein_from_proteins_model()
        self._interface_manager.status_bar_manager.show_temporary_message(
          "The protein was successfully deleted."
        )
    except Exception as e:
      logger.error(
        f"An error occurred during the protein deletion process: {e}"
      )
      self._interface_manager.status_bar_manager.show_error_message(
        "Protein delete failed!"
      )
    finally:
      self._view.ui.proteins_tree_view.selectionModel().clearSelection()
      self._interface_manager.disable_proteins_tab_buttons()
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  # </editor-fold>

  # <editor-fold desc="Save/Export protein">
  def __slot_save_selected_protein_structure_as_pdb_file(self) -> None:
    """Saves selected protein as pdb file."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Export protein' button on the 'Proteins Tab' was clicked.",
      )
      file_dialog = QtWidgets.QFileDialog()
      desktop_path = QtCore.QStandardPaths.standardLocations(
        QtCore.QStandardPaths.DesktopLocation
      )[0]
      file_dialog.setDirectory(desktop_path)
      file_path, _ = file_dialog.getSaveFileName(
        self._view,
        "Save Protein Structure",
        "",
        "Protein Data Bank File (*.pdb)",
      )
      if file_path:
        tmp_protein: "protein.Protein" = (
          self._interface_manager.get_current_protein_tree_index_object()
        )
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=protein_async.save_selected_protein_structure_as_pdb_file,
                args=(
                  tmp_protein,
                  file_path,
                  self._interface_manager.get_current_project().get_database_filepath(),
                ),
              ),
              an_await_function=self.__await_save_selected_protein_structure_as_pdb_file,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=protein_async.save_selected_protein_structure_as_pdb_file,
        #     args=(
        #         tmp_protein,
        #         file_path,
        #         self._interface_manager.get_current_project().get_database_filepath(),
        #     ),
        #     post_func=self.__await_save_selected_protein_structure_as_pdb_file,
        # )
      else:
        self._interface_manager.stop_wait_cursor()
        self._interface_manager.refresh_main_view()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui(with_wait_cursor=True)

  def __await_save_selected_protein_structure_as_pdb_file(
          self, result: tuple
  ) -> None:
    """Finishes the saving protein process.

    Args:
        result (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if result is None:
      logger.error("result is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      return
    tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(result)
    if tmp_result[0] == "":
      tmp_dialog = custom_message_box.CustomMessageBoxOk(
        "Saving the protein as .pdb file failed!",
        "Save Protein Structure",
        custom_message_box.CustomMessageBoxIcons.DANGEROUS.value,
      )
      tmp_dialog.exec_()
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.status_bar_manager.show_error_message(
        "Saving the protein as .pdb file failed!"
      )
      return

    # </editor-fold>

    try:
      self._interface_manager.status_bar_manager.show_temporary_message(
        "The protein was successfully saved as .pdb file."
      )
      self._interface_manager.refresh_protein_model()
      self._interface_manager.refresh_main_view()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()

  # </editor-fold>

  # <editor-fold desc="Clean protein">
  def __slot_clean_protein_update(self) -> None:
    """Cleans the selected protein structure."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Clean protein' context menu action was clicked.",
      )
      tmp_dialog = custom_message_box.CustomMessageBoxYesNo(
        "Are you sure you want to clean this protein?\n"
        "This will remove all organic and solvent components!",
        "Clean Protein",
        custom_message_box.CustomMessageBoxIcons.WARNING.value,
      )
      tmp_dialog.exec_()
      if tmp_dialog.response:
        tmp_main_socket, tmp_general_purpose_socket = (
          self._interface_manager.job_manager.get_general_purpose_socket_pair()
        )
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=protein_async.clean_protein_update,
                args=(
                  self._interface_manager.get_current_active_protein_object(),
                  self._interface_manager.get_current_project().get_database_filepath(),
                  tmp_main_socket,
                  tmp_general_purpose_socket,
                ),
              ),
              an_await_function=self.__await_clean_protein_update,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=protein_async.clean_protein_update,
        #     args=(
        #         self._interface_manager.get_current_active_protein_object(),
        #         self._interface_manager.get_current_project().get_database_filepath(),
        #         tmp_main_socket,
        #         tmp_general_purpose_socket,
        #     ),
        #     post_func=self.__await_clean_protein_update,
        # )
        self._interface_manager.block_gui()
        self.update_status("Cleaning protein ...")
      else:
        constants.PYSSA_LOGGER.info("No protein has been cleaned.")
        self._interface_manager.stop_wait_cursor()
        self._interface_manager.refresh_main_view()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  def __await_clean_protein_update(self, return_value: tuple) -> None:
    """Finishes the cleaning process.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return
    tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
    if tmp_result[0] == "":
      self._interface_manager.status_bar_manager.show_error_message(
        "Cleaning protein failed!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      if not self._interface_manager.pymol_session_manager.is_the_current_protein_in_session(
              self._interface_manager.get_current_active_protein_object().get_molecule_object(),
      ):
        self._interface_manager.stop_wait_cursor()
        self._interface_manager.refresh_main_view()
        return

      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=pymol_session_async.load_protein_pymol_session,
              args=(
                self._interface_manager.get_current_active_protein_object(),
                self._interface_manager.pymol_session_manager,
                False,
              ),
            ),
            an_await_function=self.__await_load_protein_pymol_session_after_cleaning,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
      # self._active_task = tasks.LegacyTask(
      #     target=pymol_session_async.load_protein_pymol_session,
      #     args=(
      #         self._interface_manager.get_current_active_protein_object(),
      #         self._interface_manager.pymol_session_manager,
      #         False,
      #     ),
      #     post_func=self.__await_load_protein_pymol_session_after_cleaning,
      # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  def __await_load_protein_pymol_session_after_cleaning(
          self,
          return_value: tuple[
            Optional["pymol_session_manager.PymolSessionManager"], bool
          ],
  ) -> None:
    """Loads the cleaned version back in PyMOL if it was in the session before cleaning.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      if tmp_result[1]:
        self._interface_manager.status_bar_manager.show_temporary_message(
          "Cleaning the protein finished."
        )
      else:
        self._interface_manager.status_bar_manager.show_error_message(
          "Cleaning the protein failed!"
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  # </editor-fold>

  # <editor-fold desc="Rename (not yet implemented properly!)">
  # TODO: This needs to be better implemented before used in production
  def __slot_rename_selected_protein_structure(self) -> None:
    """Opens a new view to rename the selected protein."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Rename protein' context menu action was clicked.",
      )
      self._external_controller = (
        rename_protein_view_controller.RenameProteinViewController(
          self._interface_manager
        )
      )
      self._external_controller.user_input.connect(
        self.post_rename_selected_protein_structure
      )
      self._external_controller.restore_ui()
      self._interface_manager.get_rename_protein_view().show()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def post_rename_selected_protein_structure(self, return_value: tuple) -> None:
    """Renames a selected protein structure."""
    try:
      if return_value[1] is True:
        self._active_task = tasks.LegacyTask(
          target=protein_async.rename_selected_protein_structure,
          args=(
            self._interface_manager.get_current_protein_tree_index_object(),
            return_value[0],
            self._interface_manager.get_current_project().get_database_filepath(),
          ),
          post_func=self.__await_post_rename_selected_protein_structure,
        )
        self._interface_manager.block_gui(with_wait_cursor=True)
        self.update_status("Renaming protein ...")
      else:
        pass
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def __await_post_rename_selected_protein_structure(
          self, result: tuple
  ) -> None:
    try:
      self._view.ui.proteins_tree_view.model().setData(
        self._interface_manager.get_current_protein_tree_index(),
        result[1],
        enums.ModelEnum.OBJECT_ROLE,
      )
      tmp_database_operation = database_operation.DatabaseOperation(
        enums.SQLQueryType.UPDATE_PYMOL_SESSION_PROTEIN,
        (
          0,
          self._view.ui.proteins_tree_view.model().data(
            self._interface_manager.get_current_protein_tree_index(),
            enums.ModelEnum.OBJECT_ROLE,
          ),
        ),
      )
      self._database_thread.put_database_operation_into_queue(
        tmp_database_operation
      )
      self._interface_manager.refresh_protein_model()
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      self.update_status("Renaming protein finished.")
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  # </editor-fold>

  def __slot_show_protein_chain_sequence(self) -> None:
    """Shows a QTextBrowser that displays the sequence of the selected chain."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Show protein sequence' context menu action was clicked.",
      )
      self.tmp_txt_browser = QtWidgets.QTextBrowser()
      try:
        tmp_chain: "chain.Chain" = (
          self._interface_manager.get_current_protein_tree_index_object()
        )
        if tmp_chain.chain_sequence.sequence == "":
          self.tmp_txt_browser.setText(
            "This chain is a non-protein chain.",
          )
        else:
          self.tmp_txt_browser.setText(
            tmp_chain.chain_sequence.sequence,
          )
      except AttributeError:
        return
      else:
        self.tmp_txt_browser.setWindowTitle("View Protein Sequence")
        self.tmp_txt_browser.setWindowIcon(
          QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH)
        )
        self.tmp_txt_browser.resize(500, 150)
        self.tmp_txt_browser.show()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  # <editor-fold desc="Update scene">
  def __slot_update_protein_scene(self) -> None:
    """Updates the current protein scene."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Update protein scene' button on the 'Proteins Tab' was clicked.",
      )
      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=pymol_session_async.update_scene,
              args=(
                self._interface_manager.pymol_session_manager, 0
              ),
            ),
            an_await_function=self.__await_update_scene_for_protein_session,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
      # self._active_task = tasks.LegacyTask(
      #     target=pymol_session_async.update_scene,
      #     args=(self._interface_manager.pymol_session_manager, 0),
      #     post_func=self.__await_update_scene_for_protein_session,
      # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __await_update_scene_for_protein_session(
          self, return_value: tuple[bool, str]
  ) -> None:
    """Finishes the update protein scene process.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      tmp_success_flag, tmp_current_scene_name = tmp_result
      if tmp_current_scene_name == "_scratch_":
        if (
                not self._interface_manager.check_if_scratch_scene_exists_in_protein_model()
        ):
          self._interface_manager.add_scratch_scene_to_protein_model()
        ui_util.set_pymol_scene_name_into_label(
          tmp_current_scene_name, self._view.ui.lbl_pymol_protein_scene
        )
        self._interface_manager.pymol_session_manager.current_scene_name = (
          tmp_current_scene_name
        )
        self._interface_manager.status_bar_manager.show_temporary_message(
          "PyMOL Scene _scratch_ updated.", a_timeout=1500
        )
      else:
        self._interface_manager.status_bar_manager.show_temporary_message(
          "PyMOL Scene updated.", a_timeout=1500
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._save_protein_pymol_session()
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  # </editor-fold>

  def _update_protein_scene_legacy(self) -> None:
    """Updates the current selected PyMOL scene.

    Notes:
        This is a legacy method that is used in many other methods!
        TODO: after discussion of workaround, this can be removed!
    """
    return
    raise NotImplementedError()
    # The code below can be useful if the scene saving feature needs to be implemented
    # if self._interface_manager.pymol_session_manager.is_the_current_pymol_scene_base is False:
    #     self._interface_manager.pymol_session_manager.user_pymol_connector.scene(
    #         a_key="auto", an_action="update"
    #     )
    #     self._interface_manager.status_bar_manager.show_temporary_message("PyMOL Scene updated.", a_timeout=1500)
    # else:
    #     if not self._interface_manager.check_if_scratch_scene_exists_in_protein_model():
    #         self._interface_manager.add_scratch_scene_to_protein_model()
    #     self._interface_manager.pymol_session_manager.user_pymol_connector.scene(
    #         a_key="_scratch_", an_action="update"
    #     )
    #     self._interface_manager.pymol_session_manager.current_scene_name = "_scratch_"
    #     ui_util.set_pymol_scene_name_into_label(self._interface_manager.pymol_session_manager.current_scene_name,
    #                                             self._view.ui.lbl_pymol_protein_scene)
    #     self._interface_manager.status_bar_manager.show_temporary_message("PyMOL Scene _scratch_ updated.", a_timeout=1500)

  # <editor-fold desc="Create/Save new scene">
  def __slot_save_scene(self) -> None:
    """Saves the current view as a new PyMOL scene."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Create pymol scene' button on the 'Proteins or Protein Pairs Tab' was clicked.",
      )
      self._external_controller = (
        add_scene_view_controller.AddSceneViewController(
          self._interface_manager
        )
      )
      self._external_controller.user_input.connect(self._post_save_scene)
      self._external_controller.restore_ui()
      self._interface_manager.get_add_scene_view().show()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def _post_save_scene(self, return_value: tuple) -> None:
    """Starts async method to create new scene in PyMOL.

    Args:
        return_value (tuple): The result data from the dialog signal.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      return

    # </editor-fold>

    try:
      tmp_scene_name, _ = return_value
      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=pymol_session_async.create_new_scene,
              args=(
                tmp_scene_name, self._interface_manager.pymol_session_manager
              ),
            ),
            an_await_function=self.__await_create_new_scene_for_protein_session,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
      # self._active_task = tasks.LegacyTask(
      #     target=pymol_session_async.create_new_scene,
      #     args=(tmp_scene_name, self._interface_manager.pymol_session_manager),
      #     post_func=self.__await_create_new_scene_for_protein_session,
      # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __await_create_new_scene_for_protein_session(
          self, return_value: tuple[bool, str]
  ) -> None:
    """Finishes the create new scene process and starts async method to save PyMOL session.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      logging.debug(
        "Returned from method 'pymol_session_async.create_new_scene'."
      )
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      tmp_success_flag, tmp_created_scene_name = tmp_result
      if self._interface_manager.current_tab_index == 1:
        # User is on Proteins tab
        logging.debug(
          "Setting up new async task 'pymol_session_async.save_protein_pymol_session_to_database'."
        )
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.save_protein_pymol_session_to_database,
                args=(
                  self._interface_manager,
                  0,
                ),
              ),
              an_await_function=self.__await_save_scene_protein,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        self._interface_manager.status_bar_manager.show_temporary_message(
          "Adding new scene to protein ...", False
        )
        logger.debug("Adding the new scene to the protein model.")
        self._interface_manager.add_scene_to_proteins_model(
          tmp_created_scene_name
        )
        self._interface_manager.pymol_session_manager.current_scene_name = (
          tmp_created_scene_name
        )
        logger.debug("Modifying the GUI for the new scene.")
        ui_util.set_pymol_scene_name_into_label(
          self._interface_manager.pymol_session_manager.current_scene_name,
          self._view.ui.lbl_pymol_protein_scene,
        )
      elif self._interface_manager.current_tab_index == 2:
        # User is on Protein Pairs tab
        # The database thread cannot be used here because the session gets loaded again
        # before the new data is in the db
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.save_protein_pair_pymol_session_to_database,
                args=(
                  self._interface_manager,
                  0,
                ),
              ),
              an_await_function=self.__await_save_scene_protein_pair,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        self._interface_manager.status_bar_manager.show_temporary_message(
          "Adding new scene to protein pair ...", False
        )
        self._interface_manager.add_scene_to_protein_pairs_model(
          tmp_created_scene_name
        )
        self._interface_manager.pymol_session_manager.current_scene_name = (
          tmp_created_scene_name
        )
        ui_util.set_pymol_scene_name_into_label(
          self._interface_manager.pymol_session_manager.current_scene_name,
          self._view.ui.lbl_pymol_protein_pair_scene,
        )
      else:
        logger.warning(
          "The current tab index is not for the proteins nor for the protein pairs tab?!"
        )
        self._interface_manager.stop_wait_cursor()
        self._interface_manager.refresh_main_view()
        self._interface_manager.status_bar_manager.show_error_message(
          "Creating scene failed! Tab index invalid!"
        )
        return
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
    else:
      logger.debug(
        "Start async task 'pymol_session_async.save_protein_pymol_session_to_database'."
      )

  def __await_save_scene_protein(self, return_value: tuple) -> None:
    """Finishes the save protein scene process.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      _, exit_flag = tmp_result
      if exit_flag:
        self._interface_manager.status_bar_manager.show_temporary_message(
          "Adding new scene to protein finished."
        )
      else:
        self._interface_manager.status_bar_manager.show_error_message(
          "Adding new scene to protein failed!"
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  def __await_save_scene_protein_pair(self, return_value: tuple) -> None:
    """Finishes the save protein pair scene process.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      _, exit_flag = tmp_result
      if exit_flag:
        self._interface_manager.status_bar_manager.show_temporary_message(
          "Adding new scene to protein pair finished."
        )
      else:
        self._interface_manager.status_bar_manager.show_error_message(
          "Adding new scene to protein pair failed!"
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  # </editor-fold>

  # <editor-fold desc="Delete scene">
  def __slot_delete_current_scene(self) -> None:
    """Starts async method to delete the selected scene."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Delete pymol scene' button on the 'Proteins or Protein Pairs Tab' was clicked.",
      )
      tmp_dialog = custom_message_box.CustomMessageBoxDelete(
        "Are you sure you want to delete this scene?",
        "Delete PyMOL Scene",
        custom_message_box.CustomMessageBoxIcons.WARNING.value,
      )
      tmp_dialog.exec_()
      if not tmp_dialog.response:
        return

      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=pymol_session_async.delete_scene,
              args=(
                self._interface_manager.pymol_session_manager, 0
              ),
            ),
            an_await_function=self.__await_delete_scene_for_protein_session,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
      # self._active_task = tasks.LegacyTask(
      #     target=pymol_session_async.delete_scene,
      #     args=(self._interface_manager.pymol_session_manager, 0),
      #     post_func=self.__await_delete_scene_for_protein_session,
      # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __await_delete_scene_for_protein_session(
          self, return_value: tuple[bool]
  ) -> None:
    """Finishes the delete scene for protein session process.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      if tmp_result[0] is False:
        self._interface_manager.stop_wait_cursor()
        self._interface_manager.refresh_main_view()
        self._interface_manager.status_bar_manager.show_error_message(
          "Deleting the selected scene failed!"
        )
        return

      if self._interface_manager.current_tab_index == 1:
        self._save_protein_pymol_session()
        self._interface_manager.remove_scene_from_proteins_model(
          self._interface_manager.get_current_protein_tree_index()
        )
        self._view.ui.lbl_pymol_protein_scene.setText(
          "PyMOL Scene: No Scene Selected"
        )
        self.__await_delete_current_scene(('generic_id', [(True, (0, True))]))
      elif self._interface_manager.current_tab_index == 2:
        # The database thread cannot be used here because the session gets loaded again
        # before the new data is in the db
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.save_protein_pair_pymol_session_to_database,
                args=(
                  self._interface_manager,
                  0,
                ),
              ),
              an_await_function=self.__await_delete_current_scene,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        self._interface_manager.status_bar_manager.show_temporary_message(
          "Deleting selected scene ...",
          False,
        )
        self._interface_manager.block_gui()
        self._interface_manager.remove_scene_from_protein_pairs_model(
          self._interface_manager.get_current_protein_pair_tree_index(),
        )
        self._view.ui.lbl_pymol_protein_pair_scene.setText(
          "PyMOL Scene: No Scene Selected"
        )
      else:
        logger.warning(
          "The current tab index is not for the proteins nor for the protein pairs tab?!"
        )
        self._interface_manager.stop_wait_cursor()
        self._interface_manager.refresh_main_view()
        return
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  def __await_delete_current_scene(self, return_value: tuple) -> None:
    """Finishes the delete current scene process.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      _, exit_flag = tmp_result
      if exit_flag:
        self._interface_manager.status_bar_manager.show_temporary_message(
          "Deleted the scene successfully."
        )
      else:
        self._interface_manager.status_bar_manager.show_error_message(
          "Deleting the scene failed!"
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  # </editor-fold>

  # </editor-fold>

  # <editor-fold desc="Protein Pairs tab methods">
  def __slot_expand_protein_pair(self) -> None:
    """Expands a protein pair branch."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "A protein pair of the tree view was expanded.",
      )
      tmp_type = (
        self._interface_manager.get_current_protein_pair_tree_index().data(
          enums.ModelEnum.TYPE_ROLE
        )
      )
      if tmp_type == "protein_pair":
        # protein pair
        self._view.ui.protein_pairs_tree_view.setExpanded(
          self._interface_manager.get_current_protein_pair_tree_index(),
          True,
        )
        # scenes
        self._view.ui.protein_pairs_tree_view.setExpanded(
          self._interface_manager.get_current_protein_pair_tree_index().child(
            0, 0
          ),
          True,
        )
        # protein 1
        self._view.ui.protein_pairs_tree_view.setExpanded(
          self._interface_manager.get_current_protein_pair_tree_index().child(
            1, 0
          ),
          True,
        )
        # chains of protein 1
        self._view.ui.protein_pairs_tree_view.expandRecursively(
          self._interface_manager.get_current_protein_pair_tree_index().child(
            1, 0
          ),
          1,
        )
        # protein 2
        self._view.ui.protein_pairs_tree_view.setExpanded(
          self._interface_manager.get_current_protein_pair_tree_index().child(
            2, 0
          ),
          True,
        )
        # chains of protein 2
        self._view.ui.protein_pairs_tree_view.expandRecursively(
          self._interface_manager.get_current_protein_pair_tree_index().child(
            2, 0
          ),
          1,
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def __slot_collapse_protein_pair(self) -> None:
    """Collapses a protein pair branch."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "A protein pair of the tree view was collapsed.",
      )
      tmp_type = (
        self._interface_manager.get_current_protein_pair_tree_index().data(
          enums.ModelEnum.TYPE_ROLE
        )
      )
      if tmp_type == "protein_pair":
        # protein pair
        self._view.ui.protein_pairs_tree_view.collapse(
          self._interface_manager.get_current_protein_pair_tree_index(),
        )
        # scenes
        self._view.ui.protein_pairs_tree_view.collapse(
          self._interface_manager.get_current_protein_pair_tree_index().child(
            0, 0
          ),
        )
        # protein 1
        self._view.ui.protein_pairs_tree_view.collapse(
          self._interface_manager.get_current_protein_pair_tree_index().child(
            1, 0
          ),
        )
        # chains of protein 1
        self._view.ui.protein_pairs_tree_view.collapse(
          self._interface_manager.get_current_protein_pair_tree_index().child(
            1, 0
          ),
        )
        # protein 2
        self._view.ui.protein_pairs_tree_view.collapse(
          self._interface_manager.get_current_protein_pair_tree_index().child(
            2, 0
          ),
        )
        # chains of protein 2
        self._view.ui.protein_pairs_tree_view.collapse(
          self._interface_manager.get_current_protein_pair_tree_index().child(
            2, 0
          ),
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  def __slot_expand_all_protein_pairs(self) -> None:
    """Expands all protein pair tree branches."""
    self._view.ui.protein_pairs_tree_view.expandAll()

  def __slot_collapse_all_protein_pairs(self) -> None:
    """Collapses all protein pair tree branches."""
    self._view.ui.protein_pairs_tree_view.collapseAll()

  def open_context_menu_for_protein_pairs(self, position) -> None:
    """Opens the context menu for the protein pairs tab."""
    try:
      tmp_protein_pair = (
        self._interface_manager.get_current_active_protein_pair_object()
      )
    except ValueError:
      tmp_is_protein_pair_in_current_session_flag = False
    else:
      tmp_is_protein_pair_in_current_session_flag = self._interface_manager.pymol_session_manager.is_the_current_protein_pair_in_session(
        tmp_protein_pair.name
      )
    tmp_is_protein_pair_expanded_flag: bool = False
    try:
      if (
              self._interface_manager.get_current_protein_pair_tree_index().data(
                enums.ModelEnum.TYPE_ROLE
              )
              == "protein_pair"
      ):
        if self._view.ui.protein_pairs_tree_view.isExpanded(
                self._interface_manager.get_current_protein_pair_tree_index()
        ):
          tmp_is_protein_pair_expanded_flag: bool = True
    except Exception as e:
      logger.error(e)
    else:
      tmp_context_menu = self._protein_pair_tree_context_menu.get_context_menu(
        self._view.ui.protein_pairs_tree_view.selectedIndexes(),
        tmp_is_protein_pair_in_current_session_flag,
        tmp_is_protein_pair_expanded_flag,
      )
      tmp_context_menu.exec_(
        self._view.ui.protein_pairs_tree_view.viewport().mapToGlobal(position)
      )
      self.__slot_get_information_about_selected_object_in_protein_pair_branch()  # fixme: This should be done in a better way than this!

  def _get_protein_name_of_a_protein_from_a_protein_pair(
          self,
          a_protein: "protein.Protein",
          a_protein_pair: "protein_pair.ProteinPair",
  ) -> str:
    """Helper function to get the correct protein name even if the protein pair consists of two identical protein names.

    Args:
        a_protein (protein.Protein): The protein to get the name from the protein pair.
        a_protein_pair (protein_pair.ProteinPair): The protein pair that contains the protein name.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_protein is None:
      logger.error("a_protein is None.")
      raise exception.IllegalArgumentError("a_protein is None.")
    if a_protein_pair is None:
      logger.error("a_protein_pair is None.")
      raise exception.IllegalArgumentError("a_protein_pair is None.")

    # </editor-fold>

    tmp_result = (
      self._interface_manager.pymol_session_manager.user_pymol_connector.get_all_object_names()
    )
    tmp_protein_name = a_protein.get_molecule_object()
    if tmp_result["success"]:
      tmp_sub_string_prot_1 = tmp_result["data"][0][
                              len(tmp_result["data"][0]) - 2: len(tmp_result["data"][0])
                              ]
      tmp_sub_string_prot_2 = tmp_result["data"][1][
                              len(tmp_result["data"][1]) - 2: len(tmp_result["data"][1])
                              ]
      if tmp_sub_string_prot_1 == "_1" and tmp_sub_string_prot_2 == "_2":
        if a_protein == a_protein_pair.protein_1:
          tmp_protein_name = f"{tmp_protein_name}_1"
        elif a_protein == a_protein_pair.protein_2:
          tmp_protein_name = f"{tmp_protein_name}_2"
        else:
          raise ValueError("Protein object not found in protein pair!")
    else:
      raise RuntimeError("PyMOL command failed!")
    return tmp_protein_name

  # <editor-fold desc="PyMOL session">
  def __slot_open_protein_pair_pymol_session(self) -> None:
    """Starts an async method that opens a protein pair pymol session."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Open protein pair pymol session' button on the 'Protein Pairs Tab' was clicked.",
      )
      self._view.tg_protein_pair_white_bg.toggle_button.setCheckState(False)
      tmp_protein_pair: "protein_pair.ProteinPair" = (
        self._interface_manager.get_current_active_protein_pair_object()
      )
      tmp_flag = False
      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=pymol_session_async.load_protein_pair_pymol_session,
              args=(
                tmp_protein_pair,
                self._interface_manager.pymol_session_manager,
                tmp_flag,
              ),
            ),
            an_await_function=self.__await_open_protein_pair_pymol_session,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
      # self._active_task = tasks.LegacyTask(
      #     target=pymol_session_async.load_protein_pair_pymol_session,
      #     args=(
      #         tmp_protein_pair,
      #         self._interface_manager.pymol_session_manager,
      #         tmp_flag,
      #     ),
      #     post_func=self.__await_open_protein_pair_pymol_session,
      # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui(with_wait_cursor=True)
      self._interface_manager.status_bar_manager.show_temporary_message(
        f"Loading PyMOL session of {tmp_protein_pair.name} ...",
        False,
      )

  def __await_open_protein_pair_pymol_session(
          self, return_value: tuple
  ) -> None:
    """Finishes the open protein pair pymol session process.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      _, exit_boolean = tmp_result
      self._view.ui.action_protein_regions.setEnabled(False)
      if exit_boolean:
        self._view.ui.action_protein_regions.setEnabled(True)
        self._view.ui.btn_create_protein_pair_scene.setEnabled(True)
        self._view.ui.btn_update_protein_pair_scene.setEnabled(True)
        self._view.ui.lbl_session_name.setText(
          f"Session Name: {self._interface_manager.pymol_session_manager.session_name}"
        )
        # self._interface_manager.pymol_session_manager.current_scene_name = self._view.ui.protein_pairs_tree_view.currentIndex().child(0, 0).child(1, 0).data(Qt.DisplayRole)
        self._interface_manager.pymol_session_manager.current_scene_name = (
          "base"
        )
        self._interface_manager.pymol_session_manager.load_current_scene()
        ui_util.set_pymol_scene_name_into_label(
          self._interface_manager.pymol_session_manager.current_scene_name,
          self._view.ui.lbl_pymol_protein_pair_scene,
        )
        logger.info("Successfully opened protein pair session.")
        self._interface_manager.status_bar_manager.show_temporary_message(
          "Loading the PyMOL session was successful."
        )
        self._view.ui.lbl_info_3.setText("Please select a chain.")
      else:
        logger.error(
          "The protein name could not be found in the object list in PyMOL!"
        )
        self._view.ui.btn_create_protein_pair_scene.setEnabled(False)
        self._view.ui.btn_update_protein_pair_scene.setEnabled(False)
        self._interface_manager.status_bar_manager.show_error_message(
          "Loading the PyMOL session failed!"
        )
        self._view.ui.lbl_info_3.setText(
          "Please load the PyMOL session of the selected protein."
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  def _save_protein_pair_pymol_session(self) -> None:
    """Saves the session as base64 string and updates the database."""
    try:
      tmp_protein_pair = (
        self._interface_manager.get_current_active_protein_pair_object()
      )
      tmp_protein_pair.pymol_session = (
        self._interface_manager.pymol_session_manager.save_current_session_as_base64()
      )
      tmp_database_operation = database_operation.DatabaseOperation(
        enums.SQLQueryType.UPDATE_PYMOL_SESSION_PROTEIN_PAIR,
        (0, tmp_protein_pair.get_id(), tmp_protein_pair),
      )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._database_thread.put_database_operation_into_queue(
        tmp_database_operation
      )

  # </editor-fold>

  def __slot_get_information_about_selected_object_in_protein_pair_branch(
          self,
  ) -> None:
    """Modifies GUI that all available options are seen for the pymol scene configuration panel."""
    try:
      tmp_type = (
        self._interface_manager.get_current_protein_pair_tree_index_type()
      )

      if tmp_type == "protein_pair":
        logger.log(
          log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
          f"The protein pair object '{self._view.ui.protein_pairs_tree_view.currentIndex().data(Qt.DisplayRole)}' on the 'Protein Pairs Tab' was clicked.",
        )
      elif tmp_type == "protein":
        logger.log(
          log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
          f"The protein object '{self._view.ui.protein_pairs_tree_view.currentIndex().data(Qt.DisplayRole)}' on the 'Protein Pairs Tab' was clicked.",
        )
      elif tmp_type == "scene":
        self._setup_protein_pair_scene()
      elif tmp_type == "chain":
        self._setup_protein_pair_pymol_scene_config()
      elif tmp_type == "header":
        logger.log(
          log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
          f"The header '{self._view.ui.protein_pairs_tree_view.currentIndex().data(Qt.DisplayRole)}' on the 'Protein Pairs Tab' was clicked.",
        )
      else:
        logger.warning("Unknown object type occurred in Protein Pairs tab.")
        return

      self._interface_manager.manage_ui_of_protein_pairs_tab(
        tmp_type,
        self._interface_manager.pymol_session_manager,
      )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  # <editor-fold desc="Setup protein pair scene of protein pairs tree">
  def _setup_protein_pair_scene(self) -> None:
    """Loads the default pymol scene."""
    try:
      tmp_scene_name = (
        self._view.ui.protein_pairs_tree_view.currentIndex().data(
          Qt.DisplayRole
        )
      )
      tmp_protein_pair_name = (
        self._view.ui.protein_pairs_tree_view.currentIndex()
        .parent()
        .parent()
        .data(Qt.DisplayRole)
      )
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        f"The scene '{tmp_scene_name}' of the protein pair '{tmp_protein_pair_name}' on the 'Protein Pairs Tab' was clicked.",
      )
      if not self._interface_manager.pymol_session_manager.is_the_current_protein_pair_in_session(
              self._interface_manager.get_current_active_protein_pair_object().name
      ):
        # The selected protein pair is not loaded into the current session, therefore nothing to do.
        return
      tmp_scene_name = (
        self._interface_manager.get_current_active_scene_name_of_protein_pair()
      )
      self._interface_manager.pymol_session_manager.current_scene_name = (
        tmp_scene_name
      )
      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=pymol_session_async.load_scene,
              args=(
                self._interface_manager.pymol_session_manager, tmp_scene_name
              ),
            ),
            an_await_function=self.__await_load_scene_protein_pair,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
      # self._active_task = tasks.LegacyTask(
      #     target=pymol_session_async.load_scene,
      #     args=(self._interface_manager.pymol_session_manager, tmp_scene_name),
      #     post_func=self.__await_load_scene_protein_pair,
      # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.status_bar_manager.show_temporary_message(
        "Loading PyMOL scene ...", a_with_timeout_flag=False
      )
      self._interface_manager.block_gui()

  def __await_load_scene_protein_pair(self, return_value: tuple[bool]) -> None:
    """Finishes the load default pymol scene process.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      ui_util.set_pymol_scene_name_into_label(
        self._interface_manager.pymol_session_manager.current_scene_name,
        self._view.ui.lbl_pymol_protein_pair_scene,
      )
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      if tmp_result[0]:
        self._interface_manager.status_bar_manager.show_temporary_message(
          "Loading PyMOL scene was successful."
        )
      else:
        self._interface_manager.status_bar_manager.show_temporary_message(
          "Loading PyMOL scene failed! Please try again."
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  # </editor-fold>

  # <editor-fold desc="Setup PyMOL scene configuration">
  def _setup_protein_pair_pymol_scene_config(self) -> None:
    """Sets up the color grid and representation section of the pymol scene configuration panel."""
    try:
      tmp_chain_letter = (
        self._view.ui.protein_pairs_tree_view.currentIndex().data(
          Qt.DisplayRole
        )
      )
      tmp_protein = (
        self._view.ui.protein_pairs_tree_view.currentIndex()
        .parent()
        .parent()
        .data(
          enums.ModelEnum.OBJECT_ROLE,
        )
      )
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        f"The chain object '{tmp_chain_letter}' of the protein '{tmp_protein.get_molecule_object}' on the 'Protein Pairs Tab' was clicked.",
      )
      if (
              self._interface_manager.pymol_session_manager.current_scene_name == ""
              and not self._interface_manager.pymol_session_manager.is_the_current_protein_pair_in_session(
        self._interface_manager.get_current_active_protein_pair_object().name
      )
      ):
        # The selected protein pair is not loaded into the current session, therefore nothing to do.
        return
      # Set icon for color grid
      self._view.color_grid_protein_pairs.reset_icon_for_selected_color()
      tmp_protein = (
        self._interface_manager.get_current_active_protein_object_of_protein_pair()
      )
      tmp_chain = (
        self._interface_manager.get_current_active_chain_object_of_protein_pair()
      )
      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=pymol_session_async.get_residue_color_config_of_a_given_protein_chain,
              args=(
                tmp_protein.get_molecule_object(),
                tmp_chain.chain_letter,
                self._interface_manager.pymol_session_manager,
              ),
            ),
            an_await_function=self.__await_get_residue_color_config_of_a_given_protein_chain_of_a_protein_pair,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
      # self._active_task = tasks.LegacyTask(
      #     target=pymol_session_async.get_residue_color_config_of_a_given_protein_chain,
      #     args=(
      #         tmp_protein.get_molecule_object(),
      #         tmp_chain.chain_letter,
      #         self._interface_manager.pymol_session_manager,
      #     ),
      #     post_func=self.__await_get_residue_color_config_of_a_given_protein_chain_of_a_protein_pair,
      # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __await_get_residue_color_config_of_a_given_protein_chain_of_a_protein_pair(
          self,
          return_value: tuple[
            bool, Optional["residue_color_config.ResidueColorConfig"]
          ],
  ) -> None:
    """Finishes the setup of the color grid and representation section process.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      logger.debug(
        "Returned from method 'pymol_session_async.get_residue_color_config_of_a_given_protein_chain'."
      )
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      tmp_success_flag, tmp_residue_color_config = tmp_result
      logger.debug(f"The return_value is: {return_value}.")
      if not tmp_success_flag:
        logger.error("Retrieving color information failed!")
        self._interface_manager.status_bar_manager.show_error_message(
          "Retrieving color information failed!"
        )
        self._interface_manager.stop_wait_cursor()
        self._interface_manager.refresh_main_view()
        return

      logger.debug("Check if chain is colored by element.")
      if tmp_residue_color_config.atoms_are_colored_by_elements():
        #self._view.tg_protein_pair_color_atoms.toggle_button.setChecked(True)
        ui_util.set_checked_async(
          self._view.tg_protein_pair_color_atoms.toggle_button, True
        )
        self._view.ui.lbl_protein_pair_current_color.setText("By Element    ")
      else:
        # self._view.tg_protein_pair_color_atoms.toggle_button.setChecked(False)
        ui_util.set_checked_async(
          self._view.tg_protein_pair_color_atoms.toggle_button, False
        )
        self._interface_manager.set_current_active_chain_color_of_protein_pair(
          tmp_residue_color_config.carbon_color
        )
        self._view.color_grid_protein_pairs.set_icon_for_selected_color(
          tmp_residue_color_config.carbon_color
        )
        self._view.ui.lbl_protein_pair_current_color.setText(
          f"{tmp_residue_color_config.carbon_color}    "
        )
      # Set representation toggle states for selected chain
      tmp_protein = (
        self._interface_manager.get_current_active_protein_object_of_protein_pair()
      )
      tmp_chain = (
        self._interface_manager.get_current_active_chain_object_of_protein_pair()
      )
      tmp_protein.pymol_selection.selection_string = f"first chain {tmp_chain.chain_letter} and {tmp_protein.get_molecule_object()}"
      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=pymol_session_async.get_representation_config_of_a_given_protein_chain,
              args=(
                tmp_protein.pymol_selection.selection_string,
                tmp_chain.chain_letter,
                self._interface_manager.pymol_session_manager,
              ),
            ),
            an_await_function=self.__await_get_representation_config_of_a_given_protein_chain_of_a_protein_pair,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
      # self._active_task = tasks.LegacyTask(
      #     target=pymol_session_async.get_representation_config_of_a_given_protein_chain,
      #     args=(
      #         tmp_protein.pymol_selection.selection_string,
      #         tmp_chain.chain_letter,
      #         self._interface_manager.pymol_session_manager,
      #     ),
      #     post_func=self.__await_get_representation_config_of_a_given_protein_chain_of_a_protein_pair,
      # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  def __await_get_representation_config_of_a_given_protein_chain_of_a_protein_pair(
          self,
          return_value: tuple[bool, Optional[dict]],
  ) -> None:
    """Fixme: This method should be obsolete with the new TEA lib is introduced!"""
    try:
      logger.debug(
        "Returned from method 'pymol_session_async.get_representation_config_of_a_given_protein_chain'."
      )
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      tmp_success_flag, tmp_representation_config = tmp_result
      logger.debug(f"The return_value is: {return_value}.")
      if tmp_success_flag:
        tmp_chain = (
          self._interface_manager.get_current_active_chain_object_of_protein_pair()
        )
        if tmp_chain.chain_type == "protein_chain":
          self._view.ui.frame_protein_pair_repr.setEnabled(True)
          self._interface_manager.manage_toggle_state_of_protein_pair_repr(
            tmp_representation_config
          )
        else:
          self._view.ui.frame_protein_pair_repr.setEnabled(False)
          ui_util.set_checked_async(
            self._view.tg_protein_pair_cartoon.toggle_button, False
          )
          ui_util.set_checked_async(
            self._view.tg_protein_pair_ribbon.toggle_button, False
          )
          ui_util.set_checked_async(
            self._view.tg_protein_pair_sticks.toggle_button, False
          )
          ui_util.set_checked_async(
            self._view.tg_protein_pair_lines.toggle_button, False
          )
          ui_util.set_checked_async(
            self._view.tg_protein_pair_spheres.toggle_button, False
          )
          ui_util.set_checked_async(
            self._view.tg_protein_pair_dots.toggle_button, False
          )
          ui_util.set_checked_async(
            self._view.tg_protein_pair_mesh.toggle_button, False
          )
          ui_util.set_checked_async(
            self._view.tg_protein_pair_surface.toggle_button, False
          )
          # self._view.tg_protein_pair_cartoon.toggle_button.setChecked(False)
          # self._view.tg_protein_pair_ribbon.toggle_button.setChecked(False)
          # self._view.tg_protein_pair_sticks.toggle_button.setChecked(False)
          # self._view.tg_protein_pair_lines.toggle_button.setChecked(False)
          # self._view.tg_protein_pair_spheres.toggle_button.setChecked(False)
          # self._view.tg_protein_pair_dots.toggle_button.setChecked(False)
          # self._view.tg_protein_pair_mesh.toggle_button.setChecked(False)
          # self._view.tg_protein_pair_surface.toggle_button.setChecked(False)
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      QtWidgets.QApplication.restoreOverrideCursor()  # Due to unknown reasons this line is needed to reset the cursor

  # </editor-fold>

  # <editor-fold desc="Change chain color of protein">
  def __slot_change_chain_color_protein_pairs(self, a_color: str) -> None:
    """Changes the 'Color' attribute of a protein chain on the 'Protein Pairs Tab' and colors the protein in User PyMOL.

    Args:
        a_color (str): The new color for the protein chain.

    Notes:
        This slot method is used if a button of the color grid is clicked.
    """
    # <editor-fold desc="Checks">
    if a_color is None or a_color == "":
      logger.error("a_color is either None or an empty string!")
      self._interface_manager.status_bar_manager.show_error_message(
        "a_color is either None or an empty string!"
      )
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      return
    if a_color not in constants.PYMOL_COLORS_WITH_INDICES.values():
      logger.error("a_color is not part of the PYMOL_COLORS_WITH_INDICES dict!")
      self._interface_manager.status_bar_manager.show_error_message(
        "a_color is not part of the PYMOL_COLORS_WITH_INDICES dict!"
      )
      return

    # </editor-fold>

    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "The 'Color' attribute of a protein chain on the 'Protein Pairs Tab' changed.",
      )
      tmp_protein_pair = (
        self._interface_manager.get_current_active_protein_pair_object()
      )
      if not self._interface_manager.pymol_session_manager.is_the_current_protein_pair_in_session(
              tmp_protein_pair.name
      ):
        logger.warning(
          f"Protein pair {tmp_protein_pair.name} is not in the current session."
        )
        self._interface_manager.status_bar_manager.show_error_message(
          f"Protein pair {tmp_protein_pair.name} is not in the current session.",
        )
        return

      tmp_protein = (
        self._interface_manager.get_current_active_protein_object_of_protein_pair()
      )
      tmp_chain = (
        self._interface_manager.get_current_active_chain_object_of_protein_pair()
      )
      tmp_protein_name = (
        self._get_protein_name_of_a_protein_from_a_protein_pair(
          tmp_protein, tmp_protein_pair
        )
      )
      tmp_protein.pymol_selection.set_custom_selection(
        f"/{tmp_protein_name}//{tmp_chain.chain_letter}"
      )
      # Color protein in User PyMOL
      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=pymol_session_async.color_pymol_selection,
              args=(
                a_color,
                tmp_protein.pymol_selection.selection_string,
                self._interface_manager.pymol_session_manager,
              ),
            ),
            an_await_function=self.__await_color_pymol_selection_for_protein_pair,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
      # self._active_task = tasks.LegacyTask(
      #     target=pymol_session_async.color_pymol_selection,
      #     args=(
      #         a_color,
      #         tmp_protein.pymol_selection.selection_string,
      #         self._interface_manager.pymol_session_manager,
      #     ),
      #     post_func=self.__await_color_pymol_selection_for_protein_pair,
      # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()  # fixme: This blocks the entire UI which leads to quick flashes. This might be a problem.

    # try:
    #     logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
    #                "The 'Color' attribute of a protein chain on the 'Protein Pairs Tab' changed.")
    #     tmp_protein_pair = self._interface_manager.get_current_active_protein_pair_object()
    #     tmp_protein = self._interface_manager.get_current_active_protein_object_of_protein_pair()
    #     tmp_chain = self._interface_manager.get_current_active_chain_object_of_protein_pair()
    #     # tmp_color = self._view.ui.lbl_protein_pair_current_color.text().strip()
    #
    #     if self._interface_manager.pymol_session_manager.session_object_type == "protein_pair" and self._interface_manager.pymol_session_manager.session_name == tmp_protein_pair.name:
    #         # Update pymol parameter in PyMOL
    #         tmp_protein_name = self._get_protein_name_of_a_protein_from_a_protein_pair(tmp_protein,
    #                                                                                    tmp_protein_pair)
    #         tmp_protein.pymol_selection.set_custom_selection(f"/{tmp_protein_name}//{tmp_chain.chain_letter}")
    #         self._interface_manager.pymol_session_manager.color_protein(
    #             a_color, tmp_protein.pymol_selection.selection_string
    #         )
    #         # Update pymol parameter in memory
    #         self._view.color_grid_protein_pairs.reset_icon_for_selected_color()
    #         if a_color != "By Element":
    #             tmp_chain.pymol_parameters["chain_color"] = a_color
    #             self._interface_manager.set_current_active_chain_color_of_protein_pair(a_color)
    #             # Update pymol parameter in database
    #             with database_manager.DatabaseManager(str(self._interface_manager.get_current_project().get_database_filepath())) as db_manager:
    #                 db_manager.update_protein_chain_color_of_protein_pair(
    #                     tmp_protein.get_id(), tmp_chain.chain_letter, tmp_protein_pair.get_id(), a_color
    #                 )
    #         self._view.color_grid_protein_pairs.set_icon_for_selected_color(a_color)
    #         self._view.ui.lbl_protein_pair_current_color.setText(f"{a_color}    ")
    #         self._update_protein_pair_scene()
    #         self._save_protein_pair_pymol_session()
    #         self._view.tg_protein_pair_color_atoms.toggle_button.setChecked(False)
    #     else:
    #         logger.warning("The color of a protein chain could not be changed. This can be due to UI setup reasons.")
    # except Exception as e:
    #     logger.error(f"An error occurred: {e}")
    #     self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

  def __await_color_pymol_selection_for_protein_pair(
          self, return_value: tuple[bool, str]
  ) -> None:
    """Updates the color of the protein chain based on the return value provided in the object and data model.

    Args:
        return_value (tuple): A tuple containing two values, a boolean flag indicating the success of the operation and a string representing the color.
    """
    # <editor-fold desc="Checks">
    if return_value is None or len(return_value) == 0:
      logger.error("return_value is either None or has a length of 0.")
      self._interface_manager.status_bar_manager.show_error_message(
        "return_value is either None or has a length of 0."
      )
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      tmp_success_flag, tmp_chain_color = tmp_result
      if tmp_success_flag:
        tmp_chain = (
          self._interface_manager.get_current_active_chain_object_of_protein_pair()
        )

        # <editor-fold desc="Setup color grid icon">
        self._view.color_grid_protein_pairs.reset_icon_for_selected_color()
        if tmp_chain_color != "By Element":
          tmp_chain.pymol_parameters["chain_color"] = tmp_chain_color
          self._interface_manager.set_current_active_chain_color_of_protein_pair(
            tmp_chain_color
          )
          self._view.color_grid_protein_pairs.set_icon_for_selected_color(
            tmp_chain_color
          )
          # self._view.tg_protein_pair_color_atoms.toggle_button.setChecked(False)
          ui_util.set_checked_async(
            self._view.tg_protein_pair_color_atoms.toggle_button, False
          )

        # </editor-fold>

        self._view.ui.lbl_protein_pair_current_color.setText(
          f"{tmp_chain_color}    "
        )
        self._update_protein_pair_scene_legacy()
        self._save_protein_pair_pymol_session()
      else:
        logger.error("The operation failed!")
        self._interface_manager.status_bar_manager.show_error_message(
          "The operation failed!"
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()

  # </editor-fold>

  # <editor-fold desc="Color protein chain atoms by element">
  def __slot_color_protein_pair_atoms_by_element(self) -> None:
    """Color protein atoms by element.

    Notes:
        This method is used to color the protein atoms based on their element.
        If the toggle button is checked, the method will color the atoms based on their element:
        grey70 for carbon atoms and atomic color for non-carbon atoms.
        It will also update the protein chain color in the database.
        If the toggle button is not checked, the method will reset the color to the chain color
        or green if the chain color is set to "By Element".

    """
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
      "Toggle 'By Elements' button on the 'Protein Pairs Tab' was clicked.",
    )
    try:
      # Create selection for atom coloring
      tmp_protein_pair = (
        self._interface_manager.get_current_active_protein_pair_object()
      )
      tmp_protein = (
        self._interface_manager.get_current_active_protein_object_of_protein_pair()
      )
      tmp_chain = (
        self._interface_manager.get_current_active_chain_object_of_protein_pair()
      )
      tmp_protein_name = (
        self._get_protein_name_of_a_protein_from_a_protein_pair(
          tmp_protein, tmp_protein_pair
        )
      )
      tmp_protein.pymol_selection.set_custom_selection(
        f"/{tmp_protein_name}//{tmp_chain.chain_letter}"
      )

      if self._view.tg_protein_pair_color_atoms.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.color_pymol_selection_atoms_by_element,
                args=(
                  tmp_protein.pymol_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_color_pymol_selection_atoms_by_element_for_protein_pair,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.color_pymol_selection_atoms_by_element,
        #     args=(
        #         tmp_protein.pymol_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_color_pymol_selection_atoms_by_element_for_protein_pair,
        # )
      else:
        if (
                self._interface_manager.get_current_active_chain_color_of_protein_pair()
                is None
        ):
          tmp_current_active_chain_color = "By Element"
        else:
          tmp_current_active_chain_color = (
            self._interface_manager.get_current_active_chain_color_of_protein_pair()
          )

        tmp_protein = (
          self._interface_manager.get_current_active_protein_object_of_protein_pair()
        )
        tmp_chain = (
          self._interface_manager.get_current_active_chain_object_of_protein_pair()
        )
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.reset_color_pymol_selection_atoms_by_element,
                args=(
                  tmp_protein.get_molecule_object(),
                  tmp_chain.chain_letter,
                  tmp_current_active_chain_color,
                  tmp_protein.pymol_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_reset_color_pymol_selection_atoms_by_element_for_protein_pair,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.reset_color_pymol_selection_atoms_by_element,
        #     args=(
        #         tmp_protein.get_molecule_object(),
        #         tmp_chain.chain_letter,
        #         tmp_current_active_chain_color,
        #         tmp_protein.pymol_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_reset_color_pymol_selection_atoms_by_element_for_protein_pair,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

    # try:
    #     tmp_protein = self._interface_manager.get_current_active_protein_object_of_protein_pair()
    #     tmp_chain = self._interface_manager.get_current_active_chain_object_of_protein_pair()
    #     tmp_protein_pair = self._interface_manager.get_current_active_protein_pair_object()
    #     if self._view.tg_protein_pair_color_atoms.toggle_button.isChecked():
    #         tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
    #         tmp_protein_name = self._get_protein_name_of_a_protein_from_a_protein_pair(tmp_protein,
    #                                                                                    tmp_protein_pair)
    #         tmp_protein.pymol_selection.set_custom_selection(f"/{tmp_protein_name}//{tmp_chain.chain_letter}")
    #         self._interface_manager.pymol_session_manager.color_protein(
    #             "atomic", f"{tmp_selection.selection_string} and not elem C"
    #         )
    #         self._interface_manager.pymol_session_manager.color_protein(
    #             "grey70", f"{tmp_selection.selection_string} and elem C"
    #         )
    #         with database_manager.DatabaseManager(
    #                 str(self._interface_manager.get_current_project().get_database_filepath())) as db_manager:
    #             db_manager.update_protein_chain_color_of_protein_pair(
    #                 tmp_protein.get_id(), tmp_chain.chain_letter, tmp_protein_pair.get_id(), "By Element"
    #             )
    #         self._view.color_grid_protein_pairs.reset_icon_for_selected_color()
    #         self._view.ui.lbl_protein_pair_current_color.setText("By Element    ")
    #     else:
    #         logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
    #                    "'Reset color atoms by element' button on the 'Protein Pairs Tab' was clicked.")
    #         tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
    #         tmp_protein_name = self._get_protein_name_of_a_protein_from_a_protein_pair(tmp_protein,
    #                                                                                    tmp_protein_pair)
    #         tmp_protein.pymol_selection.set_custom_selection(f"/{tmp_protein_name}//{tmp_chain.chain_letter}")
    #         tmp_chain_color = self._interface_manager.get_current_active_chain_color_of_protein_pair()
    #         if tmp_chain_color == "By Element":
    #             tmp_chain_color = "green"
    #         self._interface_manager.pymol_session_manager.color_protein(
    #             tmp_chain_color, f"{tmp_selection.selection_string}"
    #         )
    #         self._view.color_grid_protein_pairs.reset_icon_for_selected_color()
    #         self._view.color_grid_protein_pairs.set_icon_for_selected_color(tmp_chain_color)
    #         self._view.ui.lbl_protein_pair_current_color.setText(f"{tmp_chain_color}    ")
    #     # tmp_protein = self._interface_manager.get_current_active_protein_object_of_protein_pair()
    #     # tmp_chain = self._interface_manager.get_current_active_chain_object_of_protein_pair()
    #     # tmp_protein_pair = self._interface_manager.get_current_active_protein_pair_object()
    #     # if self._view.tg_protein_pair_color_atoms.toggle_button.isChecked():
    #     #     tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
    #     #     tmp_protein_name = self._get_protein_name_of_a_protein_from_a_protein_pair(tmp_protein,
    #     #                                                                                tmp_protein_pair)
    #     #     tmp_protein.pymol_selection.set_custom_selection(f"/{tmp_protein_name}//{tmp_chain.chain_letter}")
    #     #     self._interface_manager.pymol_session_manager.color_protein(
    #     #         "atomic", f"{tmp_selection.selection_string} and not elem C"
    #     #     )
    #     #     self._interface_manager.pymol_session_manager.color_protein(
    #     #         "grey70", f"{tmp_selection.selection_string} and elem C"
    #     #     )
    #     #     self.reset_icon_for_last_color_in_protein_pairs_tab()
    #     #     self._view.ui.lbl_protein_pair_current_color.setText("By Element    ")
    #     # else:
    #     #     logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
    #     #                "'Reset color atoms by element' button on the 'Proteins Tab' was clicked.")
    #     #     tmp_protein_name = self._get_protein_name_of_a_protein_from_a_protein_pair(tmp_protein,
    #     #                                                                                tmp_protein_pair)
    #     #     tmp_protein.pymol_selection.set_custom_selection(f"/{tmp_protein_name}//{tmp_chain.chain_letter}")
    #     #     self._interface_manager.pymol_session_manager.color_protein(
    #     #         self._view.color_grid_protein_pairs.last_clicked_color, tmp_protein.pymol_selection.selection_string
    #     #     )
    #     #     self.set_icon_for_current_color_in_protein_pairs_tab()
    # except Exception as e:
    #     logger.error(f"An error occurred: {e}")
    #     self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

  def __await_color_pymol_selection_atoms_by_element_for_protein_pair(
          self, return_value: tuple
  ) -> None:
    """Await method for the coloring of atom by their element.

    Args:
        return_value: A tuple containing the result of the operation. The first element of the tuple determines whether the operation was successful or not.

    Returns:
        None
    """
    # <editor-fold desc="Checks">
    if return_value is None or len(return_value) == 0:
      logger.error("return_value is either None or has a length of 0.")
      self._interface_manager.status_bar_manager.show_error_message(
        "return_value is either None or has a length of 0."
      )
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      if tmp_result[0]:
        self._view.color_grid_protein_pairs.reset_icon_for_selected_color()
        self._view.ui.lbl_protein_pair_current_color.setText("By Element    ")
        self._update_protein_pair_scene_legacy()
        self._save_protein_pair_pymol_session()
      else:
        logger.error("The operation failed!")
        self._interface_manager.status_bar_manager.show_error_message(
          "The operation failed!"
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()

  def __await_reset_color_pymol_selection_atoms_by_element_for_protein_pair(
          self, return_value: tuple
  ) -> None:
    """Await method for the reset coloring of atom by their element.

    Args:
        return_value: A tuple containing the result of the operation and the chain color.

    Returns:
        None
    """
    # <editor-fold desc="Checks">
    if return_value is None or len(return_value) == 0:
      logger.error("return_value is either None or has a length of 0.")
      self._interface_manager.status_bar_manager.show_error_message(
        "return_value is either None or has a length of 0."
      )
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      tmp_success_flag, tmp_chain_color = tmp_result
      if tmp_success_flag:
        self._view.color_grid_protein_pairs.reset_icon_for_selected_color()

        # <editor-fold desc="Update chain color">
        # Updates chain color in chain object and data model
        tmp_chain = (
          self._interface_manager.get_current_active_chain_object_of_protein_pair()
        )
        tmp_chain.pymol_parameters["chain_color"] = tmp_chain_color
        self._interface_manager.set_current_active_chain_color_of_protein_pair(
          tmp_chain_color
        )

        # </editor-fold>

        self._view.color_grid_protein_pairs.set_icon_for_selected_color(
          tmp_chain_color
        )
        self._view.ui.lbl_protein_pair_current_color.setText(
          f"{tmp_chain_color}    "
        )

        self._update_protein_pair_scene_legacy()
        self._save_protein_pair_pymol_session()
      else:
        logger.error("The operation failed!")
        self._interface_manager.status_bar_manager.show_error_message(
          "The operation failed!"
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.refresh_main_view()
      self._interface_manager.stop_wait_cursor()

  # </editor-fold>

  # <editor-fold desc="Color protein pair by RMSD value">
  def __slot_color_protein_pair_by_rmsd(self) -> None:
    """Colors the residues in 5 colors depending on their distance to the reference."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Color protein pair by rmsd' context menu action was clicked.",
      )
      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=protein_pair_async.color_protein_pair_by_rmsd_value,
              args=(
                self._interface_manager.get_current_active_protein_pair_object(),
                self._interface_manager.pymol_session_manager,
              ),
            ),
            an_await_function=self.__await_color_protein_pair_by_rmsd,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
      # self._active_task = tasks.LegacyTask(
      #     target=protein_pair_async.color_protein_pair_by_rmsd_value,
      #     args=(
      #         self._interface_manager.get_current_active_protein_pair_object(),
      #         self._interface_manager.pymol_session_manager,
      #     ),
      #     post_func=self.__await_color_protein_pair_by_rmsd,
      # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __await_color_protein_pair_by_rmsd(self, result: tuple) -> None:
    """Finishes the color protein pair by rmsd process.

    Args:
        result (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if result is None:
      logger.error("result is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return
    tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(result)
    if tmp_result[0] == "":
      self._interface_manager.status_bar_manager.show_error_message(
        "Coloring the protein pair failed!"
      )

    # </editor-fold>

    self._interface_manager.stop_wait_cursor()
    self._interface_manager.refresh_main_view()

  # </editor-fold>

  # <editor-fold desc="Set background color">
  def __slot_protein_pair_change_background_color(self) -> None:
    """Sets the background color for the protein pair pymol session."""
    try:
      if self._view.tg_protein_pair_white_bg.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.set_background_color,
                args=(
                  "white", self._interface_manager.pymol_session_manager
                ),
              ),
              an_await_function=self.__await_set_background_color_for_protein_pair_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.set_background_color,
        #     args=("white", self._interface_manager.pymol_session_manager),
        #     post_func=self.__await_set_background_color_for_protein_pair_session,
        # )
      else:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.set_background_color,
                args=(
                  "black", self._interface_manager.pymol_session_manager
                ),
              ),
              an_await_function=self.__await_set_background_color_for_protein_pair_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.set_background_color,
        #     args=("black", self._interface_manager.pymol_session_manager),
        #     post_func=self.__await_set_background_color_for_protein_pair_session,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __await_set_background_color_for_protein_pair_session(
          self, return_value: tuple[bool]
  ) -> None:
    """Await method for setting the background color.

    Args:
        return_value (tuple): A tuple containing the result of the operation.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      if tmp_result[0]:
        self._interface_manager.status_bar_manager.show_temporary_message(
          "Background color updated."
        )
      else:
        self._interface_manager.status_bar_manager.show_error_message(
          "Updating background color failed!"
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  # </editor-fold>

  # <editor-fold desc="Representations">
  def __await_set_representation_for_protein_pair_session(
          self, return_value: tuple[bool]
  ) -> None:
    """Saves the pymol session after changing the representation.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      if tmp_result[0]:
        self._update_protein_pair_scene_legacy()
        self._save_protein_pair_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()
        self._interface_manager.manage_hydrogen_representation_for_protein_pair_chain()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  def _create_selection_string_for_representations(self) -> str:
    """Creates a specific selection string for the representations that is also compatible with same protein protein pairs."""
    try:
      tmp_protein = (
        self._interface_manager.get_current_active_protein_object_of_protein_pair()
      )
      tmp_chain = (
        self._interface_manager.get_current_active_chain_object_of_protein_pair()
      )
      tmp_protein_pair = (
        self._interface_manager.get_current_active_protein_pair_object()
      )
      tmp_protein_name = (
        self._get_protein_name_of_a_protein_from_a_protein_pair(
          tmp_protein, tmp_protein_pair
        )
      )
    except Exception as e:
      logger.error(e)
      return ""
    return f"/{tmp_protein_name}//{tmp_chain.chain_letter}"

  def __slot_protein_pair_chain_as_cartoon(self) -> None:
    """Changes the representation based on the toggle state."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Cartoon' toggle on the 'Protein Pairs Tab' was clicked.",
      )
      tmp_selection = (
        self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
      )
      tmp_selection.set_custom_selection(
        self._create_selection_string_for_representations()
      )
      if self._view.tg_protein_pair_cartoon.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.show_specific_representation,
                args=(
                  enums.PyMOLRepresentation.CARTOON,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_pair_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.show_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.CARTOON,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_pair_session,
        # )
      else:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.hide_specific_representation,
                args=(
                  enums.PyMOLRepresentation.CARTOON,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_pair_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.hide_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.CARTOON,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_pair_session,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

    # try:
    #     logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
    #                "'Cartoon' toggle on the 'Protein Pairs Tab' was clicked.")
    #     tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
    #     tmp_selection.set_custom_selection(self._create_selection_string_for_representations())
    #     if self._view.ui.cb_protein_pair_cartoon.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 0:
    #         #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.CARTOON.value)
    #         self._interface_manager.pymol_session_manager.show_specific_representation(
    #             enums.PyMOLRepresentation.CARTOON.value, tmp_selection.selection_string
    #         )
    #     elif self._view.tg_protein_pair_cartoon.toggle_button.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 1:
    #         #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.CARTOON.value)
    #         self._interface_manager.pymol_session_manager.show_specific_representation(
    #             enums.PyMOLRepresentation.CARTOON.value, tmp_selection.selection_string
    #         )
    #     else:
    #         #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.CARTOON.value)
    #         self._interface_manager.pymol_session_manager.hide_specific_representation(
    #             enums.PyMOLRepresentation.CARTOON.value, tmp_selection.selection_string
    #         )
    #     self._update_protein_pair_scene_legacy()
    #     self._save_protein_pair_pymol_session()
    #     self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()
    #     self._interface_manager.manage_hydrogen_representation_for_protein_pair_chain()
    # except Exception as e:
    #     logger.error(f"An error occurred: {e}")
    #     self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

  def __slot_protein_pair_chain_as_sticks(self) -> None:
    """Changes the representation based on the toggle state."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Sticks' toggle on the 'Protein Pairs Tab' was clicked.",
      )
      tmp_selection = (
        self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
      )
      tmp_selection.set_custom_selection(
        self._create_selection_string_for_representations()
      )
      if self._view.tg_protein_pair_sticks.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.show_specific_representation,
                args=(
                  enums.PyMOLRepresentation.STICKS,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_pair_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.show_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.STICKS,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_pair_session,
        # )
      else:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.hide_specific_representation,
                args=(
                  enums.PyMOLRepresentation.STICKS,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_pair_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.hide_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.STICKS,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_pair_session,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __slot_protein_pair_chain_as_ribbon(self) -> None:
    """Changes the representation based on the toggle state."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Ribbon' toggle on the 'Protein Pairs Tab' was clicked.",
      )
      tmp_selection = (
        self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
      )
      tmp_selection.set_custom_selection(
        self._create_selection_string_for_representations()
      )
      if self._view.tg_protein_pair_ribbon.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.show_specific_representation,
                args=(
                  enums.PyMOLRepresentation.RIBBON,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_pair_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.show_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.RIBBON,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_pair_session,
        # )
      else:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.hide_specific_representation,
                args=(
                  enums.PyMOLRepresentation.RIBBON,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_pair_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.hide_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.RIBBON,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_pair_session,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __slot_protein_pair_chain_as_lines(self) -> None:
    """Changes the representation based on the toggle state."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Lines' toggle on the 'Protein Pairs Tab' was clicked.",
      )
      tmp_selection = (
        self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
      )
      tmp_selection.set_custom_selection(
        self._create_selection_string_for_representations()
      )
      if self._view.tg_protein_pair_lines.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.show_specific_representation,
                args=(
                  enums.PyMOLRepresentation.LINES,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_pair_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.show_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.LINES,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_pair_session,
        # )
      else:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.hide_specific_representation,
                args=(
                  enums.PyMOLRepresentation.LINES,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_pair_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.hide_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.LINES,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_pair_session,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __slot_protein_pair_chain_as_spheres(self) -> None:
    """Changes the representation based on the toggle state."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Spheres' toggle on the 'Protein Pairs Tab' was clicked.",
      )
      tmp_selection = (
        self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
      )
      tmp_selection.set_custom_selection(
        self._create_selection_string_for_representations()
      )
      if self._view.tg_protein_pair_spheres.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.show_specific_representation,
                args=(
                  enums.PyMOLRepresentation.SPHERES,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_pair_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.show_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.SPHERES,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_pair_session,
        # )
      else:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.hide_specific_representation,
                args=(
                  enums.PyMOLRepresentation.SPHERES,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_pair_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.hide_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.SPHERES,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_pair_session,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __slot_protein_pair_chain_as_dots(self) -> None:
    """Changes the representation based on the toggle state."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Dots' toggle on the 'Protein Pairs Tab' was clicked.",
      )
      tmp_selection = (
        self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
      )
      tmp_selection.set_custom_selection(
        self._create_selection_string_for_representations()
      )
      if self._view.tg_protein_pair_dots.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.show_specific_representation,
                args=(
                  enums.PyMOLRepresentation.DOTS,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_pair_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.show_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.DOTS,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_pair_session,
        # )
      else:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.hide_specific_representation,
                args=(
                  enums.PyMOLRepresentation.DOTS,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_pair_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.hide_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.DOTS,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_pair_session,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __slot_protein_pair_chain_as_mesh(self) -> None:
    """Changes the representation based on the toggle state."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Mesh' toggle on the 'Protein Pairs Tab' was clicked.",
      )
      tmp_selection = (
        self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
      )
      tmp_selection.set_custom_selection(
        self._create_selection_string_for_representations()
      )
      if self._view.tg_protein_pair_mesh.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.show_specific_representation,
                args=(
                  enums.PyMOLRepresentation.MESH,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_pair_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.show_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.MESH,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_pair_session,
        # )
      else:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.hide_specific_representation,
                args=(
                  enums.PyMOLRepresentation.MESH,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_pair_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.hide_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.MESH,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_pair_session,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __slot_protein_pair_chain_as_surface(self) -> None:
    """Changes the representation based on the toggle state."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Surface' toggle on the 'Protein Pairs Tab' was clicked.",
      )
      tmp_selection = (
        self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
      )
      tmp_selection.set_custom_selection(
        self._create_selection_string_for_representations()
      )
      if self._view.tg_protein_pair_surface.toggle_button.isChecked():
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.show_specific_representation,
                args=(
                  enums.PyMOLRepresentation.SURFACE,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_pair_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.show_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.SURFACE,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_pair_session,
        # )
      else:
        self._interface_manager.get_task_manager().append_task_result(
          task_result_factory.TaskResultFactory.run_task_result(
            a_task_result=task_result.TaskResult.from_action(
              an_action=action.Action(
                a_target=pymol_session_async.hide_specific_representation,
                args=(
                  enums.PyMOLRepresentation.SURFACE,
                  tmp_selection.selection_string,
                  self._interface_manager.pymol_session_manager,
                ),
              ),
              an_await_function=self.__await_set_representation_for_protein_pair_session,
            ),
            a_task_scheduler=self._interface_manager.get_task_scheduler(),
          )
        )
        # self._active_task = tasks.LegacyTask(
        #     target=pymol_session_async.hide_specific_representation,
        #     args=(
        #         enums.PyMOLRepresentation.SURFACE,
        #         tmp_selection.selection_string,
        #         self._interface_manager.pymol_session_manager,
        #     ),
        #     post_func=self.__await_set_representation_for_protein_pair_session,
        # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  # def __slot_protein_pair_chain_as_sticks(self) -> None:
  #     try:
  #         logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
  #                    "'Sticks' toggle on the 'Protein Pairs Tab' was clicked.")
  #         tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
  #         tmp_selection.set_custom_selection(self._create_selection_string_for_representations())
  #         if self._view.ui.cb_protein_pair_sticks.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 0:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.STICKS.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.STICKS.value, tmp_selection.selection_string
  #             )
  #         elif self._view.tg_protein_pair_sticks.toggle_button.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 1:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.STICKS.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.STICKS.value, tmp_selection.selection_string
  #             )
  #         else:
  #             #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.STICKS.value)
  #             self._interface_manager.pymol_session_manager.hide_specific_representation(
  #                 enums.PyMOLRepresentation.STICKS.value, tmp_selection.selection_string
  #             )
  #         # self.__slot_protein_pair_chain_with_hydrogens()
  #         self._update_protein_pair_scene_legacy()
  #         self._save_protein_pair_pymol_session()
  #         self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()
  #         self._interface_manager.manage_hydrogen_representation_for_protein_pair_chain()
  #     except Exception as e:
  #         logger.error(f"An error occurred: {e}")
  #         self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
  #
  # def __slot_protein_pair_chain_as_ribbon(self) -> None:
  #     try:
  #         logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
  #                    "'Ribbon' toggle on the 'Protein Pairs Tab' was clicked.")
  #         tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
  #         tmp_selection.set_custom_selection(self._create_selection_string_for_representations())
  #         if self._view.ui.cb_protein_pair_ribbon.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 0:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.RIBBON.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.RIBBON.value, tmp_selection.selection_string
  #             )
  #         elif self._view.tg_protein_pair_ribbon.toggle_button.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 1:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.RIBBON.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.RIBBON.value, tmp_selection.selection_string
  #             )
  #         else:
  #             #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.RIBBON.value)
  #             self._interface_manager.pymol_session_manager.hide_specific_representation(
  #                 enums.PyMOLRepresentation.RIBBON.value, tmp_selection.selection_string
  #             )
  #         self._update_protein_pair_scene_legacy()
  #         self._save_protein_pair_pymol_session()
  #         self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()
  #         self._interface_manager.manage_hydrogen_representation_for_protein_pair_chain()
  #     except Exception as e:
  #         logger.error(f"An error occurred: {e}")
  #         self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
  #
  # def __slot_protein_pair_chain_as_lines(self) -> None:
  #     try:
  #         logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
  #                    "'Lines' toggle on the 'Protein Pairs Tab' was clicked.")
  #         # self.__slot_protein_pair_chain_with_hydrogens()
  #         tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
  #         tmp_selection.set_custom_selection(self._create_selection_string_for_representations())
  #         if self._view.ui.cb_protein_pair_lines.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 0:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.LINES.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.LINES.value, tmp_selection.selection_string
  #             )
  #         elif self._view.tg_protein_pair_lines.toggle_button.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 1:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.LINES.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.LINES.value, tmp_selection.selection_string
  #             )
  #         else:
  #             #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.LINES.value)
  #             self._interface_manager.pymol_session_manager.hide_specific_representation(
  #                 enums.PyMOLRepresentation.LINES.value, tmp_selection.selection_string
  #             )
  #         self._update_protein_pair_scene_legacy()
  #         self._save_protein_pair_pymol_session()
  #         self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()
  #         self._interface_manager.manage_hydrogen_representation_for_protein_pair_chain()
  #     except Exception as e:
  #         logger.error(f"An error occurred: {e}")
  #         self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
  #
  # def __slot_protein_pair_chain_as_spheres(self) -> None:
  #     try:
  #         logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
  #                    "'Spheres' toggle on the 'Protein Pairs Tab' was clicked.")
  #         tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
  #         tmp_selection.set_custom_selection(self._create_selection_string_for_representations())
  #         if self._view.ui.cb_protein_pair_spheres.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 0:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SPHERES.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.SPHERES.value, tmp_selection.selection_string
  #             )
  #         elif self._view.tg_protein_pair_spheres.toggle_button.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 1:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SPHERES.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.SPHERES.value, tmp_selection.selection_string
  #             )
  #         else:
  #             #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.SPHERES.value)
  #             self._interface_manager.pymol_session_manager.hide_specific_representation(
  #                 enums.PyMOLRepresentation.SPHERES.value, tmp_selection.selection_string
  #             )
  #         # self.__slot_protein_pair_chain_with_hydrogens()
  #         self._update_protein_pair_scene_legacy()
  #         self._save_protein_pair_pymol_session()
  #         self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()
  #         self._interface_manager.manage_hydrogen_representation_for_protein_pair_chain()
  #     except Exception as e:
  #         logger.error(f"An error occurred: {e}")
  #         self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
  #
  # def __slot_protein_pair_chain_as_dots(self) -> None:
  #     try:
  #         logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
  #                    "'Dots' toggle on the 'Protein Pairs Tab' was clicked.")
  #         tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
  #         tmp_selection.set_custom_selection(self._create_selection_string_for_representations())
  #         if self._view.ui.cb_protein_pair_dots.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 0:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.DOTS.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.DOTS.value, tmp_selection.selection_string
  #             )
  #         elif self._view.tg_protein_pair_dots.toggle_button.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 1:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.DOTS.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.DOTS.value, tmp_selection.selection_string
  #             )
  #         else:
  #             #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.DOTS.value)
  #             self._interface_manager.pymol_session_manager.hide_specific_representation(
  #                 enums.PyMOLRepresentation.DOTS.value, tmp_selection.selection_string
  #             )
  #         self._update_protein_pair_scene_legacy()
  #         self._save_protein_pair_pymol_session()
  #         self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()
  #         self._interface_manager.manage_hydrogen_representation_for_protein_pair_chain()
  #     except Exception as e:
  #         logger.error(f"An error occurred: {e}")
  #         self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
  #
  # def __slot_protein_pair_chain_as_mesh(self) -> None:
  #     try:
  #         logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
  #                    "'Mesh' toggle on the 'Protein Pairs Tab' was clicked.")
  #         tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
  #         tmp_selection.set_custom_selection(self._create_selection_string_for_representations())
  #         if self._view.ui.cb_protein_pair_mesh.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 0:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.MESH.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.MESH.value, tmp_selection.selection_string
  #             )
  #         elif self._view.tg_protein_pair_mesh.toggle_button.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 1:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.MESH.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.MESH.value, tmp_selection.selection_string
  #             )
  #         else:
  #             #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.MESH.value)
  #             self._interface_manager.pymol_session_manager.hide_specific_representation(
  #                 enums.PyMOLRepresentation.MESH.value, tmp_selection.selection_string
  #             )
  #         self._update_protein_pair_scene_legacy()
  #         self._save_protein_pair_pymol_session()
  #         self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()
  #         self._interface_manager.manage_hydrogen_representation_for_protein_pair_chain()
  #     except Exception as e:
  #         logger.error(f"An error occurred: {e}")
  #         self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")
  #
  # def __slot_protein_pair_chain_as_surface(self) -> None:
  #     try:
  #         logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
  #                    "'Surface' toggle on the 'Protein Pairs Tab' was clicked.")
  #         tmp_selection = self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
  #         tmp_selection.set_custom_selection(self._create_selection_string_for_representations())
  #         if self._view.ui.cb_protein_pair_surface.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 0:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SURFACE.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.SURFACE.value, tmp_selection.selection_string
  #             )
  #         elif self._view.tg_protein_pair_surface.toggle_button.isChecked() and self._interface_manager.get_protein_pair_repr_toggle_flag() == 1:
  #             #tmp_selection.show_selection_in_a_specific_representation(enums.PyMOLRepresentation.SURFACE.value)
  #             self._interface_manager.pymol_session_manager.show_specific_representation(
  #                 enums.PyMOLRepresentation.SURFACE.value, tmp_selection.selection_string
  #             )
  #         else:
  #             #tmp_selection.hide_selection_in_a_specific_representation(enums.PyMOLRepresentation.SURFACE.value)
  #             self._interface_manager.pymol_session_manager.hide_specific_representation(
  #                 enums.PyMOLRepresentation.SURFACE.value, tmp_selection.selection_string
  #             )
  #         self._update_protein_pair_scene_legacy()
  #         self._save_protein_pair_pymol_session()
  #         self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()
  #         self._interface_manager.manage_hydrogen_representation_for_protein_pair_chain()
  #     except Exception as e:
  #         logger.error(f"An error occurred: {e}")
  #         self._interface_manager.status_bar_manager.show_error_message("An unknown error occurred!")

  def __slot_hide_protein_pair_chain_all(self) -> None:
    """Hides all representations for the selected chain."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Hide all' representations button on the 'Protein Pairs Tab' was clicked.",
      )
      tmp_selection = (
        self._interface_manager.get_current_active_protein_object_of_protein_pair().pymol_selection
      )
      tmp_selection.set_custom_selection(
        self._create_selection_string_for_representations()
      )
      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=pymol_session_async.hide_all_representations,
              args=(
                tmp_selection.selection_string,
                self._interface_manager.pymol_session_manager,
              ),
            ),
            an_await_function=self.__await_set_representation_for_protein_pair_session,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
      # self._active_task = tasks.LegacyTask(
      #     target=pymol_session_async.hide_all_representations,
      #     args=(
      #         tmp_selection.selection_string,
      #         self._interface_manager.pymol_session_manager,
      #     ),
      #     post_func=self.__await_hide_all_representations_of_protein_chain_of_a_protein_pair,
      # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __await_hide_all_representations_of_protein_chain_of_a_protein_pair(
          self, return_value: tuple[bool]
  ) -> None:
    """Finishes hide all representations for the selected chain process.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      QtWidgets.QApplication.restoreOverrideCursor()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      if tmp_result[0]:
        ui_util.set_checked_async(
          self._view.tg_protein_pair_cartoon.toggle_button, False
        )
        ui_util.set_checked_async(
          self._view.tg_protein_pair_ribbon.toggle_button, False
        )
        ui_util.set_checked_async(
          self._view.tg_protein_pair_sticks.toggle_button, False
        )
        ui_util.set_checked_async(
          self._view.tg_protein_pair_lines.toggle_button, False
        )
        ui_util.set_checked_async(
          self._view.tg_protein_pair_spheres.toggle_button, False
        )
        ui_util.set_checked_async(
          self._view.tg_protein_pair_dots.toggle_button, False
        )
        ui_util.set_checked_async(
          self._view.tg_protein_pair_mesh.toggle_button, False
        )
        ui_util.set_checked_async(
          self._view.tg_protein_pair_surface.toggle_button, False
        )
        # self._view.tg_protein_pair_cartoon.toggle_button.setChecked(False)
        # self._view.tg_protein_pair_sticks.toggle_button.setChecked(False)
        # self._view.tg_protein_pair_ribbon.toggle_button.setChecked(False)
        # self._view.tg_protein_pair_lines.toggle_button.setChecked(False)
        # self._view.tg_protein_pair_spheres.toggle_button.setChecked(False)
        # self._view.tg_protein_pair_dots.toggle_button.setChecked(False)
        # self._view.tg_protein_pair_mesh.toggle_button.setChecked(False)
        # self._view.tg_protein_pair_surface.toggle_button.setChecked(False)
        self._update_protein_pair_scene_legacy()
        self._save_protein_pair_pymol_session()
        self._interface_manager.manage_coloring_by_element_option_for_protein_pair_chain()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      QtWidgets.QApplication.restoreOverrideCursor()

  # </editor-fold>

  # <editor-fold desc="Update scene">
  def __slot_update_protein_pair_scene(self) -> None:
    """Updates the current protein pair scene."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Update protein scene' button on the 'Protein Pairs Tab' was clicked.",
      )
      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=pymol_session_async.update_scene,
              args=(
                self._interface_manager.pymol_session_manager, 0
              ),
            ),
            an_await_function=self.__await_update_scene_for_protein_pair_session,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
      # self._active_task = tasks.LegacyTask(
      #     target=pymol_session_async.update_scene,
      #     args=(self._interface_manager.pymol_session_manager, 0),
      #     post_func=self.__await_update_scene_for_protein_pair_session,
      # )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._interface_manager.block_gui()

  def __await_update_scene_for_protein_pair_session(
          self, return_value: tuple[bool, str]
  ) -> None:
    """Finishes the update protein pair scene process.

    Args:
        return_value (tuple): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>
    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      tmp_success_flag, tmp_current_scene_name = tmp_result
      if tmp_current_scene_name == "_scratch_":
        if (
                not self._interface_manager.check_if_scratch_scene_exists_in_protein_pair_model()
        ):
          self._interface_manager.add_scratch_scene_to_protein_pair_model()
        ui_util.set_pymol_scene_name_into_label(
          tmp_current_scene_name, self._view.ui.lbl_pymol_protein_pair_scene
        )
        self._interface_manager.pymol_session_manager.current_scene_name = (
          tmp_current_scene_name
        )
        self._interface_manager.status_bar_manager.show_temporary_message(
          "PyMOL Scene _scratch_ updated.", a_timeout=1500
        )
      else:
        self._interface_manager.status_bar_manager.show_temporary_message(
          "PyMOL Scene updated.", a_timeout=1500
        )
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )
    else:
      self._save_protein_pair_pymol_session()
    finally:
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  # </editor-fold>

  def _update_protein_pair_scene_legacy(self) -> None:
    """Updates the current selected PyMOL scene.

    Notes:
        This is a legacy method that is used in many other methods!
        TODO: after discussion of workaround, this can be removed!
    """
    return
    # if self._interface_manager.pymol_session_manager.is_the_current_pymol_scene_base is False:
    #     self._interface_manager.pymol_session_manager.user_pymol_connector.scene(
    #         a_key="auto", an_action="update"
    #     )
    #     self._interface_manager.status_bar_manager.show_temporary_message("PyMOL Scene updated.", a_timeout=1500)
    # else:
    #     if not self._interface_manager.check_if_scratch_scene_exists_in_protein_pair_model():
    #         self._interface_manager.add_scratch_scene_to_protein_pair_model()
    #     self._interface_manager.pymol_session_manager.user_pymol_connector.scene(
    #         a_key="_scratch_", an_action="update"
    #     )
    #     self._interface_manager.pymol_session_manager.current_scene_name = "_scratch_"
    #     ui_util.set_pymol_scene_name_into_label(self._interface_manager.pymol_session_manager.current_scene_name,
    #                                             self._view.ui.lbl_pymol_protein_pair_scene)
    #     self._interface_manager.status_bar_manager.show_temporary_message("PyMOL Scene _scratch_ updated.", a_timeout=1500)

  # <editor-fold desc="Delete protein pair">
  def __slot_delete_protein_pair_from_project(self) -> None:
    """Deletes the selected protein pair."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Delete protein pair' button on the 'Protein Pairs Tab' was clicked.",
      )
      tmp_dialog = custom_message_box.CustomMessageBoxDelete(
        "Are you sure you want to delete this protein pair?",
        "Delete Protein Pair",
        custom_message_box.CustomMessageBoxIcons.WARNING.value,
      )
      tmp_dialog.exec_()
      response: bool = tmp_dialog.response
      if not response:
        return

      tmp_protein_pair: "protein_pair.ProteinPair" = (
        self._interface_manager.get_current_active_protein_pair_object()
      )
      if not self._interface_manager.pymol_session_manager.is_the_current_protein_pair_in_session(
              tmp_protein_pair.name
      ):
        tmp_return_value = ("generic_id", [(True, (True, 0))])
        self.__await_reinitialize_session_before_delete_protein_pair(tmp_return_value)
        return
      self._interface_manager.get_task_manager().append_task_result(
        task_result_factory.TaskResultFactory.run_task_result(
          a_task_result=task_result.TaskResult.from_action(
            an_action=action.Action(
              a_target=pymol_session_async.reinitialize_session,
              args=
              (self._interface_manager.pymol_session_manager,)
            ),
            an_await_function=self.__await_reinitialize_session_before_delete_protein_pair,
          ),
          a_task_scheduler=self._interface_manager.get_task_scheduler(),
        )
      )
      # self._active_task = tasks.LegacyTask(
      #     target=pymol_session_async.reinitialize_session,
      #     args=(self._interface_manager.pymol_session_manager,),
      #     post_func=self.__await_reinitialize_session_before_delete_protein_pair,
      # )
    except Exception as e:
      logger.error(
        f"An error occurred during the protein pair deletion process: {e}"
      )
      self._interface_manager.status_bar_manager.show_error_message(
        "Protein pair delete process failed!"
      )
    else:
      self._interface_manager.block_gui()

  def __await_reinitialize_session_before_delete_protein_pair(
          self, return_value: tuple[str, list[tuple[bool, tuple]]]
  ) -> None:
    """Removes protein pair from database and refreshes the main view.

    Args:
        return_value (tuple[str, list[tuple[bool, tuple]]]): The result data from the async method.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      self._interface_manager.status_bar_manager.show_error_message(
        "No data received!"
      )
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()
      return

    # </editor-fold>

    try:
      tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
      if tmp_result[0]:
        tmp_protein_pair: "protein_pair.ProteinPair" = (
          self._interface_manager.get_current_active_protein_pair_object()
        )
        tmp_database_operation = database_operation.DatabaseOperation(
          enums.SQLQueryType.DELETE_EXISTING_PROTEIN_PAIR,
          (0, tmp_protein_pair.get_id()),
        )
        self._database_thread.put_database_operation_into_queue(
          tmp_database_operation
        )
        self._interface_manager.get_current_project().delete_specific_protein_pair(
          tmp_protein_pair.name
        )
        self._interface_manager.watcher.remove_protein_pair(
          tmp_protein_pair.name
        )
        self._interface_manager.remove_protein_pair_from_protein_pairs_model()
        self._interface_manager.status_bar_manager.show_temporary_message(
          "The protein pair was successfully deleted."
        )
    except Exception as e:
      logger.error(
        f"An error occurred during the protein pair deletion process: {e}"
      )
      self._interface_manager.status_bar_manager.show_error_message(
        "Protein pair delete process failed!"
      )
    finally:
      self._view.ui.protein_pairs_tree_view.selectionModel().clearSelection()
      self._interface_manager.disable_protein_pairs_tab_buttons()
      self._interface_manager.stop_wait_cursor()
      self._interface_manager.refresh_main_view()

  # </editor-fold>

  def __slot_check_for_results(self) -> None:
    """Checks if the results summary menu action can be enabled."""
    if (
            self._view.ui.protein_pairs_tree_view.model()
                    .data(
              self._view.ui.protein_pairs_tree_view.currentIndex(), Qt.DisplayRole
            )
                    .find("_vs_")
            != -1
    ):
      self._view.ui.action_results_summary.setEnabled(True)
    else:
      self._view.ui.action_results_summary.setEnabled(False)

  def __slot_results_summary(self) -> None:
    """Shows the results summary dialog."""
    try:
      logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "Menu entry 'Results/Summary' clicked.",
      )
      tmp_protein_pair = self._view.ui.protein_pairs_tree_view.model().data(
        self._view.ui.protein_pairs_tree_view.currentIndex(),
        enums.ModelEnum.OBJECT_ROLE,
      )
      self._external_controller = results_view_controller.ResultsViewController(
        self._interface_manager,
        tmp_protein_pair,
        self._interface_manager.pymol_session_manager,
      )
      self._interface_manager.get_results_view().show()
    except Exception as e:
      logger.error(f"An error occurred: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
        "An unknown error occurred!"
      )

  # </editor-fold>
