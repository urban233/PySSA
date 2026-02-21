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
"""Main window controller class for the PySSA frontend application.

Authors: Martin Urban, Hannah Kullik

Version: 2.0.0
"""
import collections
import copy
import logging
import os
import pathlib
import shutil
import subprocess
import time
import json
from typing import Optional, Any
from urllib import request

# import pywinctl
from Bio import SeqRecord
from Bio.Seq import Seq

from src.pyssa.gui.qt import QtWidgets
from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.qt import Qt
from src.pyssa.gui.qt import QtGui

from src.auxiliary_pymol import auxiliary_pymol_client
from src.pyssa.gui.ui import icon_resources  # this import is used for the icons! DO NOT DELETE THIS
from src.pyssa.gui.ui.custom_context_menus import protein_tree_context_menu, protein_pair_tree_context_menu, \
    sequence_list_context_menu
from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.gui.ui.custom_widgets import job_entry
from src.pyssa.internal.thread.async_pyssa import util_async, custom_signals, project_async, image_async, \
    pymol_session_async, protein_async, sequence_async, protein_pair_async
from src.pyssa.controller import results_view_controller, rename_protein_view_controller, use_project_view_controller, \
    pymol_session_manager, add_sequence_view_controller, add_scene_view_controller, add_protein_view_controller, \
    settings_view_controller, predict_protein_view_controller, import_sequence_view_controller, \
    rename_sequence_view_controller, settings_manager
from src.pyssa.internal.data_structures import chain, job
from src.pyssa.internal.data_structures.data_classes import residue_color_config
from src.pyssa.gui.ui.dialogs import dialog_settings_global, dialog_tutorial_videos, dialog_about
from src.pyssa.internal.data_structures import project, settings, protein, protein_pair
from src.pyssa.internal.data_structures.data_classes import database_operation
from src.pyssa.internal.thread import database_thread
from src.pyssa.io_pyssa import filesystem_io
from src.pyssa.logging_pyssa import log_handlers, log_levels
from src.pyssa.util import constants, enums, exit_codes, tools, ui_util, exception, main_window_util
from src.pyssa.gui import main_window, app_state
from src.pyssa.controller import interface_manager, distance_analysis_view_controller, delete_project_view_controller, \
    create_project_view_controller, open_project_view_controller, database_manager, pyssa_objects_panel_controller
from src.pyssa.util import globals
from src.pyssa_pymol import pymol_enums
from src.pyssa.internal.pymol import pml_worker
from src.pyssa.internal.thread import thread_util
from src.tea.thread import tasks, task_result_factory, task_result, action

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class MainWindowController:
    """Controller class for the main window."""

    def __init__(self, a_main_window: "main_window.MainWindow") -> None:
        """Constructor.

        Args:
          a_main_window: Main window instance
        """
        # <editor-fold desc="Checks">
        # psa_comm_api.macros.REQUIRE_NOT_NONE(a_main_window, "a_main_window is None.")
        # </editor-fold>
        # <editor-fold desc="Instance attributes">
        # <editor-fold desc="Private">
        self._main_window = a_main_window
        self._settings_manager = settings_manager.SettingsManager()
        self._app_state = app_state.AppState(self._settings_manager, self.refresh_ui)
        self._dialog_controllers = {}
        # self._interface_manager = interface_manager.InterfaceManager(a_main_window)
        # self._database_manager = database_manager.DatabaseManager("")
        # self._database_manager.set_application_settings(
        #     self._interface_manager.get_application_settings()
        # )
        # self._database_thread: "database_thread.DatabaseThread" = (
        #     database_thread.DatabaseThread("")
        # )
        self._user_pymol = self._main_window.user_pymol
        self._pyssa_objects_panel_controller = pyssa_objects_panel_controller.PySSAObjectsPanelController(
            self._app_state,
            self._main_window.pyssa_objects_panel,
            self._user_pymol
        )
        self.feedback_timer = QtCore.QTimer()
        self.custom_progress_signal = custom_signals.ProgressSignal()
        self.abort_signal = custom_signals.AbortSignal()
        self.thread_pool = QtCore.QThreadPool()
        self.thread_pool.setMaxThreadCount(os.cpu_count())
        # </editor-fold>
        # </editor-fold>
        # self._init_main_window()
        self._connect_all_signals_with_their_slots()
        # self.aux_pymol_client = aux_pymol_client.AuxPyMOLClient(self.context)
        # self.protein_model = protein_model.ProteinModel()
        # self.protein_model.add_protein(self._user_pymol.get_cmd_module().get_model())
        # self._setup_statusbar()
        # self._init_generic_help_context_menus()
        self._init_main_window()
        self._setup_application_settings()

    # <editor-fold desc="Private methods">
    def _init_main_window(self) -> None:
        """Sets up the main window."""
        self._set_models_to_views()

    def _set_models_to_views(self) -> None:
        """Sets all models to their corresponding views."""
        # self._main_window.side_panel_molecule_objects.tree_view.setModel(self.protein_model)
        pass

    def _connect_all_signals_with_their_slots(self) -> None:
        """Connects all relevant widget signals with their appropriate slots."""
        # self._main_window.dialogClosed.connect(self.__slot_close_application)

        self._main_window.action_new_project.triggered.connect(self.__slot_create_project)
        self._main_window.action_open_project.triggered.connect(self.__slot_open_project)
        self._main_window.action_use_project.triggered.connect(self.__slot_use_project)
        self._main_window.action_delete_project.triggered.connect(self.__slot_delete_project)
        # self._main_window.action_import_project.triggered.connect(self.__slot_import_project)
        # self._main_window.action_export_project.triggered.connect(self.__slot_export_current_project)
        # self._main_window.action_close_project.triggered.connect(self.__slot_close_project)
        # TODO: Add the right slot method! ;)
        # self._main_window.action_exit_application.triggered.connect(self.)



        # # <editor-fold desc="Session ribbon slots">
        # # <editor-fold desc="Session slots">
        self._main_window.viewer_toolbar_actions.get("open_session").get_action().triggered.connect(
            self.__slot_open_session
        )
        # # </editor-fold>
        # # <editor-fold desc="Scene slots">
        self._main_window.viewer_toolbar_actions.get("create_scene").get_action().triggered.connect(
            self.__slot_save_scene
        )
        self._main_window.viewer_toolbar_actions.get("save_scene").get_action().triggered.connect(
            self.__slot_update_scene
        )
        self._main_window.viewer_toolbar_actions.get("delete_scene").get_action().triggered.connect(
            self.__slot_delete_scene
        )
        # # </editor-fold>
        # # <editor-fold desc="Representation slots">
        self._main_window.viewer_toolbar_actions.get("cartoon").get_action().triggered.connect(
            self.__slot_display_as_cartoon
        )
        self._main_window.cartoon_show_action.triggered.connect(self.__slot_show_as_cartoon)
        self._main_window.cartoon_hide_action.triggered.connect(self.__slot_hide_cartoon)

        self._main_window.viewer_toolbar_actions.get("sticks").get_action().triggered.connect(
            self.__slot_display_as_sticks
        )
        self._main_window.sticks_show_action.triggered.connect(self.__slot_show_as_sticks)
        self._main_window.sticks_hide_action.triggered.connect(self.__slot_hide_sticks)

        self._main_window.viewer_toolbar_actions.get("ribbon").get_action().triggered.connect(
            self.__slot_display_as_ribbon
        )
        self._main_window.ribbon_show_action.triggered.connect(self.__slot_show_as_ribbon)
        self._main_window.ribbon_hide_action.triggered.connect(self.__slot_hide_ribbon)

        self._main_window.viewer_toolbar_actions.get("lines").get_action().triggered.connect(
            self.__slot_display_as_lines
        )
        self._main_window.lines_show_action.triggered.connect(self.__slot_show_as_lines)
        self._main_window.lines_hide_action.triggered.connect(self.__slot_hide_lines)

        self._main_window.viewer_toolbar_actions.get("spheres").get_action().triggered.connect(
            self.__slot_display_as_spheres
        )
        self._main_window.spheres_show_action.triggered.connect(self.__slot_show_as_spheres)
        self._main_window.spheres_hide_action.triggered.connect(self.__slot_hide_spheres)

        self._main_window.viewer_toolbar_actions.get("dots").get_action().triggered.connect(
            self.__slot_display_as_dots
        )
        self._main_window.dots_show_action.triggered.connect(self.__slot_show_as_dots)
        self._main_window.dots_hide_action.triggered.connect(self.__slot_hide_dots)

        self._main_window.viewer_toolbar_actions.get("mesh").get_action().triggered.connect(
            self.__slot_display_as_mesh
        )
        self._main_window.mesh_show_action.triggered.connect(self.__slot_show_as_mesh)
        self._main_window.mesh_hide_action.triggered.connect(self.__slot_hide_mesh)
        self._main_window.viewer_toolbar_actions.get("surface").get_action().triggered.connect(
            self.__slot_display_as_surface
        )
        self._main_window.surface_show_action.triggered.connect(self.__slot_show_as_surface)
        self._main_window.surface_hide_action.triggered.connect(self.__slot_hide_surface)
        # # </editor-fold>
        # # <editor-fold desc="Color slots">
        self._main_window.viewer_toolbar_actions.get("color").get_action().triggered.connect(
            self.__slot_display_color_grid
        )
        # # <editor-fold desc="Color pads">
        self._main_window.color_grid.c_red.clicked.connect(lambda: self.__slot_apply_color("red"))
        self._main_window.color_grid.c_tv_red.clicked.connect(
            lambda: self.__slot_apply_color("tv_red")
        )
        self._main_window.color_grid.c_salomon.clicked.connect(
            lambda: self.__slot_apply_color("salmon")
        )
        self._main_window.color_grid.c_raspberry.clicked.connect(
            lambda: self.__slot_apply_color("raspberry")
        )

        self._main_window.color_grid.c_green.clicked.connect(
            lambda: self.__slot_apply_color("green")
        )
        self._main_window.color_grid.c_tv_green.clicked.connect(
            lambda: self.__slot_apply_color("tv_green")
        )
        self._main_window.color_grid.c_palegreen.clicked.connect(
            lambda: self.__slot_apply_color("palegreen")
        )
        self._main_window.color_grid.c_forest.clicked.connect(
            lambda: self.__slot_apply_color("forest")
        )

        self._main_window.color_grid.c_blue.clicked.connect(lambda: self.__slot_apply_color("blue"))
        self._main_window.color_grid.c_tv_blue.clicked.connect(
            lambda: self.__slot_apply_color("tv_blue")
        )
        self._main_window.color_grid.c_lightblue.clicked.connect(
            lambda: self.__slot_apply_color("lightblue")
        )
        self._main_window.color_grid.c_skyblue.clicked.connect(
            lambda: self.__slot_apply_color("skyblue")
        )

        self._main_window.color_grid.c_yellow.clicked.connect(
            lambda: self.__slot_apply_color("yellow")
        )
        self._main_window.color_grid.c_tv_yellow.clicked.connect(
            lambda: self.__slot_apply_color("tv_yellow")
        )
        self._main_window.color_grid.c_paleyellow.clicked.connect(
            lambda: self.__slot_apply_color("paleyellow")
        )
        self._main_window.color_grid.c_sand.clicked.connect(lambda: self.__slot_apply_color("sand"))

        self._main_window.color_grid.c_magenta.clicked.connect(
            lambda: self.__slot_apply_color("magenta")
        )
        self._main_window.color_grid.c_purple.clicked.connect(
            lambda: self.__slot_apply_color("purple")
        )
        self._main_window.color_grid.c_pink.clicked.connect(lambda: self.__slot_apply_color("pink"))
        self._main_window.color_grid.c_hotpink.clicked.connect(
            lambda: self.__slot_apply_color("hotpink")
        )

        self._main_window.color_grid.c_cyan.clicked.connect(lambda: self.__slot_apply_color("cyan"))
        self._main_window.color_grid.c_aquamarine.clicked.connect(
            lambda: self.__slot_apply_color("aquamarine")
        )
        self._main_window.color_grid.c_palecyan.clicked.connect(
            lambda: self.__slot_apply_color("palecyan")
        )
        self._main_window.color_grid.c_teal.clicked.connect(lambda: self.__slot_apply_color("teal"))

        self._main_window.color_grid.c_orange.clicked.connect(
            lambda: self.__slot_apply_color("orange")
        )
        self._main_window.color_grid.c_tv_orange.clicked.connect(
            lambda: self.__slot_apply_color("tv_orange")
        )
        self._main_window.color_grid.c_lightorange.clicked.connect(
            lambda: self.__slot_apply_color("lightorange")
        )
        self._main_window.color_grid.c_olive.clicked.connect(
            lambda: self.__slot_apply_color("olive")
        )

        self._main_window.color_grid.c_white.clicked.connect(
            lambda: self.__slot_apply_color("white")
        )
        self._main_window.color_grid.c_grey_70.clicked.connect(
            lambda: self.__slot_apply_color("grey70")
        )
        self._main_window.color_grid.c_grey_30.clicked.connect(
            lambda: self.__slot_apply_color("grey30")
        )
        self._main_window.color_grid.c_black.clicked.connect(
            lambda: self.__slot_apply_color("black")
        )
        # # </editor-fold>
        self._main_window.color_config.btn_color_by_elements.clicked.connect(
            self.__slot_apply_color_by_elements
        )
        self._main_window.color_config.btn_white_bg.clicked.connect(
            lambda: self.__slot_apply_bg_color("white")
        )
        self._main_window.color_config.btn_grey_bg.clicked.connect(
            lambda: self.__slot_apply_bg_color("grey40")
        )
        self._main_window.color_config.btn_black_bg.clicked.connect(
            lambda: self.__slot_apply_bg_color("black")
        )
        # # </editor-fold>
        # # <editor-fold desc="Selection slots">
        # self._main_window.show_sele_rb_panel_item.get_action().triggered.connect(
        #     self.__slot_show_sele
        # )
        # self._main_window.hide_sele_rb_panel_item.get_action().triggered.connect(
        #     self.__slot_hide_sele
        # )
        # self._main_window.clear_sele_rb_panel_item.get_action().triggered.connect(
        #     self.__slot_clear_sele
        # )
        # # </editor-fold>
        # # </editor-fold>

        # # </editor-fold>
        self._main_window.pyssa_objects_panel.tree_view.clicked.connect(
            self.__slot_select_pymol_object
        )
        # </editor-fold>
        # self.feedback_timer.setSingleShot(True)
        # self.feedback_timer.timeout.connect(self.update_protein_structure_tree_view) # TODO: Refactor the method update_protein_structure_tree_view() first!
        # </editor-fold>

    def _setup_application_settings(self):
        # self._application_settings = settings.Settings(constants.SETTINGS_DIR, constants.SETTINGS_FILENAME)
        if not os.path.exists(constants.SETTINGS_FULL_FILEPATH):
            constants.PYSSA_LOGGER.info(
                "Settings file not found, open configuration dialog."
            )
            self._settings_manager.settings.app_launch = 1
            self._settings_manager.settings.workspace_path = constants.DEFAULT_WORKSPACE_PATH

            if not os.path.exists(pathlib.Path(f"{constants.SCRATCH_DIR}")):
                os.mkdir(pathlib.Path(f"{constants.SCRATCH_DIR}"))
            if not os.path.exists(pathlib.Path(f"{constants.CACHE_DIR}")):
                os.mkdir(pathlib.Path(f"{constants.CACHE_DIR}"))
            tools.download_file(
                constants.DEMO_PROJECT_URL,
                str(pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects.zip")),
            )
            constants.PYSSA_LOGGER.info("Demo projects are getting extracted ...")
            import zipfile

            with zipfile.ZipFile(
                    pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects.zip"), "r"
            ) as zip_ref:
                zip_ref.extractall(
                    pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects")
                )
            constants.PYSSA_LOGGER.info(
                "Demo projects are downloaded and extracted.\n Import of demo projects started ...",
            )

            path_of_demo_projects = pathlib.Path(
                f"{constants.SETTINGS_DIR}/demo-projects"
            )
            for tmp_filename in os.listdir(path_of_demo_projects):
                # Copy db file into new workspace
                tmp_project_database_filepath = str(
                    pathlib.Path(
                        f"{self.get_application_settings().workspace_path}/{tmp_filename}",
                    ),
                )
                tmp_src_filepath = str(
                    pathlib.Path(f"{path_of_demo_projects}/{tmp_filename}")
                )
                shutil.copyfile(tmp_src_filepath, tmp_project_database_filepath)
            constants.PYSSA_LOGGER.info("Import process of demo projects finished.")
            constants.PYSSA_LOGGER.info("Serialize settings ...")
            self._settings_manager.settings.serialize_settings()
            constants.PYSSA_LOGGER.info("Serialize settings finished.")

        self._settings_manager.settings = main_window_util.setup_app_settings(
            self._settings_manager.settings
        )

    # </editor-fold>

    # <editor-fold desc="Public getters">
    def get_main_window(self) -> "main_window.MainWindow":
        """Returns the main window of the controller.

        Returns:
          The main window instance of the controller
        """
        return self._main_window

    # </editor-fold>

    def refresh_ui(self) -> None:
        """Sync every piece of the main window to the current AppState.

        Called automatically by AppState whenever state changes.
        Also safe to call manually at any time.

        This method is fully declarative and idempotent: it reads
        the current ``AppState`` and unconditionally sets every relevant
        widget property based on specific object-level rules (e.g., proteins or sequences).
        UI elements that are not applicable in a given state are **disabled** (never hidden).
        """
        has_project = self._app_state.has_open_project()
        project = self._app_state.project

        # Derived booleans from detailed project data.
        has_sequences = has_project and len(project.sequences) > 0
        has_proteins = has_project and len(project.proteins) > 0
        has_protein_pairs = has_project and len(project.protein_pairs) > 0
        has_any_objects = has_sequences or has_proteins or has_protein_pairs
        has_running_jobs = len(self._app_state._cold_dbs) > 0

        # -- Project menu actions ------------------------------------------
        # Actions that replace or manage the active project context are always available.
        self._main_window.action_new_project.setEnabled(True)
        self._main_window.action_open_project.setEnabled(True)
        self._main_window.action_use_project.setEnabled(True)
        self._main_window.action_delete_project.setEnabled(True)
        self._main_window.action_import_project.setEnabled(True)
        
        # Export and close specifically act upon the *current* project.
        self._main_window.action_export_project.setEnabled(has_project)
        self._main_window.action_close_project.setEnabled(has_project)
        # action_exit_application is always enabled.

        # -- Top-level menus -----------------------------------------------
        self._main_window.menuPrediction.setEnabled(has_project)
        # Prediction needs an input sequence to execute.
        self._main_window.action_predict_monomer.setEnabled(has_sequences)
        self._main_window.action_predict_multimer.setEnabled(has_sequences)
        # Aborting requires an active background/cold database job.
        self._main_window.action_abort_prediction.setEnabled(has_running_jobs)

        self._main_window.menuAnalysis.setEnabled(has_project)
        # Distance analysis operations computationally require 3D structure models.
        self._main_window.action_distance_analysis.setEnabled(has_proteins or has_protein_pairs)

        self._main_window.menuResults.setEnabled(has_project)
        # Results summaries aggregate data from protein pair analysis/predictions.
        self._main_window.action_results_summary.setEnabled(has_protein_pairs)

        self._main_window.menuImage.setEnabled(has_project)
        # Rendering commands mathematically require actual PyMOL coordinates.
        self._main_window.action_preview_image.setEnabled(has_proteins)
        self._main_window.action_ray_tracing_image.setEnabled(has_proteins)
        self._main_window.action_simple_image.setEnabled(has_proteins)

        self._main_window.menuHotspots.setEnabled(has_project)
        # Protein region generation acts upon 3D coordinates.
        self._main_window.action_protein_regions.setEnabled(has_proteins)
        # Settings and Help menus are always enabled.

        # -- Viewer toolbar actions ----------------------------------------
        # Base scene and session commands act on the project environment.
        _project_level_toolbar_keys = ["open_session", "create_scene", "save_scene", "delete_scene"]
        for key in _project_level_toolbar_keys:
            toolbar_action = self._main_window.viewer_toolbar_actions.get(key)
            if toolbar_action is not None:
                toolbar_action.get_action().setEnabled(has_project)

        # PyMOL representation tools need a 3D structural model in the wrapper.
        _protein_level_toolbar_keys = [
            "cartoon", "sticks", "ribbon", "lines", "spheres", "dots",
            "mesh", "surface", "color",
        ]
        for key in _protein_level_toolbar_keys:
            toolbar_action = self._main_window.viewer_toolbar_actions.get(key)
            if toolbar_action is not None:
                toolbar_action.get_action().setEnabled(has_proteins)

        # General viewer state indicators are active.
        _status_level_toolbar_keys = ["running_jobs", "notifications"]
        for key in _status_level_toolbar_keys:
            toolbar_action = self._main_window.viewer_toolbar_actions.get(key)
            if toolbar_action is not None:
                toolbar_action.get_action().setEnabled(True)

        # -- PySSA Objects Panel toolbar -----------------------------------
        panel = self._main_window.pyssa_objects_panel
        # Importing sequences or structural files requires an open project.
        panel.import_file_action.get_action().setEnabled(has_project)
        panel.add_sequence_action.get_action().setEnabled(has_project)
        # Exporting or deleting explicitly requires at least one object to export/delete.
        panel.export_file_action.get_action().setEnabled(has_any_objects)
        panel.delete_object_action.get_action().setEnabled(has_any_objects)

        # -- First-pass model binding --------------------------------------
        if self._app_state.is_first_pass():
            panel.tree_view.setModel(self._app_state.pyssa_objects_model)

        # -- Window title --------------------------------------------------
        if has_project:
            self._main_window.setWindowTitle(
                f"PySSA \u2014 {project.get_project_name()}"
            )
        else:
            self._main_window.setWindowTitle("PySSA")

    # <editor-fold desc="Slot methods">
    def __slot_create_project(self):
        if not self._dialog_controllers.__contains__("create_project"):
            self._dialog_controllers["create_project"] = create_project_view_controller.CreateProjectViewController(
                self._app_state
            )
        self._dialog_controllers["create_project"].restore_default_view()
        self._dialog_controllers["create_project"].get_view().show()

    def __slot_delete_project(self):
        if not self._dialog_controllers.__contains__("delete_project"):
            from src.pyssa.controller import delete_project_view_controller
            self._dialog_controllers["delete_project"] = delete_project_view_controller.DeleteProjectViewController(
                self._app_state
            )
        self._dialog_controllers["delete_project"].restore_default_view()
        self._dialog_controllers["delete_project"].get_view().show()

    def __slot_use_project(self):
        if not self._dialog_controllers.__contains__("use_project"):
            from src.pyssa.controller import use_project_view_controller
            self._dialog_controllers["use_project"] = use_project_view_controller.UseProjectViewController(
                self._app_state
            )
        self._dialog_controllers["use_project"].restore_default_view()
        self._dialog_controllers["use_project"].get_view().show()

    def __slot_open_project(self):
        if not self._dialog_controllers.__contains__("open_project"):
            self._dialog_controllers["open_project"] = open_project_view_controller.OpenProjectViewController(
                self._app_state
            )
        self._dialog_controllers["open_project"].restore_default_view()
        self._dialog_controllers["open_project"].get_view().show()

    def __slot_close_project(self):
        pass

    # <editor-fold desc="Left side panel slots">
    def __slot_close_left_side_panel(self) -> None:
        """Closes the left side panel using the ToolWindowLayout."""
        try:
            self._main_window.tool_window_layout.set_left_panel_hidden(True)
        except Exception as e:
            logger.error(f"Failed to close left side panel: {e}")

    def __slot_show_protein_structure_side_panel(self) -> None:
        """Toggle the protein structure left side panel using ToolWindowLayout.

        - If the left panel is visible and already showing the protein structure page,
          hide the left panel.
        - Otherwise, show the left panel and switch to the protein structure page.
        """
        try:
            target_index = enums.LeftSidePanel.PROTEIN_STRUCTURE
            layout = self._main_window.tool_window_layout
            is_hidden = layout.is_left_panel_hidden
            current_index = layout.left_stack.currentIndex()
            if not is_hidden and current_index == target_index:
                layout.set_left_panel_hidden(True)
            else:
                layout.left_stack.setCurrentIndex(target_index)
                layout.set_left_panel_hidden(False)
        except Exception as e:
            logger.error(f"Failed to toggle left side panel: {e}")

    # </editor-fold>

    # <editor-fold desc="Protein structure panel">
    def __slot_select_pymol_object(self) -> None:
        selected_indexes = self._main_window.pyssa_objects_panel.tree_view.selectionModel().selectedIndexes()
        selection_strings = []
        for idx in selected_indexes:
            selection_strings.append(
                self._main_window.pyssa_objects_panel.tree_view.model().construct_selection_string(
                    idx
                )
            )
        combined_selection_string = " or ".join(selection_strings)
        self._user_pymol.get_cmd_module().select("sele", combined_selection_string, enable=1)

    # </editor-fold>

    # <editor-fold desc="Session ribbon slots">
    # <editor-fold desc="Session slots">
    def __slot_open_session(self) -> None:
        """TODO: Change this implementation to correct session opening."""
        self._interface_manager.status_bar_manager.show_temporary_message("Session opened")
        with pml_worker.PmlWorker.session(pml_worker.PmlWorker.cache_user_session(self._user_pymol, "my_test")) as worker:
            worker.do("color", ("red", "all"), sync=True)
            worker.do("draw", ("800", "600"), sync=True)
            worker.do("png", ("test.png", ), sync=True)
        # print(pml_worker.one_shot_do(
        #     self._user_pymol, "my_test", "color", ("red", "all"), True
        # ))
        # tmp_session_path = pml_worker.PmlWorker.cache_session(self._user_pymol, "my_test")
        # tmp_worker = pml_worker.PmlWorker()
        # tmp_worker.start()
        # tmp_worker.set_session_path(tmp_session_path)
        # print(
        #     tmp_worker.do("color", ("red", "all"), True)
        # )
        # tmp_worker.stop()

        # self._user_pymol.get_cmd_module().save("test.pse")
        # self._pymol_worker_connection.send(
        #     worker_command.WorkerCommand(
        #         "test.pse", "get_model", ("", ), True
        #     )
        # )
        # print(self._pymol_worker_connection.recv())

    # </editor-fold>

    # # <editor-fold desc="Scene slots">
    def _save_a_scene(self, scene_name: str):
        pass

    def __slot_save_scene(self) -> None:
      """Saves a pymol scene."""
      text, ok_pressed = QtWidgets.QInputDialog.getText(
        self._main_window,
        "Save Scene",
        "Enter scene name:",
        # QtWidgets.QLineEdit,
        # ""
      )
      if ok_pressed and text != '':
        self._user_pymol.get_cmd_module().scene(key=text, action="store")
        tmp_image_filepath = self._save_a_scene(text)
        self._main_window.side_panel_pymol_scenes.scenes_list.add_scene(text, pathlib.Path(tmp_image_filepath))

    def __slot_recall_scene(self, an_item) -> None:
      """Recalls an already created PyMOL scene."""
      tmp_widget = self._main_window.side_panel_pymol_scenes.scenes_list.list_widget.itemWidget(
        an_item)
      self._user_pymol.get_cmd_module().scene(tmp_widget.label.text(), "recall")

    def __slot_update_scene(self):
        """Update the currently selected PyMOL scene and refresh its thumbnail."""
        try:
            scene_name = self._main_window.side_panel_pymol_scenes.scenes_list.get_scene_name_from_item(
                self._main_window.side_panel_pymol_scenes.scenes_list.get_selected_item()
            )
            self._user_pymol.get_cmd_module().scene(key=scene_name, action="update")
            image_filepath = self._save_a_scene(scene_name)
            self._main_window.side_panel_pymol_scenes.scenes_list.update_scene_thumbnail(
                scene_name, pathlib.Path(image_filepath)
            )
        except Exception as error:
            logger.error(f"Failed to update scene: {error}")
            QtWidgets.QMessageBox.critical(
                self._main_window,
                "Error",
                f"Failed to update scene: {error}"
            )

    def __slot_delete_scene(self):
        """Deletes the currently selected PyMOL scene and removes it from the list."""
        try:
            scene_name = self._main_window.side_panel_pymol_scenes.scenes_list.get_scene_name_from_item(
                self._main_window.side_panel_pymol_scenes.scenes_list.get_selected_item()
            )
            self._user_pymol.get_cmd_module().scene(key=scene_name, action="clear")
            self._main_window.side_panel_pymol_scenes.scenes_list.remove_scene(scene_name)
        except Exception as error:
            logger.error(f"Failed to update scene: {error}")
            QtWidgets.QMessageBox.critical(
                self._main_window,
                "Error",
                f"Failed to update scene: {error}"
            )

    # # </editor-fold>

    # <editor-fold desc="Representation slots">
    # <editor-fold desc="Cartoon representation">
    def __slot_display_as_cartoon(self) -> None:
        """Opens a QMenu to show/hide the cartoon representation."""
        try:
            self._main_window.cartoon_show_hide_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("cartoon"))
            )
        except Exception as e:
            logger.error(e.__str__())

    def __slot_show_as_cartoon(self) -> None:
        """Shows the `pyssa_sele` selection in cartoon representation."""
        self._user_pymol.get_cmd_module().show("cartoon", "sele")

    def __slot_hide_cartoon(self) -> None:
        """Hides the cartoon representation of the `sele` selection."""
        self._user_pymol.get_cmd_module().hide("cartoon", "sele")

    # </editor-fold>

    # <editor-fold desc="Sticks representation">
    def __slot_display_as_sticks(self) -> None:
        """Opens a QMenu to show/hide the sticks representation."""
        try:
            self._main_window.sticks_show_hide_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("sticks"))
            )
        except Exception as e:
            logger.error(e.__str__())

    def __slot_show_as_sticks(self) -> None:
        """Shows the `pyssa_sele` selection in sticks representation."""
        self._user_pymol.get_cmd_module().show("sticks", "sele")

    def __slot_hide_sticks(self) -> None:
        """Hides the sticks representation of the `pyssa_sele` selection."""
        self._user_pymol.get_cmd_module().hide("sticks", "sele")

    # </editor-fold>

    # <editor-fold desc="Ribbon representation">
    def __slot_display_as_ribbon(self) -> None:
        """Opens a QMenu to show/hide the ribbon representation."""
        try:
            self._main_window.ribbon_show_hide_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("ribbon"))
            )
        except Exception as e:
            logger.error(e.__str__())

    def __slot_show_as_ribbon(self) -> None:
        """Shows the `sele` selection in ribbon representation."""
        self._user_pymol.get_cmd_module().show("ribbon", "sele")

    def __slot_hide_ribbon(self) -> None:
        """Hides the ribbon representation of the `sele` selection."""
        self._user_pymol.get_cmd_module().hide("ribbon", "sele")
    # </editor-fold>

    # <editor-fold desc="Lines representation">
    def __slot_display_as_lines(self) -> None:
        """Opens a QMenu to show/hide the lines representation."""
        try:
            self._main_window.lines_show_hide_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("lines"))
            )
        except Exception as e:
            logger.error(e.__str__())

    def __slot_show_as_lines(self) -> None:
        """Shows the `sele` selection in lines representation."""
        self._user_pymol.get_cmd_module().show("lines", "sele")

    def __slot_hide_lines(self) -> None:
        """Hides the lines representation of the `sele` selection."""
        self._user_pymol.get_cmd_module().hide("lines", "sele")
    # </editor-fold>

    # <editor-fold desc="Spheres representation">
    def __slot_display_as_spheres(self) -> None:
        """Opens a QMenu to show/hide the spheres representation."""
        try:
            self._main_window.spheres_show_hide_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("spheres"))
            )
        except Exception as e:
            logger.error(e.__str__())

    def __slot_show_as_spheres(self) -> None:
        """Shows the `pyssa_sele` selection in spheres representation."""
        self._user_pymol.get_cmd_module().show("spheres", "sele")

    def __slot_hide_spheres(self) -> None:
        """Hides the spheres representation of the `sele` selection."""
        self._user_pymol.get_cmd_module().hide("spheres", "sele")

    # </editor-fold>

    # <editor-fold desc="Dots representation">
    def __slot_display_as_dots(self) -> None:
        """Opens a QMenu to show/hide the dots representation."""
        try:
            self._main_window.dots_show_hide_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("dots"))
            )
        except Exception as e:
            logger.error(e.__str__())

    def __slot_show_as_dots(self) -> None:
        """Shows the `sele` selection in dots representation."""
        self._user_pymol.get_cmd_module().show("dots", "sele")

    def __slot_hide_dots(self) -> None:
        """Hides the dots representation of the `sele` selection."""
        self._user_pymol.get_cmd_module().hide("dots", "sele")
    # </editor-fold>

    # <editor-fold desc="Mesh representation">
    def __slot_display_as_mesh(self) -> None:
        """Opens a QMenu to show/hide the mesh representation."""
        try:
            self._main_window.mesh_show_hide_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("mesh"))
            )
        except Exception as e:
            logger.error(e.__str__())

    def __slot_show_as_mesh(self) -> None:
        """Shows the `sele` selection in mesh representation."""
        self._user_pymol.get_cmd_module().show("mesh", "sele")

    def __slot_hide_mesh(self) -> None:
        """Hides the mesh representation of the `sele` selection."""
        self._user_pymol.get_cmd_module().hide("mesh", "sele")
    # </editor-fold>

    # <editor-fold desc="Surface representation">
    def __slot_display_as_surface(self) -> None:
        """Opens a QMenu to show/hide the surface representation."""
        try:
            self._main_window.surface_show_hide_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("surface"))
            )
        except Exception as e:
            logger.error(e.__str__())

    def __slot_show_as_surface(self) -> None:
        """Shows the `pyssa_sele` selection in surface representation."""
        self._user_pymol.get_cmd_module().show("surface", "sele")

    def __slot_hide_surface(self) -> None:
        """Hides the surface representation of the `sele` selection."""
        self._user_pymol.get_cmd_module().hide("surface", "sele")

    # </editor-fold>
    # </editor-fold>

    # <editor-fold desc="Color slots">
    def __slot_display_color_grid(self) -> None:
        """Opens the color grid at the position of the 'color' quick access button."""
        try:
            # Prefer the viewer toolbar quick access button position
            self._main_window.color_grid_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("color"))
            )
            # Fallback: open without explicit position
            self._main_window.color_grid_menu.exec()
        except Exception as e:
            logger.error(e.__str__())

    def __slot_apply_color(self, a_color_name) -> None:
        """Colors the default sele selection in the given color."""
        self._user_pymol.get_cmd_module().color(a_color_name, "sele")

    def __slot_apply_color_by_elements(self) -> None:
        """Colors the default sele selection in the given color."""
        self._user_pymol.get_cmd_module().color("atomic", "sele and not elem C")
        self._user_pymol.get_cmd_module().color("grey70", "sele and elem C")

    def __slot_apply_bg_color(self, a_color_name) -> None:
        """Colors the viewer background in the given color."""
        self._user_pymol.get_cmd_module().bg_color(a_color_name)
        match a_color_name:
            case "white":
                self._main_window.tool_window_layout.apply_viewer_background("#ffffff")
            case "grey40":
                self._main_window.tool_window_layout.apply_viewer_background("#666666")
            case "black":
                self._main_window.tool_window_layout.apply_viewer_background("#000000")
            case _:
                logger.error(f"The color name {a_color_name} is not available as bg color!")

    # </editor-fold>

    # <editor-fold desc="Selection slots">
    def __slot_show_sele(self) -> None:
        """Highlights the selection in PyMOL"""
        self._user_pymol.get_cmd_module().select("sele", enable=1)

    def __slot_hide_sele(self) -> None:
        """Hides the selection highlighting in PyMOL"""
        self._user_pymol.get_cmd_module().select("sele", enable=0)

    def __slot_clear_sele(self) -> None:
        self._user_pymol.get_cmd_module().select("sele", "none", enable=0)
        self.feedback_timer.start(100)

    # </editor-fold>
    # </editor-fold>

    # <editor-fold desc="PyMOL related">
    def __slot_pymol_single_left(self) -> None:
        """PyMOL single left click event."""
        try:
            self.feedback_timer.start(100)
        except Exception as e:
            print(e.__str__())

    # TODO: Refactor below
    def update_protein_structure_tree_view(self) -> None:
        """Updates the protein structure tree view in the side panel."""
        try:
            selection_strings = []
            self._user_pymol.get_cmd_module().select(
                "sele", enable=1
            )  # Highlights the selection even if clicked on the PyMOL "canvas"
            atoms = self._user_pymol.get_cmd_module().get_model("sele")
            for at in atoms.atom:
                selection_strings.append(
                    self.parse_selection_string(
                        f"/1nb1//{at.chain}/{at.resi}+{str(at.resn)}/{at.name}"  # TODO: The hard-coded 3bmp can be fixed by storing the active protein name in some sort of pymol manager
                    )
                )  # TODO: Needs more work!
            self._main_window.pyssa_objects_panel.tree_view.selectionModel().clearSelection()
            self.select_item(
                self._main_window.pyssa_objects_panel.tree_view,
                selection_strings,
            )
        except Exception as e:
            print(e.__str__())

    @staticmethod
    def parse_selection_string(selection_string):
        """Parse a PyMOL selection string into hierarchical components."""
        components = selection_string.strip("/").split("/")
        object_name = components[0]  # Object
        chain_id = components[2]  # Chain
        residue_info = components[3].replace("+", " - ")  # Residue (e.g., "22+ASP")
        atom_name = components[4]  # Atom
        return object_name, chain_id, residue_info, atom_name

    @staticmethod
    def select_item(tree_view, pymol_selection_strings):
        """Select an item in the tree view by its text."""
        model: "protein_model.ProteinModel" = tree_view.model()

        # --- Begin insert (should only run once!)
        tmp_numbering_map = (
            model.create_counts_map()
        )  # TODO: This needs to be moved in the loading routine of a protein!
        # --- end

        selection_model = tree_view.selectionModel()
        residue_atom_map = {}
        for tmp_pymol_selection_string in pymol_selection_strings:
            resi_value, atom_value, tmp_atom_index = model.search_single_selection(
                tmp_pymol_selection_string
            )
            if resi_value not in residue_atom_map.keys():
                residue_atom_map[resi_value] = []
            residue_atom_map[resi_value].append((atom_value, tmp_atom_index))
        for tmp_key in residue_atom_map.keys():
            # Checks if all atoms were selected
            if len(residue_atom_map[tmp_key]) == tmp_numbering_map["atom_counts"][tmp_key]:
                selection_model.select(
                    residue_atom_map[tmp_key][0][1].parent(),
                    selection_model.Select | selection_model.Rows,
                )
            else:
                for _, tmp_atom in residue_atom_map[tmp_key]:
                    selection_model.select(tmp_atom, selection_model.Select | selection_model.Rows)

    # </editor-fold>

    # </editor-fold>

    def shutdown_application_processes(self) -> None:
        """Closes all threads and process as well as the application itself."""
        # if not pyssa_constants.FRONTEND_ONLY:
        #   # TODO: Add correct pyssa_core logic here
        #   raise NotImplementedError()
        # self.aux_pymol_client.shutdown_service()
        # self._pymol_worker_connection.send(worker_command.WorkerCommand("", "shutdown", ()))
        # self._pymol_worker_process.join()
        self._main_window.close()

    def _get_viewer_tool_bar_action_pos(self, an_action) -> QtCore.QPoint:
        """Return a global point beneath the toolbar button for the given action.

        This prefers the visible viewer toolbar hosted by ToolWindowLayout. If unavailable,
        it falls back to the legacy MainWindow.viewer_toolbar. If no button can be located,
        it falls back to the current cursor position to ensure a valid QPoint is always returned.
        """
        try:
            # Prefer the active viewer toolbar inside the ToolWindowLayout
            viewer_toolbar = getattr(self._main_window.tool_window_layout, "viewer_toolbar", None)
        except Exception:
            viewer_toolbar = None
        if viewer_toolbar is None:
            # Fallback to legacy attribute if present
            viewer_toolbar = getattr(self._main_window, "viewer_toolbar", None)

        if an_action is not None and viewer_toolbar is not None:
            try:
                tmp_button = viewer_toolbar.get_tool_button_for_action(an_action)
                if tmp_button is not None and tmp_button.isVisible():
                    # Position just below the button
                    return tmp_button.mapToGlobal(tmp_button.rect().bottomLeft())
            except Exception:
                pass
        # Final fallback: cursor position
        return QtGui.QCursor.pos()
