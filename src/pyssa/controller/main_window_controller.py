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
import logging
import os
import pathlib
import shutil
import subprocess
from io import BytesIO
from typing import Union

import pymol
import requests

# import pywinctl

from src.pyssa.internal.pymol.pml_worker import PmlWorker
from src.pyssa.internal.pymol.pml_enums import PmlCommand
from src.pyssa.io_pyssa.db_pyssa.write_queue import WriteOperation, OperationType
from src.pyssa.io_pyssa import bio_data, filesystem_io
from src.pyssa.internal.data_structures import protein
from src.pyssa.util import enums

from src.pyssa.gui.qt import QtWidgets
from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.qt import QtGui

from src.pyssa.controller import settings_manager, create_project_view_controller, open_project_view_controller, \
    pyssa_objects_panel_controller, welcome_screen_view_controller, help_panel_controller, \
    status_bar_manager, job_popup_controller, predict_protein_view_controller, settings_view_controller, \
    distance_analysis_view_controller, results_view_controller
from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.gui.ui.custom_filters import help_event_filter
from src.pyssa.gui.ui.dialogs import dialog_about
from src.pyssa.internal import job_definitions
from src.pyssa.internal.data_structures import protein_pair, protein
from src.pyssa.internal.data_structures.data_classes import job_descriptor
from src.pyssa.internal.pymol import pml_worker
from src.pyssa.io_pyssa.db_pyssa import WriteOperation, OperationType
from src.pyssa.logging_pyssa import log_handlers, log_levels
from src.pyssa.model import job_model, selection_snapshot
from src.pyssa.util import constants, enums, tools, main_window_util
from src.pyssa.gui import main_window, app_state
from src.pyssa.gui.ui.custom_context_menus import tree_context_menu
from src.pyssa.internal.thread.thread_api import thread_runtime

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class _PyMOLClickFilter(QtCore.QObject):
    """Event filter that detects mouse-button releases on the PyMOL widget.

    When a left-button release is detected the associated *feedback_timer*
    is (re-)started so the controller can synchronise the PyMOL selection
    back into the tree view after a short debounce.
    """

    def __init__(self, a_feedback_timer: QtCore.QTimer) -> None:
        """Constructor.

        Args:
            a_feedback_timer: The single-shot timer to start on each click.
        """
        super().__init__()
        self._feedback_timer = a_feedback_timer

    def eventFilter(self, obj: QtCore.QObject, event: QtCore.QEvent) -> bool:
        """Intercept mouse-button release events.

        Args:
            obj: The watched object (the PyMOL widget).
            event: The incoming event.

        Returns:
            Always ``False`` so the event continues to propagate to PyMOL.
        """
        if event.type() == QtCore.QEvent.Type.MouseButtonRelease:
            if event.button() == QtCore.Qt.MouseButton.LeftButton:
                self._feedback_timer.start(250)
        return False


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
        self._status_bar_manager = status_bar_manager.StatusBarManager(self._main_window)
        self._app_state = app_state.AppState(self._settings_manager, self._status_bar_manager, self.refresh_ui)
        self._user_pymol = self._main_window.user_pymol
        self._dialog_controllers = {}
        self._pyssa_objects_panel_controller = pyssa_objects_panel_controller.PySSAObjectsPanelController(
            self._app_state,
            self._main_window.pyssa_objects_panel,
            self._user_pymol
        )
        self._help_panel_controller = help_panel_controller.HelpPanelController(
            self._main_window,
            self._main_window.help_panel
        )
        self._active_jobs_controller = job_popup_controller.JobPopupController(
            self._main_window.active_jobs, self._app_state.job_model
        )
        self._active_jobs_controller.show_active_jobs()
        self._complete_jobs_controller = job_popup_controller.JobPopupController(
            self._main_window.complete_jobs, self._app_state.job_model
        )
        self._complete_jobs_controller.show_completed_jobs()
        self.feedback_timer = QtCore.QTimer()
        self.feedback_timer.setSingleShot(True)
        self._is_syncing_selection: bool = False
        
        self._auto_save_timer = QtCore.QTimer()
        self._auto_save_timer.setSingleShot(True)
        self._auto_save_timer.timeout.connect(self.save_pymol_session_to_project)
        
        self._tree_context_menu = tree_context_menu.TreeContextMenu()
        self._register_tree_context_menu_actions()
        self._current_selection_snapshot: "selection_snapshot.SelectionSnapshot | None" = None
        # self.custom_progress_signal = custom_signals.ProgressSignal()
        # self.abort_signal = custom_signals.AbortSignal()
        # self.thread_pool = QtCore.QThreadPool()
        # self.thread_pool.setMaxThreadCount(os.cpu_count())
        # </editor-fold>
        # </editor-fold>
        # self._init_main_window()
        self._connect_all_signals_with_their_slots()
        # Install an event filter on the PyMOL widget to detect mouse clicks.
        # PyMOLGLWidget does not expose a click signal, so we catch
        # MouseButtonRelease at the Qt level and start the feedback timer.
        self._pymol_click_filter = _PyMOLClickFilter(self.feedback_timer)
        self._main_window.pymolwidget.installEventFilter(self._pymol_click_filter)
        # self.aux_pymol_client = aux_pymol_client.AuxPyMOLClient(self.context)
        # self.protein_model = protein_model.ProteinModel()
        # self.protein_model.add_protein(self._user_pymol.get_cmd_module().get_model())
        # self._setup_statusbar()
        # self._init_generic_help_context_menus()
        self._setup_application_settings()

        help_map = {
            "pyssa_objects_panel": "<h3>Save Button</h3><p>Saves the current document.</p>"
        }
        self.help_filter = help_event_filter.HelpEventFilter(self._main_window.help_panel.help_text_browser, help_map)
        self._main_window.pyssa_objects_panel.installEventFilter(self.help_filter)

        self.refresh_ui()
        self.open_welcome_screen()

    # <editor-fold desc="Private methods">
    def _connect_all_signals_with_their_slots(self) -> None:
        """Connects all relevant widget signals with their appropriate slots."""
        self._main_window.dialogClosed.connect(self.__slot_exit_application)

        # <editor-fold desc="Project menu">
        self._main_window.action_new_project.triggered.connect(self.__slot_create_project)
        self._main_window.action_open_project.triggered.connect(self.__slot_open_project)
        self._main_window.action_use_project.triggered.connect(self.__slot_use_project)
        self._main_window.action_delete_project.triggered.connect(self.__slot_delete_project)
        self._main_window.action_import_project.triggered.connect(self.__slot_import_project)
        self._main_window.action_export_project.triggered.connect(self.__slot_export_current_project)
        self._main_window.action_close_project.triggered.connect(self.__slot_close_project)
        self._main_window.action_exit_application.triggered.connect(
            self.__slot_exit_application
        )
        # </editor-fold>

        # <editor-fold desc="Prediction menu">
        self._main_window.action_predict_monomer.triggered.connect(self.__slot_predict_monomer)
        self._main_window.action_predict_multimer.triggered.connect(self.__slot_predict_multimer)
        # </editor-fold>

        # <editor-fold desc="Analysis menu">
        self._main_window.action_distance_analysis.triggered.connect(self.__slot_distance_analysis)
        # </editor-fold>

        # <editor-fold desc="Results menu">
        self._main_window.action_results_summary.triggered.connect(self.__slot_results_summary)
        # </editor-fold>

        # <editor-fold desc="Image menu">
        self._main_window.action_preview_image.triggered.connect(self.__slot_preview_image)
        self._main_window.action_simple_image.triggered.connect(self.__slot_draw_image)
        self._main_window.action_ray_tracing_image.triggered.connect(self.__slot_ray_trace_image)
        # </editor-fold>

        # <editor-fold desc="Hotspots menu">
        self._main_window.action_protein_regions.triggered.connect(
            self.__slot_protein_hotspots
        )
        # </editor-fold>

        # <editor-fold desc="Settings menu">
        self._main_window.action_edit_settings.triggered.connect(
            self.__slot_open_settings_dialog
        )
        self._main_window.action_restore_settings.triggered.connect(
            self.__slot_restore_settings
        )
        # </editor-fold>

        # <editor-fold desc="Help menu">
        self._main_window.action_documentation.triggered.connect(self.__slot_toggle_help_panel)
        self._main_window.action_show_log_in_explorer.triggered.connect(self.__slot_open_logs)
        self._main_window.action_clear_logs.triggered.connect(self.__slot_clear_all_log_files)
        self._main_window.action_get_demo_projects.triggered.connect(self.__slot_get_demo_projects)
        self._main_window.action_about.triggered.connect(self.__slot_open_about)
        # </editor-fold>

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
        self._main_window.viewer_toolbar_actions.get("clean").get_action().triggered.connect(
            self.__slot_display_clean_options
        )
        self._main_window.clean_solvent_action.triggered.connect(self.__slot_clean_solvent)
        self._main_window.clean_organic_action.triggered.connect(self.__slot_clean_organic)
        # # </editor-fold>

        self._main_window.viewer_toolbar_actions.get("running_jobs").get_action().triggered.connect(
            self.__slot_open_active_jobs_popup
        )
        self._main_window.viewer_toolbar_actions.get("notifications").get_action().triggered.connect(
            self.__slot_open_completed_jobs_popup
        )

        # # </editor-fold>
        tree_view = self._main_window.pyssa_objects_panel.tree_view
        tree_view.setContextMenuPolicy(
            QtCore.Qt.ContextMenuPolicy.CustomContextMenu
        )
        tree_view.customContextMenuRequested.connect(
            self.__slot_show_tree_context_menu
        )
        self._pyssa_objects_panel_controller.selectionSnapshotUpdated.connect(
            self.__slot_on_selection_snapshot_updated
        )
        self._main_window.pyssa_objects_panel.tree_view.clicked.connect(
            self.__slot_pyssa_objects_view_clicked
        )
        # </editor-fold>
        self.feedback_timer.timeout.connect(self.__slot_sync_pymol_selection_to_tree)
        self._app_state.job_model.job_finished.connect(self._handle_job_results)

    def _register_tree_context_menu_actions(self) -> None:
        """Register all context menu actions for the tree view.

        This is the single place to add, remove, or reorganise context menu
        entries.  To add a new item, call
        ``self._tree_context_menu.register_action(section, key, label, callback)``.
        The ``configure`` method (called from ``refresh_ui``) will automatically
        show the action whenever the matching selection context is active.

        Valid sections:
            ``'sequence'``: visible when sequences are selected.
            ``'standalone_protein'``: visible when standalone proteins are selected.
            ``'protein_pair'``: visible when protein pairs are selected.
            ``'protein_pair_child'``: visible when proteins inside a pair are selected.
        """
        tmp_context_menu = self._tree_context_menu

        # -- Sequence actions --------------------------------------------------
        # No dedicated rename slot exists in this controller yet; add here when
        # a RenameSequenceViewController is integrated into MainWindowController.

        # -- Standalone protein actions ----------------------------------------
        tmp_context_menu.register_action(
            section="standalone_protein",
            key="standalone_protein__open_session",
            label="Open Session",
            callback=self.__slot_open_session,
        )
        tmp_context_menu.register_action(
            section="standalone_protein",
            key="standalone_protein__clean_solvent",
            label="Clean Solvent Molecules",
            callback=self.__slot_clean_solvent,
        )
        tmp_context_menu.register_action(
            section="standalone_protein",
            key="standalone_protein__clean_organic",
            label="Clean Organic Molecules",
            callback=self.__slot_clean_organic,
        )

        # -- Protein pair actions -----------------------------------------------
        tmp_context_menu.register_action(
            section="protein_pair",
            key="protein_pair__open_session",
            label="Open Session",
            callback=self.__slot_open_session,
        )
        tmp_context_menu.register_action(
            section="protein_pair",
            key="protein_pair__results_summary",
            label="Open Results Summary",
            callback=self.__slot_results_summary,
        )

        # -- Protein pair child actions -----------------------------------------
        tmp_context_menu.register_action(
            section="protein_pair_child",
            key="protein_pair_child__results_summary",
            label="Open Results Summary",
            callback=self.__slot_results_summary,
        )

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

    def get_application_settings(self):
        """Returns the application settings.

        Returns:
            The application settings object
        """
        return self._settings_manager.settings

    # </editor-fold>

    # <editor-fold desc="Selection snapshot helper methods">
    def _get_current_snapshot(self) -> "selection_snapshot.SelectionSnapshot | None":
        """Returns the currently cached selection snapshot.

        Returns:
            The current selection snapshot or None if no selection exists.
        """
        return self._current_selection_snapshot

    def _get_selected_protein_names(self) -> list[str]:
        """Extract protein names from the current selection snapshot.

        Returns:
            A list of protein names from the currently selected proteins.
            Returns an empty list if no snapshot or proteins are selected.
        """
        snapshot = self._get_current_snapshot()
        if not snapshot:
            return []

        protein_names = []
        for protein in snapshot.distinct_proteins:
            protein_names.append(protein.get_molecule_object())
        return protein_names

    def _get_pymol_selection_string(self) -> str:
        """Build a PyMOL selection string from the current snapshot.

        Returns:
            A PyMOL selection string targeting the selected proteins.
            Returns "sele" as fallback if no proteins are selected.
        """
        protein_names = self._get_selected_protein_names()
        if not protein_names:
            return "sele"

        # Build selection string: "protein_name_1 or protein_name_2 or ..."
        return " or ".join(protein_names)

    def _has_valid_selection(self) -> bool:
        """Check if there is any valid selection in the current snapshot.

        Returns:
            True if proteins, chains, residues, or atoms are selected, False otherwise.
        """
        snapshot = self._get_current_snapshot()
        if not snapshot:
            return False

        return (len(snapshot.distinct_proteins) > 0 or
                len(snapshot.raw_chains) > 0 or
                len(snapshot.raw_residues) > 0 or
                len(snapshot.raw_atoms) > 0)

    def _apply_pymol_command_to_selection(self, command_name: str, *args) -> None:
        """Apply a PyMOL command to the current selection.

        This method applies PyMOL commands to proteins selected in the tree view.
        Falls back to the default "sele" selection if no proteins are selected.

        Args:
            command_name: Name of the PyMOL command (e.g., "show", "hide", "color")
            *args: Additional arguments for the PyMOL command
        """
        selection_string = self._get_pymol_selection_string()
        cmd = self._user_pymol.get_cmd_module()

        # Get the command method from pymol
        pymol_command = getattr(cmd, command_name, None)
        if pymol_command is None:
            logger.error(f"PyMOL command '{command_name}' not found")
            return

        # Apply the command with the selection string and additional arguments
        try:
            pymol_command(*args, selection_string)
        except Exception as e:
            logger.error(f"Failed to apply PyMOL command '{command_name}': {e}")

    def _get_protein_context(self, protein) -> dict:
        """Get context information about a protein (standalone or part of a pair).

        Args:
            protein: The protein to get context for.

        Returns:
            A dictionary containing:
            - 'is_standalone': bool - True if protein is standalone
            - 'is_pair_child': bool - True if protein is part of a pair
            - 'protein_pair': ProteinPair | None - The pair if protein is part of one
            - 'protein': Protein - The original protein object
        """
        snapshot = self._get_current_snapshot()
        if not snapshot:
            return {
                'is_standalone': False,
                'is_pair_child': False,
                'protein_pair': None,
                'protein': protein
            }

        is_standalone = protein in snapshot.raw_standalone_proteins
        is_pair_child = protein in snapshot.raw_protein_pair_children
        protein_pair = snapshot.get_protein_pair_for_protein(protein) if is_pair_child else None

        return {
            'is_standalone': is_standalone,
            'is_pair_child': is_pair_child,
            'protein_pair': protein_pair,
            'protein': protein
        }

    def _get_all_protein_contexts(self) -> list[dict]:
        """Get context information for all currently selected proteins.

        Returns:
            A list of context dictionaries (see _get_protein_context for structure).
            Each entry indicates whether the protein is standalone or part of a pair.
        """
        snapshot = self._get_current_snapshot()
        if not snapshot:
            return []

        contexts = []
        for protein in snapshot.distinct_proteins:
            contexts.append(self._get_protein_context(protein))

        return contexts

    def _get_selected_protein_pairs(self) -> list:
        """Get all protein pairs from the current selection.

        Returns:
            A list of ProteinPair objects that are currently selected.
        """
        snapshot = self._get_current_snapshot()
        if not snapshot:
            return []

        return list(snapshot.raw_protein_pairs)

    def _group_proteins_by_context(self) -> dict:
        """Group selected proteins by their context (standalone vs pair).

        Returns:
            A dictionary with keys:
            - 'standalone_proteins': list[Protein] - Standalone proteins
            - 'pair_proteins': list[dict] - Proteins that are part of pairs, each dict contains:
                - 'protein': Protein
                - 'protein_pair': ProteinPair
            - 'protein_pairs': list[ProteinPair] - All protein pairs involved
        """
        snapshot = self._get_current_snapshot()
        if not snapshot:
            return {
                'standalone_proteins': [],
                'pair_proteins': [],
                'protein_pairs': []
            }

        pair_proteins = []
        for protein in snapshot.raw_protein_pair_children:
            protein_pair = snapshot.get_protein_pair_for_protein(protein)
            if protein_pair:
                pair_proteins.append({
                    'protein': protein,
                    'protein_pair': protein_pair
                })

        return {
            'standalone_proteins': list(snapshot.raw_standalone_proteins),
            'pair_proteins': pair_proteins,
            'protein_pairs': list(snapshot.raw_protein_pairs)
        }

    # </editor-fold>

    def _trigger_auto_save(self) -> None:
        """Trigger the auto-save timer after a PyMOL-altering interaction."""
        self._auto_save_timer.start(5000)

    def open_welcome_screen(self):
        if not self._dialog_controllers.__contains__("welcome_screen"):
            self._dialog_controllers["welcome_screen"] = welcome_screen_view_controller.WelcomeScreenViewController(
                self._main_window, self._app_state
            )
        self._dialog_controllers["welcome_screen"].restore_default_view()
        self._dialog_controllers["welcome_screen"].get_view().show()

    def refresh_ui(self, snapshot: "selection_snapshot.SelectionSnapshot | None" = None) -> None:
        """Sync every piece of the main window to the current AppState.

        Called automatically by AppState whenever state changes, natively triggered
        by full model refreshes, or by the PySSAObjectsPanelController passing a
        `SelectionSnapshot` during user interaction.

        This method is fully declarative: it reads the current ``AppState`` and the
        passed selection snapshot to unconditionally set every relevant widget
        property (e.g. enabling predicting monomer only if monomers exist and/or are selected).
        """
        has_project = self._app_state.has_open_project()
        project = self._app_state.project

        # Derived booleans from detailed project data.
        has_sequences = has_project and len(project.sequences) > 0
        has_monomer_sequences_in_project = has_sequences and any("," not in s for s in project.sequences)
        has_multimer_sequences_in_project = has_sequences and any("," in s for s in project.sequences)

        has_proteins = has_project and len(project.proteins) > 0
        has_protein_pairs = has_project and len(project.protein_pairs) > 0
        has_any_objects = has_sequences or has_proteins or has_protein_pairs
        has_running_jobs = len(self._app_state._cold_dbs) > 0
        if self._user_pymol.get_currently_loaded_object() is None:
            has_loaded_session = False
        else:
            has_loaded_session = True

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

        # -- Selection Snapshot context ------------------------------------
        has_active_monomer_selection = snapshot.has_monomer_sequences if snapshot else False
        has_active_multimer_selection = snapshot.has_multimer_sequences if snapshot else False
        has_active_protein_selection = len(snapshot.distinct_proteins) > 0 if snapshot else False
        has_active_pair_selection = len(snapshot.raw_protein_pairs) > 0 if snapshot else False
        has_active_pair_child_selection = len(snapshot.raw_protein_pair_children) > 0 if snapshot else False
        has_active_scene_selection = len(snapshot.raw_scenes) > 0 if snapshot else False

        # -- Top-level menus -----------------------------------------------
        self._main_window.menuPrediction.setEnabled(has_project)
        # Prediction needs an input sequence to execute. Prefer active selection, fallback to project existence.
        can_predict_monomer = has_active_monomer_selection or has_monomer_sequences_in_project
        can_predict_multimer = has_active_multimer_selection or has_multimer_sequences_in_project
        self._main_window.action_predict_monomer.setEnabled(can_predict_monomer)
        self._main_window.action_predict_multimer.setEnabled(can_predict_multimer)

        self._main_window.menuAnalysis.setEnabled(has_project)
        # Distance analysis operations computationally require 3D structure models.
        can_analyze = has_active_protein_selection or has_active_pair_selection or (has_proteins or has_protein_pairs)
        self._main_window.action_distance_analysis.setEnabled(can_analyze)

        self._main_window.menuResults.setEnabled(has_project)
        # Results summaries aggregate data from protein pair analysis/predictions.
        self._main_window.action_results_summary.setEnabled(has_active_pair_selection or has_active_pair_child_selection or (has_protein_pairs and not snapshot))

        self._main_window.menuImage.setEnabled(has_project)
        # Rendering commands mathematically require actual PyMOL coordinates.
        can_image = has_active_protein_selection or has_active_pair_child_selection or (has_proteins and not snapshot)
        self._main_window.action_preview_image.setEnabled(can_image)
        self._main_window.action_ray_tracing_image.setEnabled(can_image)
        self._main_window.action_simple_image.setEnabled(can_image)

        self._main_window.menuHotspots.setEnabled(has_project)
        # Protein region generation acts upon 3D coordinates.
        self._main_window.action_protein_regions.setEnabled(can_image)
        # Settings and Help menus are always enabled.

        # -- Viewer toolbar actions ----------------------------------------
        # Base scene and session commands act on the project environment.
        toolbar_action_open_session = self._main_window.viewer_toolbar_actions.get("open_session")
        if toolbar_action_open_session: toolbar_action_open_session.get_action().setEnabled(has_project)

        toolbar_action_create_scene = self._main_window.viewer_toolbar_actions.get("create_scene")
        if toolbar_action_create_scene: toolbar_action_create_scene.get_action().setEnabled(
            has_project and can_image and has_loaded_session
        )

        toolbar_action_save_scene = self._main_window.viewer_toolbar_actions.get("save_scene")
        if toolbar_action_save_scene: toolbar_action_save_scene.get_action().setEnabled(
            has_active_scene_selection and has_loaded_session
        )

        toolbar_action_delete_scene = self._main_window.viewer_toolbar_actions.get("delete_scene")
        if toolbar_action_delete_scene: toolbar_action_delete_scene.get_action().setEnabled(
            has_active_scene_selection and has_loaded_session
        )

        # PyMOL representation tools need a 3D structural model in the wrapper.
        _protein_level_toolbar_keys = [
            "cartoon", "sticks", "ribbon", "lines", "spheres", "dots",
            "mesh", "surface", "color",
        ]
        for key in _protein_level_toolbar_keys:
            toolbar_action = self._main_window.viewer_toolbar_actions.get(key)
            if toolbar_action is not None:
                toolbar_action.get_action().setEnabled(can_image)

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
            # Reconnect selection signal after model is replaced
            self._pyssa_objects_panel_controller._connect_selection_signal()

        # -- Tree context menu ---------------------------------------------
        self._tree_context_menu.configure(snapshot)

        # -- Window title --------------------------------------------------
        if has_project:
            self._main_window.setWindowTitle(
                f"PySSA \u2014 {project.get_project_name()}"
            )
        else:
            self._main_window.setWindowTitle("PySSA")

    # <editor-fold desc="Slot methods">
    # <editor-fold desc="Project menu">
    def __slot_create_project(self):
        if not self._dialog_controllers.__contains__("create_project"):
            self._dialog_controllers["create_project"] = create_project_view_controller.CreateProjectViewController(
                self._app_state
            )
        self._dialog_controllers["create_project"].restore_default_view()
        self._dialog_controllers["create_project"].get_view().show()

    def __slot_open_project(self):
        if not self._dialog_controllers.__contains__("open_project"):
            self._dialog_controllers["open_project"] = open_project_view_controller.OpenProjectViewController(
                self._app_state
            )
        self._dialog_controllers["open_project"].restore_default_view()
        self._dialog_controllers["open_project"].get_view().show()

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

    def __slot_import_project(self) -> None:
        """Imports a project into the current workspace."""
        try:
            logger.info("Menu entry 'Project/Import' clicked.")
            file_dialog = QtWidgets.QFileDialog()
            desktop_path = QtCore.QDir.homePath()
            file_dialog.setDirectory(desktop_path)
            file_path, _ = file_dialog.getOpenFileName(
                self._main_window,
                "Select a project file to import",
                "",
                "Project Database File (*.db)",
            )
            if not file_path:
                return
            tmp_import_filepath = pathlib.Path(file_path)
            tmp_project_name_input_dialog = QtWidgets.QInputDialog()
            tmp_new_project_name, ok_pressed = tmp_project_name_input_dialog.getText(
                self._main_window,
                "Project Name",
                "Enter A Project Name:",
                text=tmp_import_filepath.name.replace(".db", ""),
            )
            if not ok_pressed or not tmp_new_project_name.strip():
                return
            tmp_new_project_name = tmp_new_project_name.strip()

            db_path = str(self._app_state.workspace.construct_project_db_path(tmp_new_project_name))

            from src.pyssa.io_pyssa.db_pyssa import ProjectDatabase
            from src.pyssa.model import psa_objects_model
            from src.pyssa.internal.thread.thread_api import thread_runtime

            def import_project_task(progress_callback, is_cancelled):
                shutil.copyfile(str(tmp_import_filepath), db_path)

                if is_cancelled():
                    raise InterruptedError("Cancelled during project import.")

                tmp_db = ProjectDatabase(db_path=db_path, project_id=tmp_new_project_name)

                # Assume the new project has an id of 1 in the freshly copied database.
                tmp_db.update_project_name(tmp_new_project_name, 1)

                from src.pyssa.internal.data_structures import project
                tmp_project = project.Project(tmp_new_project_name, pathlib.Path(self._app_state.get_settings().workspace_path))
                tmp_project.set_id(1)

                if is_cancelled():
                    tmp_db.close()
                    raise InterruptedError("Cancelled after configuring project.")

                tmp_pyssa_objects_model = psa_objects_model.PSAObjectsModel()
                tmp_pyssa_objects_model.build_model(tmp_project)
                return tmp_project, tmp_db, tmp_pyssa_objects_model

            def on_success(result):
                tmp_project, tmp_db, tmp_pyssa_objects_model = result
                # self._app_state.pyssa_objects_model = tmp_pyssa_objects_model
                # self._app_state.open_project(tmp_project, tmp_db)
                self._app_state._build_workspace_model()
                self.refresh_ui()
                self._app_state.status_bar_manager.show_permanent_message("", False)
                self._app_state.status_bar_manager.show_temporary_message("Project imported.")

            def on_error(exc):
                logger.exception("Failed to import project.", exc_info=exc)
                QtWidgets.QMessageBox.critical(
                    self._main_window,
                    "Failed to import project",
                    f"Could not import the project:\n{exc}",
                )
                self._app_state.status_bar_manager.show_error_message("Failed to import project.")

            (
                thread_runtime.get_singleton_thread_runtime()
                .run(import_project_task)
                .on_success(on_success)
                .on_error(on_error)
            )
            self._app_state.status_bar_manager.show_permanent_message(
                "Importing project ...", True
            )

        except Exception as e:
            logger.error(f"An error occurred: {e}")
            QtWidgets.QMessageBox.critical(
                self._main_window,
                "Error",
                "An unknown error occurred while importing!"
            )

    def __slot_export_current_project(self) -> None:
        """Exports the current project to an importable format."""
        try:
            logger.info("Menu entry 'Project/Export' clicked.")
            file_dialog = QtWidgets.QFileDialog()
            desktop_path = QtCore.QDir.homePath()
            file_dialog.setDirectory(desktop_path)
            file_path, _ = file_dialog.getSaveFileName(
                self._main_window,
                "Export current project",
                "",
                "Project Database File (*.db)",
            )
            if file_path:
                current_project_name = self._app_state.project.get_project_name()
                db_path = str(self._app_state.workspace.construct_project_db_path(current_project_name))

                from src.pyssa.internal.thread.thread_api import thread_runtime

                def export_project_task(progress_callback, is_cancelled):
                    shutil.copyfile(db_path, file_path)

                def on_success(result):
                    logger.info("Project exported successfully to %s", file_path)
                    self._app_state.status_bar_manager.show_permanent_message("", False)
                    self._app_state.status_bar_manager.show_temporary_message("Project exported.")

                def on_error(exc):
                    logger.exception("Failed to export project.", exc_info=exc)
                    QtWidgets.QMessageBox.critical(
                        self._main_window,
                        "Failed to export project",
                        f"Could not export the project:\n{exc}",
                    )
                    self._app_state.status_bar_manager.show_error_message("Failed to export project.")

                (
                    thread_runtime.get_singleton_thread_runtime()
                    .run(export_project_task)
                    .on_success(on_success)
                    .on_error(on_error)
                )
                self._app_state.status_bar_manager.show_permanent_message(
                    "Importing project ...", True
                )

        except Exception as e:
            logger.error(f"An error occurred: {e}")
            QtWidgets.QMessageBox.critical(
                self._main_window,
                "Error",
                "An unknown error occurred while exporting!"
            )

    def __slot_close_project(self):
        if self._app_state.has_open_project():
            self._app_state.close_project()
            self._user_pymol.get_cmd_module().reinitialize()
    # </editor-fold>

    # <editor-fold desc="Prediction menu">
    def _get_sequences_for_prediction(self, is_multimer: bool) -> list:
        """Return the sequence list to pre-populate the prediction dialog.

        Selection strategy:

        1. If sequences of the requested type are currently selected in the
           tree view, use exactly those selected sequences.
        2. Otherwise fall back to *all* project sequences of the requested
           type whose names are **not** reserved in the name registry.  This
           covers both proteins that already exist in the project and names
           claimed by queued or running jobs.

        Args:
            is_multimer: ``True`` to collect multimer sequences (SeqRecord.seq
                contains a comma), ``False`` for monomers.

        Returns:
            A list of ``SeqRecord`` objects ready to pass directly to
            ``PredictProteinViewController``.
        """
        project = self._app_state.project
        if project is None:
            return []

        # Helper: decides whether a SeqRecord belongs to the requested type.
        def _is_right_type(seq_record) -> bool:
            has_comma = "," in str(seq_record.seq)
            return has_comma if is_multimer else not has_comma

        snapshot = self._get_current_snapshot()
        if snapshot and snapshot.raw_sequences:
            # Convert the name-set from the snapshot to SeqRecord objects.
            selected_names: set[str] = snapshot.raw_sequences
            selected_sequences = [
                seq for seq in project.sequences
                if seq.name in selected_names and _is_right_type(seq)
            ]
            if selected_sequences:
                return selected_sequences

        # No relevant selection — fall back to all qualifying sequences.
        registry = self._app_state.name_registry
        from src.pyssa.gui import name_registry as name_registry_module
        return [
            seq for seq in project.sequences
            if _is_right_type(seq)
            and not registry.is_reserved(name_registry_module.PROTEIN, seq.name)
        ]

    def __slot_predict_monomer(self) -> None:
        """Opens the prediction dialog pre-populated with monomer sequences."""
        sequences = self._get_sequences_for_prediction(is_multimer=False)
        if not self._dialog_controllers.__contains__("predict_monomer"):
            self._dialog_controllers["predict_monomer"] = predict_protein_view_controller.PredictProteinViewController(
                self._app_state,
                sequences,
                a_parent=self._main_window,
            )
        else:
            # Re-create when called again so it reflects the current state.
            self._dialog_controllers["predict_monomer"] = predict_protein_view_controller.PredictProteinViewController(
                self._app_state,
                sequences,
                a_parent=self._main_window,
            )
        self._dialog_controllers["predict_monomer"].get_view().show()

    def __slot_predict_multimer(self) -> None:
        """Opens the prediction dialog pre-populated with multimer sequences."""
        sequences = self._get_sequences_for_prediction(is_multimer=True)
        if not self._dialog_controllers.__contains__("predict_multimer"):
            self._dialog_controllers["predict_multimer"] = predict_protein_view_controller.PredictProteinViewController(
                self._app_state,
                sequences,
                a_parent=self._main_window,
            )
        else:
            # Re-create when called again so it reflects the current state.
            self._dialog_controllers["predict_multimer"] = predict_protein_view_controller.PredictProteinViewController(
                self._app_state,
                sequences,
                a_parent=self._main_window,
            )
        self._dialog_controllers["predict_multimer"].get_view().show()
    # </editor-fold>

    # <editor-fold desc="Analysis menu">
    def __slot_distance_analysis(self):
        if not self._dialog_controllers.__contains__("distance_analysis_dialog"):
            self._dialog_controllers["distance_analysis_dialog"] = distance_analysis_view_controller.DistanceAnalysisViewController(
                self._app_state
            )
        self._dialog_controllers["distance_analysis_dialog"].restore_default_view()
        self._dialog_controllers["distance_analysis_dialog"].get_view().show()
    # </editor-fold>

    # <editor-fold desc="Results menu">
    def __slot_results_summary(self):
        if not self._dialog_controllers.__contains__("results_summary_dialog"):
            self._dialog_controllers["results_summary_dialog"] = results_view_controller.ResultsViewController(
                self._get_selected_protein_pairs()[0], self._app_state, self._user_pymol, self._main_window
            )
        self._dialog_controllers["results_summary_dialog"].restore_default_view()
        self._dialog_controllers["results_summary_dialog"].get_view().show()
    # </editor-fold>

    # <editor-fold desc="Image menu">
    def __slot_preview_image(self):
        self._user_pymol.get_cmd_module().ray(800, 600)

    def __slot_ray_trace_image(self):
        logger.log(
            log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
            "Menu entry 'Image/Ray' clicked.",
        )
        save_dialog = QtWidgets.QFileDialog()
        full_file_name = save_dialog.getSaveFileName(
            caption="Save Image", filter="Image (*.png)"
        )
        if full_file_name == ("", ""):
            logger.info("No file has been selected.")
            return

        self._app_state.job_scheduler.submit(
            job_descriptor.JobDescriptor(
                enums.JobType.RAY_TRACING,
                self._app_state.project.get_project_name(),
                display_name=pathlib.Path(full_file_name[0]).name,
                run_fn=job_definitions.run_ray_tracing_job,
                run_args=(
                    full_file_name[0],
                    pml_worker.PmlWorker.cache_user_session(
                        self._user_pymol, "simple_image"
                    ),
                    self._app_state.get_settings().image_ray_trace_mode,
                    self._app_state.get_settings().image_ray_texture,
                    self._app_state.get_settings().image_renderer
                )
            )
        )

    def __slot_draw_image(self):
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
            return

        self._app_state.job_scheduler.submit(
            job_descriptor.JobDescriptor(
                enums.JobType.SIMPLE_IMAGE,
                self._app_state.project.get_project_name(),
                display_name=pathlib.Path(full_file_name[0]).name,
                run_fn=job_definitions.run_simple_image_job,
                run_args=(
                    full_file_name[0],
                    pml_worker.PmlWorker.cache_user_session(
                        self._user_pymol, "simple_image"
                    )
                )
            )
        )

    # </editor-fold>

    # <editor-fold desc="Hotspots menu">
    def __slot_protein_hotspots(self):
        self._user_pymol.get_cmd_module().show("sticks", "sele")
        self._user_pymol.get_cmd_module().color("atomic", "sele and not elem C")
        self._user_pymol.get_cmd_module().color("grey70", "sele and elem C")
        self._user_pymol.get_cmd_module().zoom("sele")
    # </editor-fold>

    # <editor-fold desc="Settings menu">
    def __slot_open_settings_dialog(self):
        if not self._dialog_controllers.__contains__("settings_dialog"):
            self._dialog_controllers["settings_dialog"] = settings_view_controller.SettingsViewController(
                self._app_state
            )
        self._dialog_controllers["settings_dialog"].restore_default_view()
        self._dialog_controllers["settings_dialog"].get_view().show()

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
            tmp_dialog.exec()
            if tmp_dialog.response:
                tools.restore_default_settings(self._app_state.get_settings())
                self._status_bar_manager.show_temporary_message(
                    "Settings were successfully restored."
                )
                logging.info("Settings were successfully restored.")
            else:
                self._status_bar_manager.show_temporary_message(
                    "Settings were not modified."
                )
                logging.info("Settings were not modified.")
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._status_bar_manager.show_error_message("An unknown error occurred!")
    # </editor-fold>

    # <editor-fold desc="Help menu">
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
                self._main_window, "Select a log file to open", "", "LOG File (*.log)"
            )
            if file_path:
                os.startfile(file_path)
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._status_bar_manager.show_error_message(
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
            tmp_dialog.exec()
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
                    self._status_bar_manager.show_temporary_message(
                        "All log files could be deleted."
                    )
                    constants.PYSSA_LOGGER.info("All log files were deleted.")
                else:
                    tmp_dialog = custom_message_box.CustomMessageBoxOk(
                        "Not all log files could be deleted.",
                        "Clear Log Files",
                        custom_message_box.CustomMessageBoxIcons.WARNING.value,
                    )
                    tmp_dialog.exec()
                    constants.PYSSA_LOGGER.warning("Not all log files were deleted!")
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._status_bar_manager.show_error_message(
                "An unknown error occurred!"
            )

    def __slot_get_demo_projects(self) -> None:
        """Downloads, extracts, and integrates demo projects into the workspace.

        This method:
        1. Checks for internet connectivity
        2. Downloads the demo projects ZIP file from the configured URL
        3. Extracts the ZIP file to the settings directory
        4. Copies all project database files to the user's workspace
        5. Refreshes the workspace model to display the new projects

        All operations run asynchronously with proper error handling and user feedback.
        """
        try:
            logger.log(
                log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                "Menu entry 'Help/Get Demo Projects' clicked.",
            )

            # Check internet connectivity first
            if not tools.check_internet_connectivity():
                tmp_dialog = custom_message_box.CustomMessageBoxOk(
                    "You do not have a working internet connection\nbut that is necessary for this operation!",
                    "Internet Connection",
                    custom_message_box.CustomMessageBoxIcons.ERROR.value,
                )
                tmp_dialog.exec()
                return

            workspace_path = self._app_state.get_settings().workspace_path

            def download_demo_projects_task(progress_callback, is_cancelled):
                """Async task to download, extract and import demo projects."""
                import zipfile

                download_dest = pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects.zip")
                extract_dest = pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects")

                try:
                    # Step 1: Download the ZIP file
                    constants.PYSSA_LOGGER.info("Starting download of demo projects...")

                    if is_cancelled():
                        raise InterruptedError("Download cancelled by user.")

                    # Remove old files if they exist
                    if os.path.exists(download_dest):
                        os.remove(download_dest)
                    if os.path.exists(extract_dest):
                        shutil.rmtree(extract_dest)

                    # Download with streaming
                    response = requests.get(constants.DEMO_PROJECT_URL, stream=True, timeout=60)
                    response.raise_for_status()

                    with open(download_dest, 'wb') as file:
                        for chunk in response.iter_content(chunk_size=8192):
                            if is_cancelled():
                                raise InterruptedError("Download cancelled by user.")
                            if chunk:
                                file.write(chunk)

                    constants.PYSSA_LOGGER.info("Demo projects downloaded successfully.")

                    if is_cancelled():
                        raise InterruptedError("Operation cancelled by user.")

                    # Step 2: Extract the ZIP file
                    constants.PYSSA_LOGGER.info("Extracting demo projects...")

                    with zipfile.ZipFile(download_dest, "r") as zip_ref:
                        zip_ref.extractall(extract_dest)

                    constants.PYSSA_LOGGER.info("Demo projects extracted successfully.")

                    if is_cancelled():
                        raise InterruptedError("Operation cancelled by user.")

                    # Step 3: Import projects into workspace
                    constants.PYSSA_LOGGER.info("Importing demo projects into workspace...")

                    imported_count = 0
                    for tmp_filename in os.listdir(extract_dest):
                        if is_cancelled():
                            raise InterruptedError("Operation cancelled by user.")

                        # Only process .db files
                        if not tmp_filename.endswith('.db'):
                            continue

                        tmp_src_filepath = pathlib.Path(extract_dest / tmp_filename)
                        tmp_dest_filepath = pathlib.Path(workspace_path) / tmp_filename

                        # Copy database file to workspace
                        shutil.copyfile(str(tmp_src_filepath), str(tmp_dest_filepath))
                        imported_count += 1
                        constants.PYSSA_LOGGER.info(f"Imported project: {tmp_filename}")

                    constants.PYSSA_LOGGER.info(
                        f"Import process finished. {imported_count} demo project(s) imported."
                    )

                    # Clean up downloaded files
                    if os.path.exists(download_dest):
                        os.remove(download_dest)
                    if os.path.exists(extract_dest):
                        shutil.rmtree(extract_dest)

                    return (True, imported_count)

                except requests.exceptions.HTTPError as e:
                    constants.PYSSA_LOGGER.error(f"HTTP Error during download: {e}")
                    return (False, f"HTTP Error: {e}")
                except requests.exceptions.ConnectionError as e:
                    constants.PYSSA_LOGGER.error(f"Connection Error during download: {e}")
                    return (False, f"Connection Error: {e}")
                except requests.exceptions.Timeout as e:
                    constants.PYSSA_LOGGER.error(f"Timeout Error during download: {e}")
                    return (False, f"Timeout Error: {e}")
                except requests.exceptions.RequestException as e:
                    constants.PYSSA_LOGGER.error(f"Request Error during download: {e}")
                    return (False, f"Request Error: {e}")
                except zipfile.BadZipFile as e:
                    constants.PYSSA_LOGGER.error(f"Invalid ZIP file: {e}")
                    return (False, f"Invalid ZIP file: {e}")
                except InterruptedError as e:
                    constants.PYSSA_LOGGER.info(f"Operation cancelled: {e}")
                    return (False, "Operation cancelled by user.")
                except Exception as e:
                    constants.PYSSA_LOGGER.error(f"Unexpected error: {e}")
                    return (False, f"Unexpected error: {e}")

            def on_success(result):
                """Callback when download completes successfully."""
                success, data = result

                if success:
                    # Rebuild workspace model to show new projects
                    self._app_state._build_workspace_model()
                    self.refresh_ui()

                    self._app_state.status_bar_manager.show_permanent_message("", False)
                    self._app_state.status_bar_manager.show_temporary_message(
                        f"Demo projects downloaded and imported successfully. {data} project(s) added."
                    )
                else:
                    # Show error message
                    error_msg = str(data)
                    tmp_dialog = custom_message_box.CustomMessageBoxOk(
                        f"The download of the demo projects failed.\n\n{error_msg}",
                        "Get Demo Projects",
                        custom_message_box.CustomMessageBoxIcons.ERROR.value,
                    )
                    tmp_dialog.exec()

                    self._app_state.status_bar_manager.show_permanent_message("", False)
                    self._app_state.status_bar_manager.show_error_message(
                        "Failed to download demo projects."
                    )

            def on_error(exc):
                """Callback when download encounters an error."""
                logger.exception("Failed to download demo projects.", exc_info=exc)

                tmp_dialog = custom_message_box.CustomMessageBoxOk(
                    f"An error occurred while downloading demo projects:\n\n{exc}",
                    "Get Demo Projects",
                    custom_message_box.CustomMessageBoxIcons.ERROR.value,
                )
                tmp_dialog.exec()

                self._app_state.status_bar_manager.show_permanent_message("", False)
                self._app_state.status_bar_manager.show_error_message(
                    "Failed to download demo projects."
                )

            # Start the async operation
            (
                thread_runtime.get_singleton_thread_runtime()
                .run(download_demo_projects_task)
                .on_success(on_success)
                .on_error(on_error)
            )

            # Show progress message
            self._app_state.status_bar_manager.show_permanent_message(
                "Downloading demo projects ...", True
            )

        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._app_state.status_bar_manager.show_error_message(
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
            dialog.exec()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._status_bar_manager.show_error_message(
                "An unknown error occurred!"
            )
    # </editor-fold>

    def __slot_toggle_help_panel(self):
        layout = self._main_window.tool_window_layout
        layout.set_right_panel_hidden(not layout.is_right_panel_hidden)

    # <editor-fold desc="Session ribbon slots">
    # <editor-fold desc="Session slots">
    def __slot_open_session(self) -> None:
        """Opens a PyMOL session based on the current selection.

        Behavior depends on selection context:
        - If protein pair(s) selected: Opens session bound to the protein pair
        - If standalone protein(s) selected: Opens session bound directly to the protein
        - If no selection: Operates on all objects
        """
        # Get detailed protein context information
        grouped_context = self._group_proteins_by_context()
        protein_pairs = grouped_context['protein_pairs']
        standalone_proteins = grouped_context['standalone_proteins']
        pair_proteins = grouped_context['pair_proteins']

        # Determine session context and log appropriate information
        if protein_pairs:
            # Protein pairs are selected - session is bound to the pair
            for pair in protein_pairs:
                logger.info(f"Opening session for protein pair: {pair.name}")
                logger.info(f"  - Protein 1: {pair.protein_1.get_molecule_object()}")
                logger.info(f"  - Protein 2: {pair.protein_2.get_molecule_object()}")
                self._user_pymol.load_session(pair.pymol_session, pair)

        elif pair_proteins:
            # Individual proteins from pairs are selected - identify their parent pairs
            for protein_info in pair_proteins:
                protein = protein_info['protein']
                protein_pair = protein_info['protein_pair']
                logger.info(f"Opening session for protein {protein.get_molecule_object()} "
                           f"(part of pair: {protein_pair.name})")
                self._user_pymol.load_session(protein_pair.pymol_session, protein_pair)

        elif standalone_proteins:
            # Standalone proteins are selected - session is bound directly to protein
            for protein in standalone_proteins:
                logger.info(f"Opening session for standalone protein: {protein.get_molecule_object()}")
                self._user_pymol.load_session(protein.pymol_session, protein)

        else:
            # No selection - operate on all objects
            logger.info("Opening session for all objects (no specific selection)")

        # Example implementation (replace with actual logic)
        # with pml_worker.PmlWorker.session(pml_worker.PmlWorker.cache_user_session(self._user_pymol, "my_test")) as worker:
        #     worker.do("color", ("red", "all"), sync=True)
        #     worker.do("draw", ("800", "600"), sync=True)
        #     worker.do("png", ("test.png", ), sync=True)
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
    def __slot_save_scene(self) -> None:
        """Saves a PyMOL scene.

        The scene includes the current view and all visible objects.
        Selection context is available via _get_current_snapshot() if needed.
        """
        tmp_input_dialog = QtWidgets.QInputDialog()
        tmp_name, ok_pressed = tmp_input_dialog.getText(
            self._main_window,
            "Scene Name",
            "Enter A Scene Name:",
            text="",
        )
        if not ok_pressed or not tmp_name.strip():
            return
        tmp_scene_name = tmp_name.strip()

        self._user_pymol.get_cmd_module().scene(key=tmp_scene_name, action="append")

        self._app_state.pyssa_objects_model.add_scene(
            tmp_scene_name, self._user_pymol.get_currently_loaded_object()
        )

        # # Log selection context for debugging
        # snapshot = self._get_current_snapshot()
        # if snapshot and len(snapshot.distinct_proteins) > 0:
        #     logger.info(f"Scene saved with {len(snapshot.distinct_proteins)} protein(s) selected")

    def __slot_recall_scene(self) -> None:
        """Recalls an already created PyMOL scene."""
        snapshot = self._get_current_snapshot()
        if snapshot:
            for tmp_raw_scene in snapshot.raw_scenes:
                self._user_pymol.get_cmd_module().scene(tmp_raw_scene, "recall")
                # Log scene recall
                logger.info(f"Recalled scene: {tmp_raw_scene}")

    def __slot_update_scene(self):
        """Update the currently selected PyMOL scene and refresh its thumbnail."""
        snapshot = self._get_current_snapshot()
        if snapshot:
            for tmp_raw_scene in snapshot.raw_scenes:
                self._user_pymol.get_cmd_module().scene(tmp_raw_scene, "update")
                # Log scene recall
                logger.info(f"Updated scene: {tmp_raw_scene}")

    def __slot_delete_scene(self):
        """Deletes the currently selected PyMOL scene and removes it from the list."""
        snapshot = self._get_current_snapshot()
        if snapshot:
            for tmp_raw_scene in snapshot.raw_scenes:
                self._user_pymol.get_cmd_module().scene(tmp_raw_scene, "clear")
                # Log scene recall
                logger.info(f"Updated scene: {tmp_raw_scene}")
                self._app_state.pyssa_objects_model.remove_scene(
                    tmp_raw_scene, self._user_pymol.get_currently_loaded_object()
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
        self._trigger_auto_save()

    def __slot_hide_cartoon(self) -> None:
        """Hides the cartoon representation of the `sele` selection."""
        self._user_pymol.get_cmd_module().hide("cartoon", "sele")
        self._trigger_auto_save()

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
        self._trigger_auto_save()

    def __slot_hide_sticks(self) -> None:
        """Hides the sticks representation of the `pyssa_sele` selection."""
        self._user_pymol.get_cmd_module().hide("sticks", "sele")
        self._trigger_auto_save()

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
        self._trigger_auto_save()

    def __slot_hide_ribbon(self) -> None:
        """Hides the ribbon representation of the `sele` selection."""
        self._user_pymol.get_cmd_module().hide("ribbon", "sele")
        self._trigger_auto_save()
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
        self._trigger_auto_save()

    def __slot_hide_lines(self) -> None:
        """Hides the lines representation of the `sele` selection."""
        self._user_pymol.get_cmd_module().hide("lines", "sele")
        self._trigger_auto_save()
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
        self._trigger_auto_save()

    def __slot_hide_spheres(self) -> None:
        """Hides the spheres representation of the `sele` selection."""
        self._user_pymol.get_cmd_module().hide("spheres", "sele")
        self._trigger_auto_save()

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
        self._trigger_auto_save()

    def __slot_hide_dots(self) -> None:
        """Hides the dots representation of the `sele` selection."""
        self._user_pymol.get_cmd_module().hide("dots", "sele")
        self._trigger_auto_save()
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
        self._trigger_auto_save()

    def __slot_apply_color_by_elements(self) -> None:
        """Colors the default sele selection in the given color."""
        self._user_pymol.get_cmd_module().color("atomic", "sele and not elem C")
        self._user_pymol.get_cmd_module().color("grey70", "sele and elem C")
        self._trigger_auto_save()

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
        self._trigger_auto_save()

    # </editor-fold>

    def __slot_display_clean_options(self):
        try:
            self._main_window.clean_solvent_organic_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("clean"))
            )
        except Exception as e:
            logger.error(e.__str__())


    def __slot_clean_solvent(self):
        active_object = self._user_pymol.get_currently_loaded_object()
        self._user_pymol.get_cmd_module().remove("solvent")
        self._run_protein_structure_update_async(active_object)

    def __slot_clean_organic(self):
        active_object = self._user_pymol.get_currently_loaded_object()
        self._user_pymol.get_cmd_module().remove("organic")
        self._run_protein_structure_update_async(active_object)

    def _run_protein_structure_update_async(self, active_object):
        if not active_object or not isinstance(active_object, protein.Protein):
            self._trigger_auto_save()
            return

        # Retrieve current session string to pass to the worker
        session_str = self._user_pymol.save_session()
        
        # Show loading indicator in status bar
        self._app_state.status_bar_manager.show_permanent_message("Synchronizing structure data...", True)
        
        def background_task(progress_callback, is_cancelled):
            # 1. Update structure via PmlWorker
            tmp_reply_data = PmlWorker.one_shot_do(
                PmlCommand.CLEAN_PROTEIN_UPDATE_STRUCTURE,
                args=(session_str, active_object.get_molecule_object())
            )
            if not tmp_reply_data or tmp_reply_data[0] == "":
                raise ValueError("Clean protein failed inside PyMOL.")
                
            new_session, tmp_pdb_filepath = tmp_reply_data
            
            # 2. Parse new PDB data
            tmp_pdb_data, tmp_more_than_one_ca = bio_data.parse_pdb_file(tmp_pdb_filepath)
            
            # 3. Update the protein object in memory
            active_object.pymol_session = new_session
            active_object.set_pdb_data(tmp_pdb_data)
            
            # 4. Use new Database API (ProjectWriteQueue) to persist changes
            hot_db = self._app_state.hot_db
            if hot_db:
                hot_db.write_queue.submit(
                    WriteOperation(OperationType.UPDATE_PROTEIN_SESSION, active_object)
                )
                
                # Update PDB atoms
                hot_db.write_queue.submit(
                    WriteOperation(OperationType.UPDATE_PROTEIN_PDB_DATA, active_object) # We need to verify if this operation type exists, else we write a custom operation
                )

                # Find and remove non-protein chains
                for tmp_chain in active_object.chains:
                    if tmp_chain.chain_type == enums.ChainTypeEnum.NON_PROTEIN_CHAIN.value:
                        hot_db.write_queue.submit(
                            WriteOperation(OperationType.DELETE_CHAIN, (active_object.get_id(), tmp_chain.get_id()))
                        )
                
            return "Success"

        def on_success(result):
            # Trigger a UI refresh to rebuild the tree view
            if isinstance(active_object, protein.Protein):
                self._app_state.pyssa_objects_model.update_protein(active_object)
            self._trigger_auto_save()
            self._app_state.status_bar_manager.show_permanent_message("", False)
            
        def on_error(exc):
            logger.exception("Failed to update structure data.", exc_info=exc)
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                f"An error occurred while updating the structure:\n\n{exc}",
                "Update Structure",
                custom_message_box.CustomMessageBoxIcons.ERROR.value,
            )
            tmp_dialog.exec()
            self._app_state.status_bar_manager.show_permanent_message("", False)
            
        (
            thread_runtime.get_singleton_thread_runtime()
            .run(background_task)
            .on_success(on_success)
            .on_error(on_error)
        )

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

    # <editor-fold desc="Job popups">
    def __slot_open_active_jobs_popup(self):
        self._main_window.active_jobs_menu.exec(
            self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("running_jobs"))
        )

    def __slot_open_completed_jobs_popup(self):
        self._main_window.complete_jobs_menu.exec(
            self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("notifications"))
        )
    # </editor-fold>

    def __slot_on_selection_snapshot_updated(self, snapshot: "selection_snapshot.SelectionSnapshot") -> None:
        """Slot that receives the SelectionSnapshot when the project tree selection changes.

        This snapshot simplifies evaluating what the user currently has selected in the PySSA Objects Panel.
        Instead of traversing tree nodes, UI logic can query the snapshot directly.

        Args:
            snapshot: An immutable snapshot of the current tree selection.
        """
        # Cache the snapshot for use in other slot methods
        self._current_selection_snapshot = snapshot
        # Skip PyMOL sync when the tree is being updated from the PyMOL side
        # to prevent an infinite feedback loop.
        if not self._is_syncing_selection:
            if snapshot.pymol_selection_string:
                try:
                    self._user_pymol.get_cmd_module().select(
                        "sele",
                        selection=snapshot.pymol_selection_string,
                        enable=1,
                    )
                except pymol.CmdException:
                    logger.warning(f"The PyMOL selection string '{snapshot.pymol_selection_string}' is invalid!")
            else:
                # No 3-D selectable items chosen; clear PyMOL selection
                self._user_pymol.get_cmd_module().select("sele", "none", enable=0)

        # We defer all UI enabling/disabling logic to refresh_ui, allowing a single
        # source of truth for the entire application state.
        self.refresh_ui(snapshot)

    def __slot_pyssa_objects_view_clicked(self):
        snapshot = self._get_current_snapshot()
        if snapshot and snapshot.raw_scenes:
            self.__slot_recall_scene()

    def __slot_show_tree_context_menu(self, pos: QtCore.QPoint) -> None:
        """Show the tree view context menu at the right-click position.

        The menu items are already configured by the most recent ``refresh_ui``
        call, so this slot only needs to translate the local widget position to
        global screen coordinates and delegate display to the menu object.

        Args:
            pos: Local-widget-space position of the right-click, provided
                 automatically by the ``customContextMenuRequested`` signal.
        """
        tree_view = self._main_window.pyssa_objects_panel.tree_view
        self._tree_context_menu.show_at(tree_view.mapToGlobal(pos))

    def __slot_sync_pymol_selection_to_tree(self) -> None:
        """Sync the current PyMOL 'sele' selection into the tree view.

        Triggered by the ``feedback_timer`` after a PyMOL left-click event.
        Reads the atoms in ``sele`` via ``cmd.get_model``, finds the
        corresponding tree-view indexes through
        :meth:`PSAObjectsModel.find_indexes_for_pymol_atoms`, and updates
        the tree selection accordingly.

        The panel controller's selection signal is suppressed during the
        update to prevent a feedback loop (PyMOL → tree → PyMOL).
        """
        try:
            self._is_syncing_selection = True
            cmd = self._user_pymol.get_cmd_module()
            chempy_model = cmd.get_model("sele")

            model = self._app_state.pyssa_objects_model
            tree_view = self._main_window.pyssa_objects_panel.tree_view
            selection_model = tree_view.selectionModel()
            if selection_model is None:
                return

            # Suppress tree → PyMOL sync while we modify the tree selection.
            self._pyssa_objects_panel_controller.suppress_selection_signal()
            try:
                selection_model.clearSelection()
                indexes = model.find_indexes_for_pymol_atoms(chempy_model)
                for index in indexes:
                    selection_model.select(
                        index,
                        selection_model.SelectionFlag.Select | selection_model.SelectionFlag.Rows,
                    )
            finally:
                self._pyssa_objects_panel_controller.restore_selection_signal()

            # Because we suppressed the panel signal, the normal
            # snapshot → refresh_ui path was skipped.  Resolve a snapshot
            # from the tree's current selection and refresh the UI manually.
            current_indexes = list(selection_model.selectedIndexes())
            snapshot = model.resolve_selection(current_indexes)
            self._current_selection_snapshot = snapshot
            self.refresh_ui(snapshot)
        except Exception as e:
            logger.error(f"Failed to sync PyMOL selection to tree: {e}")
        finally:
            self._is_syncing_selection = False

    # </editor-fold>

    # <editor-fold desc="Handle job results">
    def _handle_job_results(self, descriptor: "job_descriptor.JobDescriptor", result: dict):
        match descriptor.job_type:
            case enums.JobType.PREDICTION:
                self._handle_prediction_job_result(descriptor, result)
            case enums.JobType.DISTANCE_ANALYSIS:
                self._handle_distance_analysis_job_result(descriptor, result)
            case enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS:
                self._handle_prediction_and_distance_analysis_job_result(descriptor, result)
            case enums.JobType.RAY_TRACING:
                self._handle_ray_tracing_job_result(descriptor, result)
                self._status_bar_manager.show_temporary_message(
                    "Raytracing job completed."
                )
            case _:
                logger.warning(f"Unhandled job type: {descriptor.job_type}")

    # <editor-fold desc="Handle specific job results">
    def _handle_distance_analysis_job_result(
            self,
            descriptor: "job_descriptor.JobDescriptor",
            result: dict[str, Union[list["protein_pair.ProteinPair"], bool]]
    ):
        if descriptor.is_hot:
            # Results belong to the hot project
            for tmp_protein_pair in result["protein_pairs"]:
                self._app_state.project.add_protein_pair(tmp_protein_pair)
                self._app_state.pyssa_objects_model.add_protein_pair(tmp_protein_pair)
                self._app_state.hot_db.write_queue.submit(
                    WriteOperation(OperationType.INSERT_PROTEIN_PAIR, tmp_protein_pair)
                )
        else:
            # Results belong to a cold project
            for tmp_protein_pair in result["protein_pairs"]:
                descriptor.cold_handle.submit(
                    WriteOperation(OperationType.INSERT_PROTEIN_PAIR, tmp_protein_pair)
                )

    def _handle_prediction_job_result(
            self,
            descriptor: "job_descriptor.JobDescriptor",
            result: dict[str, Union[list["protein.Protein"], bool]]
    ):
        if descriptor.is_hot:
            current_project_id = self._app_state.project.get_id()
            for tmp_protein in result["predicted_proteins"]:
                tmp_protein.db_project_id = current_project_id
                self._app_state.project.add_existing_protein(tmp_protein)
                self._app_state.pyssa_objects_model.add_protein(tmp_protein)
                self._app_state.hot_db.write_queue.submit(
                    WriteOperation(OperationType.INSERT_PROTEIN, tmp_protein)
                )
        else:
            project_id = descriptor.cold_handle._db.get_project_id(descriptor.project_name)
            for tmp_protein in result["predicted_proteins"]:
                tmp_protein.db_project_id = project_id
                descriptor.cold_handle.submit(
                    WriteOperation(OperationType.INSERT_PROTEIN, tmp_protein)
                )

    def _handle_prediction_and_distance_analysis_job_result(
            self,
            descriptor: "job_descriptor.JobDescriptor",
            result: dict[str, Union[bool, list["protein.Protein"], list["protein_pair.ProteinPair"]]]
    ):
        if descriptor.is_hot:
            for tmp_protein in result["predicted_proteins"]:
                self._app_state.project.add_existing_protein(tmp_protein)
                self._app_state.pyssa_objects_model.add_protein(tmp_protein)
                self._app_state.hot_db.write_queue.submit(
                    WriteOperation(OperationType.INSERT_PROTEIN, tmp_protein)
                )
            for tmp_protein_pair in result["protein_pairs"]:
                self._app_state.project.add_protein_pair(tmp_protein_pair)
                self._app_state.pyssa_objects_model.add_protein_pair(tmp_protein_pair)
                self._app_state.hot_db.write_queue.submit(
                    WriteOperation(OperationType.INSERT_PROTEIN_PAIR, tmp_protein_pair)
                )
        else:
            # Results belong to a cold project
            for tmp_protein in result["predicted_proteins"]:
                descriptor.cold_handle.submit(
                    WriteOperation(OperationType.INSERT_PROTEIN, tmp_protein)
                )
            for tmp_protein_pair in result["protein_pairs"]:
                descriptor.cold_handle.submit(
                    WriteOperation(OperationType.INSERT_PROTEIN_PAIR, tmp_protein_pair)
                )

    def _handle_ray_tracing_job_result(
            self,
            descriptor: "job_descriptor.JobDescriptor",
            result: dict
    ):
        # Until now, there is no implementation needed after the ray-tracing job.
        # However, it might be helpful to have such a method scaffold for
        # later use.
        pass
    # </editor-fold>
    # </editor-fold>

    def save_pymol_session_to_project(self) -> None:
        """Saves the current PyMOL session to the active project object and database.

        Retrieves the current PyMOL session as a base64 string, assigns it to the
        currently loaded Protein or ProteinPair object, and asynchronously updates
        the corresponding record in the project database.
        """
        current_object = self._user_pymol.get_currently_loaded_object()
        if not current_object:
            return

        try:
            session_str = self._user_pymol.save_session()
            current_object.pymol_session = session_str

            from src.pyssa.io_pyssa.db_pyssa.write_queue import WriteOperation, OperationType
            from src.pyssa.internal.data_structures.protein import Protein
            from src.pyssa.internal.data_structures.protein_pair import ProteinPair

            hot_db = self._app_state.hot_db
            if hot_db:
                if isinstance(current_object, Protein):
                    hot_db.write_queue.submit(
                        WriteOperation(OperationType.UPDATE_PROTEIN_SESSION, current_object)
                    )
                elif isinstance(current_object, ProteinPair):
                    hot_db.write_queue.submit(
                        WriteOperation(OperationType.UPDATE_PAIR_SESSION, current_object)
                    )

            logger.info("Saved PyMOL session to project database successfully.")
            self._app_state.status_bar_manager.show_temporary_message("PyMOL session saved.")

        except Exception as e:
            logger.error(f"Failed to capture or trigger PyMOL session save: {e}")

    def __slot_exit_application(self) -> None:
        """Closes all threads and process as well as the application itself."""
        tmp_message = "Are you sure you want to close PySSA?"
        tmp_jobs_are_running = self._app_state.job_scheduler.has_running_jobs()
        if tmp_jobs_are_running:
            tmp_message = "There are still jobs running.\nAre you sure you want to close PySSA?\n\n The progress of the running job(s) are lost!"
        tmp_dialog = custom_message_box.CustomMessageBoxYesNo(
            tmp_message,
            "Close PySSA",
            custom_message_box.CustomMessageBoxIcons.WARNING.value,
        )
        tmp_dialog.exec()
        if tmp_dialog.response:
            if tmp_jobs_are_running:
                subprocess.run(["wsl", "--terminate", "almaColabfold9"], creationflags=subprocess.CREATE_NO_WINDOW)
                filesystem_io.FilesystemCleaner.clean_prediction_scratch_folder()
                constants.PYSSA_LOGGER.info("Shutdown of wsl environment.")
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
