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
from io import BytesIO
from typing import Union

import pymol
import requests

# import pywinctl

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
            "pyssa_objects_panel": """
                <p>PySSA Objects Panel:</p>
                <p>Displays all PySSA objects in the current project:
                <ul>
                    <li>Sequences</li>
                    <li>Proteins (with chains, residues and atoms)</li>
                    <li>Protein Pairs (with proteins including chains, residues and atoms)</li>
                </ul>
            """,

          # <editor-fold desc="PyMOL Viewer">
            "pymolwidget": """
                <p>PyMOL Viewer:</p>
                <p>This is the interactive 3D molecular viewer powered by PyMOL.
                It displays protein structures of a protein or protein pair in 3D.</p>
                <p>Navigation in the viewer:</p>
                <ul>
                    <li>Left-click and drag down to zoom in; drag up to zoom out.</li>
                    <li>Right-click and drag to rotate the protein or protein pair.</li>
                    <li>Hold the middle mouse button and drag to move the
                     protein or protein pair.</li>
                    <li>Scroll the mouse wheel up to reveal more of the structure
                    (reduce clipping); scroll down to hide more.</li>
                </ul>
                <p>Loading structures:</p>
                <ol>
                    <li>Select a protein or protein pair in the PySSA Objects Panel.</li>
                    <li>Click Open Session in the toolbar, or right-click and choose
                    Open Session from the context menu.</li>
                </ol>
            """,

          # <editor-fold desc="Viewer toolbar buttons">
          # <editor-fold desc="Session/Scene">
          "viewer_open_session": """
                <p>Open Session:</p>
                <p>Load the stored PyMOL session for the currently selected protein or
                protein pair into the PyMOL viewer.</p>
                <ol>
                    <li>Select a protein or protein pair in the PySSA Objects Panel.</li>
                    <li>Click Open Session in the toolbar, or right-click and choose
                    Open Session from the context menu.</li>
                    <li>This leads to load its 3D structure into the PyMOL viewer.</li>
                </ol>
            """,

            "viewer_create_scene": """
                <p>Create Scene:</p>
                <p>Save the current viewport state as a new named scene.</p>
                <ol>
                    <li>Set up the desired camera angle, zoom, and representation in the PyMOL viewer.</li>
                    <li>Click Create Scene.</li>
                    <li>Enter a name for the new scene.</li>
                    <li>The scene appears under Scenes under the protein or
                    protein pair in the PySSA Objects Panel.</li>
                </ol>
            """,

            "viewer_save_scene": """
                <p>Save Scene:</p>
                <p>Overwrite the currently active scene with the present viewport state.</p>
                <ol>
                    <li>Select the scene you want to update under Scenes under
                    the protein or protein pair in the PySSA Objects Panel.</li>
                    <li>Adjust the viewport to the desired state.</li>
                    <li>Click Save Scene to overwrite the existing scene data.</li>
                </ol>
            """,

            "viewer_delete_scene": """
                <p>Delete Scene:</p>
                <p>Remove the currently selected scene.</p>
                <ol>
                    <li>Select the scene to delete under Scenes under the
                    protein or protein pair in the PySSA Objects Panel.</li>
                    <li>Click Delete Scene.</li>
                    <li>The scene is permanently removed.</li>
                </ol>
                <p style="color: red;"><b>&#9888; CAUTION:</b></p>
                <p>This action cannot be undone.</p>
            """,
          # </editor-fold>

          # <editor-fold desc="Representations">
          "viewer_cartoon": """
                <p>Cartoon Representation:</p>
                <p>Switch the selected protein to cartoon representation.</p>
                <p>Cartoon shows secondary structure elements (helices as ribbons,
                sheets as flat arrows, loops as thin tubes).
                This is the most common representation for protein overview.</p>
                <ol>
                    <li>Select a protein or a chain in the PySSA Objects Panel.</li>
                    <li>Click Cartoon to apply the representation.</li>
                    <li>Use Show or Hide from the dropdown to toggle visibility.</li>
                </ol>
            """,

            "viewer_sticks": """
                <p>Sticks Representation:</p>
                <p>Switch the selected structure to sticks representation.</p>
                <p>Sticks draws every bond as a cylinder, making individual bonds,
                side chains, and small molecules clearly visible.</p>
                <ol>
                    <li>Select a protein, residue, or atom in the PySSA Objects Panel.</li>
                    <li>Click Sticks to apply the representation.</li>
                    <li>Use Show or Hide from the dropdown to toggle visibility.</li>
                </ol>
            """,

            "viewer_ribbon": """
                <p>Ribbon Representation:</p>
                <p>Switch the selected protein to ribbon representation.</p>
                <p>Ribbon shows the protein backbone as a smooth continuous ribbon,
                giving an overview of the fold without showing individual atoms.</p>
                <ol>
                    <li>Select a protein or chain.</li>
                    <li>Click Ribbon to apply the representation.</li>
                    <li>Use Show or Hide from the dropdown to toggle visibility.</li>
                </ol>
            """,

            "viewer_lines": """
                <p>Lines Representation:</p>
                <p>Switch the selected structure to lines representation.</p>
                <p>Lines draws bonds as thin lines. Useful for a lightweight
                overview when rendering many atoms at once.</p>
                <ol>
                    <li>Select a protein or chain.</li>
                    <li>Click Lines to apply the representation.</li>
                    <li>Use Show or Hide from the dropdown to toggle visibility.</li>
                </ol>
            """,

            "viewer_spheres": """
                <p>Spheres Representation:</p>
                <p>Switch the selected structure to Van der Waals spheres representation.</p>
                <p>Each atom is drawn as a sphere scaled to its Van der Waals radius.
                Useful for visualising molecular surface and packing.</p>
                <ol>
                    <li>Select a protein or residue.</li>
                    <li>Click Spheres to apply the representation.</li>
                    <li>Use Show or Hide from the dropdown to toggle visibility.</li>
                </ol>
            """,

            "viewer_dots": """
                <p>Dots Representation:</p>
                <p>Switch the selected structure to dots representation.</p>
                <p>Each atom is shown as a dot cloud at its Van der Waals radius.
                Gives a lightweight surface-like appearance.</p>
                <ol>
                    <li>Select a protein or residue.</li>
                    <li>Click Dots to apply the representation.</li>
                    <li>Use Show or Hide from the dropdown to toggle visibility.</li>
                </ol>
            """,

            "viewer_mesh": """
                <p>Mesh Representation:</p>
                <p>Switch the selected structure to mesh (wireframe surface) representation.</p>
                <p>A semi-transparent 3D mesh is drawn around the molecular surface,
                useful for showing the overall shape while keeping interior visible.</p>
                <ol>
                    <li>Select a protein or chain.</li>
                    <li>Click Mesh to apply the representation.</li>
                    <li>Use Show or Hide from the dropdown to toggle visibility.</li>
                </ol>
            """,

            "viewer_surface": """
                <p>Surface Representation:</p>
                <p>Switch the selected structure to solid molecular surface representation.</p>
                <p>The solvent-accessible surface is rendered as a solid skin.
                Ideal for visualising binding pockets and surface features.</p>
                <ol>
                    <li>Select a protein or chain.</li>
                    <li>Click Surface to apply the representation.</li>
                    <li>Use Show or Hide from the dropdown to toggle visibility.</li>
                </ol>
            """,
          # </editor-fold>

            "viewer_color": """
                <p>Color:</p>
                <p>Apply a colour scheme to the selected protein or structure in the viewport.</p>
                <ol>
                    <li>Select a protein or chain in the PySSA Objects Panel.</li>
                    <li>Click Color to open the colour palette popup.</li>
                    <li>Choose a colour or colouring scheme (e.g. by element, chain, spectrum).</li>
                    <li>The selected structure is recoloured immediately.</li>
                </ol>
            """,

            "viewer_clean": """
                <p>Clean:</p>
                <p>Remove solvent molecules or organic ligands from the viewport display.</p>
                <ol>
                    <li>Click Clean to open the cleanup popup.</li>
                    <li>Select Solvent Molecules to hide water and ions.</li>
                    <li>Select Organic Molecules to hide small-molecule ligands.</li>
                </ol>
                <p style="color: red;"><b>&#9888; CAUTION:</b></p>
                <p>This hides molecules from the viewport only. No data is deleted from the project.</p>
            """,

            "viewer_running_jobs": """
                <p>Active Jobs:</p>
                <p>View and monitor currently running background jobs
                (predictions, analyses, image exports).</p>
                <ol>
                    <li>Click Active Jobs to open the running jobs popup.</li>
                    <li>Each active job is listed with its type and current status.</li>
                    <li>When a job finishes, its result appears automatically in the
                    PySSA Objects Panel and a notification is created.</li>
                </ol>
            """,

            "viewer_notifications": """
                <p>Completed Jobs:</p>
                <p>View the list of all finished background jobs and their results.</p>
                <ol>
                    <li>Click Completed Jobs to open the finished jobs popup.</li>
                    <li>Each completed job is listed with its name and outcome.</li>
                    <li>Results are also reflected automatically in the
                    PySSA Objects Panel.</li>
                </ol>
            """,

          # <editor-fold desc="Popup menus — representation show/hide">
            "popup_cartoon_show": """
                <p>Show Cartoon:</p>
                <p>Make the cartoon representation visible for the currently selected
                protein or chain in the PyMOL Viewer.</p>
            """,
            "popup_cartoon_hide": """
                <p>Hide Cartoon:</p>
                <p>Hide the cartoon representation from the PyMOL Viewer.
                The data is not deleted; click Show to restore it.</p>
            """,

            "popup_sticks_show": """
                <p>Show Sticks:</p>
                <p>Make the sticks representation visible for the currently selected
                protein, residue, or chain.</p>
            """,
            "popup_sticks_hide": """
                <p>Hide Sticks:</p>
                <p>Hide the sticks representation from the PyMOL Viewer.
                The data is not deleted; click Show to restore it.</p>
            """,

            "popup_ribbon_show": """
                <p>Show Ribbon:</p>
                <p>Make the ribbon representation visible for the currently selected
                protein or chain.</p>
            """,
            "popup_ribbon_hide": """
                <p>Hide Ribbon:</p>
                <p>Hide the ribbon representation from the PyMOL Viewer.
                The data is not deleted; click Show to restore it.</p>
            """,

            "popup_lines_show": """
                <p>Show Lines:</p>
                <p>Make the lines representation visible for the currently selected
                protein or chain.</p>
            """,
            "popup_lines_hide": """
                <p>Hide Lines:</p>
                <p>Hide the lines representation from the PyMOL Viewer.
                The data is not deleted; click Show to restore it.</p>
            """,

            "popup_spheres_show": """
                <p>Show Spheres:</p>
                <p>Make the van der Waals spheres representation visible for the
                currently selected protein or residue.</p>
            """,
            "popup_spheres_hide": """
                <p>Hide Spheres:</p>
                <p>Hide the spheres representation from the PyMOL Viewer.
                The data is not deleted; click Show to restore it.</p>
            """,

            "popup_dots_show": """
                <p>Show Dots:</p>
                <p>Make the dots representation visible for the currently selected
                protein or residue.</p>
            """,
            "popup_dots_hide": """
                <p>Hide Dots:</p>
                <p>Hide the dots representation from the PyMOL Viewer.
                The data is not deleted; click Show to restore it.</p>
            """,

            "popup_mesh_show": """
                <p>Show Mesh:</p>
                <p>Make the wireframe mesh representation visible for the currently
                selected protein or chain.</p>
            """,
            "popup_mesh_hide": """
                <p>Hide Mesh:</p>
                <p>Hide the mesh representation from the PyMOL Viewer.
                The data is not deleted; click Show to restore it.</p>
            """,

            "popup_surface_show": """
                <p>Show Surface:</p>
                <p>Make the solid molecular surface visible for the currently selected
                protein or chain.</p>
            """,
            "popup_surface_hide": """
                <p>Hide Surface:</p>
                <p>Hide the molecular surface from the PyMOL Viewer.
                The data is not deleted; click Show to restore it.</p>
            """,
          # </editor-fold>

          # <editor-fold desc="Popup menus — Clean">
            "popup_clean_solvent": """
                <p>Clean: Solvent Molecules:</p>
                <p>Hide all solvent molecules (water, ions) from the PyMOL Viewer.</p>
                <p style="color: red;"><b>&#9888; CAUTION:</b></p>
                <p>This hides solvent molecules from the viewport only.
                No data is deleted from the project.</p>
            """,

            "popup_clean_organic": """
                <p>Clean: Organic Molecules:</p>
                <p>Hide all organic small-molecule ligands from the PyMOL Viewer.</p>
                <p style="color: red;"><b>&#9888; CAUTION:</b></p>
                <p>This hides organic molecules from the viewport only.
                No data is deleted from the project.</p>
            """,
          # </editor-fold>

          # <editor-fold desc="Popup panels — Color">
            "popup_color_bg_white": """
                <p>Background Color: White:</p>
                <p>Set the PyMOL Viewer background to white.
                Useful for publication-ready images.</p>
            """,

            "popup_color_bg_grey": """
                <p>Background Color: Grey:</p>
                <p>Set the PyMOL Viewer background to a medium grey (grey40).</p>
            """,

            "popup_color_bg_black": """
                <p>Background Color: Black:</p>
                <p>Set the PyMOL Viewer background to black.
                The default background; provides high contrast for protein colours.</p>
            """,

            "popup_color_by_elements": """
                <p>Color Atoms by Element:</p>
                <p>Apply standard CPK element colouring to all atoms of the
                selected protein or structure.
                Carbon = grey, Nitrogen = blue, Oxygen = red, Sulphur = yellow.</p>
            """,

          # <editor-fold desc="Color grid buttons">
            # --- Reds ---
            "color_red":        "<p>Color: <b>red</b> (#ff0000) &mdash; Click to apply pure red to the selected protein or structure.</p>",
            "color_tv_red":     "<p>Color: <b>tv_red</b> (#ff3333) &mdash; Click to apply TV red to the selected protein or structure.</p>",
            "color_salmon":     "<p>Color: <b>salmon</b> (#ff9999) &mdash; Click to apply salmon to the selected protein or structure.</p>",
            "color_raspberry":  "<p>Color: <b>raspberry</b> (#b24c66) &mdash; Click to apply raspberry to the selected protein or structure.</p>",

            # --- Greens ---
            "color_green":      "<p>Color: <b>green</b> (#00ff00) &mdash; Click to apply pure green to the selected protein or structure.</p>",
            "color_tv_green":   "<p>Color: <b>tv_green</b> (#33ff33) &mdash; Click to apply TV green to the selected protein or structure.</p>",
            "color_palegreen":  "<p>Color: <b>palegreen</b> (#a5e5a5) &mdash; Click to apply pale green to the selected protein or structure.</p>",
            "color_forest":     "<p>Color: <b>forest</b> (#339933) &mdash; Click to apply forest green to the selected protein or structure.</p>",

            # --- Blues ---
            "color_blue":       "<p>Color: <b>blue</b> (#0000ff) &mdash; Click to apply pure blue to the selected protein or structure.</p>",
            "color_tv_blue":    "<p>Color: <b>tv_blue</b> (#4c4cff) &mdash; Click to apply TV blue to the selected protein or structure.</p>",
            "color_lightblue":  "<p>Color: <b>lightblue</b> (#bfbfff) &mdash; Click to apply light blue to the selected protein or structure.</p>",
            "color_skyblue":    "<p>Color: <b>skyblue</b> (#337fcc) &mdash; Click to apply sky blue to the selected protein or structure.</p>",

            # --- Yellows ---
            "color_yellow":     "<p>Color: <b>yellow</b> (#ffff00) &mdash; Click to apply pure yellow to the selected protein or structure.</p>",
            "color_tv_yellow":  "<p>Color: <b>tv_yellow</b> (#ffff33) &mdash; Click to apply TV yellow to the selected protein or structure.</p>",
            "color_paleyellow": "<p>Color: <b>paleyellow</b> (#ffff7f) &mdash; Click to apply pale yellow to the selected protein or structure.</p>",
            "color_sand":       "<p>Color: <b>sand</b> (#b78c4c) &mdash; Click to apply sand to the selected protein or structure.</p>",

            # --- Magentas ---
            "color_magenta":    "<p>Color: <b>magenta</b> (#ff00ff) &mdash; Click to apply magenta to the selected protein or structure.</p>",
            "color_purple":     "<p>Color: <b>purple</b> (#bf00bf) &mdash; Click to apply purple to the selected protein or structure.</p>",
            "color_pink":       "<p>Color: <b>pink</b> (#ffa5d8) &mdash; Click to apply pink to the selected protein or structure.</p>",
            "color_hotpink":    "<p>Color: <b>hotpink</b> (#ff007f) &mdash; Click to apply hot pink to the selected protein or structure.</p>",

            # --- Cyans ---
            "color_cyan":       "<p>Color: <b>cyan</b> (#00ffff) &mdash; Click to apply cyan to the selected protein or structure.</p>",
            "color_aquamarine": "<p>Color: <b>aquamarine</b> (#7fffff) &mdash; Click to apply aquamarine to the selected protein or structure.</p>",
            "color_palecyan":   "<p>Color: <b>palecyan</b> (#ccffff) &mdash; Click to apply pale cyan to the selected protein or structure.</p>",
            "color_teal":       "<p>Color: <b>teal</b> (#00bfbf) &mdash; Click to apply teal to the selected protein or structure.</p>",

            # --- Oranges ---
            "color_orange":      "<p>Color: <b>orange</b> (#ff7f00) &mdash; Click to apply orange to the selected protein or structure.</p>",
            "color_tv_orange":   "<p>Color: <b>tv_orange</b> (#ff8c26) &mdash; Click to apply TV orange to the selected protein or structure.</p>",
            "color_lightorange": "<p>Color: <b>lightorange</b> (#ffcc7f) &mdash; Click to apply light orange to the selected protein or structure.</p>",
            "color_olive":       "<p>Color: <b>olive</b> (#c4b200) &mdash; Click to apply olive to the selected protein or structure.</p>",

            # --- Greys / Black / White ---
            "color_white":  "<p>Color: <b>white</b> (#ffffff) &mdash; Click to apply white to the selected protein or structure.</p>",
            "color_grey70": "<p>Color: <b>grey70</b> (#b2b2b2) &mdash; Click to apply grey70 to the selected protein or structure.</p>",
            "color_grey30": "<p>Color: <b>grey30</b> (#4c4c4c) &mdash; Click to apply grey30 to the selected protein or structure.</p>",
            "color_black":  "<p>Color: <b>black</b> (#000000) &mdash; Click to apply black to the selected protein or structure.</p>",
          # </editor-fold>

          # </editor-fold>

          # <editor-fold desc="Popup panels — Jobs">
            "popup_active_jobs": """
                <p>Active Jobs panel:</p>
                <p>This panel lists all currently running background jobs,
                including protein structure predictions and distance analyses.</p>
                <p>Each entry shows the job type and its current status.
                Results appear automatically in the PySSA Objects Panel when
                a job completes.</p>
            """,

            "popup_complete_jobs": """
                <p>Completed Jobs panel:</p>
                <p>This panel lists all background jobs that have finished.</p>
                <p>Each entry shows the job type and its final outcome.
                Successful results are already reflected in the
                PySSA Objects Panel.</p>
            """,
          # </editor-fold>

          # </editor-fold>
          # </editor-fold>

          # <editor-fold desc="PySSA Objects Panel toolbar actions">
            "expand_all": """
                <p>Expand all:</p>
                <p>Click to expand the entire tree in the PySSA Objects Panel,
                showing all sequences, proteins (with scenes, chains, residues, atoms)
                and protein pairs (with scenes, proteins, chains, residues, atoms) at once.</p>
            """,

            "collapse_all": """
                <p>Collapse all:</p>
                <p>Click to collapse the entire tree in the PySSA Objects Panel,
                hiding all child items and showing only the top-level entries.</p>
            """,

            "import_file": """
                <p>Import a file into the current project:</p>
                <ol>
                    <li>Click on the Import File button to open the import menu.</li>
                    <li>Select Sequence to import a FASTA sequence file or
                    select Protein to import a PDB protein structure file.</li>
                    <li>Choose the file from your computer in the file dialog.</li>
                    <li>Click Open. The imported sequence or protein appears in
                    the PySSA Objects Panel.</li>
                </ol>
            """,

            "add_sequence": """
                <p>Add a new sequence manually to the current project:</p>
                <ol>
                    <li>Click on the Add Sequence button.</li>
                    <li>Enter the sequence name in the dialog.</li>
                    <li>Enter the amino acid sequence (single-letter code) in the text field.</li>
                    <li>Click on Add to save the sequence to the project.</li>
                </ol>
                <p style="color: red;"><b>&#9888; CAUTION:</b></p>
                <p>The only characters for the sequence name that can be used are
                0-9, a-z, A-Z, -, _. Moreover, multimer sequences must contain
                a comma between the individual chain sequences (e.g. MKABC,MKLMN).
                </p>
            """,

            "export_file": """
                <p>Export selected sequence or protein to a file on your computer:</p>
                <ol>
                    <li>Select one or more sequences or standalone proteins
                    in the PySSA Objects Panel. The selection of >1 sequence or
                    protein can be done via Ctrl+Click, Shift+Click,
                    Ctrl+arrow keys up and down and Shift+arrow keys up and down.</li>
                    <li>Click on the Export File button.</li>
                    <li>Choose an output directory in the file dialog.</li>
                    <li>Click Select Folder.</li>
                    <li>Sequences are exported as .fasta files;
                    proteins are exported as .pdb files.</li>
                </ol>
                <p style="color: red;"><b>&#9888; CAUTION:</b></p>
                <p>Protein pairs cannot be exported with this button.
                Export individual proteins from the pair instead.</p>
            """,

            "delete_object": """
                <p>Delete selected sequence, protein or protein pair from
                the current project:</p>
                <ol>
                    <li>Select one sequence, protein or protein pair in the
                    PySSA Objects Panel.</li>
                    <li>Click on the Delete Object button.</li>
                    <li>Confirm the deletion in the confirmation dialog.</li>
                    <li>The selected objects are permanently removed from the project.</li>
                </ol>
                <p style="color: red;"><b>&#9888; CAUTION:</b></p>
                <p>This action cannot be undone.
                A standalone protein cannot be deleted while it is still part
                of a protein pair. Delete the protein pair first, then delete
                the protein.</p>
            """,
          # </editor-fold>

          # <editor-fold desc="PySSA Objects Panel tree node types">
            "sequence_item": """
                <p>Sequence:</p>
                <p>A sequence item stores an amino acid sequence (single-letter code)
                in the project.</p>
                <p>What you can do with a sequence:</p>
                <ul>
                    <li>Select it and click Export File to save it as a .fasta file.</li>
                    <li>Select it and run Prediction to generate a 3D protein structure.</li>
                    <li>Select it and click Delete Object to remove it from the project.</li>
                </ul>
            """,

            "protein_item": """
                <p>Protein:</p>
                <p>A protein item stores a 3D protein structure (PDB data) in the project.</p>
                <p>What you can do with a protein:</p>
                <ul>
                    <li>Select it and click Export File to save it as a .pdb file.</li>
                    <li>Select it and run Analysis to compare it with another protein.</li>
                    <li>Select it and load it into the PyMOL Viewer.
                    After that, modify the protein and use Image to render a viewport image.</li>
                    <li>Select it and click Delete Object to remove it from the project.</li>
                </ul>
                <p style="color: red;"><b>&#9888; CAUTION:</b></p>
                <p>A protein that is part of a protein pair cannot be deleted independently.
                Delete the protein pair first.</p>
            """,

            "protein_pair_item": """
                <p>Protein pair:</p>
                <p>A protein pair links two proteins that have been compared via distance analysis.</p>
                <p>What you can do with a protein pair:</p>
                <ul>
                    <li>Select it and click Results to view the analysis summary.</li>
                    <li>Select it and load it into the PyMOL Viewer.
                    After that, modify the protein pair and use Image to render
                    a viewport image.</li>
                    <li>Select it and click Delete Object to remove the pair
                    (the underlying proteins are kept).</li>
                </ul>
            """,

            "scene_item": """
                <p>Scene:</p>
                <p>A scene stores a saved PyMOL viewport state (camera position,
                representation, colors) for a specific protein or protein pair.</p>
                <p>What you can do with a scene:</p>
                <ul>
                    <li>Select it to load the saved viewport state in PyMOL.</li>
                    <li>Use the viewer toolbar to resave the current viewport as a new scene.</li>
                    <li>Modify the scene and add a new scene.</li>
                    <li>Select a scene and delete it.</li>
                </ul>
            """,

          # <editor-fold desc="PySSA Objects Panel tree — section nodes">
            "section_sequences": """
                <p>Sequences section:</p>
                <p>This section groups all amino acid sequences stored in the current project.</p>
                <p>From here you can:</p>
                <ol>
                    <li>Expand the section to see all individual sequences.</li>
                    <li>Select one or more sequences and run Prediction to generate 3D structures.</li>
                    <li>Select a sequence and click Export File to save it as a .fasta file.</li>
                    <li>Select a sequence and click Delete Object to remove it from the project.</li>
                    <li>Click Add Sequence in the toolbar to create a new sequence manually.</li>
                    <li>Click Import File in the toolbar to import a .fasta file.</li>
                </ol>
            """,

            "section_proteins": """
                <p>Proteins section:</p>
                <p>This section groups all standalone protein structures stored in the current project.</p>
                <p>From here you can:</p>
                <ol>
                    <li>Expand the section to see all individual proteins.</li>
                    <li>Select a protein and run Analysis to compare it with another protein.</li>
                    <li>Select a protein and use Image to render a viewport image.</li>
                    <li>Select a protein and click Export File to save it as a .pdb file.</li>
                    <li>Select a protein and click Delete Object to remove it from the project.</li>
                    <li>Click Import File in the toolbar to import a .pdb file.</li>
                </ol>
                <p style="color: red;"><b>&#9888; CAUTION:</b></p>
                <p>A protein that is still part of a protein pair cannot be deleted.
                Delete the protein pair first.</p>
            """,

            "section_protein_pairs": """
                <p>Protein Pairs section:</p>
                <p>This section groups all protein pairs that have been analysed
                via distance analysis.</p>
                <p>From here you can:</p>
                <ol>
                    <li>Expand the section to see all protein pairs.</li>
                    <li>Select a protein pair and click Results to view the analysis summary.</li>
                    <li>Expand a protein pair to inspect its two child proteins and their scenes.</li>
                    <li>Select a protein pair and click Delete Object to remove it
                    (the underlying proteins are kept).</li>
                </ol>
            """,
          # </editor-fold>

          # <editor-fold desc="PySSA Objects Panel tree — individual items">
            "protein_item_in_pair": """
                <p>Protein (inside a protein pair):</p>
                <p>This protein is one of the two members of the parent protein pair.</p>
                <p>From here you can:</p>
                <ol>
                    <li>Expand it to see its Chains (and their residues and atoms).</li>
                    <li>Select a chain, residue, or atom under this protein to zoom
                    to that region in the PyMOL viewport.</li>
                    <li>Select it and use Image to render a viewport image of this protein alone.</li>
                </ol>
                <p style="color: red;"><b>&#9888; CAUTION:</b></p>
                <p>This protein cannot be deleted on its own.
                Delete the protein pair first to free both proteins.</p>
            """,

            "header_scenes": """
                <p>Scenes sub-section:</p>
                <p>Lists all saved PyMOL viewport states for the parent protein or protein pair.</p>
                <p>Each scene remembers camera angle, representation style, and colouring.</p>
                <p>From here you can:</p>
                <ol>
                    <li>Expand it to see all scene entries.</li>
                    <li>Select a scene to restore its viewport state in PyMOL.</li>
                </ol>
            """,

            "header_chains": """
                <p>Chains sub-section:</p>
                <p>Lists all polypeptide chains contained in the parent protein structure.</p>
                <p>From here you can:</p>
                <ol>
                    <li>Expand it to see all chain identifiers.</li>
                    <li>Select a chain to highlight it in the PyMOL viewport.</li>
                    <li>Expand a chain to browse its individual residues and atoms.</li>
                </ol>
            """,

            "chain_item": """
                <p>Chain:</p>
                <p>A chain is a single polypeptide chain within a protein structure.
                Each chain is identified by a one-letter ID (e.g. A, B, C).</p>
                <p>From here you can:</p>
                <ol>
                    <li>Select it to highlight the entire chain in the PyMOL viewport.</li>
                    <li>Expand it to see all residues belonging to this chain.</li>
                </ol>
            """,

            "residue_item": """
                <p>Residue:</p>
                <p>A residue is a single amino acid within a chain.
                It is shown as its sequence number followed by the three-letter residue name
                (e.g. 42 - ALA).</p>
                <p>From here you can:</p>
                <ol>
                    <li>Select it to highlight and zoom to this residue in the PyMOL viewport.</li>
                    <li>Expand it to see the individual atoms of this residue.</li>
                </ol>
            """,

            "atom_item": """
                <p>Atom:</p>
                <p>An atom is a single atom within a residue (e.g. CA = alpha carbon, N, O, CB).</p>
                <p>From here you can:</p>
                <ol>
                    <li>Select it to highlight and zoom to this exact atom in the PyMOL viewport.</li>
                </ol>
            """,
          # </editor-fold>

          # <editor-fold desc="Dialogs">
          # <editor-fold desc="Project">
          # Create Project Dialog
            "Dialog": """
                <p>Create a new project:</p>
                <ol>
                    <li>Enter a new project name.</li>
                    <li>Click on Create.</li>
                </ol>
                <p style="color: red;"><b>&#9888; CAUTION:</b></p>
                <p>The project name must have 20 characters or fewer.
                Moreover, the only characters that can be used are 0-9, a-z, A-Z, -, _.</p>
            """,

            # Open Project Dialog
            "OpenProjectDialog": """
                <p>Open a project from your workspace:</p>
                <ol>
                    <li>Enter a new project name.</li>
                    <li><i>Optional:</i> Check the text box under Selected Project to see if the selected name is displayed.</li>
                    <li>Click on Open.</li>
                </ol>
                <p>Search for a project name in your workspace:</p>
                <ol>
                    <li>Type part of the project name in the first textbox.</li>
                    <li>Verify if the text box under 'Selected Project' displays the correct project.</li>
                    <li>Click on Open.</li>
                </ol>
            """,

            # Use Project Dialog
            "UseProjectDialog": """
                        <p>Use a project from your workspace:</p>
                        <ol>
                            <li>Select a project from the list.</li>
                            <li>Choose an existing project.</li>
                            <li><i>Optional:</i> Choose another existing project
                            in combobox or going back to select another existing project.</li>
                            <li>Select the protein under Available Proteins that you want to have in your new project.</li>
                            <li>Click on Add.</li>
                            <li><i>Optional:</i> Choose another protein(s) under
                            Available Proteins that you want to have in your new project. Click on Add.
                            Repeat this until you have all proteins.</li>
                            <li><i>Optional:</i> If you added a wrong protein to Proteins In New Project
                            you can select this protein and click on Remove.
                            Repeat this until you have all the proteins you really want.</li>
                            <li>Click on Create.</li>
                        </ol>
                        <p style="color: red;"><b>&#9888; CAUTION:</b></p>
                        <p>The project name must have 20 characters or fewer.
                        Moreover, the only characters that can be used are 0-9, a-z, A-Z, -, _.</p>
                    """,

            # Delete Project Dialog
            "DeleteProjectDialog": """
                <p>Delete a project from your workspace:</p>
                <ol>
                    <li>Click on one of the projects from the list.</li>
                    <li><i>Optional:</i> Check the text box under 'Selected Project' to confirm that the selected name is displayed.</li>
                    <li>Click on Delete.</li>
                </ol>
            """,
          # </editor-fold>

          # <editor-fold desc="Prediction">
          # Predict Protein Dialog
            "PredictProteinDialog": """
                <p>Protein Structure Prediction dialog:</p>
                <ol>
                    <li>Review the pre-selected sequences in the sequence list on the left.</li>
                    <li><i>Optional:</i> Uncheck sequences you do not want to include
                    in this prediction run.</li>
                    <li>Click on Start Prediction to submit the job.</li>
                    <li>Wait until the prediction finishes.
                    You can monitor progress in the toolbar (running jobs indicator).</li>
                    <li>Once finished, the predicted structure appears in the
                    PySSA Objects Panel automatically.</li>
                </ol>
            """,
          # </editor-fold>

          # <editor-fold desc="Analysis">
          # Distance Analysis Dialog
            "DistanceAnalysisDialog": """
                <p>Distance Analysis dialog:</p>
                <ol>
                    <li>Select the protein pair you want to analyse from the drop-down list.</li>
                    <li>Set the number of Cycles for the structural alignment
                    (higher = more accurate, slower).</li>
                    <li>Set the Cutoff value (&Aring;) to define which C-alpha distances
                    are considered significant.</li>
                    <li>Click on Start Analysis to submit the job.</li>
                    <li>Wait until the progress indicator disappears.</li>
                    <li>Open Results &rarr; Summary to inspect the results.</li>
                </ol>
            """,
          # </editor-fold>

          # <editor-fold desc="Results">
          # Results Summary Dialog
            "ResultsSummaryDialog": """
                <p>Results Summary dialog:</p>
                <ol>
                    <li>Review the RMSD value and the number of aligned residues at the top.</li>
                    <li><i>Optional:</i> Click View Plots to open the distance histogram
                    in a separate window.</li>
                    <li><i>Optional:</i> Click Export Data to save all distance values
                    to a CSV file on your computer.</li>
                </ol>
            """,
          # </editor-fold>

          # <editor-fold desc="Hotspots">
          # Hotspots Protein Regions Dialog
            "HotspotsProteinRegionsDialog": """
                <p>Protein Regions dialog:</p>
                <ol>
                    <li>Select a protein or residues in the PyMOL viewport.</li>
                    <li>Use the controls in the dialog to define the region of interest.</li>
                    <li>Click Apply to highlight the region as sticks with atomic colours.</li>
                    <li>Click Cancel to close the dialog without applying changes.</li>
                </ol>
            """,
          # </editor-fold>

          # <editor-fold desc="Settings">
          # Settings Dialog
            "SettingsDialog": """
                <p>Settings dialog:</p>
                <ol>
                    <li>Enter the ColabFold server address to connect to the
                    prediction server (e.g. http://colabfold-server:8080).</li>
                    <li>Set the desired Cycles and Cutoff values that will be
                    used as the default for new distance analyses.</li>
                    <li>Choose the preferred image renderer, ray trace mode,
                    and ray texture for image exports.</li>
                    <li>Click OK to save all changes.</li>
                </ol>
            """,
          # </editor-fold>


          # </editor-fold>

          # <editor-fold desc="Menu items">
          # <editor-fold desc="Project">
          "action_new_project": """
                <p>Create a new project in your workspace:</p>
                <p>Click to open the Create Project dialog.</p>
            """,

            "action_open_project": """
                <p>Open an existing project from your workspace:</p>
                <p>Click to open the Open Project dialog.</p>
            """,

            "action_use_project": """
                        <p>Use an existing project from your workspace:</p>
                        <p>Click to open the Use Project dialog.</p>
            """,

            "action_delete_project": """
                <p>Delete a project from your workspace:</p>
                <p>Click to open the Delete Project dialog.</p>
            """,

            "action_import_project": """
                <p>Import a project into your workspace:</p>
                <ol>
                    <li>Select Import from the Project menu.</li>
                    <li>Select a Project Database File (.db) from your
                    computer in the file dialog.</li>
                    <li>Click on Open.</li>
                    <li>Either accept the project name or enter
                    a new project name.</li>
                    <li>Wait until the project is imported.</li>
                </ol>
            """,

            "action_export_project": """
                <p>Export the currently active project:</p>
                <ol>
                    <li>Under the Project menu, click on Export.</li>
                    <li>Choose a location to save your Project Database File on
                    your computer in the file dialog.</li>
                    <li>Click on Save.</li>
                    <li>Wait for the export process to finish.</li>
                </ol>
                <p>Share your work with others: </p>
                <ol>
                    <li>Export your project.</li>
                    <li>Send the Project Database File to others.</li>
                </ol>
            """,

            "action_close_project": """
                <p>Close the currently active project:</p>
                <p>Click to close the active project.</p>
                <p>The project is automatically saved.</p>
            """
          # </editor-fold>

          # <editor-fold desc="Prediction">
            ,
            "action_predict_monomer": """
                <p>Run a monomer protein structure prediction:</p>
                <ol>
                    <li>Select one or more monomer sequences in the PySSA Objects Panel.</li>
                    <li>Click on Prediction &rarr; Monomer to open the prediction dialog.</li>
                    <li>Verify the pre-filled sequence list in the dialog.</li>
                    <li>Click on Start Prediction to begin the calculation.</li>
                    <li>Wait for the job to finish. Progress is shown in the toolbar.</li>
                </ol>
                <p style="color: red;"><b>&#9888; CAUTION:</b></p>
                <p>An internet connection and a valid ColabFold server address
                (configured in Settings) are required for prediction.</p>
            """,

            "action_predict_multimer": """
                <p>Run a multimer protein structure prediction:</p>
                <ol>
                    <li>Select one or more multimer sequences (comma-separated chains)
                    in the PySSA Objects Panel.</li>
                    <li>Click on Prediction &rarr; Multimer to open the prediction dialog.</li>
                    <li>Verify the pre-filled sequence list in the dialog.</li>
                    <li>Click on Start Prediction to begin the calculation.</li>
                    <li>Wait for the job to finish. Progress is shown in the toolbar.</li>
                </ol>
                <p style="color: red;"><b>&#9888; CAUTION:</b></p>
                <p>Multimer sequences must contain a comma between chain sequences
                (e.g. <i>MKABC,MKLMN</i>). An internet connection is required.</p>
            """,
          # </editor-fold>

          # <editor-fold desc="Analysis">
            "action_distance_analysis": """
                <p>Run a distance analysis between two protein structures:</p>
                <ol>
                    <li>Select at least one protein pair in the PySSA Objects Panel,
                    or ensure the project contains protein pairs.</li>
                    <li>Click on Analysis &rarr; Distance to open the dialog.</li>
                    <li>Configure the protein pair, cycles, and cutoff value.</li>
                    <li>Click on Start Analysis.</li>
                    <li>Wait for the job to finish. Results appear automatically.</li>
                </ol>
            """,
          # </editor-fold>

          # <editor-fold desc="Results">
            "action_results_summary": """
                <p>View a summary of a completed distance analysis:</p>
                <ol>
                    <li>Select a protein pair (or a protein inside a pair) in the
                    PySSA Objects Panel.</li>
                    <li>Click on Results &rarr; Summary to open the dialog.</li>
                    <li>Inspect the RMSD value and the number of aligned residues.</li>
                    <li><i>Optional:</i> Click on View Plots to see distance histogram plots.</li>
                    <li><i>Optional:</i> Click on Export Data to save results as a CSV file.</li>
                </ol>
            """,
          # </editor-fold>

          # <editor-fold desc="Image">
            "action_preview_image": """
                <p>Preview a ray-traced image in the PyMOL viewport:</p>
                <ol>
                    <li>Load a protein structure session (open a session or select an object).</li>
                    <li>Arrange the view in the PyMOL viewport as desired.</li>
                    <li>Click on Image &rarr; Preview.</li>
                    <li>A preview render (800 &times; 600 px) appears directly in the viewport.
                    No file is saved.</li>
                </ol>
            """,

            "action_ray_tracing_image": """
                <p>Save a high-quality ray-traced image to disk:</p>
                <ol>
                    <li>Load a protein structure session and arrange the viewport view.</li>
                    <li>Click on Image &rarr; Ray-Tracing.</li>
                    <li>Choose a save location and filename in the file dialog.</li>
                    <li>Click Save.</li>
                    <li>Wait for the rendering job to finish.
                    Progress is shown in the toolbar.</li>
                </ol>
                <p style="color: red;"><b>&#9888; CAUTION:</b></p>
                <p>Ray-tracing can take several minutes depending on scene complexity
                and your hardware.</p>
            """,

            "action_simple_image": """
                <p>Save a fast, non-ray-traced image to disk:</p>
                <ol>
                    <li>Load a protein structure session and arrange the viewport view.</li>
                    <li>Click on Image &rarr; Simple.</li>
                    <li>Choose a save location and filename in the file dialog.</li>
                    <li>Click Save.</li>
                    <li>The image is saved immediately without ray-tracing.</li>
                </ol>
            """,
          # </editor-fold>

          # <editor-fold desc="Hotspots">
            "action_protein_regions": """
                <p>Highlight selected protein regions directly in the PyMOL viewport:</p>
                <ol>
                    <li>Select atoms or residues of interest in the PyMOL viewport
                    (they will be stored as the <i>sele</i> selection).</li>
                    <li>Click on Hotspots &rarr; Protein Regions.</li>
                    <li>The selected region is displayed as sticks with atomic colours
                    (non-carbon atoms coloured by element, carbon atoms in grey).</li>
                    <li>The viewport automatically zooms to the selection.</li>
                </ol>
            """,
          # </editor-fold>

          # <editor-fold desc="Settings">
            "action_edit_settings": """
                <p>Edit the application settings:</p>
                <ol>
                    <li>Click on Settings &rarr; Edit to open the Settings dialog.</li>
                    <li>Adjust the ColabFold server address, image options, or analysis
                    parameters as needed.</li>
                    <li>Click OK to save your changes.</li>
                </ol>
            """,

            "action_restore_settings": """
                <p>Reset all settings to their factory defaults:</p>
                <ol>
                    <li>Click on Settings &rarr; Restore.</li>
                    <li>Confirm the reset in the confirmation dialog that appears.</li>
                    <li>All application settings are restored to their default values.</li>
                </ol>
                <p style="color: red;"><b>&#9888; CAUTION:</b></p>
                <p>This action cannot be undone. Any customised settings
                (server address, image quality, cutoff values) will be lost.</p>
            """,
          # </editor-fold>

          # <editor-fold desc="Help menu">
            "action_documentation": """
                <p>Toggle the Help panel:</p>
                <ol>
                    <li>Click on Help &rarr; Documentation.</li>
                    <li>The Help panel on the right side opens or closes.</li>
                    <li>Hover over any menu item or dialog to see context-sensitive help
                    text appear in this panel.</li>
                </ol>
            """,

            "action_get_demo_projects": """
                <p>Download and install demo projects into your workspace:</p>
                <ol>
                    <li>Click on Help &rarr; Get Demo Projects.</li>
                    <li>PySSA downloads a set of pre-built demo projects from the internet.</li>
                    <li>The demo projects are automatically imported into your workspace.</li>
                    <li>Open them via Project &rarr; Open to explore example results.</li>
                </ol>
                <p style="color: red;"><b>&#9888; CAUTION:</b></p>
                <p>An active internet connection is required to download the demo projects.</p>
            """,

            "action_show_log_in_explorer": """
                <p>Open a log file for inspection:</p>
                <ol>
                    <li>Click on Help &rarr; Show Logs in Explorer.</li>
                    <li>A file dialog opens showing all available log files.</li>
                    <li>Select the log file you want to read.</li>
                    <li>Click Open to view the file in the default text application.</li>
                </ol>
            """,

            "action_clear_logs": """
                <p>Delete all generated log files:</p>
                <ol>
                    <li>Click on Help &rarr; Clear Logs.</li>
                    <li>All log files stored in the <i>.pyssa/logs</i> folder are deleted.</li>
                </ol>
                <p style="color: red;"><b>&#9888; CAUTION:</b></p>
                <p>This action cannot be undone. Deleted log files cannot be recovered.</p>
            """,

            "action_about": """
                <p>View information about PySSA:</p>
                <p>Click on Help &rarr; About to open the About dialog.</p>
                <p>The dialog shows the current version, authors, and licence information.</p>
            """,
          # </editor-fold>

          # </editor-fold>
          # </editor-fold>
        }
        self.help_filter = help_event_filter.HelpEventFilter(
            self._main_window.help_panel.help_text_browser,
            help_map,
            self._main_window.help_panel
        )
        self._main_window.pyssa_objects_panel.installEventFilter(self.help_filter)

        # Install hover help on Objects Panel toolbar buttons
        _panel = self._main_window.pyssa_objects_panel
        _toolbar = _panel.get_toolbar()
        _btn_map = [
            (_panel.expand_all, "expand_all"),
            (_panel.collapse_all, "collapse_all"),
            (_panel.import_file_action, "import_file"),
            (_panel.add_sequence_action, "add_sequence"),
            (_panel.export_file_action, "export_file"),
            (_panel.delete_object_action, "delete_object"),
        ]
        for action_wrapper, obj_name in _btn_map:
            btn = _toolbar.get_tool_button_for_action(action_wrapper)
            if btn is not None:
                btn.setObjectName(obj_name)
                btn.installEventFilter(self.help_filter)
                logger.info(f"Installed hover help on Objects Panel button: {obj_name}")

        # Connect the tree view hover signals for item-level help text.
        # setMouseTracking ensures Enter events fire when the cursor moves over items
        # even without clicking.
        _tree = _panel.tree_view
        _tree.setMouseTracking(True)
        _tree.entered.connect(self.help_filter.handle_tree_item_entered)
        _tree.viewport().installEventFilter(self.help_filter)
        logger.info("Connected tree view hover help signals.")

        # Install hover help on the PyMOL viewport widget.
        _main_window = self._main_window
        _main_window.pymolwidget.setObjectName("pymolwidget")
        _main_window.pymolwidget.installEventFilter(self.help_filter)

        # Install hover help on all viewer toolbar buttons.
        # NOTE: The actual visible toolbar is inside tool_window_layout, not
        # main_window.viewer_toolbar (which is never added to any layout).
        _viewer_toolbar = _main_window.tool_window_layout.viewer_toolbar
        _viewer_btn_map = [
            ("open_session",    "viewer_open_session"),
            ("create_scene",    "viewer_create_scene"),
            ("save_scene",      "viewer_save_scene"),
            ("delete_scene",    "viewer_delete_scene"),
            ("cartoon",         "viewer_cartoon"),
            ("sticks",          "viewer_sticks"),
            ("ribbon",          "viewer_ribbon"),
            ("lines",           "viewer_lines"),
            ("spheres",         "viewer_spheres"),
            ("dots",            "viewer_dots"),
            ("mesh",            "viewer_mesh"),
            ("surface",         "viewer_surface"),
            ("color",           "viewer_color"),
            ("clean",           "viewer_clean"),
            ("running_jobs",    "viewer_running_jobs"),
            ("notifications",   "viewer_notifications"),
        ]
        for action_key, obj_name in _viewer_btn_map:
            action_wrapper = _main_window.viewer_toolbar_actions.get(action_key)
            if action_wrapper is not None:
                btn = _viewer_toolbar.get_tool_button_for_action(action_wrapper)
                if btn is not None:
                    btn.setObjectName(obj_name)
                    btn.installEventFilter(self.help_filter)
                    logger.info(f"Installed hover help on viewer toolbar button: {obj_name}")

        # Connect hovered signals for all 8 representation show/hide popup menus.
        _repr_menus = [
            _main_window.cartoon_show_hide_menu,
            _main_window.sticks_show_hide_menu,
            _main_window.ribbon_show_hide_menu,
            _main_window.lines_show_hide_menu,
            _main_window.spheres_show_hide_menu,
            _main_window.dots_show_hide_menu,
            _main_window.mesh_show_hide_menu,
            _main_window.surface_show_hide_menu,
        ]
        for _menu in _repr_menus:
            _menu.hovered.connect(self.help_filter.handle_menu_action_hovered)
            _menu.aboutToHide.connect(self.help_filter.handle_menu_about_to_hide)
        logger.info("Connected hover help for all representation show/hide menus.")

        # Connect the Clean popup menu.
        _main_window.clean_solvent_organic_menu.hovered.connect(self.help_filter.handle_menu_action_hovered)
        _main_window.clean_solvent_organic_menu.aboutToHide.connect(self.help_filter.handle_menu_about_to_hide)

        # Install event filters on the Color popup panel's interactive widgets.
        _color_config = _main_window.color_config
        for _btn in (
            _color_config.btn_white_bg,
            _color_config.btn_grey_bg,
            _color_config.btn_black_bg,
            _color_config.btn_color_by_elements,
        ):
            _btn.installEventFilter(self.help_filter)

        # Install event filters on all 32 color grid buttons.
        for _btn in _main_window.color_grid.get_all_color_buttons().values():
            _btn.installEventFilter(self.help_filter)
        logger.info("Installed hover help on all color config and color grid buttons.")

        # Install event filters on the Active Jobs and Completed Jobs panels.
        _main_window.active_jobs.installEventFilter(self.help_filter)
        _main_window.complete_jobs.installEventFilter(self.help_filter)
        logger.info("Installed hover help on job panels.")

        # Connect menu hovered signal to show help for menu items
        self._main_window.menuProject.hovered.connect(self.help_filter.handle_menu_action_hovered)
        self._main_window.menuProject.aboutToHide.connect(self.help_filter.handle_menu_about_to_hide)
        self._main_window.menuPrediction.hovered.connect(self.help_filter.handle_menu_action_hovered)
        self._main_window.menuPrediction.aboutToHide.connect(self.help_filter.handle_menu_about_to_hide)
        self._main_window.menuAnalysis.hovered.connect(self.help_filter.handle_menu_action_hovered)
        self._main_window.menuAnalysis.aboutToHide.connect(self.help_filter.handle_menu_about_to_hide)
        self._main_window.menuResults.hovered.connect(self.help_filter.handle_menu_action_hovered)
        self._main_window.menuResults.aboutToHide.connect(self.help_filter.handle_menu_about_to_hide)
        self._main_window.menuImage.hovered.connect(self.help_filter.handle_menu_action_hovered)
        self._main_window.menuImage.aboutToHide.connect(self.help_filter.handle_menu_about_to_hide)
        self._main_window.menuHotspots.hovered.connect(self.help_filter.handle_menu_action_hovered)
        self._main_window.menuHotspots.aboutToHide.connect(self.help_filter.handle_menu_about_to_hide)
        self._main_window.menuSettings.hovered.connect(self.help_filter.handle_menu_action_hovered)
        self._main_window.menuSettings.aboutToHide.connect(self.help_filter.handle_menu_about_to_hide)
        self._main_window.menuAbout.hovered.connect(self.help_filter.handle_menu_action_hovered)
        self._main_window.menuAbout.aboutToHide.connect(self.help_filter.handle_menu_about_to_hide)
        logger.info("Connected hover help for all menu items")

        self.refresh_ui()
        self.open_welcome_screen()

    # <editor-fold desc="Private methods">
    def _connect_all_signals_with_their_slots(self) -> None:
        """Connects all relevant widget signals with their appropriate slots."""
        # self._main_window.dialogClosed.connect(self.__slot_close_application)

        # <editor-fold desc="Project menu">
        self._main_window.action_new_project.triggered.connect(self.__slot_create_project)
        self._main_window.action_open_project.triggered.connect(self.__slot_open_project)
        self._main_window.action_use_project.triggered.connect(self.__slot_use_project)
        self._main_window.action_delete_project.triggered.connect(self.__slot_delete_project)
        self._main_window.action_import_project.triggered.connect(self.__slot_import_project)
        self._main_window.action_export_project.triggered.connect(self.__slot_export_current_project)
        self._main_window.action_close_project.triggered.connect(self.__slot_close_project)
        # </editor-fold>
        # TODO: Add the right slot method! ;)
        # self._main_window.action_exit_application.triggered.connect(self.)

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
            # Install hover help event filter on the dialog
            dialog_view = self._dialog_controllers["create_project"].get_view()
            dialog_view.installEventFilter(self.help_filter)
            logger.info(f"Installed hover help on Create Project dialog (objectName: {dialog_view.objectName()})")
        self._dialog_controllers["create_project"].restore_default_view()
        self._dialog_controllers["create_project"].get_view().show()

    def __slot_open_project(self):
        if not self._dialog_controllers.__contains__("open_project"):
            self._dialog_controllers["open_project"] = open_project_view_controller.OpenProjectViewController(
                self._app_state
            )
            # Install hover help event filter on the dialog
            dialog_view = self._dialog_controllers["open_project"].get_view()
            dialog_view.installEventFilter(self.help_filter)
            logger.info(f"Installed hover help on Open Project dialog (objectName: {dialog_view.objectName()})")
        self._dialog_controllers["open_project"].restore_default_view()
        self._dialog_controllers["open_project"].get_view().show()

    def __slot_delete_project(self):
        if not self._dialog_controllers.__contains__("delete_project"):
            from src.pyssa.controller import delete_project_view_controller
            self._dialog_controllers["delete_project"] = delete_project_view_controller.DeleteProjectViewController(
                self._app_state
            )
            # Install hover help event filter on the dialog
            dialog_view = self._dialog_controllers["delete_project"].get_view()
            dialog_view.installEventFilter(self.help_filter)
            logger.info(f"Installed hover help on Delete Project dialog (objectName: {dialog_view.objectName()})")
        self._dialog_controllers["delete_project"].restore_default_view()
        self._dialog_controllers["delete_project"].get_view().show()

    def __slot_use_project(self):
        if not self._dialog_controllers.__contains__("use_project"):
            from src.pyssa.controller import use_project_view_controller
            self._dialog_controllers["use_project"] = use_project_view_controller.UseProjectViewController(
                self._app_state
            )
            # Install hover help event filter on the dialog
            dialog_view = self._dialog_controllers["use_project"].get_view()
            dialog_view.installEventFilter(self.help_filter)
            logger.info(f"Installed hover help on Use Project dialog (objectName: {dialog_view.objectName()})")
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

        # No relevant selection â€” fall back to all qualifying sequences.
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
            # Install hover help event filter on the dialog
            dialog_view = self._dialog_controllers["predict_monomer"].get_view()
            dialog_view.installEventFilter(self.help_filter)
            logger.info(f"Installed hover help on Predict Protein dialog (objectName: {dialog_view.objectName()})")
        else:
            # Re-create when called again so it reflects the current state.
            self._dialog_controllers["predict_monomer"] = predict_protein_view_controller.PredictProteinViewController(
                self._app_state,
                sequences,
                a_parent=self._main_window,
            )
            self._dialog_controllers["predict_monomer"].get_view().installEventFilter(self.help_filter)
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
            # Install hover help event filter on the dialog
            dialog_view = self._dialog_controllers["predict_multimer"].get_view()
            dialog_view.installEventFilter(self.help_filter)
            logger.info(f"Installed hover help on Predict Multimer dialog (objectName: {dialog_view.objectName()})")
        else:
            # Re-create when called again so it reflects the current state.
            self._dialog_controllers["predict_multimer"] = predict_protein_view_controller.PredictProteinViewController(
                self._app_state,
                sequences,
                a_parent=self._main_window,
            )
            self._dialog_controllers["predict_multimer"].get_view().installEventFilter(self.help_filter)
        self._dialog_controllers["predict_multimer"].get_view().show()
    # </editor-fold>

    # <editor-fold desc="Analysis menu">
    def __slot_distance_analysis(self):
        if not self._dialog_controllers.__contains__("distance_analysis_dialog"):
            self._dialog_controllers["distance_analysis_dialog"] = distance_analysis_view_controller.DistanceAnalysisViewController(
                self._app_state
            )
            # Install hover help event filter on the dialog
            dialog_view = self._dialog_controllers["distance_analysis_dialog"].get_view()
            dialog_view.installEventFilter(self.help_filter)
            logger.info(f"Installed hover help on Distance Analysis dialog (objectName: {dialog_view.objectName()})")
        self._dialog_controllers["distance_analysis_dialog"].restore_default_view()
        self._dialog_controllers["distance_analysis_dialog"].get_view().show()
    # </editor-fold>

    # <editor-fold desc="Results menu">
    def __slot_results_summary(self):
        if not self._dialog_controllers.__contains__("results_summary_dialog"):
            self._dialog_controllers["results_summary_dialog"] = results_view_controller.ResultsViewController(
                self._get_selected_protein_pairs()[0], self._app_state, self._user_pymol, self._main_window
            )
            # Install hover help event filter on the dialog
            dialog_view = self._dialog_controllers["results_summary_dialog"].get_view()
            dialog_view.installEventFilter(self.help_filter)
            logger.info(f"Installed hover help on Results Summary dialog (objectName: {dialog_view.objectName()})")
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
            # Install hover help event filter on the dialog
            dialog_view = self._dialog_controllers["settings_dialog"].get_view()
            dialog_view.installEventFilter(self.help_filter)
            logger.info(f"Installed hover help on Settings dialog (objectName: {dialog_view.objectName()})")
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

    def __slot_display_clean_options(self):
        try:
            self._main_window.clean_solvent_organic_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("clean"))
            )
        except Exception as e:
            logger.error(e.__str__())


    def __slot_clean_solvent(self):
        self._user_pymol.get_cmd_module().remove("solvent")

    def __slot_clean_organic(self):
        self._user_pymol.get_cmd_module().remove("organic")

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
        update to prevent a feedback loop (PyMOL â†’ tree â†’ PyMOL).
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

            # Suppress tree â†’ PyMOL sync while we modify the tree selection.
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
            # snapshot â†’ refresh_ui path was skipped.  Resolve a snapshot
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
            for tmp_protein in result["predicted_proteins"]:
                self._app_state.project.add_existing_protein(tmp_protein)
                self._app_state.pyssa_objects_model.add_protein(tmp_protein)
                self._app_state.hot_db.write_queue.submit(
                    WriteOperation(OperationType.INSERT_PROTEIN, tmp_protein)
                )
        else:
            for tmp_protein in result["predicted_proteins"]:
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
        print("Hi")
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
