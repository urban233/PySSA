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
"""Module for the interface manager."""
import glob
import os
import pathlib
import shutil
import subprocess
import sys

from PyQt5 import QtGui, QtCore
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt

from application_process import application_process_manager
from pyssa.gui.ui import icon_resources  # this import is used for the icons! DO NOT DELETE THIS
from pyssa.controller import pymol_session_manager, settings_manager, main_tasks_manager, \
    status_bar_manager, job_manager, watcher
from pyssa.controller.database_manager import logger
from pyssa.gui.ui.custom_widgets import job_entry
from pyssa.gui.ui.dialogs import dialog_startup
from pyssa.gui.ui.views import rename_protein_view, use_project_view, \
    predict_multimer_view, add_sequence_view, add_scene_view, settings_view, predict_protein_view, \
    fasta_file_import_preview_view, rename_sequence_view, add_protein_pair_view, advanced_prediction_configurations
from pyssa.gui.ui.views import create_project_view, open_project_view, delete_project_view, import_sequence_view
from pyssa.gui.ui.views import main_view, predict_monomer_view, distance_analysis_view, results_view, add_protein_view
from pyssa.gui.ui.views import hotspots_protein_regions_view
from pyssa.gui.ui.styles import styles
from pyssa.internal.data_structures import project, settings, chain, protein, protein_pair
from pyssa.internal.data_structures.data_classes import current_session
from pyssa.internal.thread import tasks
from pyssa.internal.thread.async_pyssa import locks, custom_signals
from pyssa.io_pyssa import filesystem_io
from pyssa.model import proteins_model, protein_pairs_model
from pyssa.util import enums, constants, main_window_util
from pyssa.util.void import rvoid


class InterfaceManager:

    """A manager for all views."""
    string_model = QtCore.QStringListModel()

    _main_view: "main_view.MainView"
    _settings_view: "settings_view.SettingsView"
    _predict_monomer_view: "predict_monomer_view.PredictMonomerView"
    _predict_protein_view: "predict_protein_view.PredictProteinView"
    _distance_analysis_view: "distance_analysis_view.DistanceAnalysisView"
    _create_project_view: "create_project_view.CreateProjectView"
    _open_project_view: "open_project_view.OpenProjectView"
    _delete_project_view: "delete_project_view.DeleteProjectView"
    _results_view: "results_view.ResultsView"
    _add_protein_view: "add_protein_view.AddProteinView"
    _import_sequence_view: "import_sequence_view.ImportSequenceView"
    _fasta_file_import_preview_view: "fasta_file_import_preview_view.FastaFileImportPreviewView"
    _add_sequence_view: "add_sequence_view.AddSequenceView"
    _rename_protein_view: "rename_protein_view.RenameProteinView"
    _use_project_view: "use_project_view.UseProjectView"
    _hotspots_protein_regions_view: "hotspots_protein_regions_view.HotspotsProteinRegionsView"
    _add_scene_view: "add_scene_view.AddSceneView"
    _add_protein_pair_view: "add_protein_pair_view.AddProteinPairView"

    _current_workspace: pathlib.Path
    _current_project: "project.Project"
    _current_pymol_session: "current_session.CurrentPymolSession"

    project_lock: QtCore.QMutex
    pymol_lock: "locks.PyMOL_LOCK"

    # _application_settings: "settings.Settings"
    current_tab_index: int = 0

    _workspace_model: QtGui.QStandardItemModel
    _sequence_model: QtGui.QStandardItemModel
    _protein_model: QtGui.QStandardItemModel
    _protein_pair_model: QtGui.QStandardItemModel

    def __init__(self) -> None:
        # View definitions
        self._main_view = main_view.MainView()
        self._settings_view = settings_view.SettingsView()
        self._predict_monomer_view = predict_monomer_view.PredictMonomerView()
        self._predict_multimer_view = predict_multimer_view.PredictMultimerView()
        self._predict_protein_view = predict_protein_view.PredictProteinView()
        self._distance_analysis_view = distance_analysis_view.DistanceAnalysisView()
        self._create_project_view = create_project_view.CreateProjectView()
        self._open_project_view = open_project_view.OpenProjectView()
        self._delete_project_view = delete_project_view.DeleteProjectView()
        self._hotspots_protein_regions_view = hotspots_protein_regions_view.HotspotsProteinRegionsView()
        self._results_view = results_view.ResultsView()
        self._add_protein_view = add_protein_view.AddProteinView()
        self._import_sequence_view: "import_sequence_view.ImportSequenceView" = import_sequence_view.ImportSequenceView()
        self._fasta_file_import_preview_view = fasta_file_import_preview_view.FastaFileImportPreviewView()
        self._add_sequence_view = add_sequence_view.AddSequenceView()
        self._rename_protein_view = rename_protein_view.RenameProteinView()
        self._rename_sequence_view = rename_sequence_view.RenameSequenceView()
        self._use_project_view = use_project_view.UseProjectView()
        self._add_scene_view = add_scene_view.AddSceneView()
        self._add_protein_pair_view = add_protein_pair_view.AddProteinPairView()
        self._advanced_prediction_configurations = advanced_prediction_configurations.AdvancedPredictionConfigurationsView()

        self.main_tasks_manager = main_tasks_manager.MainTasksManager()
        self.app_process_manager = application_process_manager.ApplicationProcessManager(self._reset_pymol_session)
        self.pymol_session_manager = pymol_session_manager.PymolSessionManager(self.app_process_manager)
        self.job_manager = job_manager.JobManager()
        self.status_bar_manager = status_bar_manager.StatusBarManager(self._main_view,
                                                                      self.main_tasks_manager)

        self.watcher = watcher.Watcher()

        self.documentation_window = None
        self.job_entry_widgets = []

        self.refresh_after_job_finished_signal = custom_signals.RefreshAfterJobFinishedSignal()

        self._settings_manager = settings_manager.SettingsManager()

        # <editor-fold desc="Setup App Settings">
        # self._application_settings = settings.Settings(constants.SETTINGS_DIR, constants.SETTINGS_FILENAME)
        if not os.path.exists(constants.SETTINGS_FULL_FILEPATH):
            constants.PYSSA_LOGGER.info("Settings file not found, open configuration dialog.")
            # Configuration dialog to setup setting file
            dialog = dialog_startup.DialogStartup()
            dialog.exec_()

            # checks if the cancel button was pressed
            if dialog_startup.global_var_terminate_app == 1:
                os.remove(constants.SETTINGS_FULL_FILEPATH)
                constants.PYSSA_LOGGER.info("Configuration dialog closed, and removed new settings file.")
                sys.exit()

            self._settings_manager.settings.app_launch = 1
            self._settings_manager.settings.workspace_path = pathlib.Path(dialog_startup.global_var_startup_workspace)
            # self._application_settings.app_launch = 1
            # self._application_settings.workspace_path = pathlib.Path(dialog_startup.global_var_startup_workspace)

            constants.PYSSA_LOGGER.info("Demo projects are getting extracted ...")
            import zipfile

            with zipfile.ZipFile(pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects.zip"), "r") as zip_ref:
                zip_ref.extractall(pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects"))
            constants.PYSSA_LOGGER.info(
                "Demo projects are downloaded and extracted.\n Import of demo projects started ...",
            )

            path_of_demo_projects = pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects")
            for tmp_filename in os.listdir(path_of_demo_projects):
                # Copy db file into new workspace
                tmp_project_database_filepath = str(
                    pathlib.Path(
                        f"{self.get_application_settings().workspace_path}/{tmp_filename}"
                    )
                )
                tmp_src_filepath = str(pathlib.Path(f"{path_of_demo_projects}/{tmp_filename}"))
                shutil.copyfile(tmp_src_filepath, tmp_project_database_filepath)
            constants.PYSSA_LOGGER.info("Import process of demo projects finished.")
            constants.PYSSA_LOGGER.info("Serialize settings ...")
            self._settings_manager.settings.serialize_settings()
            # self._application_settings.serialize_settings()
            constants.PYSSA_LOGGER.info("Serialize settings finished.")

            QtWidgets.QApplication.restoreOverrideCursor()
        self._settings_manager.settings = main_window_util.setup_app_settings(self._settings_manager.settings)

        # </editor-fold>

        # General attributes definitions
        self._current_project = project.Project()
        self._current_pymol_session = current_session.CurrentPymolSession("", "")

        self.project_lock = QtCore.QMutex()
        self.pymol_lock: "locks.PyMOL_LOCK" = locks.PyMOL_LOCK()

        self.job_manager.start_auxiliary_pymol()
        self.start_app_process_manager()
        #self.start_pymol()

        # Model definitions
        self._workspace_model = QtGui.QStandardItemModel()
        self._sequence_model = QtGui.QStandardItemModel()
        self._protein_model: "proteins_model.ProteinsModel" = proteins_model.ProteinsModel()
        self._protein_pair_model: "protein_pairs_model.ProteinPairsModel" = protein_pairs_model.ProteinPairsModel()
        self._build_workspace_model()

    def start_app_process_manager(self):
        self._app_process_manager_thread = tasks.Task(
            target=self.app_process_manager.check_process,
            args=(0, 0),
            post_func=self._closed_app_process_manager
        )
        self._app_process_manager_thread.start()

    def _closed_app_process_manager(self):
        logger.info("Application process manager closed.")

    def start_pymol(self):
        process = subprocess.Popen([f"{constants.PLUGIN_PATH}\\scripts\\batch\\start_pymol.bat"],  # r"C:\Users\martin\github_repos\PySSA\scripts\batch\start_pymol.bat"
                                   creationflags=subprocess.CREATE_NO_WINDOW)
        if process.poll() is None:
            print("PyMOL from interface manager started correctly.")
        else:
            print("PyMOL failed to start.")

    def _reset_pymol_session(self):
        """Resets the pymol session like it was before the User PyMOL crash.

        Notes:
            This function gets executed in the check_process() function of the application process manager!
        """
        if self.pymol_session_manager.session_object_type == "protein":
            self.pymol_session_manager.load_protein_session(self.pymol_session_manager.session_objects[0])
        elif self.pymol_session_manager.session_object_type == "protein_pair":
            self.pymol_session_manager.load_protein_pair_session(self.pymol_session_manager.session_objects[0])
        else:
            return
        self.pymol_session_manager.load_current_scene()

    # <editor-fold desc="Getter Methods">
    # <editor-fold desc="Settings">
    def get_settings_manager(self):
        return self._settings_manager

    def get_application_settings(self) -> "settings.Settings":
        return self._settings_manager.settings
    # </editor-fold>

    # <editor-fold desc="Workspace">
    def get_workspace_path(self):
        return self._settings_manager.settings.workspace_path

    def get_workspace_model(self) -> QtGui.QStandardItemModel:
        """Returns the current workspace model"""
        return self._workspace_model

    def get_workspace_projects(self):
        """Returns the workspace projects."""
        db_pattern = os.path.join(self._settings_manager.settings.get_workspace_path(), '*.db')
        self.string_model.setStringList(
            # Filters the workspace for all project files based on the xml extension
            [os.path.basename(file).replace(".db", "") for file in glob.glob(db_pattern)]
        )
        return self.string_model

    def get_workspace_projects_as_list(self) -> list:
        db_pattern = os.path.join(self._settings_manager.settings.get_workspace_path(), '*.db')
        return [os.path.basename(file).replace(".db", "") for file in glob.glob(db_pattern)]
    # </editor-fold>

    # <editor-fold desc="Getter Methods for view">
    def get_main_view(self) -> "main_view.MainView":
        return self._main_view

    def get_settings_view(self):
        return self._settings_view

    def get_open_view(self) -> "open_project_view.OpenProjectView":
        return self._open_project_view

    def get_create_view(self) -> "create_project_view.CreateProjectView":
        return self._create_project_view

    def get_delete_view(self) -> "delete_project_view.DeleteProjectView":
        return self._delete_project_view

    def get_use_project_view(self):
        return self._use_project_view

    def get_fasta_file_import_preview_view(self):
        return self._fasta_file_import_preview_view

    def get_add_sequence_view(self):
        return self._add_sequence_view

    def get_rename_sequence_view(self):
        return self._rename_sequence_view

    # <editor-fold desc="Prediction">
    def get_predict_monomer_view(self) -> "predict_monomer_view.PredictMonomerView":
        return self._predict_monomer_view

    def get_predict_multimer_view(self) -> "predict_multimer_view.PredictMultimerView":
        return self._predict_multimer_view

    def get_predict_protein_view(self):
        return self._predict_protein_view

    # </editor-fold>

    def get_distance_analysis_view(self) -> "distance_analysis_view.DistanceAnalysisView":
        return self._distance_analysis_view

    def get_results_view(self):
        return self._results_view

    def get_add_scene_view(self):
        return self._add_scene_view

    def get_add_protein_pair_view(self):
        return self._add_protein_pair_view

    def get_advanced_prediction_configurations_view(self):
        return self._advanced_prediction_configurations

    # <editor-fold desc="Getter methods for Sequence Tab in main view">
    def get_current_sequence_list_index(self):
        return self._main_view.ui.seqs_list_view.currentIndex()

    def get_current_sequence_list_index_object(self):
        """Returns the selected seq record object from the list view."""
        return self.get_current_sequence_list_index().data(enums.ModelEnum.OBJECT_ROLE)

    def get_import_sequence_view(self):
        return self._import_sequence_view

    def get_add_sequence_view(self):
        return self._add_sequence_view

    # </editor-fold>

    # <editor-fold desc="Getter methods for Protein Tab in main view">
    def get_current_protein_tree_index(self):
        return self._main_view.ui.proteins_tree_view.currentIndex()

    def get_child_index_of_get_current_protein_tree_index(self):
        return self._main_view.ui.proteins_tree_view.currentIndex().child(0, 0)

    def get_current_protein_tree_index_type(self):
        return self._main_view.ui.proteins_tree_view.model().data(
            self.get_current_protein_tree_index(), enums.ModelEnum.TYPE_ROLE
        )

    def get_current_protein_tree_index_object(self):
        return self.get_current_protein_tree_index().data(enums.ModelEnum.OBJECT_ROLE)

    def get_parent_index_object_of_current_protein_tree_index(self):
        return self.get_current_protein_tree_index().parent().data(enums.ModelEnum.OBJECT_ROLE)

    def get_current_active_protein_object(self) -> "protein.Protein":
        """Returns the protein object of the current active branch.

        Note:
            This function is handles also the header type and does not throw an error!

        Raises:
            ValueError: if the type is unknown

        Returns:
            a protein object
        """
        tmp_type = self._main_view.ui.proteins_tree_view.currentIndex().data(enums.ModelEnum.TYPE_ROLE)
        if tmp_type == "protein":
            return self._main_view.ui.proteins_tree_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE)
        elif tmp_type == "header":
            return self._main_view.ui.proteins_tree_view.currentIndex().parent().data(enums.ModelEnum.OBJECT_ROLE)
        elif tmp_type == "scene":
            return self._main_view.ui.proteins_tree_view.currentIndex().parent().parent().data(
                enums.ModelEnum.OBJECT_ROLE
            )
        elif tmp_type == "chain":
            return self._main_view.ui.proteins_tree_view.currentIndex().parent().parent().data(
                enums.ModelEnum.OBJECT_ROLE
            )
        else:
            raise ValueError("Unknown type!")

    def get_current_header_name(self) -> str:
        """Returns the name of the current header."""
        tmp_type = self._main_view.ui.proteins_tree_view.currentIndex().data(enums.ModelEnum.TYPE_ROLE)
        if tmp_type == "protein":
            raise ValueError(f"Cannot get a header object if the type is: {tmp_type}!")
        elif tmp_type == "header":
            return self._main_view.ui.proteins_tree_view.currentIndex().parent().data(enums.ModelEnum.OBJECT_ROLE)
        elif tmp_type == "scene":
            raise ValueError(f"Cannot get a header object if the type is: {tmp_type}!")
        elif tmp_type == "chain":
            raise ValueError(f"Cannot get a header object if the type is: {tmp_type}!")
        else:
            raise ValueError("Unknown type!")

    def get_current_active_scene_name(self) -> str:
        """Returns the scene name of the current active branch.

        Returns:
            a scene name
        """
        tmp_type = self._main_view.ui.proteins_tree_view.currentIndex().data(enums.ModelEnum.TYPE_ROLE)
        if tmp_type == "protein":
            raise ValueError(f"Cannot get a scene name if the type is: {tmp_type}!")
        elif tmp_type == "header":
            raise ValueError(f"Cannot get a scene name if the type is: {tmp_type}!")
        elif tmp_type == "scene":
            return self._main_view.ui.proteins_tree_view.currentIndex().data(Qt.DisplayRole)
        elif tmp_type == "chain":
            raise ValueError(f"Cannot get a scene name if the type is: {tmp_type}!")
        else:
            raise ValueError("Unknown type!")

    def get_current_active_chain_object(self) -> "chain.Chain":
        """Returns the chain object of the current active branch.

        Returns:
            a chain object
        """
        tmp_type = self._main_view.ui.proteins_tree_view.currentIndex().data(enums.ModelEnum.TYPE_ROLE)
        if tmp_type == "protein":
            raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
        elif tmp_type == "header":
            raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
        elif tmp_type == "scene":
            raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
        elif tmp_type == "chain":
            return self._main_view.ui.proteins_tree_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE)
        else:
            raise ValueError("Unknown type!")

    def get_protein_repr_toggle_flag(self) -> int:
        return self._settings_manager.settings.proteins_tab_use_toggle

    def get_add_protein_view(self):
        return self._add_protein_view

    def get_rename_protein_view(self):
        return self._rename_protein_view

    def get_hotspots_protein_regions_view(self) -> "hotspots_protein_regions_view.HotspotsProteinRegionsView":
        return self._hotspots_protein_regions_view

    # </editor-fold>

    # <editor-fold desc="Getter methods for Protein Pairs Tab in main view">
    def get_current_protein_pair_tree_index(self):
        return self._main_view.ui.protein_pairs_tree_view.currentIndex()

    def get_child_index_of_get_current_protein_pair_tree_index(self):
        return self._main_view.ui.protein_pairs_tree_view.currentIndex().child(0, 0)

    def get_current_protein_pair_tree_index_type(self):
        return self._main_view.ui.protein_pairs_tree_view.model().data(
            self.get_current_protein_pair_tree_index(), enums.ModelEnum.TYPE_ROLE
        )

    def get_current_protein_pair_tree_index_object(self):
        return self.get_current_protein_pair_tree_index().data(enums.ModelEnum.OBJECT_ROLE)

    def get_parent_index_object_of_current_protein_pair_tree_index(self):
        return self.get_current_protein_pair_tree_index().parent().data(enums.ModelEnum.OBJECT_ROLE)

    def get_grand_parent_index_object_of_current_protein_pair_tree_index(self):
        return self.get_current_protein_pair_tree_index().parent().parent().data(enums.ModelEnum.OBJECT_ROLE)

    def get_current_active_protein_pair_object(self) -> "protein_pair.ProteinPair":
        """Returns the protein pair object of the current active branch.

        Note:
            This function is handles also the header type and does not throw an error!

        Raises:
            ValueError: if the type is unknown

        Returns:
            a protein pair object
        """
        tmp_type = self._main_view.ui.protein_pairs_tree_view.currentIndex().data(enums.ModelEnum.TYPE_ROLE)
        tmp_display_role = self._main_view.ui.protein_pairs_tree_view.currentIndex().data(Qt.DisplayRole)
        if tmp_type == "protein_pair":
            return self._main_view.ui.protein_pairs_tree_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE)
        elif tmp_type == "protein":
            return self._main_view.ui.protein_pairs_tree_view.currentIndex().parent().data(
                enums.ModelEnum.OBJECT_ROLE)
        elif tmp_type == "scene":
            return self._main_view.ui.protein_pairs_tree_view.currentIndex().parent().parent().data(
                enums.ModelEnum.OBJECT_ROLE
            )
        elif tmp_type == "chain":
            return self._main_view.ui.protein_pairs_tree_view.currentIndex().parent().parent().parent().data(
                enums.ModelEnum.OBJECT_ROLE
            )
        elif tmp_type == "header" and tmp_display_role == "Scenes":
            return self._main_view.ui.protein_pairs_tree_view.currentIndex().parent().data(
                enums.ModelEnum.OBJECT_ROLE)
        elif tmp_type == "header" and tmp_display_role == "Chains":
            return self._main_view.ui.protein_pairs_tree_view.currentIndex().parent().parent().data(
                enums.ModelEnum.OBJECT_ROLE
            )
        elif tmp_type == "":
            raise ValueError("Nothing there!")
        else:
            raise ValueError("Unknown type!")

    def get_current_active_protein_object_of_protein_pair(self) -> "protein.Protein":
        """Returns the protein object of the current active protein pair branch.

        Note:
            This function is handles also the header type and does not throw an error!

        Raises:
            ValueError: if the type is unknown

        Returns:
            a protein object
        """
        tmp_type = self._main_view.ui.protein_pairs_tree_view.currentIndex().data(enums.ModelEnum.TYPE_ROLE)
        tmp_display_role = self._main_view.ui.protein_pairs_tree_view.currentIndex().data(Qt.DisplayRole)
        if tmp_type == "protein":
            return self._main_view.ui.protein_pairs_tree_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE)
        elif tmp_type == "scene":
            raise ValueError(f"Cannot get a protein object with the type: {tmp_type}")
        elif tmp_type == "chain":
            return self._main_view.ui.protein_pairs_tree_view.currentIndex().parent().parent().data(
                enums.ModelEnum.OBJECT_ROLE
            )
        elif tmp_type == "header" and tmp_display_role == "Scenes":
            raise ValueError(f"Cannot get a protein object with the type: {tmp_type}")
        elif tmp_type == "header" and tmp_display_role == "Chains":
            return self._main_view.ui.protein_pairs_tree_view.currentIndex().parent().data(
                enums.ModelEnum.OBJECT_ROLE)
        elif tmp_type == "protein_pair":
            raise ValueError(f"Cannot get a protein object if the type is: {tmp_type}!")
        else:
            raise ValueError("Unknown type!")

    def get_current_active_scene_name_of_protein_pair(self) -> str:
        """Returns the scene name of the current active protein pair branch.

        Note:
            This function is handles also the header type and does not throw an error!

        Raises:
            ValueError: if the type is unknown

        Returns:
            a scene name
        """
        tmp_type = self._main_view.ui.protein_pairs_tree_view.currentIndex().data(enums.ModelEnum.TYPE_ROLE)
        tmp_display_role = self._main_view.ui.protein_pairs_tree_view.currentIndex().data(Qt.DisplayRole)
        if tmp_type == "protein":
            raise ValueError(f"Cannot get a scene name if the type is: {tmp_type}!")
        elif tmp_type == "scene":
            return tmp_display_role
        elif tmp_type == "chain":
            raise ValueError(f"Cannot get a scene name if the type is: {tmp_type}!")
        elif tmp_type == "header" and tmp_display_role == "Scenes":
            raise ValueError(f"Cannot get a scene name if the type is: {tmp_type}!")
        elif tmp_type == "header" and tmp_display_role == "Chains":
            raise ValueError(f"Cannot get a scene name if the type is: {tmp_type}!")
        elif tmp_type == "protein_pair":
            raise ValueError(f"Cannot get a scene name if the type is: {tmp_type}!")
        else:
            raise ValueError("Unknown type!")

    def get_current_active_chain_object_of_protein_pair(self) -> "chain.Chain":
        """Returns the chain object of the current active protein pair branch.

        Note:
            This function is handles also the header type and does not throw an error!

        Raises:
            ValueError: if the type is unknown

        Returns:
            a chain object
        """
        tmp_type = self._main_view.ui.protein_pairs_tree_view.currentIndex().data(enums.ModelEnum.TYPE_ROLE)
        tmp_display_role = self._main_view.ui.protein_pairs_tree_view.currentIndex().data(Qt.DisplayRole)
        if tmp_type == "protein":
            raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
        elif tmp_type == "scene":
            raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
        elif tmp_type == "chain":
            return self._main_view.ui.protein_pairs_tree_view.currentIndex().data(enums.ModelEnum.OBJECT_ROLE)
        elif tmp_type == "header" and tmp_display_role == "Scenes":
            raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
        elif tmp_type == "header" and tmp_display_role == "Chains":
            raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
        elif tmp_type == "protein_pair":
            raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
        else:
            raise ValueError("Unknown type!")

    def get_protein_pair_repr_toggle_flag(self) -> int:
        return self._settings_manager.settings.protein_pairs_tab_use_toggle

    # </editor-fold>

    # </editor-fold>

    def get_current_project(self) -> "project.Project":
        """Returns the current project."""
        return self._current_project

    def get_information_about_current_session(self):
        return self._current_pymol_session.session_name, self._current_pymol_session.object_type

    def get_protein_model(self):
        return self._protein_model

    # </editor-fold>

    def add_protein_to_current_project(self, a_protein: "protein.Protein"):
        self._current_project.add_existing_protein(a_protein)

    def add_protein_pair_to_current_project(self, a_protein_pair: "protein_pair.ProteinPair"):
        self._current_project.add_protein_pair(a_protein_pair)

    # <editor-fold desc="Setter Methods">
    def set_new_project(self, the_current_project: "project.Project") -> None:
        """Sets the new current project into the interface manager."""
        self._current_project = the_current_project
        self._sequence_model.clear()
        self._build_sequences_model()
        self._protein_model.clear()
        self._build_proteins_model()
        self._protein_pair_model.clear()
        self._build_protein_pairs_model()

    def set_new_workspace(self, the_current_workspace) -> None:
        """Sets the new current workspace into the interface manager."""
        self._settings_manager.settings.workspace_path = the_current_workspace
        self._workspace_model.clear()
        self._build_workspace_model()

    def set_new_session_information(self, a_session_name: str, an_object_name: str) -> None:
        self._current_pymol_session.session_name = a_session_name
        self._current_pymol_session.object_type = an_object_name

    # <editor-fold desc="Protein">
    def set_current_chain_color_for_ui_for_proteins(
            self,
            the_pymol_session_manager: "pymol_session_manager.PymolSessionManager"
    ) -> None:
        """Checks which color the protein chain has and sets the index of the combobox accordingly."""
        tmp_protein = self.get_current_active_protein_object()
        tmp_chain = self.get_current_active_chain_object()
        if the_pymol_session_manager.is_the_current_protein_in_session(self.get_current_active_protein_object().get_molecule_object()):
            # fixme: This can easily be bypassed by a power user if the first residue color is changed
            if tmp_chain.chain_type == "protein_chain":
                tmp_protein.pymol_selection.selection_string = f"first chain {tmp_chain.chain_letter} and name CA"
            else:
                tmp_protein.pymol_selection.selection_string = f"first chain {tmp_chain.chain_letter}"
            if self.pymol_session_manager.get_residue_color_config(f"(first elem N, first elem C, first elem O) and chain {tmp_chain.chain_letter}",
                                                                   tmp_chain.chain_letter).atoms_are_colored_by_elements():
                self._main_view.ui.lbl_protein_current_color.setText("By Element    ")
                self._main_view.tg_protein_color_atoms.toggle_button.setChecked(True)
            else:
                rvoid(tmp_chain.get_color(tmp_protein.pymol_selection.selection_string, self.pymol_session_manager))
                self._main_view.ui.lbl_protein_current_color.setText(f"{tmp_chain.pymol_parameters['chain_color']}    ")
                self._main_view.tg_protein_color_atoms.toggle_button.setChecked(False)

    def set_repr_state_in_ui_for_protein_chain(self, the_pymol_session_manager: "pymol_session_manager.PymolSessionManager"):
        tmp_protein = self.get_current_active_protein_object()
        tmp_chain = self.get_current_active_chain_object()
        if the_pymol_session_manager.is_the_current_protein_in_session(self.get_current_active_protein_object().get_molecule_object()):
            # fixme: This can easily be bypassed by a power user if the first residue color is changed
            tmp_protein.pymol_selection.selection_string = f"first chain {tmp_chain.chain_letter}"
            print(f"This is a chain type: {tmp_chain.chain_type}")
            if tmp_chain.chain_type == "protein_chain":
                self._main_view.ui.frame_protein_repr.setEnabled(True)
                tmp_repr_state = self.pymol_session_manager.get_chain_repr_state(tmp_protein.pymol_selection.selection_string,
                                                                                 tmp_chain.chain_letter)
                if self._settings_manager.settings.proteins_tab_use_toggle == 1:
                    self.manage_toggle_state_of_protein_repr(tmp_repr_state)
                else:
                    self.manage_check_state_of_protein_repr(tmp_repr_state)
            else:
                self._main_view.ui.frame_protein_repr.setEnabled(False)
                self._main_view.tg_protein_cartoon.toggle_button.setChecked(False)
                self._main_view.tg_protein_ribbon.toggle_button.setChecked(False)
                self._main_view.tg_protein_sticks.toggle_button.setChecked(False)
                self._main_view.tg_protein_lines.toggle_button.setChecked(False)
                self._main_view.tg_protein_spheres.toggle_button.setChecked(False)
                self._main_view.tg_protein_dots.toggle_button.setChecked(False)
                self._main_view.tg_protein_mesh.toggle_button.setChecked(False)
                self._main_view.tg_protein_surface.toggle_button.setChecked(False)
    # </editor-fold>

    def set_repr_state_in_ui_for_protein_pair_chain(self,
                                                    a_protein_name: str,
                                                    the_pymol_session_manager: "pymol_session_manager.PymolSessionManager"):
        tmp_protein = self.get_current_active_protein_object_of_protein_pair()
        tmp_chain = self.get_current_active_chain_object_of_protein_pair()
        if the_pymol_session_manager.is_the_current_protein_pair_in_session(self.get_current_active_protein_pair_object().name):
            # fixme: This can easily be bypassed by a power user if the first residue color is changed
            tmp_protein.pymol_selection.selection_string = f"first chain {tmp_chain.chain_letter} and {a_protein_name}"
            if tmp_chain.chain_type == "protein_chain":
                self._main_view.ui.frame_protein_pair_repr.setEnabled(True)
                tmp_repr_state = self.pymol_session_manager.get_chain_repr_state(tmp_protein.pymol_selection.selection_string,
                                                                                 tmp_chain.chain_letter)
                if self._settings_manager.settings.protein_pairs_tab_use_toggle == 1:
                    self.manage_toggle_state_of_protein_pair_repr(tmp_repr_state)
                else:
                    self.manage_check_state_of_protein_pair(tmp_repr_state)
            else:
                self._main_view.ui.frame_protein_pair_repr.setEnabled(False)
                self._main_view.tg_protein_pair_cartoon.toggle_button.setChecked(False)
                self._main_view.tg_protein_pair_ribbon.toggle_button.setChecked(False)
                self._main_view.tg_protein_pair_sticks.toggle_button.setChecked(False)
                self._main_view.tg_protein_pair_lines.toggle_button.setChecked(False)
                self._main_view.tg_protein_pair_spheres.toggle_button.setChecked(False)
                self._main_view.tg_protein_pair_dots.toggle_button.setChecked(False)
                self._main_view.tg_protein_pair_mesh.toggle_button.setChecked(False)
                self._main_view.tg_protein_pair_surface.toggle_button.setChecked(False)

    def manage_check_state_of_protein_pair(self, tmp_repr_state):
        if tmp_repr_state[enums.PyMOLRepresentation.CARTOON.value] == 0:
            self._main_view.ui.cb_protein_pair_cartoon.setChecked(False)
        else:
            self._main_view.ui.cb_protein_pair_cartoon.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.STICKS.value] == 0:
            self._main_view.ui.cb_protein_pair_sticks.setChecked(False)
        else:
            self._main_view.ui.cb_protein_pair_sticks.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.RIBBON.value] == 0:
            self._main_view.ui.cb_protein_pair_ribbon.setChecked(False)
        else:
            self._main_view.ui.cb_protein_pair_ribbon.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.LINES.value] == 0:
            self._main_view.ui.cb_protein_pair_lines.setChecked(False)
        else:
            self._main_view.ui.cb_protein_pair_lines.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.SPHERES.value] == 0:
            self._main_view.ui.cb_protein_pair_spheres.setChecked(False)
        else:
            self._main_view.ui.cb_protein_pair_spheres.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.DOTS.value] == 0:
            self._main_view.ui.cb_protein_pair_dots.setChecked(False)
        else:
            self._main_view.ui.cb_protein_pair_dots.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.MESH.value] == 0:
            self._main_view.ui.cb_protein_pair_mesh.setChecked(False)
        else:
            self._main_view.ui.cb_protein_pair_mesh.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.SURFACE.value] == 0:
            self._main_view.ui.cb_protein_pair_surface.setChecked(False)
        else:
            self._main_view.ui.cb_protein_pair_surface.setChecked(True)

    def manage_toggle_state_of_protein_pair_repr(self, tmp_repr_state):
        if tmp_repr_state[enums.PyMOLRepresentation.CARTOON.value] == 0:
            self._main_view.tg_protein_pair_cartoon.toggle_button.setChecked(False)
        else:
            self._main_view.tg_protein_pair_cartoon.toggle_button.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.STICKS.value] == 0:
            self._main_view.tg_protein_pair_sticks.toggle_button.setChecked(False)
        else:
            self._main_view.tg_protein_pair_sticks.toggle_button.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.RIBBON.value] == 0:
            self._main_view.tg_protein_pair_ribbon.toggle_button.setChecked(False)
        else:
            self._main_view.tg_protein_pair_ribbon.toggle_button.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.LINES.value] == 0:
            self._main_view.tg_protein_pair_lines.toggle_button.setChecked(False)
        else:
            self._main_view.tg_protein_pair_lines.toggle_button.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.SPHERES.value] == 0:
            self._main_view.tg_protein_pair_spheres.toggle_button.setChecked(False)
        else:
            self._main_view.tg_protein_pair_spheres.toggle_button.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.DOTS.value] == 0:
            self._main_view.tg_protein_pair_dots.toggle_button.setChecked(False)
        else:
            self._main_view.tg_protein_pair_dots.toggle_button.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.MESH.value] == 0:
            self._main_view.tg_protein_pair_mesh.toggle_button.setChecked(False)
        else:
            self._main_view.tg_protein_pair_mesh.toggle_button.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.SURFACE.value] == 0:
            self._main_view.tg_protein_pair_surface.toggle_button.setChecked(False)
        else:
            self._main_view.tg_protein_pair_surface.toggle_button.setChecked(True)

    # </editor-fold>

    # <editor-fold desc="Build Methods">
    def _build_sequences_model(self):
        if len(self._current_project.sequences) > 0:
            tmp_root_item = self._sequence_model.invisibleRootItem()
            for tmp_sequence in self._current_project.sequences:
                tmp_sequence_item = QtGui.QStandardItem(tmp_sequence.name)
                tmp_sequence_item.setData(tmp_sequence, enums.ModelEnum.OBJECT_ROLE)
                if "," in tmp_sequence.seq:
                    tmp_sequence_item.setData(enums.ModelTypeEnum.MULTIMER_SEQ, enums.ModelEnum.TYPE_ROLE)
                else:
                    tmp_sequence_item.setData(enums.ModelTypeEnum.MONOMER_SEQ, enums.ModelEnum.TYPE_ROLE)
                tmp_root_item.appendRow(tmp_sequence_item)

    def _build_proteins_model(self) -> None:
        if len(self._current_project.proteins) > 0:
            tmp_main_socket, tmp_general_purpose_socket = self.job_manager.get_general_purpose_socket_pair()
            self._protein_model.build_model_from_scratch(self._current_project.proteins,
                                                         tmp_main_socket,
                                                         tmp_general_purpose_socket)
            # tmp_root_item = self._protein_model.invisibleRootItem()
            # for tmp_protein in self._current_project.proteins:
            #     tmp_protein_item = QtGui.QStandardItem(tmp_protein.get_molecule_object())
            #     tmp_protein_item.setData(tmp_protein, enums.ModelEnum.OBJECT_ROLE)
            #     tmp_protein_item.setData("protein", enums.ModelEnum.TYPE_ROLE)
            #     tmp_root_item.appendRow(tmp_protein_item)
            #     tmp_scenes_item = QtGui.QStandardItem("Scenes")
            #     tmp_scenes_item.setData("header", enums.ModelEnum.TYPE_ROLE)
            #     tmp_protein_item.appendRow(tmp_scenes_item)
            #     tmp_protein.load_protein_pymol_session()
            #     for tmp_scene in pymol_io.get_all_scenes_from_pymol_session():
            #         tmp_scene_item = QtGui.QStandardItem(tmp_scene)
            #         tmp_scene_item.setData("scene", enums.ModelEnum.TYPE_ROLE)
            #         tmp_scenes_item.appendRow(tmp_scene_item)
            #     tmp_chains_item = QtGui.QStandardItem("Chains")
            #     tmp_chains_item.setData("header", enums.ModelEnum.TYPE_ROLE)
            #     tmp_protein_item.appendRow(tmp_chains_item)
            #     for tmp_chain in tmp_protein.chains:
            #         tmp_chain_item = QtGui.QStandardItem(tmp_chain.chain_letter)
            #         tmp_chain_item.setData(tmp_chain, enums.ModelEnum.OBJECT_ROLE)
            #         tmp_chain_item.setData("chain", enums.ModelEnum.TYPE_ROLE)
            #         tmp_chains_item.appendRow(tmp_chain_item)

    def _build_protein_pairs_model(self) -> None:
        if len(self._current_project.protein_pairs) > 0:
            tmp_main_socket, tmp_general_purpose_socket = self.job_manager.get_general_purpose_socket_pair()
            self._protein_pair_model.build_model_from_scratch(self._current_project.protein_pairs,
                                                              tmp_main_socket,
                                                              tmp_general_purpose_socket)
            # tmp_root_item = self._protein_pair_model.invisibleRootItem()
            # for tmp_protein_pair in self._current_project.protein_pairs:
            #     tmp_protein_pair_item = QtGui.QStandardItem(tmp_protein_pair.name)
            #     tmp_protein_pair_item.setData(tmp_protein_pair, enums.ModelEnum.OBJECT_ROLE)
            #     tmp_protein_pair_item.setData("protein_pair", enums.ModelEnum.TYPE_ROLE)
            #     tmp_scenes_item = QtGui.QStandardItem("Scenes")
            #     tmp_scenes_item.setData("header", enums.ModelEnum.TYPE_ROLE)
            #     tmp_protein_pair_item.appendRow(tmp_scenes_item)
            #     tmp_protein_pair.load_pymol_session()
            #     for tmp_scene in pymol_io.get_all_scenes_from_pymol_session():
            #         tmp_scene_item = QtGui.QStandardItem(tmp_scene)
            #         tmp_scene_item.setData("scene", enums.ModelEnum.TYPE_ROLE)
            #         tmp_scenes_item.appendRow(tmp_scene_item)
            #     # Create protein 1 item
            #     tmp_protein_item_1 = QtGui.QStandardItem(tmp_protein_pair.protein_1.get_molecule_object())
            #     tmp_protein_item_1.setData(tmp_protein_pair.protein_1, enums.ModelEnum.OBJECT_ROLE)
            #     tmp_protein_item_1.setData("protein", enums.ModelEnum.TYPE_ROLE)
            #     tmp_chains_item_1 = QtGui.QStandardItem("Chains")
            #     tmp_chains_item_1.setData("header", enums.ModelEnum.TYPE_ROLE)
            #     tmp_protein_item_1.appendRow(tmp_chains_item_1)
            #     for tmp_chain in tmp_protein_pair.protein_1.chains:
            #         tmp_chain_item = QtGui.QStandardItem(tmp_chain.chain_letter)
            #         tmp_chain_item.setData(tmp_chain, enums.ModelEnum.OBJECT_ROLE)
            #         tmp_chain_item.setData("chain", enums.ModelEnum.TYPE_ROLE)
            #         tmp_chains_item_1.appendRow(tmp_chain_item)
            #     # Create protein 2 item
            #     tmp_protein_item_2 = QtGui.QStandardItem(tmp_protein_pair.protein_2.get_molecule_object())
            #     tmp_protein_item_2.setData(tmp_protein_pair.protein_2, enums.ModelEnum.OBJECT_ROLE)
            #     tmp_protein_item_2.setData("protein", enums.ModelEnum.TYPE_ROLE)
            #     tmp_chains_item_2 = QtGui.QStandardItem("Chains")
            #     tmp_chains_item_2.setData("header", enums.ModelEnum.TYPE_ROLE)
            #     tmp_protein_item_2.appendRow(tmp_chains_item_2)
            #     for tmp_chain in tmp_protein_pair.protein_2.chains:
            #         tmp_chain_item = QtGui.QStandardItem(tmp_chain.chain_letter)
            #         tmp_chain_item.setData(tmp_chain, enums.ModelEnum.OBJECT_ROLE)
            #         tmp_chain_item.setData("chain", enums.ModelEnum.TYPE_ROLE)
            #         tmp_chains_item_2.appendRow(tmp_chain_item)
            #     tmp_protein_pair_item.appendRow(tmp_protein_item_1)
            #     tmp_protein_pair_item.appendRow(tmp_protein_item_2)
            #     tmp_root_item.appendRow(tmp_protein_pair_item)
    # </editor-fold>

    # <editor-fold desc="Refresh Methods">
    def refresh_main_view(self):
        """Modifies the UI of the main view based on an app model."""
        self._main_view.ui.project_tab_widget.setEnabled(True)
        self._main_view.ui.lbl_logo.hide()
        # Settings
        self._main_view.ui.menuSettings.setEnabled(True)
        self._main_view.ui.action_edit_settings.setEnabled(True)
        self._main_view.ui.action_restore_settings.setEnabled(True)
        # Help
        self._main_view.ui.menuAbout.setEnabled(True)
        self._main_view.ui.action_documentation.setEnabled(True)
        self._main_view.ui.action_get_demo_projects.setEnabled(True)
        self._main_view.ui.action_show_log_in_explorer.setEnabled(True)
        self._main_view.ui.action_clear_logs.setEnabled(True)
        self._main_view.ui.action_about.setEnabled(True)

        # Protein hotspots
        self._main_view.ui.lbl_protein_protein_regions.hide()
        self._main_view.ui.frame_protein_protein_regions.hide()
        self._main_view.ui.lbl_protein_pair_protein_regions.hide()
        self._main_view.ui.frame_protein_pair_protein_regions.hide()

        self._main_view.ui.menuProject.setEnabled(True)
        # A project is open
        if self._current_project.get_project_name() != "":
            styles.set_stylesheet(self._main_view)
            self._main_view.ui.lbl_project_name.show()
            self._main_view.ui.lbl_project_name.setText(f"Project Name: {self._current_project.get_project_name()}")
            self._main_view.ui.lbl_session_name.show()
            self._main_view.ui.project_tab_widget.show()
            # Project options
            self._main_view.ui.action_new_project.setEnabled(False)
            self._main_view.ui.action_open_project.setEnabled(False)
            self._main_view.ui.action_use_project.setEnabled(True)
            self._main_view.ui.action_delete_project.setEnabled(False)
            self._main_view.ui.action_import_project.setEnabled(False)
            self._main_view.ui.action_export_project.setEnabled(True)
            self._main_view.ui.action_close_project.setEnabled(True)
            # Demo Projects
            self._main_view.ui.action_get_demo_projects.setEnabled(False)
            # Sequence objects
            if len(self._current_project.sequences) > 0:
                # A project has sequence(s)
                self._main_view.ui.seqs_list_view.setModel(self._sequence_model)
                # It is possible to do a prediction
                self._main_view.ui.menuPrediction.setEnabled(True)

                # <editor-fold desc="Checks type(s) of sequences">
                tmp_sequence_model_state = self._check_sequence_model_state()
                if tmp_sequence_model_state == "monomer":
                    self._main_view.ui.action_predict_monomer.setEnabled(True)
                    self._main_view.ui.action_predict_multimer.setEnabled(False)
                elif tmp_sequence_model_state == "multimer":
                    self._main_view.ui.action_predict_monomer.setEnabled(False)
                    self._main_view.ui.action_predict_multimer.setEnabled(True)
                elif tmp_sequence_model_state == "both":
                    self._main_view.ui.action_predict_monomer.setEnabled(True)
                    self._main_view.ui.action_predict_multimer.setEnabled(True)
                elif tmp_sequence_model_state == "nothing":
                    self._main_view.ui.action_predict_monomer.setEnabled(False)
                    self._main_view.ui.action_predict_multimer.setEnabled(False)

                tmp_are_all_monomer_sequences_predicted = self._check_if_sequences_are_already_predicted()[0]
                tmp_are_all_multimer_sequences_predicted = self._check_if_sequences_are_already_predicted()[1]
                if tmp_are_all_monomer_sequences_predicted:
                    self._main_view.ui.action_predict_monomer.setEnabled(False)
                if tmp_are_all_multimer_sequences_predicted:
                    self._main_view.ui.action_predict_multimer.setEnabled(False)
                if tmp_are_all_monomer_sequences_predicted and tmp_are_all_multimer_sequences_predicted:
                    self._main_view.ui.menuPrediction.setEnabled(
                        False)  # the entire menu can be disabled because all sequences are predicted
                # </editor-fold>
            else:
                # A project has no sequence(s)
                self._sequence_model = QtGui.QStandardItemModel()
                self._main_view.ui.seqs_list_view.setModel(self._sequence_model)
                # It isn't possible to do a prediction
                self._main_view.ui.menuPrediction.setEnabled(False)
            # Protein objects
            if len(self._current_project.proteins) > 0:
                # A project has protein(s)
                self._main_view.ui.proteins_tree_view.setModel(self._protein_model)
                self._main_view.ui.proteins_tree_view.setHeaderHidden(True)
                # It is possible to do an analysis, image and hotspots
                # Analysis
                self._main_view.ui.menuAnalysis.setEnabled(True)
                self._main_view.ui.action_distance_analysis.setEnabled(True)
                # Image/ Hotspots
                if self.pymol_session_manager.is_the_current_session_empty():
                    self._main_view.ui.menuImage.setEnabled(False)
                    self._main_view.ui.menuHotspots.setEnabled(False)
                else:
                    self._main_view.ui.menuImage.setEnabled(True)
                    self._main_view.ui.action_preview_image.setEnabled(True)
                    self._main_view.ui.action_ray_tracing_image.setEnabled(True)
                    self._main_view.ui.action_simple_image.setEnabled(True)
                    # Hotspots
                    self._main_view.ui.menuHotspots.setEnabled(True)
                    self._main_view.ui.action_protein_regions.setEnabled(True)

                if len(self._current_project.protein_pairs) > 0:
                    # A project has protein pair(s)
                    self._main_view.ui.protein_pairs_tree_view.setModel(self._protein_pair_model)
                    self._main_view.ui.protein_pairs_tree_view.setHeaderHidden(True)
                    # It is possible to view results
                    self._main_view.ui.menuResults.setEnabled(True)
                else:
                    # A project has no protein pair(s)
                    # Protein Pairs tab
                    self._main_view.ui.btn_delete_protein_pair.setEnabled(False)
                    self._main_view.ui.btn_open_protein_pair_session.setEnabled(False)
                    self._main_view.ui.btn_create_protein_pair_scene.setEnabled(False)
                    self._main_view.ui.btn_update_protein_pair_scene.setEnabled(False)
                    self._main_view.ui.btn_delete_protein_pair_scene.setEnabled(False)
                    # It isn't possible to view results
                    self._main_view.ui.menuResults.setEnabled(False)
            else:
                # A project has no protein(s)
                # Proteins tab
                self._main_view.ui.btn_delete_protein.setEnabled(False)
                self._main_view.ui.btn_save_protein.setEnabled(False)
                self._main_view.ui.btn_open_protein_session.setEnabled(False)
                self._main_view.ui.btn_create_protein_scene.setEnabled(False)
                self._main_view.ui.btn_update_protein_scene.setEnabled(False)
                self._main_view.ui.btn_delete_protein_scene.setEnabled(False)
                # It isn't possible to do an analysis, image and hotspots
                # Analysis
                self._main_view.ui.menuAnalysis.setEnabled(False)
                # Image
                self._main_view.ui.menuImage.setEnabled(False)
                # Hotspots
                self._main_view.ui.menuHotspots.setEnabled(False)
        # No project open
        else:
            # Homepage view
            styles.set_stylesheet_homepage(self._main_view)
            # No project is open
            # No project(s) available
            self._main_view.ui.lbl_logo.show()
            self._main_view.ui.lbl_project_name.hide()
            self._main_view.ui.project_tab_widget.hide()
            if len(self.get_workspace_projects_as_list()) == 0:
                # Project
                self._main_view.ui.action_new_project.setEnabled(True)
                self._main_view.ui.action_open_project.setEnabled(False)
                self._main_view.ui.action_use_project.setEnabled(False)
                self._main_view.ui.action_delete_project.setEnabled(False)
                self._main_view.ui.action_import_project.setEnabled(True)
                self._main_view.ui.action_export_project.setEnabled(False)
                self._main_view.ui.action_close_project.setEnabled(False)
            else:
                # Project
                self._main_view.ui.action_new_project.setEnabled(True)
                self._main_view.ui.action_open_project.setEnabled(True)
                self._main_view.ui.action_use_project.setEnabled(False)
                self._main_view.ui.action_delete_project.setEnabled(True)
                self._main_view.ui.action_import_project.setEnabled(True)
                self._main_view.ui.action_export_project.setEnabled(False)
                self._main_view.ui.action_close_project.setEnabled(False)
            # Prediction
            self._main_view.ui.menuPrediction.setEnabled(False)
            self._main_view.ui.action_predict_monomer.setEnabled(False)
            self._main_view.ui.action_predict_multimer.setEnabled(False)
            # Analysis
            self._main_view.ui.menuAnalysis.setEnabled(False)
            self._main_view.ui.action_distance_analysis.setEnabled(True)
            # Results
            self._main_view.ui.menuResults.setEnabled(False)
            self._main_view.ui.action_results_summary.setEnabled(False)
            # Image
            self._main_view.ui.menuImage.setEnabled(False)
            self._main_view.ui.action_preview_image.setEnabled(False)
            self._main_view.ui.action_ray_tracing_image.setEnabled(False)
            self._main_view.ui.action_simple_image.setEnabled(False)
            # Hotspots
            self._main_view.ui.menuHotspots.setEnabled(False)
            self._main_view.ui.action_protein_regions.setEnabled(False)

        # # Menu bar view for prediction
        # if self.main_tasks_manager.prediction_task is not None:
        #     if not self.main_tasks_manager.check_if_prediction_task_is_finished():
        #         # A prediction is running.
        #         self._main_view.ui.menuProject.setEnabled(False)
        #         self._main_view.ui.action_predict_monomer.setEnabled(False)
        #         self._main_view.ui.action_predict_multimer.setEnabled(False)
        #         self._main_view.ui.action_abort_prediction.setEnabled(True)
        #         self._main_view.ui.menuAnalysis.setEnabled(False)
        #         self._main_view.ui.menuImage.setEnabled(False)
        #         # Check for existing protein(s)
        #         if len(self._current_project.proteins) > 0:
        #             self._main_view.ui.menuHotspots.setEnabled(True)
        #         else:
        #             self._main_view.ui.menuHotspots.setEnabled(False)
        #         # Check for existing results
        #         if len(self._current_project.protein_pairs) > 0:
        #             self._main_view.ui.menuResults.setEnabled(True)
        #             self._main_view.ui.action_results_summary.setEnabled(True)
        #         else:
        #             self._main_view.ui.menuResults.setEnabled(False)
        #     else:
        #         # A prediction is finished or not running.
        #         pass
        #
        # # Menu bar view for analysis
        # if self.main_tasks_manager.distance_analysis_task is not None:
        #     if not self.main_tasks_manager.check_if_distance_analysis_task_is_finished():
        #         # An analysis is running.
        #         self._main_view.ui.menuProject.setEnabled(False)
        #         self._main_view.ui.menuPrediction.setEnabled(False)
        #         self._main_view.ui.menuAnalysis.setEnabled(False)
        #         self._main_view.ui.menuImage.setEnabled(False)
        #         self._main_view.ui.menuHotspots.setEnabled(False)
        #         # Check for existing results
        #         if len(self._current_project.protein_pairs) > 0:
        #             self._main_view.ui.menuResults.setEnabled(True)
        #             self._main_view.ui.action_results_summary.setEnabled(True)
        #         else:
        #             self._main_view.ui.menuResults.setEnabled(False)
        #     else:
        #         # An analysis is finished or not running.
        #         pass
        #
        # # TODO: this could be not needed
        # #self._main_view.ui.project_tab_widget.setCurrentIndex(0)
        #
        # # # Fixme: The following code below is for clicking on sequence, protein or protein pair with available options.
        # # # Fixme: Is this important with the code above?
        # # self._main_view.ui.project_tab_widget.setCurrentIndex(self.current_tab_index)
        # # # Sequence
        # # if len(self._current_project.sequences) > 0 and self._main_view.ui.seqs_list_view.currentIndex().data(Qt.DisplayRole) is not None:
        # #     # There are sequences in the project and there is also one selected
        # #     self._main_view.ui.seqs_list_view.setModel(self._sequence_model)
        # #     self.show_menu_options_with_seq()
        # # elif len(self._current_project.sequences) > 0 and self._main_view.ui.seqs_list_view.currentIndex().data(Qt.DisplayRole) is None:
        # #     # There are sequences in the project and there is NO one selected
        # #     self._main_view.ui.seqs_list_view.setModel(self._sequence_model)
        # #     self.show_menu_options_without_seq()
        # # else:
        # #     # There are no sequences in the project
        # #     self.show_menu_options_without_seq()
        # #     self._sequence_model = QtGui.QStandardItemModel()
        # #     self._main_view.ui.seqs_list_view.setModel(self._sequence_model)
        # #
        # # # Proteins Tab
        # # if len(self._current_project.proteins) > 0:
        # #     self.show_menu_options_with_protein()
        # # if len(self._current_project.proteins) == 0 or self._main_view.ui.proteins_tree_view.currentIndex().data(Qt.DisplayRole) is None:
        # #     self.show_menu_options_without_protein()
        # #
        # # # Protein Pairs Tab
        # # if len(self._current_project.protein_pairs) > 0:
        # #     self.show_menu_options_with_protein_pair()
        # # if len(self._current_project.protein_pairs) == 0 or self._main_view.ui.protein_pairs_tree_view.currentIndex().data(Qt.DisplayRole) is None:
        # #     self.show_menu_options_without_protein_pair()
        #
        # # tmp_projects = self.get_workspace_projects_as_list()
        # # if len(tmp_projects) > 0 and not self._main_view.ui.lbl_logo.isHidden():
        # #     self._main_view.ui.action_open_project.setEnabled(True)
        # #     self._main_view.ui.action_delete_project.setEnabled(True)
        # # else:
        # #     self._main_view.ui.action_open_project.setEnabled(False)
        # #     self._main_view.ui.action_delete_project.setEnabled(False)
        #
        # # Fixme: Is the text in statusbar and logger not enough?
        # if self.main_tasks_manager.prediction_task is not None:
        #     if not self.main_tasks_manager.check_if_prediction_task_is_finished():
        #         logger.info("Running prediction in the background ...")
        #         self._main_view.ui.action_abort_prediction.setEnabled(True)
        #         self._main_view.ui.action_predict_monomer.setEnabled(False)
        #         self._main_view.ui.action_predict_multimer.setEnabled(False)
        #     else:
        #         logger.info("Prediction finished.")
        #         self._main_view.status_bar.setStyleSheet("""
        #             QStatusBar {
        #                 background-color: white;
        #                 border-style: solid;
        #                 border-width: 2px;
        #                 border-radius: 4px;
        #                 border-color: #DCDBE3;
        #             }
        #         """)
        #         self._main_view.ui.action_abort_prediction.setEnabled(False)
        # else:
        #     self._main_view.ui.action_abort_prediction.setEnabled(False)

    def refresh_workspace_model(self):
        self._workspace_model.clear()
        self._build_workspace_model()

    def refresh_sequence_model(self):
        self._sequence_model.clear()
        self._build_sequences_model()

    def refresh_protein_model(self):
        self._protein_model.clear()
        self._build_proteins_model()

    def refresh_protein_pair_model(self):
        self._protein_pair_model.clear()
        self._build_protein_pairs_model()

    def _build_workspace_model(self) -> None:
        tmp_workspace = self._settings_manager.settings.workspace_path
        db_pattern = os.path.join(tmp_workspace, '*.db')
        tmp_root_item = self._workspace_model.invisibleRootItem()
        for tmp_filename in [os.path.basename(file).replace(".db", "") for file in glob.glob(db_pattern)]:
            tmp_project_item = QtGui.QStandardItem(tmp_filename)
            tmp_filepath = pathlib.Path(f"{tmp_workspace}/{tmp_filename}.db")
            tmp_project_item.setData(tmp_filepath, enums.ModelEnum.FILEPATH_ROLE)
            tmp_root_item.appendRow(tmp_project_item)

    def add_project_to_workspace_model(self, a_filename):
        tmp_root_item = self._workspace_model.invisibleRootItem()
        tmp_project_item = QtGui.QStandardItem(a_filename)
        tmp_filepath = pathlib.Path(f"{self._settings_manager.settings.workspace_path}/{a_filename}.db")
        tmp_project_item.setData(tmp_filepath, enums.ModelEnum.FILEPATH_ROLE)
        tmp_root_item.appendRow(tmp_project_item)

    def disable_proteins_tab_buttons(self):
        self._main_view.ui.btn_save_protein.setEnabled(False)
        self._main_view.ui.btn_delete_protein.setEnabled(False)
        self._main_view.ui.btn_open_protein_session.setEnabled(False)
        self._main_view.ui.btn_create_protein_scene.setEnabled(False)
        self._main_view.ui.btn_update_protein_scene.setEnabled(False)
        self._main_view.ui.btn_delete_protein_scene.setEnabled(False)

    def disable_protein_pairs_tab_buttons(self):
        self._main_view.ui.btn_save_protein_pair.setEnabled(False)
        self._main_view.ui.btn_delete_protein_pair.setEnabled(False)
        self._main_view.ui.btn_open_protein_pair_session.setEnabled(False)
        self._main_view.ui.btn_create_protein_pair_scene.setEnabled(False)
        self._main_view.ui.btn_update_protein_pair_scene.setEnabled(False)
        self._main_view.ui.btn_delete_protein_pair_scene.setEnabled(False)

    # </editor-fold>

    # <editor-fold desc="Progress bar methods">
    def update_progress_bar(self, value: int, message: str):
        if value < 0 or value > 100:
            raise ValueError("Value for progress bar must be between 0 and 100!")
        self._main_view.progress_bar.show()
        self._main_view.progress_bar.setFormat(message)
        self._main_view.progress_bar.setValue(value)

    def hide_progress_bar(self):
        self._main_view.progress_bar.hide()

    # </editor-fold>

    # <editor-fold desc="Sequences">
    def show_sequence_parameters(self, a_sequence_item: QtGui.QStandardItem):
        self._main_view.setup_sequences_table(2)
        tmp_sequence = a_sequence_item.data(enums.ModelEnum.OBJECT_ROLE)
        # Table label items
        tmp_name_label_item = QtWidgets.QTableWidgetItem("Name")
        self._main_view.ui.seqs_table_widget.setItem(0, 0, tmp_name_label_item)
        tmp_sequence_label_item = QtWidgets.QTableWidgetItem("Sequence")
        self._main_view.ui.seqs_table_widget.setItem(1, 0, tmp_sequence_label_item)
        # Table value items
        tmp_seq_name_item = QtWidgets.QTableWidgetItem(tmp_sequence.name)
        tmp_seq_name_item.setToolTip("Click to edit name")
        self._main_view.ui.seqs_table_widget.setItem(0, 1, tmp_seq_name_item)
        tmp_sequence_item = QtWidgets.QTableWidgetItem(f"{tmp_sequence.seq[:15]} ...")
        tmp_sequence_item.setToolTip("Click to view complete sequence")
        tmp_sequence_item.setData(enums.ModelEnum.OBJECT_ROLE, tmp_sequence)
        self._main_view.ui.seqs_table_widget.setItem(1, 1, tmp_sequence_item)
        # Table item flags
        tmp_seq_name_item.setFlags(tmp_seq_name_item.flags() & ~Qt.ItemIsEditable)
        tmp_sequence_item.setFlags(tmp_sequence_item.flags() & ~Qt.ItemIsEditable)
        tmp_name_label_item.setFlags(tmp_name_label_item.flags() & ~Qt.ItemIsEditable)
        tmp_sequence_label_item.setFlags(tmp_sequence_label_item.flags() & ~Qt.ItemIsEditable)
        self._main_view.ui.seqs_table_widget.resizeColumnsToContents()

    def show_menu_options_with_seq(self):
        self._main_view.ui.menuAnalysis.setEnabled(False)
        self._main_view.ui.menuResults.setEnabled(False)
        self._main_view.ui.menuImage.setEnabled(False)
        self._main_view.ui.menuHotspots.setEnabled(False)

    def show_menu_options_without_seq(self):
        self._main_view.ui.btn_save_sequence.setEnabled(False)
        self._main_view.ui.btn_delete_sequence.setEnabled(False)

        # <editor-fold desc="Checks type(s) of sequences">
        tmp_sequence_model_state = self._check_sequence_model_state()
        if tmp_sequence_model_state == "monomer":
            self._main_view.ui.action_predict_monomer.setEnabled(True)
            self._main_view.ui.action_predict_multimer.setEnabled(False)
        elif tmp_sequence_model_state == "multimer":
            self._main_view.ui.action_predict_monomer.setEnabled(False)
            self._main_view.ui.action_predict_multimer.setEnabled(True)
        elif tmp_sequence_model_state == "both":
            self._main_view.ui.action_predict_monomer.setEnabled(True)
            self._main_view.ui.action_predict_multimer.setEnabled(True)
        elif tmp_sequence_model_state == "nothing":
            self._main_view.ui.action_predict_monomer.setEnabled(False)
            self._main_view.ui.action_predict_multimer.setEnabled(False)

        # </editor-fold>

    def _check_sequence_model_state(self) -> str:
        """Checks what type(s) of sequences are in the sequence model.

        Returns:
            a string representing the values: "both", "monomer", "multimer", "nothing"
        """
        tmp_contains_monomer = False
        tmp_contains_multimer = False
        for tmp_row in range(self._sequence_model.rowCount()):
            tmp_item = self._sequence_model.item(tmp_row, 0)
            if tmp_item.data(enums.ModelEnum.TYPE_ROLE) == enums.ModelTypeEnum.MONOMER_SEQ:
                tmp_contains_monomer = True
            elif tmp_item.data(enums.ModelEnum.TYPE_ROLE) == enums.ModelTypeEnum.MULTIMER_SEQ:
                tmp_contains_multimer = True

        if tmp_contains_monomer and tmp_contains_multimer:
            return "both"
        elif tmp_contains_monomer and not tmp_contains_multimer:
            return "monomer"
        elif not tmp_contains_monomer and tmp_contains_multimer:
            return "multimer"
        else:
            return "nothing"

    def _check_if_sequences_are_already_predicted(self) -> tuple[bool, bool]:
        """Checks if sequences are already predicted or not.

        Returns:
            a tuple of booleans where True means all sequences of the type are already predicted.
        """
        tmp_number_of_monomer_sequences = 0
        tmp_number_of_multimer_sequences = 0
        tmp_number_of_predicted_monomer_sequences = 0
        tmp_number_of_predicted_multimer_sequences = 0
        for tmp_row in range(self._sequence_model.rowCount()):
            tmp_item = self._sequence_model.item(tmp_row, 0)
            if tmp_item.data(enums.ModelEnum.TYPE_ROLE) == enums.ModelTypeEnum.MONOMER_SEQ:
                tmp_number_of_monomer_sequences += 1
            elif tmp_item.data(enums.ModelEnum.TYPE_ROLE) == enums.ModelTypeEnum.MULTIMER_SEQ:
                tmp_number_of_multimer_sequences += 1
            else:
                logger.warning("Invalid TYPE_ROLE!")

            if self._current_project.is_sequence_as_protein_in_project(tmp_item.data(enums.ModelEnum.OBJECT_ROLE).name) and tmp_item.data(enums.ModelEnum.TYPE_ROLE) == enums.ModelTypeEnum.MONOMER_SEQ:
                tmp_number_of_predicted_monomer_sequences += 1
            elif self._current_project.is_sequence_as_protein_in_project(tmp_item.data(enums.ModelEnum.OBJECT_ROLE).name) and tmp_item.data(enums.ModelEnum.TYPE_ROLE) == enums.ModelTypeEnum.MULTIMER_SEQ:
                tmp_number_of_predicted_multimer_sequences += 1

        if tmp_number_of_predicted_monomer_sequences == tmp_number_of_monomer_sequences:
            tmp_monomer_sequence_flag = True
        else:
            tmp_monomer_sequence_flag = False
        if tmp_number_of_predicted_multimer_sequences == tmp_number_of_multimer_sequences:
            tmp_multimer_sequence_flag = True
        else:
            tmp_multimer_sequence_flag = False
        return tmp_monomer_sequence_flag, tmp_multimer_sequence_flag

    # </editor-fold>

    # <editor-fold desc="Proteins">
    def add_protein_to_proteins_model(self, a_protein):
        self._protein_model.add_protein(a_protein)

    def remove_protein_from_proteins_model(self):
        self._protein_model.remove_protein(self.get_current_protein_tree_index())

    # <editor-fold desc="Menu Options">
    def show_menu_options_with_protein(self):
        self._main_view.ui.menuAnalysis.setEnabled(True)
        self._main_view.ui.menuResults.setEnabled(True)
        self._main_view.ui.menuImage.setEnabled(True)
        self._main_view.ui.menuHotspots.setEnabled(True)

        self._main_view.ui.proteins_tree_view.setModel(self._protein_model)
        self._main_view.ui.proteins_tree_view.setHeaderHidden(True)

    def show_menu_options_without_protein(self):
        self._main_view.ui.btn_save_protein.setEnabled(False)
        self._main_view.ui.btn_delete_protein.setEnabled(False)
        self._main_view.ui.btn_open_protein_session.setEnabled(False)
        self._main_view.ui.btn_create_protein_scene.setEnabled(False)
        self._main_view.ui.btn_update_protein_scene.setEnabled(False)

        self._main_view.ui.menuAnalysis.setEnabled(False)
        self._main_view.ui.menuResults.setEnabled(False)
        self._main_view.ui.menuImage.setEnabled(False)
        self._main_view.ui.menuHotspots.setEnabled(False)

    # </editor-fold>

    def manage_ui_of_protein_tab(
            self,
            an_object_type: str,
            is_protein_in_pair: bool,
            the_pymol_session_manager: "pymol_session_manager.PymolSessionManager"
    ) -> None:
        self._main_view.ui.lbl_info_2.hide()

        tmp_is_protein_in_session_flag: bool = the_pymol_session_manager.is_the_current_protein_in_session(self.get_current_active_protein_object().get_molecule_object())
        tmp_current_scene_name: str = the_pymol_session_manager.current_scene_name

        self._main_view.ui.btn_delete_protein.setEnabled(False)
        self._main_view.ui.btn_delete_protein_scene.setEnabled(False)

        if an_object_type == "protein":
            if not is_protein_in_pair:
                self._main_view.ui.btn_delete_protein.setEnabled(True)
            self._main_view.ui.btn_save_protein.setEnabled(True)
            self._main_view.ui.btn_open_protein_session.setEnabled(True)
            self.hide_protein_pymol_scene_configuration()

            if tmp_is_protein_in_session_flag and tmp_current_scene_name == "":
                self._main_view.ui.lbl_info.setText("Please select a scene.")
            elif tmp_is_protein_in_session_flag and tmp_current_scene_name != "":
                self._main_view.ui.lbl_info.setText("Please select a chain.")
            else:
                self._main_view.ui.lbl_info.setText("Please load the PyMOL session of the selected protein.")

        elif an_object_type == "scene":
            self._main_view.ui.btn_save_protein.setEnabled(False)
            self._main_view.ui.btn_open_protein_session.setEnabled(False)
            self.hide_protein_pymol_scene_configuration()

            if tmp_is_protein_in_session_flag and the_pymol_session_manager.current_scene_name == "base":
                self._main_view.ui.btn_delete_protein_scene.setEnabled(False)
                self._main_view.ui.lbl_info.setText("Please select a chain.")
            elif tmp_is_protein_in_session_flag and the_pymol_session_manager.current_scene_name != "base":
                self._main_view.ui.btn_delete_protein_scene.setEnabled(True)
                self._main_view.ui.lbl_info.setText("Please select a chain.")
            else:
                self._main_view.ui.lbl_info.setText("Please load the PyMOL session of the selected protein.")

        elif an_object_type == "chain":
            self._main_view.ui.btn_save_protein.setEnabled(False)
            self._main_view.ui.btn_open_protein_session.setEnabled(False)

            if tmp_is_protein_in_session_flag:
                self.show_protein_pymol_scene_configuration()
            else:
                self.hide_protein_pymol_scene_configuration()
            self.manage_coloring_by_element_option_for_protein_chain()

        elif an_object_type == "header":
            self._main_view.ui.btn_save_protein.setEnabled(False)
            self._main_view.ui.btn_open_protein_session.setEnabled(False)
            self.hide_protein_pymol_scene_configuration()

            if tmp_is_protein_in_session_flag and tmp_current_scene_name == "":
                self._main_view.ui.lbl_info.setText("Please select a scene.")
            elif tmp_is_protein_in_session_flag and tmp_current_scene_name != "":
                self._main_view.ui.lbl_info.setText("Please select a chain.")
            else:
                self._main_view.ui.lbl_info.setText("Please load the PyMOL session of the selected protein.")

        else:
            constants.PYSSA_LOGGER.warning("Unknown object type on proteins tab selected.")

        if tmp_is_protein_in_session_flag:
            self._main_view.ui.btn_create_protein_scene.setEnabled(True)
            self._main_view.ui.btn_update_protein_scene.setEnabled(True)
            self._main_view.ui.action_protein_regions.setEnabled(True)
        else:
            self._main_view.ui.btn_create_protein_scene.setEnabled(False)
            self._main_view.ui.btn_update_protein_scene.setEnabled(False)
            self._main_view.ui.action_protein_regions.setEnabled(False)

    def manage_coloring_by_element_option_for_protein_chain(self):
        if self.get_protein_repr_toggle_flag() == 1:
            if (self._main_view.tg_protein_sticks.toggle_button.isChecked()
                    or self._main_view.tg_protein_lines.toggle_button.isChecked()
                    or self._main_view.tg_protein_spheres.toggle_button.isChecked()
                    or self._main_view.tg_protein_dots.toggle_button.isChecked()
                    or self._main_view.tg_protein_mesh.toggle_button.isChecked()
                    or self._main_view.tg_protein_surface.toggle_button.isChecked()):
                # self._main_view.ui.btn_protein_color_atoms.setEnabled(True)
                # self._main_view.ui.btn_protein_reset_atoms.setEnabled(True)
                self._main_view.tg_protein_color_atoms.setEnabled(True)
            else:
                # self._main_view.ui.btn_protein_color_atoms.setEnabled(False)
                # self._main_view.ui.btn_protein_reset_atoms.setEnabled(False)
                self._main_view.tg_protein_color_atoms.setEnabled(False)
        else:
            if (self._main_view.ui.cb_protein_sticks.isChecked()
                    or self._main_view.ui.cb_protein_lines.isChecked()
                    or self._main_view.ui.cb_protein_spheres.isChecked()
                    or self._main_view.ui.cb_protein_dots.isChecked()
                    or self._main_view.ui.cb_protein_mesh.isChecked()
                    or self._main_view.ui.cb_protein_surface.isChecked()):
                self._main_view.ui.btn_protein_color_atoms.setEnabled(True)
                self._main_view.ui.btn_protein_reset_atoms.setEnabled(True)
            else:
                self._main_view.ui.btn_protein_color_atoms.setEnabled(False)
                self._main_view.ui.btn_protein_reset_atoms.setEnabled(False)

    def manage_hydrogen_representation_for_protein_chain(self):
        if (self._main_view.tg_protein_sticks.toggle_button.isChecked()
                    or self._main_view.tg_protein_lines.toggle_button.isChecked()
                    or self._main_view.tg_protein_spheres.toggle_button.isChecked()):
            # self._main_view.ui.btn_protein_show_hydrogens.setEnabled(True)
            # self._main_view.ui.btn_protein_hide_hydrogens.setEnabled(True)
            self._main_view.tg_protein_hydrogen_atoms.setEnabled(True)
        else:
            # self._main_view.ui.btn_protein_show_hydrogens.setEnabled(False)
            # self._main_view.ui.btn_protein_hide_hydrogens.setEnabled(False)
            self._main_view.tg_protein_hydrogen_atoms.setEnabled(False)

    def manage_toggle_state_of_protein_repr(self, tmp_repr_state):
        if tmp_repr_state[enums.PyMOLRepresentation.CARTOON.value] == 0:
            self._main_view.tg_protein_cartoon.toggle_button.setChecked(False)
        else:
            self._main_view.tg_protein_cartoon.toggle_button.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.STICKS.value] == 0:
            self._main_view.tg_protein_sticks.toggle_button.setChecked(False)
        else:
            self._main_view.tg_protein_sticks.toggle_button.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.RIBBON.value] == 0:
            self._main_view.tg_protein_ribbon.toggle_button.setChecked(False)
        else:
            self._main_view.tg_protein_ribbon.toggle_button.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.LINES.value] == 0:
            self._main_view.tg_protein_lines.toggle_button.setChecked(False)
        else:
            self._main_view.tg_protein_lines.toggle_button.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.SPHERES.value] == 0:
            self._main_view.tg_protein_spheres.toggle_button.setChecked(False)
        else:
            self._main_view.tg_protein_spheres.toggle_button.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.DOTS.value] == 0:
            self._main_view.tg_protein_dots.toggle_button.setChecked(False)
        else:
            self._main_view.tg_protein_dots.toggle_button.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.MESH.value] == 0:
            self._main_view.tg_protein_mesh.toggle_button.setChecked(False)
        else:
            self._main_view.tg_protein_mesh.toggle_button.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.SURFACE.value] == 0:
            self._main_view.tg_protein_surface.toggle_button.setChecked(False)
        else:
            self._main_view.tg_protein_surface.toggle_button.setChecked(True)

    def manage_check_state_of_protein_repr(self, tmp_repr_state):
        if tmp_repr_state[enums.PyMOLRepresentation.CARTOON.value] == 0:
            self._main_view.ui.cb_protein_cartoon.setChecked(False)
        else:
            self._main_view.ui.cb_protein_cartoon.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.STICKS.value] == 0:
            self._main_view.ui.cb_protein_sticks.setChecked(False)
        else:
            self._main_view.ui.cb_protein_sticks.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.RIBBON.value] == 0:
            self._main_view.ui.cb_protein_ribbon.setChecked(False)
        else:
            self._main_view.ui.cb_protein_ribbon.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.LINES.value] == 0:
            self._main_view.ui.cb_protein_lines.setChecked(False)
        else:
            self._main_view.ui.cb_protein_lines.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.SPHERES.value] == 0:
            self._main_view.ui.cb_protein_spheres.setChecked(False)
        else:
            self._main_view.ui.cb_protein_spheres.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.DOTS.value] == 0:
            self._main_view.ui.cb_protein_dots.setChecked(False)
        else:
            self._main_view.ui.cb_protein_dots.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.MESH.value] == 0:
            self._main_view.ui.cb_protein_mesh.setChecked(False)
        else:
            self._main_view.ui.cb_protein_mesh.setChecked(True)
        if tmp_repr_state[enums.PyMOLRepresentation.SURFACE.value] == 0:
            self._main_view.ui.cb_protein_surface.setChecked(False)
        else:
            self._main_view.ui.cb_protein_surface.setChecked(True)

    # <editor-fold desc="Scene">
    def add_scene_to_proteins_model(
            self,
            a_scene_name
    ):
        tmp_scene_item = QtGui.QStandardItem(a_scene_name)
        tmp_scene_item.setData("scene", enums.ModelEnum.TYPE_ROLE)
        self._protein_model.add_scene(
            self.get_current_protein_tree_index(),
            tmp_scene_item
        )

    def show_protein_pymol_scene_configuration(self):
        self._main_view.ui.frame_protein_color.show()
        self._main_view.ui.frame_protein_repr.show()
        if self._settings_manager.settings.proteins_tab_use_combobox_for_colors == 1:
            self._main_view.ui.box_protein_color.show()
            self._main_view.ui.lbl_protein_current_color.hide()
            self._main_view.ui.lbl_protein_pymol_colors.hide()
            self._main_view.color_grid_proteins.hide()
        else:
            self._main_view.ui.box_protein_color.hide()
            self._main_view.ui.lbl_protein_current_color.show()
            self._main_view.ui.lbl_protein_pymol_colors.show()
            self._main_view.color_grid_proteins.show()

        if self._settings_manager.settings.proteins_tab_use_toggle == 1:
            self._main_view.ui.verticalLayout_15.setSpacing(0)  # layout for the representation section
            """IMPORTANT: 
                The layout for the hide all repr frame must have a top and bottom margin of 6, 
                set in the QDesigner's settings
            """
            # toggles should be used
            self._main_view.ui.lbl_protein_current_color.show()
            self._main_view.ui.lbl_protein_atoms.show()
            self._main_view.ui.lbl_protein_cartoon.show()
            self._main_view.ui.lbl_protein_sticks.show()
            self._main_view.ui.lbl_protein_ribbon.show()
            self._main_view.ui.lbl_protein_lines.show()
            self._main_view.ui.lbl_protein_spheres.show()
            self._main_view.ui.lbl_protein_dots.show()
            self._main_view.ui.lbl_protein_mesh.show()
            self._main_view.ui.lbl_protein_surface.show()

            self._main_view.tg_protein_cartoon.show()
            self._main_view.tg_protein_sticks.show()
            self._main_view.tg_protein_ribbon.show()
            self._main_view.tg_protein_lines.show()
            self._main_view.tg_protein_spheres.show()
            self._main_view.tg_protein_dots.show()
            self._main_view.tg_protein_mesh.show()
            self._main_view.tg_protein_surface.show()
            # hide ui elements from checkbox options
            self._main_view.ui.cb_protein_cartoon.hide()
            self._main_view.ui.cb_protein_sticks.hide()
            self._main_view.ui.cb_protein_ribbon.hide()
            self._main_view.ui.cb_protein_lines.hide()
            self._main_view.ui.cb_protein_spheres.hide()
            self._main_view.ui.cb_protein_dots.hide()
            self._main_view.ui.cb_protein_mesh.hide()
            self._main_view.ui.cb_protein_surface.hide()
        else:
            self._main_view.ui.cb_protein_cartoon.show()
            self._main_view.ui.cb_protein_sticks.show()
            self._main_view.ui.cb_protein_ribbon.show()
            self._main_view.ui.cb_protein_lines.show()
            self._main_view.ui.cb_protein_spheres.show()
            self._main_view.ui.cb_protein_dots.show()
            self._main_view.ui.cb_protein_mesh.show()
            self._main_view.ui.cb_protein_surface.show()
            # hide ui elements from toggle options
            self._main_view.ui.lbl_protein_atoms.show()
            self._main_view.ui.lbl_protein_cartoon.hide()
            self._main_view.ui.lbl_protein_sticks.hide()
            self._main_view.ui.lbl_protein_ribbon.hide()
            self._main_view.ui.lbl_protein_lines.hide()
            self._main_view.ui.lbl_protein_spheres.hide()
            self._main_view.ui.lbl_protein_dots.hide()
            self._main_view.ui.lbl_protein_mesh.hide()
            self._main_view.ui.lbl_protein_surface.hide()
            self._main_view.tg_protein_cartoon.hide()
            self._main_view.tg_protein_sticks.hide()
            self._main_view.tg_protein_ribbon.hide()
            self._main_view.tg_protein_lines.hide()
            self._main_view.tg_protein_spheres.hide()
            self._main_view.tg_protein_dots.hide()
            self._main_view.tg_protein_mesh.hide()
            self._main_view.tg_protein_surface.hide()

        self._main_view.ui.lbl_protein_color.show()
        self._main_view.ui.lbl_protein_all_representations.show()
        self._main_view.ui.btn_protein_hide_all_representations.show()
        self._main_view.ui.lbl_info.hide()
        self._main_view.ui.lbl_info_2.hide()

        # self._main_view.ui.lbl_protein_color.show()
        # self._main_view.ui.lbl_protein_atoms.show()
        # self._main_view.ui.lbl_protein_cartoon.show()
        # self._main_view.ui.lbl_protein_sticks.show()
        # self._main_view.ui.lbl_protein_ribbon.show()
        # self._main_view.ui.lbl_protein_lines.show()
        # self._main_view.ui.lbl_protein_spheres.show()
        # self._main_view.ui.lbl_protein_dots.show()
        # self._main_view.ui.lbl_protein_mesh.show()
        # self._main_view.ui.lbl_protein_surface.show()
        #
        # self._main_view.ui.lbl_protein_all_representations.show()
        # #self._main_view.ui.box_protein_color.show()
        # #self._main_view.ui.btn_protein_color_atoms.show()
        # #self._main_view.ui.btn_protein_reset_atoms.show()
        # # self._main_view.ui.btn_protein_show_cartoon.show()
        # # self._main_view.ui.btn_protein_hide_cartoon.show()
        # # self._main_view.ui.btn_protein_show_sticks.show()
        # # self._main_view.ui.btn_protein_hide_sticks.show()
        # # self._main_view.ui.btn_protein_show_ribbon.show()
        # # self._main_view.ui.btn_protein_hide_ribbon.show()
        # # self._main_view.ui.cb_protein_cartoon.show()
        # # self._main_view.ui.cb_protein_sticks.show()
        # # self._main_view.ui.cb_protein_ribbon.show()
        # # self._main_view.ui.cb_protein_lines.show()
        # # self._main_view.ui.cb_protein_spheres.show()
        # # self._main_view.ui.cb_protein_dots.show()
        # # self._main_view.ui.cb_protein_mesh.show()
        # # self._main_view.ui.cb_protein_surface.show()
        # self._main_view.tg_protein_cartoon.show()
        # self._main_view.tg_protein_sticks.show()
        # self._main_view.tg_protein_ribbon.show()
        # self._main_view.tg_protein_lines.show()
        # self._main_view.tg_protein_spheres.show()
        # self._main_view.tg_protein_dots.show()
        # self._main_view.tg_protein_mesh.show()
        # self._main_view.tg_protein_surface.show()
        # self._main_view.ui.btn_protein_hide_all_representations.show()
        # self._main_view.ui.lbl_info.hide()
        # self._main_view.ui.lbl_info_2.hide()

    def hide_protein_pymol_scene_configuration(self):
        self._main_view.ui.frame_protein_color.hide()
        self._main_view.ui.frame_protein_repr.hide()
        if self._settings_manager.settings.proteins_tab_use_toggle == 1:
            # toggles should be used
            self._main_view.ui.lbl_protein_atoms.hide()
            self._main_view.ui.lbl_protein_cartoon.hide()
            self._main_view.ui.lbl_protein_sticks.hide()
            self._main_view.ui.lbl_protein_ribbon.hide()
            self._main_view.ui.lbl_protein_lines.hide()
            self._main_view.ui.lbl_protein_spheres.hide()
            self._main_view.ui.lbl_protein_dots.hide()
            self._main_view.ui.lbl_protein_mesh.hide()
            self._main_view.ui.lbl_protein_surface.hide()
            self._main_view.tg_protein_cartoon.hide()
            self._main_view.tg_protein_sticks.hide()
            self._main_view.tg_protein_ribbon.hide()
            self._main_view.tg_protein_lines.hide()
            self._main_view.tg_protein_spheres.hide()
            self._main_view.tg_protein_dots.hide()
            self._main_view.tg_protein_mesh.hide()
            self._main_view.tg_protein_surface.hide()
        else:
            self._main_view.ui.cb_protein_cartoon.hide()
            self._main_view.ui.cb_protein_sticks.hide()
            self._main_view.ui.cb_protein_ribbon.hide()
            self._main_view.ui.cb_protein_lines.hide()
            self._main_view.ui.cb_protein_spheres.hide()
            self._main_view.ui.cb_protein_dots.hide()
            self._main_view.ui.cb_protein_mesh.hide()
            self._main_view.ui.cb_protein_surface.hide()

        self._main_view.ui.lbl_protein_color.hide()
        self._main_view.ui.box_protein_color.hide()
        self._main_view.ui.lbl_protein_all_representations.hide()
        self._main_view.ui.btn_protein_hide_all_representations.hide()
        # self._main_view.ui.btn_protein_show_cartoon.hide()
        # self._main_view.ui.btn_protein_hide_cartoon.hide()
        # self._main_view.ui.btn_protein_show_sticks.hide()
        # self._main_view.ui.btn_protein_hide_sticks.hide()
        # self._main_view.ui.btn_protein_show_ribbon.hide()
        # self._main_view.ui.btn_protein_hide_ribbon.hide()
        self._main_view.ui.btn_protein_color_atoms.hide()
        self._main_view.ui.btn_protein_reset_atoms.hide()
        self._main_view.ui.lbl_info.show()

    def remove_scene_from_proteins_model(self, the_model_index_of_the_scene: QtCore.QModelIndex):
        self._protein_model.remove_scene(the_model_index_of_the_scene)
    # </editor-fold>

    # </editor-fold>

    # <editor-fold desc="Protein Pairs">
    def add_protein_pair_to_protein_pairs_model(self, a_protein_pair: "protein_pair.ProteinPair"):
        tmp_main_socket, the_general_purpose_socket = self.job_manager.get_general_purpose_socket_pair()
        self._protein_pair_model.add_protein_pair(a_protein_pair, tmp_main_socket, the_general_purpose_socket)

    def remove_protein_pair_from_protein_pairs_model(self):
        self._protein_pair_model.remove_protein_pair(self.get_current_protein_pair_tree_index())

    # <editor-fold desc="Menu Options">
    def show_menu_options_with_protein_pair(self):
        self._main_view.ui.protein_pairs_tree_view.setModel(self._protein_pair_model)
        self._main_view.ui.protein_pairs_tree_view.setHeaderHidden(True)

    def show_menu_options_without_protein_pair(self):
        self._main_view.ui.btn_delete_protein_pair.setEnabled(False)
        self._main_view.ui.btn_open_protein_pair_session.setEnabled(False)
        self._main_view.ui.btn_create_protein_pair_scene.setEnabled(False)
        self._main_view.ui.btn_update_protein_pair_scene.setEnabled(False)
    # </editor-fold>

    def manage_ui_of_protein_pairs_tab(self,
                                       an_object_type: str,
                                       the_pymol_session_manager: pymol_session_manager.PymolSessionManager):
        self._main_view.ui.lbl_info_4.hide()

        tmp_is_protein_pair_in_session_flag: bool = the_pymol_session_manager.is_the_current_protein_pair_in_session(self.get_current_active_protein_pair_object().name)
        tmp_current_scene_name: str = the_pymol_session_manager.current_scene_name
        self._main_view.ui.btn_delete_protein_pair_scene.setEnabled(False)

        if an_object_type == "protein_pair":
            self._main_view.ui.btn_delete_protein_pair.setEnabled(True)
            self._main_view.ui.btn_open_protein_pair_session.setEnabled(True)
            self.hide_protein_pair_pymol_scene_configuration()

            if tmp_is_protein_pair_in_session_flag and tmp_current_scene_name == "":
                self._main_view.ui.lbl_info_3.setText("Please select a scene.")
            elif tmp_is_protein_pair_in_session_flag and tmp_current_scene_name != "":
                self._main_view.ui.lbl_info_3.setText("Please select a chain.")
            else:
                self._main_view.ui.lbl_info_3.setText("Please load the PyMOL session of the\nselected protein pair.")

        elif an_object_type == "protein":
            self._main_view.ui.btn_delete_protein_pair.setEnabled(False)
            self._main_view.ui.btn_open_protein_pair_session.setEnabled(False)
            self.hide_protein_pair_pymol_scene_configuration()

            if tmp_is_protein_pair_in_session_flag and tmp_current_scene_name == "":
                self._main_view.ui.lbl_info_3.setText("Please select a scene.")
            elif tmp_is_protein_pair_in_session_flag and tmp_current_scene_name != "":
                self._main_view.ui.lbl_info_3.setText("Please select a chain.")
            else:
                self._main_view.ui.lbl_info_3.setText("Please load the PyMOL session of the\nselected protein pair.")

        elif an_object_type == "scene":
            self._main_view.ui.btn_delete_protein_pair.setEnabled(False)
            self._main_view.ui.btn_open_protein_pair_session.setEnabled(False)
            self.hide_protein_pair_pymol_scene_configuration()

            if tmp_is_protein_pair_in_session_flag and tmp_current_scene_name == "base":
                self._main_view.ui.btn_delete_protein_pair_scene.setEnabled(False)
                self._main_view.ui.lbl_info_3.setText("Please select a chain.")
            elif tmp_is_protein_pair_in_session_flag and tmp_current_scene_name != "base":
                self._main_view.ui.btn_delete_protein_pair_scene.setEnabled(True)
                self._main_view.ui.lbl_info_3.setText("Please select a chain.")
            else:
                self._main_view.ui.lbl_info_3.setText("Please load the PyMOL session of the \nselected protein pair.")

            if tmp_is_protein_pair_in_session_flag:
                self._main_view.ui.lbl_info_3.setText("Please select a chain.")
            else:
                self._main_view.ui.lbl_info_3.setText("Please load the PyMOL session of the \nselected protein pair.")

        elif an_object_type == "chain":
            self._main_view.ui.btn_delete_protein_pair.setEnabled(False)
            self._main_view.ui.btn_open_protein_pair_session.setEnabled(False)
            if tmp_is_protein_pair_in_session_flag:
                self.show_protein_pair_pymol_scene_configuration()
            else:
                self.hide_protein_pair_pymol_scene_configuration()

        elif an_object_type == "header":
            self._main_view.ui.btn_delete_protein_pair.setEnabled(False)
            self._main_view.ui.btn_open_protein_pair_session.setEnabled(False)
            self.hide_protein_pair_pymol_scene_configuration()
            if tmp_is_protein_pair_in_session_flag and tmp_current_scene_name == "":
                self._main_view.ui.lbl_info_3.setText("Please select a scene.")
            elif tmp_is_protein_pair_in_session_flag and tmp_current_scene_name != "":
                self._main_view.ui.lbl_info_3.setText("Please select a chain.")
            else:
                self._main_view.ui.lbl_info_3.setText("Please load the PyMOL session of the \nselected protein pair.")
        else:
            constants.PYSSA_LOGGER.warning("Unknown object type on proteins tab selected.")

        # if tmp_is_protein_pair_in_session_flag and tmp_current_scene_name == "base":
        #     self._main_view.ui.btn_delete_protein_scene.setEnabled(False)
        #     self._main_view.ui.lbl_info.setText("Please select a chain.")
        # elif tmp_is_protein_pair_in_session_flag and tmp_current_scene_name != "base":
        #     self._main_view.ui.btn_delete_protein_pair_scene.setEnabled(True)
        #     self._main_view.ui.lbl_info_3.setText("Please select a chain.")
        # else:
        #     self._main_view.ui.lbl_info_3.setText("Please load the PyMOL session of the selected protein.")

        if tmp_is_protein_pair_in_session_flag:
            self._main_view.ui.btn_create_protein_pair_scene.setEnabled(True)
            self._main_view.ui.btn_update_protein_pair_scene.setEnabled(True)
            self._main_view.ui.action_protein_regions.setEnabled(True)
        else:
            self._main_view.ui.btn_create_protein_pair_scene.setEnabled(False)
            self._main_view.ui.btn_update_protein_pair_scene.setEnabled(False)
            self._main_view.ui.btn_delete_protein_pair_scene.setEnabled(False)
            self._main_view.ui.action_protein_regions.setEnabled(False)

    def manage_coloring_by_element_option_for_protein_pair_chain(self):
        if self.get_protein_pair_repr_toggle_flag() == 1:
            if (self._main_view.tg_protein_pair_sticks.toggle_button.isChecked()
                    or self._main_view.tg_protein_pair_lines.toggle_button.isChecked()
                    or self._main_view.tg_protein_pair_spheres.toggle_button.isChecked()
                    or self._main_view.tg_protein_pair_dots.toggle_button.isChecked()
                    or self._main_view.tg_protein_pair_mesh.toggle_button.isChecked()
                    or self._main_view.tg_protein_pair_surface.toggle_button.isChecked()):
                # self._main_view.ui.btn_protein_pair_color_atoms.setEnabled(True)
                # self._main_view.ui.btn_protein_pair_reset_atoms.setEnabled(True)
                self._main_view.tg_protein_pair_color_atoms.setEnabled(True)
            else:
                # self._main_view.ui.btn_protein_pair_color_atoms.setEnabled(False)
                # self._main_view.ui.btn_protein_pair_reset_atoms.setEnabled(False)
                self._main_view.tg_protein_pair_color_atoms.setEnabled(False)
        else:
            if (self._main_view.ui.cb_protein_pair_sticks.isChecked()
                    or self._main_view.ui.cb_protein_pair_lines.isChecked()
                    or self._main_view.ui.cb_protein_pair_spheres.isChecked()
                    or self._main_view.ui.cb_protein_pair_dots.isChecked()
                    or self._main_view.ui.cb_protein_pair_mesh.isChecked()
                    or self._main_view.ui.cb_protein_pair_surface.isChecked()):
                self._main_view.ui.btn_protein_pair_color_atoms.setEnabled(True)
                self._main_view.ui.btn_protein_pair_reset_atoms.setEnabled(True)
            else:
                self._main_view.ui.btn_protein_pair_color_atoms.setEnabled(False)
                self._main_view.ui.btn_protein_pair_reset_atoms.setEnabled(False)

    def manage_hydrogen_representation_for_protein_pair_chain(self):
        if (self._main_view.tg_protein_pair_sticks.toggle_button.isChecked()
                    or self._main_view.tg_protein_pair_lines.toggle_button.isChecked()
                    or self._main_view.tg_protein_pair_spheres.toggle_button.isChecked()):
            # self._main_view.ui.btn_protein_pair_show_hydrogens.setEnabled(True)
            # self._main_view.ui.btn_protein_pair_hide_hydrogens.setEnabled(True)
            self._main_view.tg_protein_pair_hydrogen_atoms.setEnabled(True)
        else:
            # self._main_view.ui.btn_protein_pair_show_hydrogens.setEnabled(False)
            # self._main_view.ui.btn_protein_pair_hide_hydrogens.setEnabled(False)
            self._main_view.tg_protein_pair_hydrogen_atoms.setEnabled(False)

    def set_current_chain_color_for_ui_for_protein_pairs(
            self,
            a_protein_name: str,
            the_pymol_session_manager: "pymol_session_manager.PymolSessionManager"
    ):
        tmp_protein = self.get_current_active_protein_object_of_protein_pair()
        tmp_chain = self.get_current_active_chain_object_of_protein_pair()
        if the_pymol_session_manager.is_the_current_protein_pair_in_session(self.get_current_active_protein_pair_object().name):
            if tmp_chain.chain_type == "protein_chain":
                tmp_protein.pymol_selection.selection_string = f"first chain {tmp_chain.chain_letter} and name CA and {a_protein_name}"
            else:
                tmp_protein.pymol_selection.selection_string = f"first chain {tmp_chain.chain_letter}"
            if self.pymol_session_manager.get_residue_color_config(f"(first elem N, first elem C, first elem O) and chain {tmp_chain.chain_letter} and {a_protein_name}",
                                                                   tmp_chain.chain_letter).atoms_are_colored_by_elements():
                self._main_view.ui.lbl_protein_pair_current_color.setText("By Element    ")
                self._main_view.tg_protein_pair_color_atoms.toggle_button.setChecked(True)
            else:
                rvoid(tmp_chain.get_color(tmp_protein.pymol_selection.selection_string, self.pymol_session_manager))
                self._main_view.ui.lbl_protein_pair_current_color.setText(f"{tmp_chain.pymol_parameters['chain_color']}    ")
                self._main_view.tg_protein_pair_color_atoms.toggle_button.setChecked(False)



        #     # fixme: This can easily be bypassed by a power user if the first residue color is changed
        #     if tmp_chain.chain_type == "protein_chain":
        #         tmp_protein.pymol_selection.selection_string = f"first chain {tmp_chain.chain_letter} and {a_protein_name} and name CA"
        #     else:
        #         tmp_protein.pymol_selection.selection_string = f"first chain {tmp_chain.chain_letter} and {a_protein_name}"
        #     rvoid(tmp_chain.get_color(tmp_protein.pymol_selection.selection_string, self.pymol_session_manager))
        # self._main_view.ui.lbl_protein_pair_current_color.setText(tmp_chain.pymol_parameters["chain_color"])

    # <editor-fold desc="Scene">
    def add_scene_to_protein_pairs_model(
            self,
            a_scene_name
    ):
        tmp_scene_item = QtGui.QStandardItem(a_scene_name)
        tmp_scene_item.setData("scene", enums.ModelEnum.TYPE_ROLE)
        self._protein_pair_model.add_scene(
            self.get_current_protein_pair_tree_index(),
            tmp_scene_item
        )

    def show_protein_pair_pymol_scene_configuration(self):
        self._main_view.ui.frame_protein_pair_color.show()
        self._main_view.ui.frame_protein_pair_repr.show()
        if self._settings_manager.settings.protein_pairs_tab_use_combobox_for_colors == 1:
            self._main_view.ui.box_protein_pair_color.show()
            self._main_view.ui.lbl_protein_pair_current_color.hide()
            self._main_view.ui.lbl_protein_pair_pymol_colors.hide()
            self._main_view.color_grid_protein_pairs.hide()
        else:
            self._main_view.ui.box_protein_pair_color.hide()
            self._main_view.ui.lbl_protein_pair_current_color.show()
            self._main_view.ui.lbl_protein_pair_pymol_colors.show()
            self._main_view.color_grid_protein_pairs.show()

        if self._settings_manager.settings.protein_pairs_tab_use_toggle == 1:
            self._main_view.ui.verticalLayout_18.setSpacing(0)  # layout for the representation section
            """IMPORTANT: 
                The layout for the hide all repr frame must have a top and bottom margin of 6, 
                set in the QDesigner's settings
            """
            # toggles should be used
            self._main_view.ui.lbl_protein_pair_atoms.show()
            self._main_view.ui.lbl_protein_pair_cartoon.show()
            self._main_view.ui.lbl_protein_pair_sticks.show()
            self._main_view.ui.lbl_protein_pair_ribbon.show()
            self._main_view.ui.lbl_protein_pair_lines.show()
            self._main_view.ui.lbl_protein_pair_spheres.show()
            self._main_view.ui.lbl_protein_pair_dots.show()
            self._main_view.ui.lbl_protein_pair_mesh.show()
            self._main_view.ui.lbl_protein_pair_surface.show()

            self._main_view.tg_protein_pair_cartoon.show()
            self._main_view.tg_protein_pair_sticks.show()
            self._main_view.tg_protein_pair_ribbon.show()
            self._main_view.tg_protein_pair_lines.show()
            self._main_view.tg_protein_pair_spheres.show()
            self._main_view.tg_protein_pair_dots.show()
            self._main_view.tg_protein_pair_mesh.show()
            self._main_view.tg_protein_pair_surface.show()
            # hide ui elements from checkbox options
            self._main_view.ui.cb_protein_pair_cartoon.hide()
            self._main_view.ui.cb_protein_pair_sticks.hide()
            self._main_view.ui.cb_protein_pair_ribbon.hide()
            self._main_view.ui.cb_protein_pair_lines.hide()
            self._main_view.ui.cb_protein_pair_spheres.hide()
            self._main_view.ui.cb_protein_pair_dots.hide()
            self._main_view.ui.cb_protein_pair_mesh.hide()
            self._main_view.ui.cb_protein_pair_surface.hide()
        else:
            self._main_view.ui.lbl_protein_pair_atoms.show()
            self._main_view.ui.cb_protein_pair_cartoon.show()
            self._main_view.ui.cb_protein_pair_sticks.show()
            self._main_view.ui.cb_protein_pair_ribbon.show()
            self._main_view.ui.cb_protein_pair_lines.show()
            self._main_view.ui.cb_protein_pair_spheres.show()
            self._main_view.ui.cb_protein_pair_dots.show()
            self._main_view.ui.cb_protein_pair_mesh.show()
            self._main_view.ui.cb_protein_pair_surface.show()
            # hide ui elements from toggle options
            self._main_view.ui.lbl_protein_pair_cartoon.hide()
            self._main_view.ui.lbl_protein_pair_sticks.hide()
            self._main_view.ui.lbl_protein_pair_ribbon.hide()
            self._main_view.ui.lbl_protein_pair_lines.hide()
            self._main_view.ui.lbl_protein_pair_spheres.hide()
            self._main_view.ui.lbl_protein_pair_dots.hide()
            self._main_view.ui.lbl_protein_pair_mesh.hide()
            self._main_view.ui.lbl_protein_pair_surface.hide()

            self._main_view.tg_protein_pair_color_atoms.hide()
            self._main_view.tg_protein_pair_cartoon.hide()
            self._main_view.tg_protein_pair_sticks.hide()
            self._main_view.tg_protein_pair_ribbon.hide()
            self._main_view.tg_protein_pair_lines.hide()
            self._main_view.tg_protein_pair_spheres.hide()
            self._main_view.tg_protein_pair_dots.hide()
            self._main_view.tg_protein_pair_mesh.hide()
            self._main_view.tg_protein_pair_surface.hide()

        self._main_view.ui.lbl_protein_pair_color.show()
        self._main_view.ui.lbl_protein_pair_all_representations.show()
        self._main_view.ui.btn_protein_pair_hide_all_representations.show()
        self._main_view.ui.lbl_info_3.hide()
        self._main_view.ui.lbl_info_4.hide()

    def hide_protein_pair_pymol_scene_configuration(self):
        self._main_view.ui.frame_protein_pair_color.hide()
        self._main_view.ui.frame_protein_pair_repr.hide()
        self._main_view.ui.lbl_info_3.show()
        # self._main_view.ui.lbl_protein_pair_color.hide()
        # self._main_view.ui.lbl_protein_pair_atoms.hide()
        # # self._main_view.ui.lbl_protein_pair_cartoon.hide()
        # # self._main_view.ui.lbl_protein_pair_sticks.hide()
        # # self._main_view.ui.lbl_protein_pair_ribbon.hide()
        # self._main_view.ui.lbl_protein_pair_all_representations.hide()
        # self._main_view.ui.box_protein_pair_color.hide()
        # self._main_view.ui.btn_protein_pair_color_atoms.hide()
        # self._main_view.ui.btn_protein_pair_reset_atoms.hide()
        # # self._main_view.ui.btn_protein_pair_show_cartoon.hide()
        # # self._main_view.ui.btn_protein_pair_hide_cartoon.hide()
        # # self._main_view.ui.btn_protein_pair_show_sticks.hide()
        # # self._main_view.ui.btn_protein_pair_hide_sticks.hide()
        # # self._main_view.ui.btn_protein_pair_show_ribbon.hide()
        # # self._main_view.ui.btn_protein_pair_hide_ribbon.hide()
        # self._main_view.ui.cb_protein_pair_cartoon.hide()
        # self._main_view.ui.cb_protein_pair_sticks.hide()
        # self._main_view.ui.cb_protein_pair_ribbon.hide()
        # self._main_view.ui.cb_protein_pair_lines.hide()
        # self._main_view.ui.cb_protein_pair_spheres.hide()
        # self._main_view.ui.cb_protein_pair_dots.hide()
        # self._main_view.ui.cb_protein_pair_mesh.hide()
        # self._main_view.ui.cb_protein_pair_surface.hide()
        # # toggles
        # self._main_view.tg_protein_pair_cartoon.hide()
        # self._main_view.tg_protein_pair_sticks.hide()
        # self._main_view.tg_protein_pair_ribbon.hide()
        # self._main_view.tg_protein_pair_lines.hide()
        # self._main_view.tg_protein_pair_spheres.hide()
        # self._main_view.tg_protein_pair_dots.hide()
        # self._main_view.tg_protein_pair_mesh.hide()
        # self._main_view.tg_protein_pair_surface.hide()
        #
        # self._main_view.ui.btn_protein_pair_hide_all_representations.hide()
        #
        # self._main_view.ui.lbl_info_3.show()

    def remove_scene_from_protein_pairs_model(self, the_model_index_of_the_scene: QtCore.QModelIndex):
        self._protein_pair_model.remove_scene(the_model_index_of_the_scene)
    # </editor-fold>

    # </editor-fold>

    def update_settings(self):
        """Deserializes the settings json file."""
        self._settings_manager.settings = self._settings_manager.settings.deserialize_settings()

    def restore_default_main_view(self):
        # Restore sequences table
        logger.info("Restoring default main view at seq table")
        self._main_view.ui.seqs_table_widget.clear()
        logger.info("seq table cleared")
        self._main_view.ui.seqs_table_widget.setRowCount(0)
        logger.info("seq table set row count to 0")
        # # Restore proteins table
        # logger.info("Restoring default main view at protein table")
        # self._main_view.ui.proteins_table_widget.clear()
        # logger.info("protein table cleared")
        # self._main_view.ui.proteins_table_widget.setRowCount(0)
        # logger.info("protein table set row count to 0")
        # # Restore protein pairs table
        # logger.info("Restoring default main view at protein pair table")
        # self._main_view.ui.protein_pairs_table_widget.clear()
        # logger.info("protein pair table cleared")
        # self._main_view.ui.protein_pairs_table_widget.setRowCount(0)
        # logger.info("protein pair table set row count to 0")
        # Initialize UI
        self._main_view.initialize_ui()

    def start_wait_cursor(self) -> None:
        """Starts the wait cursor."""
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        self._main_view.disable_menu_bar()
        self._main_view.disable_tab_widget()
        self._main_view.disable_job_panels()

    def stop_wait_cursor(self) -> None:
        """Stops the spinner."""
        QtWidgets.QApplication.restoreOverrideCursor()
        self._main_view.enable_job_panels()

    # <editor-fold desc="Job related methods">
    def add_job_entry_to_job_overview_layout(self, a_job_entry_widget):
        self.job_entry_widgets.append(a_job_entry_widget)
        self._main_view.ui.job_overview_layout.insertWidget(self._main_view.ui.job_overview_layout.count() - 1,
                                                            a_job_entry_widget)
        self._main_view.lbl_job_overview.hide()
        self._main_view.btn_open_job_overview.setIcon(self._main_view.icon_jobs_running)

    def update_job_entry(self, update_job_entry_signal_values):
        """

        Notes:
            Gets signal from the classes of job.py with the signal "update_job_entry_signal"

        """
        _, a_description, a_value = update_job_entry_signal_values
        a_job_entry_widget: "job_entry.JobEntryWidget" = update_job_entry_signal_values[0]  # Unpack of 0 because of type annotation need
        a_job_entry_widget.ui.progress_bar_job.setValue(a_value)
        if a_value == 100:
            tmp_type = a_job_entry_widget.job_base_information.job_type
            tmp_progress = a_job_entry_widget.job_base_information.job_progress
            if tmp_type == enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS and tmp_progress == enums.JobProgress.RUNNING:
                return
            self._create_job_notification(a_job_entry_widget, a_description)
        elif a_value == 0 or a_value == 25:
            a_job_entry_widget.ui.btn_cancel_job.setEnabled(True)
        else:
            a_job_entry_widget.ui.btn_cancel_job.setEnabled(False)

    def _create_job_notification(self, a_job_entry_widget: "job_entry.JobEntryWidget", a_description):
        """Helper function for setting up a job notification after a job finished."""
        if a_job_entry_widget.job_base_information.project_name == self._current_project.get_project_name():
            tmp_job_is_from_current_project = True
        else:
            tmp_job_is_from_current_project = False

        tmp_job_notification_widget = job_entry.JobNotificationWidget(a_description,
                                                                      a_job_entry_widget.job_base_information,
                                                                      tmp_job_is_from_current_project,
                                                                      self.refresh_after_job_finished_signal)
        # remove job entry widget
        if a_job_entry_widget.parent():
            a_job_entry_widget.setParent(None)
        a_job_entry_widget.deleteLater()
        # Add new notification widget to notification panel
        self._main_view.ui.job_notification_layout.insertWidget(self._main_view.ui.job_notification_layout.count() - 1,
                                                                tmp_job_notification_widget)
        self._main_view.btn_open_job_notification.setIcon(self._main_view.icon_notify_unread)
        self._main_view.lbl_job_notification.hide()
        if self._main_view.ui.job_overview_layout.count() == 2:
            self._main_view.lbl_job_overview.show()
            self._main_view.btn_open_job_overview.setIcon(self._main_view.icon_jobs)

    def remove_job_notification_widget(self, a_job_notification_widget):
        """Removes a job notification after pressing open or refresh.

        Notes:
            This function is called from the main_view_controller.
        """
        print(a_job_notification_widget)
        self._main_view.ui.job_notification_layout.removeWidget(a_job_notification_widget)
        if self._main_view.ui.job_notification_layout.count() == 2:
            self._main_view.lbl_job_notification.show()
            self._main_view.btn_open_job_notification.setIcon(self._main_view.icon_notify)

    def close_job_notification_panel(self):
        """Closes the job notification panel after pressing open or refresh.

        Notes:
            This function is called from the main_view_controller.
        """
        self._main_view.ui.frame_job_notification.hide()

    def close_job_overview_panel(self):
        self._main_view.ui.frame_job_overview.hide()

    def cancel_job(self, signal_tuple):
        if signal_tuple[0] == enums.JobType.PREDICTION:
            # signal_tuple = (job_type, job_entry_widget, job object, protein_prediction_infos)
            for tmp_protein_info in signal_tuple[3]:
                self.watcher.remove_protein(tmp_protein_info.name)
            if signal_tuple[1].parent():
                signal_tuple[1].setParent(None)
            signal_tuple[1].deleteLater()
            if self.job_manager.is_prediction_job_currently_running(signal_tuple[2]):
                constants.PYSSA_LOGGER.info("Structure prediction process was aborted manually.")
                subprocess.run(["wsl", "--shutdown"])
                constants.PYSSA_LOGGER.info("Shutdown of wsl environment.")
                filesystem_io.FilesystemCleaner.clean_prediction_scratch_folder()
                constants.PYSSA_LOGGER.info("Cleaned scratch directory.")
                self.job_manager.stop_prediction_queue()
            else:
                self.job_manager.pop_job_from_queue(signal_tuple[2])
    # </editor-fold>
