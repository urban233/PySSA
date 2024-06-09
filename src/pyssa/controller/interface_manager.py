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
import collections
import glob
import logging
import os
import pathlib
import shutil
import subprocess
import sys
import time

from Bio import SeqRecord
from PyQt5 import QtGui, QtCore
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt

from src.application_process import application_process_manager
from src.pyssa.gui.ui import icon_resources  # this import is used for the icons! DO NOT DELETE THIS  # noqa: F401
from src.pyssa.controller import pymol_session_manager, settings_manager, status_bar_manager, job_manager, watcher, \
  help_manager
from src.pyssa.gui.ui.custom_widgets import job_entry
from src.pyssa.gui.ui.dialogs import dialog_startup
from src.pyssa.gui.ui.views import rename_protein_view, use_project_view, add_sequence_view, add_scene_view, settings_view, predict_protein_view, fasta_file_import_preview_view, rename_sequence_view, add_protein_pair_view, advanced_prediction_configurations, restart_user_pymol_view
from src.pyssa.gui.ui.views import create_project_view, open_project_view, delete_project_view, import_sequence_view
from src.pyssa.gui.ui.views import main_view, distance_analysis_view, results_view, add_protein_view
from src.pyssa.gui.ui.views import hotspots_protein_regions_view
from src.pyssa.gui.ui.styles import styles
from src.pyssa.internal.data_structures import project, protein_pair
from src.pyssa.internal.data_structures.data_classes import current_session
from src.pyssa.internal.thread import tasks
from src.pyssa.internal.thread.async_pyssa import custom_signals
from src.pyssa.io_pyssa import filesystem_io
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.model import proteins_model, protein_pairs_model
from src.pyssa.util import enums, constants, main_window_util, ui_util, exception
from src.pyssa.util.void import rvoid
from src.pyssa.internal.thread import thread_util
from src.tea.thread import task_manager, task_scheduler

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class InterfaceManager:
  """A manager for all views."""

  # <editor-fold desc="Class attributes">
  string_model = QtCore.QStringListModel()  # TODO: this should be removed!

  _main_view: "main_view.MainView"
  """The main view window."""

  _settings_view: "settings_view.SettingsView"
  """The settings view window."""

  _predict_protein_view: "predict_protein_view.PredictProteinView"
  """The predict protein view window."""

  _distance_analysis_view: "distance_analysis_view.DistanceAnalysisView"
  """The distance analysis view window."""

  _create_project_view: "create_project_view.CreateProjectView"
  """The create project view window."""

  _open_project_view: "open_project_view.OpenProjectView"
  """The open project view window."""

  _delete_project_view: "delete_project_view.DeleteProjectView"
  """The delete project view window."""

  _results_view: "results_view.ResultsView"
  """The results view window."""

  _add_protein_view: "add_protein_view.AddProteinView"
  """The add protein view window."""

  _import_sequence_view: "import_sequence_view.ImportSequenceView"
  """The import sequence view window."""

  _fasta_file_import_preview_view: (
      "fasta_file_import_preview_view.FastaFileImportPreviewView"
  )
  """The fasta file import preview view window."""

  _add_sequence_view: "add_sequence_view.AddSequenceView"
  """The add sequence view window."""

  _rename_protein_view: "rename_protein_view.RenameProteinView"
  """The rename protein view window."""

  _use_project_view: "use_project_view.UseProjectView"
  """The use project view window."""

  _hotspots_protein_regions_view: (
      "hotspots_protein_regions_view.HotspotsProteinRegionsView"
  )
  """The hotspots protein regions view window."""

  _add_scene_view: "add_scene_view.AddSceneView"
  """The add scene view window."""

  _add_protein_pair_view: "add_protein_pair_view.AddProteinPairView"
  """The add protein pair view window."""

  _current_workspace: pathlib.Path
  """The current workspace path."""

  _current_project: "project.Project"
  """The current opened project."""

  _current_pymol_session: "current_session.CurrentPymolSession"
  """The current active pymol session."""

  project_lock: QtCore.QMutex
  """A QMutex to lock the current project."""
  # pymol_lock: "locks.PyMOL_LOCK"
  # _application_settings: "settings.Settings"
  current_tab_index: int = 0
  """The current tab index."""

  _workspace_model: QtGui.QStandardItemModel
  """The model for the workspace."""

  _sequence_model: QtGui.QStandardItemModel
  """The model for the sequences."""

  _protein_model: QtGui.QStandardItemModel
  """The model for the proteins."""

  _protein_pair_model: QtGui.QStandardItemModel
  """The model for the protein pairs."""

  _task_manager: "task_manager.TaskManager"
  """The global manager of tasks."""

  _task_scheduler: "task_scheduler.TaskScheduler"
  """The global scheduler of tasks."""

  # </editor-fold>

  def __init__(self) -> None:
    """Constructor."""
    # View definitions
    self._main_view = main_view.MainView()
    self._settings_view = settings_view.SettingsView()
    # self._predict_monomer_view = predict_monomer_view.PredictMonomerView()
    # self._predict_multimer_view = predict_multimer_view.PredictMultimerView()
    self._predict_protein_view = predict_protein_view.PredictProteinView()
    self._distance_analysis_view = distance_analysis_view.DistanceAnalysisView()
    self._create_project_view = create_project_view.CreateProjectView()
    self._open_project_view = open_project_view.OpenProjectView()
    self._delete_project_view = delete_project_view.DeleteProjectView()
    self._hotspots_protein_regions_view = (
        hotspots_protein_regions_view.HotspotsProteinRegionsView()
    )
    self._results_view = results_view.ResultsView()
    self._add_protein_view = add_protein_view.AddProteinView()
    self._import_sequence_view: "import_sequence_view.ImportSequenceView" = (
        import_sequence_view.ImportSequenceView()
    )
    self._fasta_file_import_preview_view = (
        fasta_file_import_preview_view.FastaFileImportPreviewView()
    )
    self._add_sequence_view = add_sequence_view.AddSequenceView(self._main_view)
    self._rename_protein_view = rename_protein_view.RenameProteinView()
    self._rename_sequence_view = rename_sequence_view.RenameSequenceView()
    self._use_project_view = use_project_view.UseProjectView()
    self._add_scene_view = add_scene_view.AddSceneView()
    self._add_protein_pair_view = add_protein_pair_view.AddProteinPairView()
    self._restart_user_pymol_view = (
        restart_user_pymol_view.RestartUserPyMOLView()
    )
    self._advanced_prediction_configurations = (
        advanced_prediction_configurations.AdvancedPredictionConfigurationsView()
    )
    self.help_manager = help_manager.HelpManager()
    self.app_process_manager = (
        application_process_manager.ApplicationProcessManager(
            self._reset_pymol_session
        )
    )
    self.pymol_session_manager = pymol_session_manager.PymolSessionManager(
        self.app_process_manager
    )
    self.job_manager = job_manager.JobManager()
    self.status_bar_manager = status_bar_manager.StatusBarManager(
        self._main_view
    )
    self._settings_manager = settings_manager.SettingsManager()
    self._task_manager = task_manager.TaskManager()

    self._task_scheduler = task_scheduler.TaskScheduler()

    self.watcher = watcher.Watcher()

    self.documentation_window = None
    self.job_entry_widgets = []

    self.refresh_after_job_finished_signal = (
        custom_signals.RefreshAfterJobFinishedSignal()
    )

    # <editor-fold desc="Setup App Settings">
    # self._application_settings = settings.Settings(constants.SETTINGS_DIR, constants.SETTINGS_FILENAME)
    if not os.path.exists(constants.SETTINGS_FULL_FILEPATH):
      constants.PYSSA_LOGGER.info(
          "Settings file not found, open configuration dialog."
      )
      # Configuration dialog to setup setting file
      dialog = dialog_startup.DialogStartup()
      dialog.exec_()

      # checks if the cancel button was pressed
      if dialog_startup.global_var_terminate_app == 1:
        os.remove(constants.SETTINGS_FULL_FILEPATH)
        constants.PYSSA_LOGGER.info(
            "Configuration dialog closed, and removed new settings file."
        )
        sys.exit()

      self._settings_manager.settings.app_launch = 1
      self._settings_manager.settings.workspace_path = pathlib.Path(
          dialog_startup.global_var_startup_workspace
      )
      # self._application_settings.app_launch = 1
      # self._application_settings.workspace_path = pathlib.Path(dialog_startup.global_var_startup_workspace)

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
      # self._application_settings.serialize_settings()
      constants.PYSSA_LOGGER.info("Serialize settings finished.")

      QtWidgets.QApplication.restoreOverrideCursor()
    self._settings_manager.settings = main_window_util.setup_app_settings(
        self._settings_manager.settings
    )

    # </editor-fold>

    # General attributes definitions
    self._current_project = project.Project()
    self._current_pymol_session = current_session.CurrentPymolSession("", "")

    self.project_lock = QtCore.QMutex()
    # self.pymol_lock: "locks.PyMOL_LOCK" = locks.PyMOL_LOCK()

    self.job_manager.start_auxiliary_pymol()
    self.start_app_process_manager()
    # self.start_pymol()

    # Model definitions
    self._workspace_model = QtGui.QStandardItemModel()
    self._sequence_model = QtGui.QStandardItemModel()
    self._protein_model: "proteins_model.ProteinsModel" = (
        proteins_model.ProteinsModel()
    )
    self._protein_pair_model: "protein_pairs_model.ProteinPairsModel" = (
        protein_pairs_model.ProteinPairsModel()
    )
    self._build_workspace_model()

  # <editor-fold desc="Application process manager related methods">
  def start_app_process_manager(self) -> None:
    """Starts a LegacyTasks for the application process manager."""
    self._app_process_manager_thread = tasks.LegacyTask(
        target=self.app_process_manager.check_process,
        args=(0, 0),
        post_func=self._closed_app_process_manager,
    )
    self._app_process_manager_thread.start()

  def _closed_app_process_manager(self) -> None:
    """Await method after the app process manager closed."""
    if self.app_process_manager.pymol_closed():
      logger.info("PyMOL did not crash. The user requested to close PySSA.")
      return
    logger.warning(
        "Check process method of application process manager closed (likely due to a User PyMOL crash)."
    )
    self._restart_user_pymol_view.show()
    self._restart_user_pymol_view.move(30, 300)
    self._app_process_manager_thread = tasks.LegacyTask(
        target=self._recover_user_pymol,
        args=(0, 0),
        post_func=self.__await_recover_user_pymol,
    )
    self._app_process_manager_thread.start()

  def _recover_user_pymol(
      self, a_placeholder_1: int, a_placeholder_2: int
  ) -> tuple[str, str]:
    """Starts recover process of User PyMOL."""
    try:
      logger.info("Starting recovery process ...")
      self.app_process_manager.start_pymol()
      self.app_process_manager.arrange_windows()
      time.sleep(10)
      self.pymol_session_manager.user_pymol_connector.reset_connection()
      self._reset_pymol_session()
    except Exception as e:
      logger.error(
          "The following error occurred during the PyMOL recovery process: (see line below)"
      )
      logger.error(e)
    return "", ""  # These two empty strings are needed for the task class

  def __await_recover_user_pymol(self) -> None:
    """Await method that runs after the recovery process."""
    try:
      logger.info("Finished recovery process.")
      logger.info(
          "Restarting check process routine of application process manager."
      )
      self.start_app_process_manager()
    except Exception as e:
      logger.error(
          "The following error occurred during the PyMOL startup after recovery: (see line below)"
      )
      logger.error(e)
    finally:
      self._restart_user_pymol_view.close()

  def _reset_pymol_session(self) -> None:
    """Resets the pymol session like it was before the User PyMOL crash.

    Notes:
        This function gets executed in the check_process() function of the application process manager!
    """
    if self.pymol_session_manager.session_object_type == "protein":
      self.pymol_session_manager.load_protein_session(
          self.pymol_session_manager.session_objects[0]
      )
    elif self.pymol_session_manager.session_object_type == "protein_pair":
      self.pymol_session_manager.load_protein_pair_session(
          self.pymol_session_manager.session_objects[0]
      )
    else:
      return
    self.pymol_session_manager.load_current_scene()

  # </editor-fold>

  # <editor-fold desc="Getter Methods">
  # <editor-fold desc="Settings">
  def get_settings_manager(self) -> "settings_manager.SettingsManager":
    """Gets the settings manager instance.

    Returns:
        The settings manager instance.
    """
    return self._settings_manager

  def get_application_settings(self) -> "settings.Settings":
    """Gets the application settings.

    Returns:
        The current application settings.
    """
    return self._settings_manager.settings

  # </editor-fold>

  # <editor-fold desc="Workspace">
  def get_workspace_path(self) -> pathlib.Path:
    """Returns the workspace path from the settings manager.

    Returns:
        The workspace path.
    """
    return self._settings_manager.settings.workspace_path

  def get_workspace_model(self) -> QtGui.QStandardItemModel:
    """Gets the workspace model.

    Returns:
        The workspace model.
    """
    return self._workspace_model

  def get_workspace_projects(self) -> QtCore.QStringListModel:
    """Returns the workspace projects.

    Returns:
        A QStringListModel containing the names of all the project files in the workspace.
    """
    db_pattern = os.path.join(
        self._settings_manager.settings.get_workspace_path(), "*.db"
    )
    self.string_model.setStringList(
        # Filters the workspace for all project files based on the xml extension
        [
            os.path.basename(file).replace(".db", "")
            for file in glob.glob(db_pattern)
        ],
    )
    return self.string_model

  def get_workspace_projects_as_list(self) -> list:
    """Returns a list of project names present in the workspace.

    Returns:
        A list of project names without the '.db' extension.
    """
    db_pattern = os.path.join(
        self._settings_manager.settings.get_workspace_path(), "*.db"
    )
    return [
        os.path.basename(file).replace(".db", "")
        for file in glob.glob(db_pattern)
    ]

  # </editor-fold>

  # <editor-fold desc="Getter Methods for view">
  def get_main_view(self) -> "main_view.MainView":
    """Gets the main view of the application.

    Returns:
        The main view of the application.
    """
    return self._main_view

  def get_restart_pymol_view(
      self,
  ) -> "restart_user_pymol_view.RestartUserPyMOLView":
    """Gets the restart pymol view.

    Returns:
        The restart pymol view.
    """
    return self._restart_user_pymol_view

  def get_settings_view(self) -> "settings_view.SettingsView":
    """Gets the settings view.

    Returns:
        The settings view.
    """
    return self._settings_view

  def get_open_view(self) -> "open_project_view.OpenProjectView":
    """Gets the open project view.

    Returns:
        The open project view.
    """
    return self._open_project_view

  def get_create_view(self) -> "create_project_view.CreateProjectView":
    """Gets the create project view.

    Returns:
        The create project view.
    """
    return self._create_project_view

  def get_delete_view(self) -> "delete_project_view.DeleteProjectView":
    """Gets the delete project view.

    Returns:
        The delete project view.
    """
    return self._delete_project_view

  def get_use_project_view(self) -> "use_project_view.UseProjectView":
    """Gets the use project view.

    Returns:
        The use project view.
    """
    return self._use_project_view

  def get_fasta_file_import_preview_view(
      self,
  ) -> "fasta_file_import_preview_view.FastaFileImportPreviewView":
    """Gets the fasta file import preview view.

    Returns:
        The fasta file import preview view.
    """
    return self._fasta_file_import_preview_view

  def get_add_sequence_view(self) -> "add_sequence_view.AddSequenceView":
    """Gets the add sequence view.

    Returns:
        The add sequence view.
    """
    return self._add_sequence_view

  def get_rename_sequence_view(
      self,
  ) -> "rename_sequence_view.RenameSequenceView":
    """Gets the rename sequence view.

    Returns:
        The rename sequence view.
    """
    return self._rename_sequence_view

  def get_predict_protein_view(
      self,
  ) -> "predict_protein_view.PredictProteinView":
    """Gets the predict protein view.

    Returns:
        The predict protein view.
    """
    return self._predict_protein_view

  def get_distance_analysis_view(
      self,
  ) -> "distance_analysis_view.DistanceAnalysisView":
    """Gets the distance analysis view.

    Returns:
        The distance analysis view.
    """
    return self._distance_analysis_view

  def get_results_view(self) -> "results_view.ResultsView":
    """Gets the results view.

    Returns:
        The results view.
    """
    return self._results_view

  def get_add_scene_view(self) -> "add_scene_view.AddSceneView":
    """Gets the add scene view.

    Returns:
        The add scene view.
    """
    return self._add_scene_view

  def get_add_protein_pair_view(
      self,
  ) -> "add_protein_pair_view.AddProteinPairView":
    """Gets the add protein pair view.

    Returns:
        The add protein pair view.
    """
    return self._add_protein_pair_view

  def get_advanced_prediction_configurations_view(
      self,
  ) -> (
      "advanced_prediction_configurations.AdvancedPredictionConfigurationsView"
  ):
    """Gets the advanced configurations view.

    Returns:
        The advanced configurations view.
    """
    return self._advanced_prediction_configurations

  # <editor-fold desc="Getter methods for Sequence Tab in main view">
  def get_current_sequence_list_index(self) -> QtCore.QModelIndex:
    """Gets the current index of the sequence list view.

    Returns:
        The current index of the sequence list view.
    """
    return self._main_view.ui.seqs_list_view.currentIndex()

  def get_current_sequence_list_index_object(self) -> SeqRecord.SeqRecord:
    """Retrieves the SeqRecord object associated with the current index in the sequence list.

    Returns:
        The SeqRecord object from the current index in the sequence list.
    """
    return self.get_current_sequence_list_index().data(
        enums.ModelEnum.OBJECT_ROLE
    )

  def get_import_sequence_view(
      self,
  ) -> "import_sequence_view.ImportSequenceView":
    """Gets the import sequence view.

    Returns:
        The import sequence view.
    """
    return self._import_sequence_view

  # </editor-fold>

  # <editor-fold desc="Getter methods for Protein Tab in main view">

  def get_current_protein_tree_index(self) -> QtCore.QModelIndex:
    """Gets the index of the current protein tree item in the proteins_tree_view.

    Returns:
        The index of the current protein tree item.
    """
    return self._main_view.ui.proteins_tree_view.currentIndex()

  def get_child_index_of_get_current_protein_tree_index(
      self,
  ) -> QtCore.QModelIndex:
    """Get the child index of the current protein tree index.

    Returns:
        The child index of the current protein tree index.
    """
    return self._main_view.ui.proteins_tree_view.currentIndex().child(0, 0)

  def get_current_protein_tree_index_type(self) -> enums.ModelEnum.TYPE_ROLE:
    """Get the current protein tree index type.

    Returns:
        The protein tree index type.
    """
    return self._main_view.ui.proteins_tree_view.model().data(
        self.get_current_protein_tree_index(),
        enums.ModelEnum.TYPE_ROLE,
    )

  def get_current_protein_tree_index_object(self) -> "protein.Protein":
    """Retrieve the current protein tree index object.

    Returns:
        The protein object of the current tree index.
    """
    return self.get_current_protein_tree_index().data(
        enums.ModelEnum.OBJECT_ROLE
    )

  def get_parent_index_object_of_current_protein_tree_index(
      self,
  ) -> "protein.Protein":
    """Gets the parent object of the current protein tree index.

    Returns:
        The parent object of the current protein tree index.
    """
    return (
        self.get_current_protein_tree_index()
        .parent()
        .data(enums.ModelEnum.OBJECT_ROLE)
    )

  def get_current_active_protein_object(self) -> "protein.Protein":
    """Returns the protein object of the current active branch.

    Note:
        This function is handles also the header type and does not throw an error!

    Raises:
        ValueError: if the type is unknown

    Returns:
        a protein object
    """
    tmp_type = self._main_view.ui.proteins_tree_view.currentIndex().data(
        enums.ModelEnum.TYPE_ROLE
    )
    if tmp_type == "protein":
      return self._main_view.ui.proteins_tree_view.currentIndex().data(
          enums.ModelEnum.OBJECT_ROLE
      )
    elif tmp_type == "header":
      return (
          self._main_view.ui.proteins_tree_view.currentIndex()
          .parent()
          .data(enums.ModelEnum.OBJECT_ROLE)
      )
    elif tmp_type == "scene":
      return (
          self._main_view.ui.proteins_tree_view.currentIndex()
          .parent()
          .parent()
          .data(
              enums.ModelEnum.OBJECT_ROLE,
          )
      )
    elif tmp_type == "chain":
      return (
          self._main_view.ui.proteins_tree_view.currentIndex()
          .parent()
          .parent()
          .data(
              enums.ModelEnum.OBJECT_ROLE,
          )
      )
    else:
      raise ValueError("Unknown type!")

  def get_current_header_name(self) -> str:
    """Returns the name of the current header."""
    tmp_type = self._main_view.ui.proteins_tree_view.currentIndex().data(
        enums.ModelEnum.TYPE_ROLE
    )
    if tmp_type == "protein":
      raise ValueError(
          f"Cannot get a header object if the type is: {tmp_type}!"
      )
    elif tmp_type == "header":
      return (
          self._main_view.ui.proteins_tree_view.currentIndex()
          .parent()
          .data(enums.ModelEnum.OBJECT_ROLE)
      )
    elif tmp_type == "scene":
      raise ValueError(
          f"Cannot get a header object if the type is: {tmp_type}!"
      )
    elif tmp_type == "chain":
      raise ValueError(
          f"Cannot get a header object if the type is: {tmp_type}!"
      )
    else:
      raise ValueError("Unknown type!")

  def get_current_active_scene_name(self) -> str:
    """Returns the scene name of the current active branch.

    Returns:
        A scene name.
    """
    tmp_type = self._main_view.ui.proteins_tree_view.currentIndex().data(
        enums.ModelEnum.TYPE_ROLE
    )
    if tmp_type == "protein":
      raise ValueError(f"Cannot get a scene name if the type is: {tmp_type}!")
    elif tmp_type == "header":
      raise ValueError(f"Cannot get a scene name if the type is: {tmp_type}!")
    elif tmp_type == "scene":
      return self._main_view.ui.proteins_tree_view.currentIndex().data(
          Qt.DisplayRole
      )
    elif tmp_type == "chain":
      raise ValueError(f"Cannot get a scene name if the type is: {tmp_type}!")
    else:
      raise ValueError("Unknown type!")

  def get_current_active_chain_object(self) -> "chain.Chain":
    """Returns the chain object of the current active branch.

    Returns:
        A chain object.
    """
    tmp_type = self._main_view.ui.proteins_tree_view.currentIndex().data(
        enums.ModelEnum.TYPE_ROLE
    )
    if tmp_type == "protein":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "header":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "scene":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "chain":
      return self._main_view.ui.proteins_tree_view.currentIndex().data(
          enums.ModelEnum.OBJECT_ROLE
      )
    else:
      raise ValueError("Unknown type!")

  def get_current_active_chain_color_of_protein(self) -> str:
    """Returns the chain color of the current active branch.

    Returns:
        A chain color.
    """
    tmp_type = self._main_view.ui.proteins_tree_view.currentIndex().data(
        enums.ModelEnum.TYPE_ROLE
    )
    if tmp_type == "protein":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "header":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "scene":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "chain":
      return self._main_view.ui.proteins_tree_view.currentIndex().data(
          enums.ModelEnum.CHAIN_COLOR_ROLE
      )
    else:
      raise ValueError("Unknown type!")

  def set_current_active_chain_color_of_protein(self, a_color: str) -> None:
    """Sets the color for the current active chain.

    Args:
        a_color (str): The color to set as the current active chain color of the protein.

    Raises:
        exception.IllegalArgumentError: If `a_color` is either None or an empty string.
        ValueError: If the type of the current index in the proteins tree view is "protein", "header", "scene", or if it is an unknown type.
    """
    # <editor-fold desc="Checks">
    if a_color is None or a_color == "":
      logger.error("a_color is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_color is either None or an empty string."
      )

    # </editor-fold>

    tmp_type = self._main_view.ui.proteins_tree_view.currentIndex().data(
        enums.ModelEnum.TYPE_ROLE
    )
    if tmp_type == "protein":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "header":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "scene":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "chain":
      self._main_view.ui.proteins_tree_view.model().setData(
          self._main_view.ui.proteins_tree_view.currentIndex(),
          a_color,
          enums.ModelEnum.CHAIN_COLOR_ROLE,
      )
    else:
      raise ValueError("Unknown type!")

  def get_protein_repr_toggle_flag(self) -> int:
    """Gets the toggle flag value for protein representation.

    Returns:
        The toggle flag value for protein representation.
    """
    return self._settings_manager.settings.proteins_tab_use_toggle

  def get_add_protein_view(self) -> "add_protein_view.AddProteinView":
    """Gets the instance of the AddProteinView class.

    Returns:
        An instance of the AddProteinView class.
    """
    return self._add_protein_view

  def get_rename_protein_view(self) -> "rename_protein_view.RenameProteinView":
    """Gets the RenameProteinView instance associated with this object.

    Returns:
        The RenameProteinView instance.
    """
    return self._rename_protein_view

  def get_hotspots_protein_regions_view(
      self,
  ) -> "hotspots_protein_regions_view.HotspotsProteinRegionsView":
    """Get the hotspots protein regions view.

    Returns:
        The hotspots protein regions view.
    """
    return self._hotspots_protein_regions_view

  # </editor-fold>

  # <editor-fold desc="Getter methods for Protein Pairs Tab in main view">
  def get_current_protein_pair_tree_index(self) -> QtCore.QModelIndex:
    """Returns the current index of the protein pair tree view.

    Returns:
        The current index of the protein pair tree view.
    """
    return self._main_view.ui.protein_pairs_tree_view.currentIndex()

  def get_child_index_of_get_current_protein_pair_tree_index(
      self,
  ) -> QtCore.QModelIndex:
    """Gets the child index of the current protein pair tree index.

    Returns:
        The child index of the current protein pair tree index.
    """
    return self._main_view.ui.protein_pairs_tree_view.currentIndex().child(0, 0)

  def get_current_protein_pair_tree_index_type(
      self,
  ) -> enums.ModelEnum.TYPE_ROLE:
    """Gets the protein pair tree index type for the current protein pair in the UI.

    Returns:
        enums.ModelEnum.TYPE_ROLE value representing the protein pair tree index type.
    """
    return self._main_view.ui.protein_pairs_tree_view.model().data(
        self.get_current_protein_pair_tree_index(),
        enums.ModelEnum.TYPE_ROLE,
    )

  def get_current_protein_pair_tree_index_object(
      self,
  ) -> "protein_pair.ProteinPair":
    """Gets the ProteinPair object associated with the current protein pair tree index.

    Returns:
        The ProteinPair object associated with the current protein pair tree index.
    """
    return self.get_current_protein_pair_tree_index().data(
        enums.ModelEnum.OBJECT_ROLE
    )

  def get_parent_index_object_of_current_protein_pair_tree_index(
      self,
  ) -> "protein_pair.ProteinPair":
    """Gets the parent object of the current protein pair tree index.

    Returns:
        protein_pair.ProteinPair: The parent object of the current protein pair tree index.
    """
    return (
        self.get_current_protein_pair_tree_index()
        .parent()
        .data(enums.ModelEnum.OBJECT_ROLE)
    )

  def get_grand_parent_index_object_of_current_protein_pair_tree_index(
      self,
  ) -> "protein_pair.ProteinPair":
    """Gets the grand-parent index object of the current protein pair tree index.

    Returns:
        The grand-parent index object of the current protein pair tree index.
    """
    return (
        self.get_current_protein_pair_tree_index()
        .parent()
        .parent()
        .data(enums.ModelEnum.OBJECT_ROLE)
    )

  def get_current_active_protein_pair_object(
      self,
  ) -> "protein_pair.ProteinPair":
    """Returns the protein pair object of the current active branch.

    Note:
        This function is handles also the header type and does not throw an error!

    Raises:
        ValueError: if the type is unknown

    Returns:
        A protein pair object.
    """
    tmp_type = self._main_view.ui.protein_pairs_tree_view.currentIndex().data(
        enums.ModelEnum.TYPE_ROLE
    )
    tmp_display_role = (
        self._main_view.ui.protein_pairs_tree_view.currentIndex().data(
            Qt.DisplayRole
        )
    )
    if tmp_type == "protein_pair":
      return self._main_view.ui.protein_pairs_tree_view.currentIndex().data(
          enums.ModelEnum.OBJECT_ROLE
      )
    elif tmp_type == "protein":
      return (
          self._main_view.ui.protein_pairs_tree_view.currentIndex()
          .parent()
          .data(enums.ModelEnum.OBJECT_ROLE)
      )
    elif tmp_type == "scene":
      return (
          self._main_view.ui.protein_pairs_tree_view.currentIndex()
          .parent()
          .parent()
          .data(
              enums.ModelEnum.OBJECT_ROLE,
          )
      )
    elif tmp_type == "chain":
      return (
          self._main_view.ui.protein_pairs_tree_view.currentIndex()
          .parent()
          .parent()
          .parent()
          .data(
              enums.ModelEnum.OBJECT_ROLE,
          )
      )
    elif tmp_type == "header" and tmp_display_role == "Scenes":
      return (
          self._main_view.ui.protein_pairs_tree_view.currentIndex()
          .parent()
          .data(enums.ModelEnum.OBJECT_ROLE)
      )
    elif tmp_type == "header" and tmp_display_role == "Chains":
      return (
          self._main_view.ui.protein_pairs_tree_view.currentIndex()
          .parent()
          .parent()
          .data(
              enums.ModelEnum.OBJECT_ROLE,
          )
      )
    elif tmp_type == "":
      raise ValueError("Nothing there!")
    else:
      raise ValueError("Unknown type!")

  def get_current_active_protein_object_of_protein_pair(
      self,
  ) -> "protein.Protein":
    """Returns the protein object of the current active protein pair branch.

    Note:
        This function is handles also the header type and does not throw an error!

    Raises:
        ValueError: if the type is unknown

    Returns:
        A protein object.
    """
    tmp_type = self._main_view.ui.protein_pairs_tree_view.currentIndex().data(
        enums.ModelEnum.TYPE_ROLE
    )
    tmp_display_role = (
        self._main_view.ui.protein_pairs_tree_view.currentIndex().data(
            Qt.DisplayRole
        )
    )
    if tmp_type == "protein":
      return self._main_view.ui.protein_pairs_tree_view.currentIndex().data(
          enums.ModelEnum.OBJECT_ROLE
      )
    elif tmp_type == "scene":
      raise ValueError(f"Cannot get a protein object with the type: {tmp_type}")
    elif tmp_type == "chain":
      return (
          self._main_view.ui.protein_pairs_tree_view.currentIndex()
          .parent()
          .parent()
          .data(
              enums.ModelEnum.OBJECT_ROLE,
          )
      )
    elif tmp_type == "header" and tmp_display_role == "Scenes":
      raise ValueError(f"Cannot get a protein object with the type: {tmp_type}")
    elif tmp_type == "header" and tmp_display_role == "Chains":
      return (
          self._main_view.ui.protein_pairs_tree_view.currentIndex()
          .parent()
          .data(enums.ModelEnum.OBJECT_ROLE)
      )
    elif tmp_type == "protein_pair":
      raise ValueError(
          f"Cannot get a protein object if the type is: {tmp_type}!"
      )
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
    tmp_type = self._main_view.ui.protein_pairs_tree_view.currentIndex().data(
        enums.ModelEnum.TYPE_ROLE
    )
    tmp_display_role = (
        self._main_view.ui.protein_pairs_tree_view.currentIndex().data(
            Qt.DisplayRole
        )
    )
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
        A chain object.
    """
    tmp_type = self._main_view.ui.protein_pairs_tree_view.currentIndex().data(
        enums.ModelEnum.TYPE_ROLE
    )
    tmp_display_role = (
        self._main_view.ui.protein_pairs_tree_view.currentIndex().data(
            Qt.DisplayRole
        )
    )
    if tmp_type == "protein":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "scene":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "chain":
      return self._main_view.ui.protein_pairs_tree_view.currentIndex().data(
          enums.ModelEnum.OBJECT_ROLE
      )
    elif tmp_type == "header" and tmp_display_role == "Scenes":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "header" and tmp_display_role == "Chains":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "protein_pair":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    else:
      raise ValueError("Unknown type!")

  def get_protein_pair_repr_toggle_flag(self) -> int:
    """Gets the toggle flag for displaying protein pair representations in the protein pairs tab.

    Returns:
        The toggle flag value.
    """
    return self._settings_manager.settings.protein_pairs_tab_use_toggle

  # </editor-fold>

  # </editor-fold>

  def get_current_project(self) -> "project.Project":
    """Retrieve the current project.

    Returns:
        A project object.
    """
    return self._current_project

  def get_information_about_current_session(self) -> tuple:
    """Returns information about the current session.

    Returns a tuple containing the session name and object type of the current session.

    Returns:
        A tuple of two elements - the session name (str) and object type (str) of the current session.
    """
    return (
        self._current_pymol_session.session_name,
        self._current_pymol_session.object_type,
    )

  def get_protein_model(self) -> QtGui.QStandardItemModel:
    """Gets the protein model.

    Returns:
        The protein model used for displaying protein data.
    """
    return self._protein_model

  def get_task_manager(self) -> "task_manager.TaskManager":
    """Gets the task manager.

    Returns:
        The task manager instance.
    """
    return self._task_manager

  def get_task_scheduler(self) -> "task_scheduler.TaskScheduler":
    return self._task_scheduler

  # </editor-fold>

  def add_protein_to_current_project(
      self, a_protein: "protein.Protein"
  ) -> None:
    """Adds a given Protein object to the current project.

    The Protein object must already exist.

    Args:
        a_protein (protein.Protein): A Protein object that will be added to the current project.

    Raises:
        exception.IllegalArgumentError: If `a_protein` is None.
    """
    # <editor-fold desc="Checks">
    if a_protein is None:
      logger.error("a_protein is None.")
      raise exception.IllegalArgumentError("a_protein is None.")

    # </editor-fold>

    self._current_project.add_existing_protein(a_protein)

  def add_protein_pair_to_current_project(
      self, a_protein_pair: "protein_pair.ProteinPair"
  ) -> None:
    """Adds a protein pair to the current project.

    Args:
        a_protein_pair (protein_pair.ProteinPair): A ProteinPair object that will be added to the current project.

    Raises:
        exception.IllegalArgumentError: If `a_protein_pair` is None.
    """
    # <editor-fold desc="Checks">
    if a_protein_pair is None:
      logger.error("a_protein_pair is None.")
      raise exception.IllegalArgumentError("a_protein_pair is None.")

    # </editor-fold>

    self._current_project.add_protein_pair(a_protein_pair)

  # <editor-fold desc="Setter Methods">
  def set_new_project(self, the_current_project: "project.Project") -> None:
    """Sets the new current project into the interface manager.

    Args:
        the_current_project (project.Project): The new project to be set as the current project.

    Raises:
        exception.IllegalArgumentError: If `the_current_project` is None.
    """
    # <editor-fold desc="Checks">
    if the_current_project is None:
      logger.error("the_current_project is None.")
      raise exception.IllegalArgumentError("the_current_project is None.")

    # </editor-fold>

    self._current_project = the_current_project
    self._sequence_model.clear()
    self._build_sequences_model()
    self._protein_model.clear()
    self._build_proteins_model()
    self._protein_pair_model.clear()
    self._build_protein_pairs_model()

  # fixme: This function does not used anymore!!!
  def set_new_workspace(self, the_current_workspace: str) -> None:
    """Sets the new current workspace into the interface manager.

    Args:
        the_current_workspace (str): The path to the new workspace.

    Raises:
        exception.IllegalArgumentError: If `the_current_workspace` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if the_current_workspace is None or the_current_workspace == "":
      logger.error("the_current_workspace is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "the_current_workspace is either None or an empty string."
      )

    # </editor-fold>

    self._settings_manager.settings.workspace_path = the_current_workspace
    self._workspace_model.clear()
    self._build_workspace_model()

  def set_new_session_information(
      self, a_session_name: str, an_object_name: str
  ) -> None:
    """Sets the new session information for the current session.

    Args:
        a_session_name (str): The name of the session to be set. Cannot be None or an empty string.
        an_object_name (str): The name of the object to be set. Cannot be None or an empty string.

    Raises:
        exception.IllegalArgumentError: If either `a_session_name` or `an_object_name` is None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_session_name is None or a_session_name == "":
      logger.error("a_session_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_session_name is either None or an empty string."
      )
    if an_object_name is None or an_object_name == "":
      logger.error("an_object_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "an_object_name is either None or an empty string."
      )

    # </editor-fold>

    self._current_pymol_session.session_name = a_session_name
    self._current_pymol_session.object_type = an_object_name

  # <editor-fold desc="Protein">
  def set_current_chain_color_for_ui_for_proteins(
      self,
      the_pymol_session_manager: "pymol_session_manager.PymolSessionManager",
  ) -> None:
    """Checks which color the protein chain has and sets the index of the combobox accordingly.

    Args:
        the_pymol_session_manager (pymol_session_manager.PymolSessionManager): The PymolSessionManager object.
    """
    # <editor-fold desc="Checks">
    if the_pymol_session_manager is None:
      logger.error("the_pymol_session_manager is None.")
      raise exception.IllegalArgumentError("the_pymol_session_manager is None.")

    # </editor-fold>

    tmp_protein = self.get_current_active_protein_object()
    tmp_chain = self.get_current_active_chain_object()
    if the_pymol_session_manager.is_the_current_protein_in_session(
        self.get_current_active_protein_object().get_molecule_object()
    ):
      # fixme: This can easily be bypassed by a power user if the first residue color is changed
      if tmp_chain.chain_type == "protein_chain":
        tmp_protein.pymol_selection.selection_string = (
            f"first chain {tmp_chain.chain_letter} and name CA"
        )
      else:
        tmp_protein.pymol_selection.selection_string = (
            f"first chain {tmp_chain.chain_letter}"
        )
      if self.pymol_session_manager.get_residue_color_config_of_a_given_selection(
          f"(first elem N, first elem C, first elem O) and chain {tmp_chain.chain_letter}",
          tmp_chain.chain_letter,
      ).atoms_are_colored_by_elements():
        self._main_view.ui.lbl_protein_current_color.setText("By Element    ")
        self._main_view.tg_protein_color_atoms.toggle_button.setChecked(True)
      else:
        rvoid(
            tmp_chain.get_color(
                tmp_protein.pymol_selection.selection_string,
                self.pymol_session_manager,
            )
        )
        self._main_view.ui.lbl_protein_current_color.setText(
            f"{tmp_chain.pymol_parameters['chain_color']}    "
        )
        self._main_view.tg_protein_color_atoms.toggle_button.setChecked(False)

  def set_repr_state_in_ui_for_protein_chain(
      self,
      the_pymol_session_manager: "pymol_session_manager.PymolSessionManager",
  ) -> None:
    """Sets the representation state in the user interface for a protein chain.

    Args:
        the_pymol_session_manager (pymol_session_manager.PymolSessionManager): The PymolSessionManager object.

    Raises:
        exception.IllegalArgumentError: If `the_pymol_session_manager` is None.
    """
    # <editor-fold desc="Checks">
    if the_pymol_session_manager is None:
      logger.error("the_pymol_session_manager is None.")
      raise exception.IllegalArgumentError("the_pymol_session_manager is None.")

    # </editor-fold>

    tmp_protein = self.get_current_active_protein_object()
    tmp_chain = self.get_current_active_chain_object()
    if the_pymol_session_manager.is_the_current_protein_in_session(
        self.get_current_active_protein_object().get_molecule_object()
    ):
      # fixme: This can easily be bypassed by a power user if the first residue color is changed
      tmp_protein.pymol_selection.selection_string = (
          f"first chain {tmp_chain.chain_letter}"
      )
      print(f"This is a chain type: {tmp_chain.chain_type}")
      if tmp_chain.chain_type == "protein_chain":
        self._main_view.ui.frame_protein_repr.setEnabled(True)
        tmp_repr_state = self.pymol_session_manager.get_chain_repr_state(
            tmp_protein.pymol_selection.selection_string, tmp_chain.chain_letter
        )
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

  def set_repr_state_in_ui_for_protein_pair_chain(
      self,
      a_protein_name: str,
      the_pymol_session_manager: "pymol_session_manager.PymolSessionManager",
  ) -> None:
    """Updates the representation state in the user interface for a protein pair chain.

    Args:
        a_protein_name (str): The name of the protein.
        the_pymol_session_manager (pymol_session_manager.PymolSessionManager): The instance of the PymolSessionManager class.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None or if `a_protein_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_protein_name is None or a_protein_name == "":
      logger.error("a_protein_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_protein_name is either None or an empty string."
      )
    if the_pymol_session_manager is None:
      logger.error("the_pymol_session_manager is None.")
      raise exception.IllegalArgumentError("the_pymol_session_manager is None.")

    # </editor-fold>

    tmp_protein = self.get_current_active_protein_object_of_protein_pair()
    tmp_chain = self.get_current_active_chain_object_of_protein_pair()
    if the_pymol_session_manager.is_the_current_protein_pair_in_session(
        self.get_current_active_protein_pair_object().name
    ):
      # fixme: This can easily be bypassed by a power user if the first residue color is changed
      tmp_protein.pymol_selection.selection_string = (
          f"first chain {tmp_chain.chain_letter} and {a_protein_name}"
      )
      if tmp_chain.chain_type == "protein_chain":
        self._main_view.ui.frame_protein_pair_repr.setEnabled(True)
        tmp_repr_state = self.pymol_session_manager.get_chain_repr_state(
            tmp_protein.pymol_selection.selection_string, tmp_chain.chain_letter
        )
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

  def manage_check_state_of_protein_pair(self, tmp_repr_state: dict) -> None:
    """Updates the state of the checkboxes for protein pair representations in the user interface based on the values in the given `tmp_repr_state` dictionary.

    Args:
        tmp_repr_state (dict): A dictionary representing the state of protein pair representations.

    Raises:
        exception.IllegalArgumentError: If `tmp_repr_state` is None.
    """
    # <editor-fold desc="Checks">
    if tmp_repr_state is None:
      logger.error("tmp_repr_state is None.")
      raise exception.IllegalArgumentError("tmp_repr_state is None.")

    # </editor-fold>

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

  def manage_toggle_state_of_protein_pair_repr(
      self, tmp_repr_state: dict
  ) -> None:
    """Updates the state of the toggles for protein pair representations in the user interface based on the values in the given `tmp_repr_state` dictionary.

    Args:
        tmp_repr_state (dict): A dictionary representing the state of protein pair representations.

    Raises:
        exception.IllegalArgumentError: If `tmp_repr_state` is None.
    """
    # <editor-fold desc="Checks">
    if tmp_repr_state is None:
      logger.error("tmp_repr_state is None.")
      raise exception.IllegalArgumentError("tmp_repr_state is None.")

    # </editor-fold>

    if tmp_repr_state[enums.PyMOLRepresentation.CARTOON.value] == 0:
      ui_util.set_checked_async(
        self._main_view.tg_protein_pair_cartoon.toggle_button, False
      )
    else:
      ui_util.set_checked_async(
        self._main_view.tg_protein_pair_cartoon.toggle_button, True
      )
    if tmp_repr_state[enums.PyMOLRepresentation.STICKS.value] == 0:
      ui_util.set_checked_async(
        self._main_view.tg_protein_pair_sticks.toggle_button, False
      )
    else:
      ui_util.set_checked_async(
        self._main_view.tg_protein_pair_sticks.toggle_button, True
      )
    if tmp_repr_state[enums.PyMOLRepresentation.RIBBON.value] == 0:
      ui_util.set_checked_async(
        self._main_view.tg_protein_pair_ribbon.toggle_button, False
      )
    else:
      ui_util.set_checked_async(
        self._main_view.tg_protein_pair_ribbon.toggle_button, True
      )
    if tmp_repr_state[enums.PyMOLRepresentation.LINES.value] == 0:
      ui_util.set_checked_async(
        self._main_view.tg_protein_pair_lines.toggle_button, False
      )
    else:
      ui_util.set_checked_async(
        self._main_view.tg_protein_pair_lines.toggle_button, True
      )
    if tmp_repr_state[enums.PyMOLRepresentation.SPHERES.value] == 0:
      ui_util.set_checked_async(
        self._main_view.tg_protein_pair_spheres.toggle_button, False
      )
    else:
      ui_util.set_checked_async(
        self._main_view.tg_protein_pair_spheres.toggle_button, True
      )
    if tmp_repr_state[enums.PyMOLRepresentation.DOTS.value] == 0:
      ui_util.set_checked_async(
        self._main_view.tg_protein_pair_dots.toggle_button, False
      )
    else:
      ui_util.set_checked_async(
        self._main_view.tg_protein_pair_dots.toggle_button, True
      )
    if tmp_repr_state[enums.PyMOLRepresentation.MESH.value] == 0:
      ui_util.set_checked_async(
        self._main_view.tg_protein_pair_mesh.toggle_button, False
      )
    else:
      ui_util.set_checked_async(
        self._main_view.tg_protein_pair_mesh.toggle_button, True
      )
    if tmp_repr_state[enums.PyMOLRepresentation.SURFACE.value] == 0:
      ui_util.set_checked_async(
        self._main_view.tg_protein_pair_surface.toggle_button, False
      )
    else:
      ui_util.set_checked_async(
        self._main_view.tg_protein_pair_surface.toggle_button, True
      )

  # </editor-fold>

  # <editor-fold desc="Build Methods">
  def _build_sequences_model(self) -> None:
    """Builds the sequences model for the current project."""
    if len(self._current_project.sequences) > 0:
      tmp_root_item = self._sequence_model.invisibleRootItem()
      for tmp_sequence in self._current_project.sequences:
        tmp_sequence_item = QtGui.QStandardItem(tmp_sequence.name)
        tmp_sequence_item.setData(tmp_sequence, enums.ModelEnum.OBJECT_ROLE)
        if "," in tmp_sequence.seq:
          tmp_sequence_item.setData(
              enums.ModelTypeEnum.MULTIMER_SEQ, enums.ModelEnum.TYPE_ROLE
          )
        else:
          tmp_sequence_item.setData(
              enums.ModelTypeEnum.MONOMER_SEQ, enums.ModelEnum.TYPE_ROLE
          )
        tmp_root_item.appendRow(tmp_sequence_item)

  def _build_proteins_model(self) -> None:
    """Builds the proteins model for the current project."""
    if len(self._current_project.proteins) > 0:
      tmp_main_socket, tmp_general_purpose_socket = (
          self.job_manager.get_general_purpose_socket_pair()
      )
      self._protein_model.build_model_from_scratch(
          self._current_project.proteins,
          tmp_main_socket,
          tmp_general_purpose_socket,
      )
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
    """Builds the protein pairs model for the current project."""
    if len(self._current_project.protein_pairs) > 0:
      tmp_main_socket, tmp_general_purpose_socket = (
          self.job_manager.get_general_purpose_socket_pair()
      )
      self._protein_pair_model.build_model_from_scratch(
          self._current_project.protein_pairs,
          tmp_main_socket,
          tmp_general_purpose_socket,
      )
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
  def refresh_main_view(self) -> None:
    """Modifies the UI of the main view based on an app model."""
    # <editor-fold desc="Thread check">
    if thread_util.is_main_thread() is False:
      logger.warning(
          "Method 'refresh_main_view' was called from a separate thread. "
          "Cannot run this method because it has to be called from the main thread!",
      )
      raise exception.NotMainThreadError()

    # </editor-fold>

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

    self._main_view.ui.menuProject.setEnabled(True)
    # A project is open
    if self._current_project.get_project_name() != "":
      styles.set_stylesheet(self._main_view)
      self._main_view.ui.lbl_project_name.show()
      self._main_view.ui.lbl_project_name.setText(
          f"Project Name: {self._current_project.get_project_name()}"
      )
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

        tmp_are_all_monomer_sequences_predicted = (
            self._check_if_sequences_are_already_predicted()[0]
        )
        tmp_are_all_multimer_sequences_predicted = (
            self._check_if_sequences_are_already_predicted()[1]
        )
        if tmp_are_all_monomer_sequences_predicted:
          self._main_view.ui.action_predict_monomer.setEnabled(False)
        if tmp_are_all_multimer_sequences_predicted:
          self._main_view.ui.action_predict_multimer.setEnabled(False)
        if (
            tmp_are_all_monomer_sequences_predicted
            and tmp_are_all_multimer_sequences_predicted
        ):
          self._main_view.ui.menuPrediction.setEnabled(
              False
          )  # the entire menu can be disabled because all sequences are predicted
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
        self._main_view.ui.btn_protein_tree_view_expand.setEnabled(True)
        self._main_view.ui.btn_protein_tree_view_collapse.setEnabled(True)
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
          if self._main_view.ui.proteins_tree_view.model().data(
                  self._main_view.ui.proteins_tree_view.currentIndex(), enums.ModelEnum.TYPE_ROLE
          ) == "chain":
            self._main_view.ui.action_protein_regions.setEnabled(True)
          else:
            self._main_view.ui.action_protein_regions.setEnabled(False)

        try:
          if self.get_current_active_scene_name() == "base":
            self._main_view.ui.btn_delete_protein_scene.setEnabled(False)
        except ValueError:
          pass  # Is necessary because tree view selection might be something else than a scene.

        if len(self._current_project.protein_pairs) > 0:
          # A project has protein pair(s)
          self._main_view.ui.protein_pairs_tree_view.setModel(
              self._protein_pair_model
          )
          self._main_view.ui.protein_pairs_tree_view.setHeaderHidden(True)
          self._main_view.ui.btn_protein_pair_tree_view_expand.setEnabled(True)
          self._main_view.ui.btn_protein_pair_tree_view_collapse.setEnabled(
              True
          )
          # It is possible to view results
          self._main_view.ui.menuResults.setEnabled(True)
          try:
            if self.get_current_active_scene_name_of_protein_pair() == "base":
              self._main_view.ui.btn_delete_protein_pair_scene.setEnabled(False)
          except ValueError:
            pass  # Is necessary because tree view selection might be something else than a scene.
          # Check for any possible protein region menu activation
          if self._main_view.ui.project_tab_widget.currentIndex() == 2:
            # User is on the protein pairs tab
            if self.pymol_session_manager.is_the_current_session_empty():
              self._main_view.ui.action_protein_regions.setEnabled(False)
            else:
              tmp_current_protein_pairs_tree_type = self._main_view.ui.protein_pairs_tree_view.model().data(
                      self._main_view.ui.protein_pairs_tree_view.currentIndex(), enums.ModelEnum.TYPE_ROLE
              )
              tmp_current_protein_pairs_tree_object: "protein_pair.ProteinPair" = self.get_current_active_protein_pair_object()
              if tmp_current_protein_pairs_tree_type == "chain" and self.pymol_session_manager.is_the_current_protein_pair_in_session(tmp_current_protein_pairs_tree_object.name):
                self._main_view.ui.action_protein_regions.setEnabled(True)
              else:
                self._main_view.ui.action_protein_regions.setEnabled(False)
        else:
          # A project has no protein pair(s)
          # Protein Pairs tab
          self._main_view.ui.btn_protein_pair_tree_view_expand.setEnabled(False)
          self._main_view.ui.btn_protein_pair_tree_view_collapse.setEnabled(
              False
          )
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
        self._main_view.ui.btn_protein_tree_view_expand.setEnabled(False)
        self._main_view.ui.btn_protein_tree_view_collapse.setEnabled(False)
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

  def refresh_workspace_model(self) -> None:
    """Clears the workspace model and rebuilds it.

    This method is responsible for refreshing the workspace model by clearing it and then rebuilding it using the
    _build_workspace_model method.
    """
    self._workspace_model.clear()
    self._build_workspace_model()
    self.status_bar_manager.show_temporary_message(
      "Workspace changed successfully."
    )

  def refresh_sequence_model(self) -> None:
    """Refreshes the sequence model.

    This method clears the sequence model and then builds it again.
    """
    self._sequence_model.clear()
    self._build_sequences_model()

  def refresh_protein_model(self) -> None:
    """Refreshes the protein model.

    This method clears the protein pair model and then builds it again.
    """
    self._protein_model.clear()
    self._build_proteins_model()

  def refresh_protein_pair_model(self) -> None:
    """Refreshes the protein pair model.

    This method clears the protein pair model and then builds it again.
    """
    self._protein_pair_model.clear()
    self._build_protein_pairs_model()

  def _build_workspace_model(self) -> None:
    """Builds the workspace model.

    This method populates the workspace model with project items based on the database files found in the workspace directory.
    """
    tmp_workspace = self._settings_manager.settings.workspace_path
    db_pattern = os.path.join(tmp_workspace, "*.db")
    tmp_root_item = self._workspace_model.invisibleRootItem()
    for tmp_filename in [
        os.path.basename(file).replace(".db", "")
        for file in glob.glob(db_pattern)
    ]:
      tmp_project_item = QtGui.QStandardItem(tmp_filename)
      tmp_filepath = pathlib.Path(f"{tmp_workspace}/{tmp_filename}.db")
      tmp_project_item.setData(tmp_filepath, enums.ModelEnum.FILEPATH_ROLE)
      tmp_root_item.appendRow(tmp_project_item)

  def add_project_to_workspace_model(self, a_filename: str) -> None:
    """Adds a project to the workspace model.

    Args:
        a_filename (str): The filename of the project to be added.

    Raises:
        IllegalArgumentError: If `a_filename` is None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_filename is None or a_filename == "":
      logger.error("a_filename is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_filename is either None or an empty string."
      )

    # </editor-fold>

    tmp_root_item = self._workspace_model.invisibleRootItem()
    tmp_project_item = QtGui.QStandardItem(a_filename)
    tmp_filepath = pathlib.Path(
        f"{self._settings_manager.settings.workspace_path}/{a_filename}.db"
    )
    tmp_project_item.setData(tmp_filepath, enums.ModelEnum.FILEPATH_ROLE)
    tmp_root_item.appendRow(tmp_project_item)

  def disable_proteins_tab_buttons(self) -> None:
    """Disable Proteins Tab buttons.

    Raises:
        NotMainThreadError: If the method is called from a separate thread instead of the main thread.
    """
    # <editor-fold desc="Thread check">
    if thread_util.is_main_thread() is False:
      logger.warning(
          "Method 'disable_proteins_tab_buttons' was called from a separate thread. "
          "Cannot run this method because it has to be called from the main thread!",
      )
      raise exception.NotMainThreadError()

    # </editor-fold>

    self._main_view.ui.btn_save_protein.setEnabled(False)
    self._main_view.ui.btn_delete_protein.setEnabled(False)
    self._main_view.ui.btn_open_protein_session.setEnabled(False)
    self._main_view.ui.btn_create_protein_scene.setEnabled(False)
    self._main_view.ui.btn_update_protein_scene.setEnabled(False)
    self._main_view.ui.btn_delete_protein_scene.setEnabled(False)

  def disable_protein_pairs_tab_buttons(self) -> None:
    """Disable Protein Pairs Tab buttons.

    Raises:
        NotMainThreadError: If the method is called from a separate thread instead of the main thread.
    """
    # <editor-fold desc="Thread check">
    if thread_util.is_main_thread() is False:
      logger.warning(
          "Method 'disable_protein_pairs_tab_buttons' was called from a separate thread. "
          "Cannot run this method because it has to be called from the main thread!",
      )
      raise exception.NotMainThreadError()

    # </editor-fold>

    self._main_view.ui.btn_delete_protein_pair.setEnabled(False)
    self._main_view.ui.btn_open_protein_pair_session.setEnabled(False)
    self._main_view.ui.btn_create_protein_pair_scene.setEnabled(False)
    self._main_view.ui.btn_update_protein_pair_scene.setEnabled(False)
    self._main_view.ui.btn_delete_protein_pair_scene.setEnabled(False)

  # </editor-fold>

  # <editor-fold desc="Progress bar methods">
  def update_progress_bar(self, value: int, message: str) -> None:
    """Updates the progress bar with the given value and message.

    Args:
        value (int): The value to update the progress bar with. Must be between 0 and 100.
        message (str): The message to display on the progress bar.

    Raises:
        exception.IllegalArgumentError: If `value` is None.
        exception.IllegalArgumentError: If `message` is either None or an empty string.
        ValueError: If the value is not between 0 and 100.
    """
    # <editor-fold desc="Checks">
    if value is None:
      logger.error("value is None.")
      raise exception.IllegalArgumentError("value is None.")
    if message is None or message == "":
      logger.error("message is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "message is either None or an empty string."
      )
    if value < 0 or value > 100:
      raise ValueError("Value for progress bar must be between 0 and 100!")

    # </editor-fold>

    self._main_view.progress_bar.show()
    self._main_view.progress_bar.setFormat(message)
    self._main_view.progress_bar.setValue(value)

  def hide_progress_bar(self) -> None:
    """Hide the progress bar."""
    self._main_view.progress_bar.hide()

  # </editor-fold>

  # <editor-fold desc="Sequences">
  def show_sequence_parameters(
      self, a_sequence_item: QtGui.QStandardItem
  ) -> None:
    """Sets up the sequences table in the main view with the provided sequence item.

    Args:
        a_sequence_item (QtGui.QStandardItem): The QStandardItem representing a sequence.

    Raises:
        exception.IllegalArgumentError: If `a_sequence_item` is None.
    """
    # <editor-fold desc="Checks">
    if a_sequence_item is None:
      logger.error("a_sequence_item is None.")
      raise exception.IllegalArgumentError("a_sequence_item is None.")

    # </editor-fold>

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
    tmp_sequence_item = QtWidgets.QTableWidgetItem(
        f"{tmp_sequence.seq[:15]} ..."
    )
    tmp_sequence_item.setToolTip("Click to view complete sequence")
    tmp_sequence_item.setData(enums.ModelEnum.OBJECT_ROLE, tmp_sequence)
    self._main_view.ui.seqs_table_widget.setItem(1, 1, tmp_sequence_item)
    # Table item flags
    tmp_seq_name_item.setFlags(tmp_seq_name_item.flags() & ~Qt.ItemIsEditable)
    tmp_sequence_item.setFlags(tmp_sequence_item.flags() & ~Qt.ItemIsEditable)
    tmp_name_label_item.setFlags(
        tmp_name_label_item.flags() & ~Qt.ItemIsEditable
    )
    tmp_sequence_label_item.setFlags(
        tmp_sequence_label_item.flags() & ~Qt.ItemIsEditable
    )
    self._main_view.ui.seqs_table_widget.resizeColumnsToContents()

  def show_menu_options_with_seq(self) -> None:
    """Disables specific menu options in the UI."""
    self._main_view.ui.menuAnalysis.setEnabled(False)
    self._main_view.ui.menuResults.setEnabled(False)
    self._main_view.ui.menuImage.setEnabled(False)
    self._main_view.ui.menuHotspots.setEnabled(False)

  def show_menu_options_without_seq(self) -> None:
    """Show menu options without a sequence in a project."""
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
      if (
          tmp_item.data(enums.ModelEnum.TYPE_ROLE)
          == enums.ModelTypeEnum.MONOMER_SEQ
      ):
        tmp_contains_monomer = True
      elif (
          tmp_item.data(enums.ModelEnum.TYPE_ROLE)
          == enums.ModelTypeEnum.MULTIMER_SEQ
      ):
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
        A tuple of booleans where True means all sequences of the type are already predicted.
    """
    tmp_number_of_monomer_sequences = 0
    tmp_number_of_multimer_sequences = 0
    tmp_number_of_predicted_monomer_sequences = 0
    tmp_number_of_predicted_multimer_sequences = 0
    for tmp_row in range(self._sequence_model.rowCount()):
      tmp_item = self._sequence_model.item(tmp_row, 0)
      if (
          tmp_item.data(enums.ModelEnum.TYPE_ROLE)
          == enums.ModelTypeEnum.MONOMER_SEQ
      ):
        tmp_number_of_monomer_sequences += 1
      elif (
          tmp_item.data(enums.ModelEnum.TYPE_ROLE)
          == enums.ModelTypeEnum.MULTIMER_SEQ
      ):
        tmp_number_of_multimer_sequences += 1
      else:
        logger.warning("Invalid TYPE_ROLE!")

      if (
          self._current_project.is_sequence_as_protein_in_project(
              tmp_item.data(enums.ModelEnum.OBJECT_ROLE).name
          )
          and tmp_item.data(enums.ModelEnum.TYPE_ROLE)
          == enums.ModelTypeEnum.MONOMER_SEQ
      ):
        tmp_number_of_predicted_monomer_sequences += 1
      elif (
          self._current_project.is_sequence_as_protein_in_project(
              tmp_item.data(enums.ModelEnum.OBJECT_ROLE).name
          )
          and tmp_item.data(enums.ModelEnum.TYPE_ROLE)
          == enums.ModelTypeEnum.MULTIMER_SEQ
      ):
        tmp_number_of_predicted_multimer_sequences += 1

    if (
        tmp_number_of_predicted_monomer_sequences
        == tmp_number_of_monomer_sequences
    ):
      tmp_monomer_sequence_flag = True
    else:
      tmp_monomer_sequence_flag = False
    if (
        tmp_number_of_predicted_multimer_sequences
        == tmp_number_of_multimer_sequences
    ):
      tmp_multimer_sequence_flag = True
    else:
      tmp_multimer_sequence_flag = False
    return tmp_monomer_sequence_flag, tmp_multimer_sequence_flag

  # </editor-fold>

  # <editor-fold desc="Proteins">
  def check_if_scratch_scene_exists_in_protein_model(self) -> bool:
    """Checks if a scratch scene exists in the protein model.

    Returns:
        True if a scratch scene exists, False otherwise.
    """
    return self._protein_model.check_if_scratch_scene_exists(
        self.get_current_protein_tree_index()
    )

  def add_scratch_scene_to_protein_model(self) -> None:
    """Adds a scratch scene to the protein model."""
    self.add_scene_to_proteins_model("_scratch_")

  def add_protein_to_proteins_model(self, a_protein: "protein.Protein") -> None:
    """Add a protein to the proteins model.

    Args:
        a_protein (protein.Protein): An instance of the Protein class representing the protein to be added.

    Raises:
        exception.IllegalArgumentError: If `a_protein` is None.
    """
    # <editor-fold desc="Checks">
    if a_protein is None:
      logger.error("a_protein is None.")
      raise exception.IllegalArgumentError("a_protein is None.")

    # </editor-fold>

    self._protein_model.add_protein(a_protein)

  def remove_protein_from_proteins_model(self) -> None:
    """Removes the current protein from the protein model."""
    self._protein_model.remove_protein(self.get_current_protein_tree_index())

  def remove_non_protein_chain_from_protein(self) -> None:
    """Removes the non-protein chains from protein."""
    tmp_chains_index = self.get_current_protein_tree_index().child(1, 0)
    tmp_non_protein_chain_indexes: collections.deque = collections.deque()
    for i in range(self._protein_model.rowCount(tmp_chains_index)):
      if tmp_chains_index.child(i, 0).data(enums.ModelEnum.OBJECT_ROLE).chain_type == enums.ChainTypeEnum.NON_PROTEIN_CHAIN.value:
        tmp_non_protein_chain_indexes.append(tmp_chains_index.child(i, 0))
    tmp_non_protein_chain_indexes.reverse()
    for tmp_non_protein_chain_index in tmp_non_protein_chain_indexes:
      self._protein_model.remove_chain_of_protein(tmp_non_protein_chain_index)

  # <editor-fold desc="Menu Options">
  def show_menu_options_with_protein(self) -> None:
    """Enable menu options related to protein analysis, results, image, and hotspots."""
    self._main_view.ui.menuAnalysis.setEnabled(True)
    self._main_view.ui.menuResults.setEnabled(True)
    self._main_view.ui.menuImage.setEnabled(True)
    self._main_view.ui.menuHotspots.setEnabled(True)

    self._main_view.ui.proteins_tree_view.setModel(self._protein_model)
    self._main_view.ui.proteins_tree_view.setHeaderHidden(True)

  def show_menu_options_without_protein(self) -> None:
    """Disables menu options related to protein tasks."""
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
      the_pymol_session_manager: "pymol_session_manager.PymolSessionManager",
  ) -> None:
    """Manages the ui protein tab.

    Args:
        an_object_type (str): A string indicating the type of object ("protein", "scene", "chain", "header").
        is_protein_in_pair (bool): A boolean indicating whether the protein is in a pair.
        the_pymol_session_manager (pymol_session_manager.PymolSessionManager): An instance of the PymolSessionManager class.
    """
    self._main_view.ui.lbl_info_2.hide()

    tmp_is_protein_in_session_flag: bool = (
        the_pymol_session_manager.is_the_current_protein_in_session(
            self.get_current_active_protein_object().get_molecule_object()
        )
    )
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
        self._main_view.ui.lbl_info.setText(
            "Please load the PyMOL session of the selected protein."
        )
      self._main_view.ui.action_protein_regions.setEnabled(False)

    elif an_object_type == "scene":
      self._main_view.ui.btn_save_protein.setEnabled(False)
      self._main_view.ui.btn_open_protein_session.setEnabled(False)
      self.hide_protein_pymol_scene_configuration()

      if (
          tmp_is_protein_in_session_flag
          and the_pymol_session_manager.current_scene_name == "base"
      ):
        self._main_view.ui.btn_delete_protein_scene.setEnabled(False)
        self._main_view.ui.lbl_info.setText("Please select a chain.")
      elif (
          tmp_is_protein_in_session_flag
          and the_pymol_session_manager.current_scene_name != "base"
      ):
        self._main_view.ui.btn_delete_protein_scene.setEnabled(True)
        self._main_view.ui.lbl_info.setText("Please select a chain.")
      else:
        self._main_view.ui.lbl_info.setText(
            "Please load the PyMOL session of the selected protein."
        )
      self._main_view.ui.action_protein_regions.setEnabled(False)

    elif an_object_type == "chain":
      self._main_view.ui.btn_save_protein.setEnabled(False)
      self._main_view.ui.btn_open_protein_session.setEnabled(False)

      if tmp_is_protein_in_session_flag:
        self.show_protein_pymol_scene_configuration()
        self._main_view.ui.action_protein_regions.setEnabled(True)
      else:
        self.hide_protein_pymol_scene_configuration()
        self._main_view.ui.action_protein_regions.setEnabled(False)
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
        self._main_view.ui.lbl_info.setText(
            "Please load the PyMOL session of the selected protein."
        )
      self._main_view.ui.action_protein_regions.setEnabled(False)

    else:
      constants.PYSSA_LOGGER.warning(
          "Unknown object type on proteins tab selected."
      )

    if tmp_is_protein_in_session_flag:
      self._main_view.ui.btn_create_protein_scene.setEnabled(True)
      self._main_view.ui.btn_update_protein_scene.setEnabled(True)
    else:
      self._main_view.ui.btn_create_protein_scene.setEnabled(False)
      self._main_view.ui.btn_update_protein_scene.setEnabled(False)

  def manage_coloring_by_element_option_for_protein_chain(self) -> None:
    """Manages coloring options for protein chain based on element selection."""
    if self.get_protein_repr_toggle_flag() == 1:
      if (
          self._main_view.tg_protein_sticks.toggle_button.isChecked()
          or self._main_view.tg_protein_lines.toggle_button.isChecked()
          or self._main_view.tg_protein_spheres.toggle_button.isChecked()
          or self._main_view.tg_protein_dots.toggle_button.isChecked()
          or self._main_view.tg_protein_mesh.toggle_button.isChecked()
          or self._main_view.tg_protein_surface.toggle_button.isChecked()
      ):
        # self._main_view.ui.btn_protein_color_atoms.setEnabled(True)
        # self._main_view.ui.btn_protein_reset_atoms.setEnabled(True)
        self._main_view.tg_protein_color_atoms.setEnabled(True)
      else:
        # self._main_view.ui.btn_protein_color_atoms.setEnabled(False)
        # self._main_view.ui.btn_protein_reset_atoms.setEnabled(False)
        self._main_view.tg_protein_color_atoms.setEnabled(False)
    else:
      if (
          self._main_view.ui.cb_protein_sticks.isChecked()
          or self._main_view.ui.cb_protein_lines.isChecked()
          or self._main_view.ui.cb_protein_spheres.isChecked()
          or self._main_view.ui.cb_protein_dots.isChecked()
          or self._main_view.ui.cb_protein_mesh.isChecked()
          or self._main_view.ui.cb_protein_surface.isChecked()
      ):
        self._main_view.ui.btn_protein_color_atoms.setEnabled(True)
        self._main_view.ui.btn_protein_reset_atoms.setEnabled(True)
      else:
        self._main_view.ui.btn_protein_color_atoms.setEnabled(False)
        self._main_view.ui.btn_protein_reset_atoms.setEnabled(False)

  def manage_hydrogen_representation_for_protein_chain(self) -> None:
    """Manages hydrogen representation for a protein chain.

    This method is responsible for managing the hydrogen representation for a protein chain in the user interface. It checks the state of the toggle buttons for protein sticks, lines, and spheres. If any of these toggle buttons are checked, the method enables certain buttons and elements in the user interface.
    Otherwise, it disables them.
    """
    if (
        self._main_view.tg_protein_sticks.toggle_button.isChecked()
        or self._main_view.tg_protein_lines.toggle_button.isChecked()
        or self._main_view.tg_protein_spheres.toggle_button.isChecked()
    ):
      # self._main_view.ui.btn_protein_show_hydrogens.setEnabled(True)
      # self._main_view.ui.btn_protein_hide_hydrogens.setEnabled(True)
      # self._main_view.tg_protein_hydrogen_atoms.setEnabled(True)
      pass
    else:
      # self._main_view.ui.btn_protein_show_hydrogens.setEnabled(False)
      # self._main_view.ui.btn_protein_hide_hydrogens.setEnabled(False)
      # self._main_view.tg_protein_hydrogen_atoms.setEnabled(False)
      pass

  def manage_toggle_state_of_protein_repr(self, tmp_repr_state: dict) -> None:
    """Manages the toggle state of protein representations.

    Args:
        tmp_repr_state (dict): A dictionary representing the state of protein representations.

    Raises:
        exception.IllegalArgumentError: If `tmp_repr_state` is None.
    """
    # <editor-fold desc="Checks">
    if tmp_repr_state is None:
      logger.error("tmp_repr_state is None.")
      raise exception.IllegalArgumentError("tmp_repr_state is None.")

    # </editor-fold>

    if tmp_repr_state[enums.PyMOLRepresentation.CARTOON.value] == 0:
      ui_util.set_checked_async(
          self._main_view.tg_protein_cartoon.toggle_button, False
      )
      # self._main_view.tg_protein_cartoon.toggle_button.setChecked(False)
    else:
      ui_util.set_checked_async(
          self._main_view.tg_protein_cartoon.toggle_button, True
      )
      # self._main_view.tg_protein_cartoon.toggle_button.setChecked(True)
    if tmp_repr_state[enums.PyMOLRepresentation.STICKS.value] == 0:
      ui_util.set_checked_async(
          self._main_view.tg_protein_sticks.toggle_button, False
      )
      # self._main_view.tg_protein_sticks.toggle_button.setChecked(False)
    else:
      ui_util.set_checked_async(
          self._main_view.tg_protein_sticks.toggle_button, True
      )
      # self._main_view.tg_protein_sticks.toggle_button.setChecked(True)
    if tmp_repr_state[enums.PyMOLRepresentation.RIBBON.value] == 0:
      ui_util.set_checked_async(
          self._main_view.tg_protein_ribbon.toggle_button, False
      )
      # self._main_view.tg_protein_ribbon.toggle_button.setChecked(False)
    else:
      ui_util.set_checked_async(
          self._main_view.tg_protein_ribbon.toggle_button, True
      )
      # self._main_view.tg_protein_ribbon.toggle_button.setChecked(True)
    if tmp_repr_state[enums.PyMOLRepresentation.LINES.value] == 0:
      ui_util.set_checked_async(
          self._main_view.tg_protein_lines.toggle_button, False
      )
      # self._main_view.tg_protein_lines.toggle_button.setChecked(False)
    else:
      ui_util.set_checked_async(
          self._main_view.tg_protein_lines.toggle_button, True
      )
      # self._main_view.tg_protein_lines.toggle_button.setChecked(True)
    if tmp_repr_state[enums.PyMOLRepresentation.SPHERES.value] == 0:
      ui_util.set_checked_async(
          self._main_view.tg_protein_spheres.toggle_button, False
      )
      # self._main_view.tg_protein_spheres.toggle_button.setChecked(False)
    else:
      ui_util.set_checked_async(
          self._main_view.tg_protein_spheres.toggle_button, True
      )
      # self._main_view.tg_protein_spheres.toggle_button.setChecked(True)
    if tmp_repr_state[enums.PyMOLRepresentation.DOTS.value] == 0:
      ui_util.set_checked_async(
          self._main_view.tg_protein_dots.toggle_button, False
      )
      # self._main_view.tg_protein_dots.toggle_button.setChecked(False)
    else:
      ui_util.set_checked_async(
          self._main_view.tg_protein_dots.toggle_button, True
      )
      # self._main_view.tg_protein_dots.toggle_button.setChecked(True)
    if tmp_repr_state[enums.PyMOLRepresentation.MESH.value] == 0:
      ui_util.set_checked_async(
          self._main_view.tg_protein_mesh.toggle_button, False
      )
      # self._main_view.tg_protein_mesh.toggle_button.setChecked(False)
    else:
      ui_util.set_checked_async(
          self._main_view.tg_protein_mesh.toggle_button, True
      )
      # self._main_view.tg_protein_mesh.toggle_button.setChecked(True)
    if tmp_repr_state[enums.PyMOLRepresentation.SURFACE.value] == 0:
      ui_util.set_checked_async(
          self._main_view.tg_protein_surface.toggle_button, False
      )
      # self._main_view.tg_protein_surface.toggle_button.setChecked(False)
    else:
      ui_util.set_checked_async(
          self._main_view.tg_protein_surface.toggle_button, True
      )
      # self._main_view.tg_protein_surface.toggle_button.setChecked(True)

  def manage_check_state_of_protein_repr(self, tmp_repr_state: dict) -> None:
    """Manages the check state of the protein representations in the user interface.

    Args:
        tmp_repr_state: A dictionary containing the state of protein representations.

    Raises:
        exception.IllegalArgumentError: If `tmp_repr_state` is None.
    """
    # <editor-fold desc="Checks">
    if tmp_repr_state is None:
      logger.error("tmp_repr_state is None.")
      raise exception.IllegalArgumentError("tmp_repr_state is None.")

    # </editor-fold>

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

  def get_current_protein_representation_states(
      self,
  ) -> list[tuple[enums.PyMOLRepresentation, bool]]:
    """Gets the representation toggle states of a protein chain on the Proteins tab."""
    tmp_representation_states: list[tuple[enums.PyMOLRepresentation, bool]] = []
    if self._main_view.tg_protein_cartoon.toggle_button.isChecked():
      tmp_representation_states.append(
          (enums.PyMOLRepresentation.CARTOON, True)
      )
    else:
      tmp_representation_states.append(
          (enums.PyMOLRepresentation.CARTOON, False)
      )
    if self._main_view.tg_protein_sticks.toggle_button.isChecked():
      tmp_representation_states.append((enums.PyMOLRepresentation.STICKS, True))
    else:
      tmp_representation_states.append(
          (enums.PyMOLRepresentation.STICKS, False)
      )
    if self._main_view.tg_protein_ribbon.toggle_button.isChecked():
      tmp_representation_states.append((enums.PyMOLRepresentation.RIBBON, True))
    else:
      tmp_representation_states.append(
          (enums.PyMOLRepresentation.RIBBON, False)
      )
    if self._main_view.tg_protein_lines.toggle_button.isChecked():
      tmp_representation_states.append((enums.PyMOLRepresentation.LINES, True))
    else:
      tmp_representation_states.append((enums.PyMOLRepresentation.LINES, False))
    if self._main_view.tg_protein_spheres.toggle_button.isChecked():
      tmp_representation_states.append(
          (enums.PyMOLRepresentation.SPHERES, True)
      )
    else:
      tmp_representation_states.append(
          (enums.PyMOLRepresentation.SPHERES, False)
      )
    if self._main_view.tg_protein_dots.toggle_button.isChecked():
      tmp_representation_states.append((enums.PyMOLRepresentation.DOTS, True))
    else:
      tmp_representation_states.append((enums.PyMOLRepresentation.DOTS, False))
    if self._main_view.tg_protein_mesh.toggle_button.isChecked():
      tmp_representation_states.append((enums.PyMOLRepresentation.MESH, True))
    else:
      tmp_representation_states.append((enums.PyMOLRepresentation.MESH, False))
    if self._main_view.tg_protein_surface.toggle_button.isChecked():
      tmp_representation_states.append(
          (enums.PyMOLRepresentation.SURFACE, True)
      )
    else:
      tmp_representation_states.append(
          (enums.PyMOLRepresentation.SURFACE, False)
      )
    return tmp_representation_states

  # <editor-fold desc="Scene">
  def add_scene_to_proteins_model(self, a_scene_name: str) -> None:
    """Adds a scene to the proteins model.

    Args:
        a_scene_name (str): The name of the scene to be added to the proteins model.

    Raises:
        exception.IllegalArgumentError: If `a_scene_name` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_scene_name is None or a_scene_name == "":
      logger.error("a_scene_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_scene_name is either None or an empty string."
      )

    # </editor-fold>

    tmp_scene_item = QtGui.QStandardItem(a_scene_name)
    tmp_scene_item.setData("scene", enums.ModelEnum.TYPE_ROLE)
    self._protein_model.add_scene(
        self.get_current_protein_tree_index(),
        tmp_scene_item,
    )

  def show_protein_pymol_scene_configuration(self) -> None:
    """Shows the protein scene configuration in PyMOL."""
    self._main_view.ui.frame_protein_color.show()
    self._main_view.ui.frame_protein_repr.show()
    if (
        self._settings_manager.settings.proteins_tab_use_combobox_for_colors
        == 1
    ):
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
      self._main_view.ui.verticalLayout_15.setSpacing(
          0
      )  # layout for the representation section
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

  def hide_protein_pymol_scene_configuration(self) -> None:
    """Hides the protein scene configuration in PyMOL."""
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

  def remove_scene_from_proteins_model(
      self, the_model_index_of_the_scene: QtCore.QModelIndex
  ) -> None:
    """Removes a scene from the proteins model.

    Args:
        the_model_index_of_the_scene (QtCore.QModelIndex): The index of the scene to be removed.

    Raises:
        exception.IllegalArgumentError: If `the_model_index_of_the_scene` is None.
    """
    # <editor-fold desc="Checks">
    if the_model_index_of_the_scene is None:
      logger.error("the_model_index_of_the_scene is None.")
      raise exception.IllegalArgumentError(
          "the_model_index_of_the_scene is None."
      )

    # </editor-fold>

    self._protein_model.remove_scene(the_model_index_of_the_scene)

  # </editor-fold>

  # </editor-fold>

  # <editor-fold desc="Protein Pairs">
  def check_if_scratch_scene_exists_in_protein_pair_model(self) -> bool:
    """Checks if a scratch scene exists in the protein pair model.

    Returns:
        True if a scratch scene exists, False otherwise.
    """
    return self._protein_pair_model.check_if_scratch_scene_exists(
        self.get_current_protein_pair_tree_index()
    )

  def add_scratch_scene_to_protein_pair_model(self) -> None:
    """Adds a scratch scene to the protein pair model."""
    self.add_scene_to_protein_pairs_model("_scratch_")

  def add_protein_pair_to_protein_pairs_model(
      self, a_protein_pair: "protein_pair.ProteinPair"
  ) -> None:
    """Adds a protein pair to the protein pair model.

    Args:
        a_protein_pair: The protein pair object to be added to the protein pairs model.

    Raises:
        exception.IllegalArgumentError: If `a_protein_pair` is None.
    """
    # <editor-fold desc="Checks">
    if a_protein_pair is None:
      logger.error("a_protein_pair is None.")
      raise exception.IllegalArgumentError("a_protein_pair is None.")

    # </editor-fold>

    tmp_main_socket, the_general_purpose_socket = (
        self.job_manager.get_general_purpose_socket_pair()
    )
    self._protein_pair_model.add_protein_pair(
        a_protein_pair, tmp_main_socket, the_general_purpose_socket
    )

  def remove_protein_pair_from_protein_pairs_model(self) -> None:
    """Removes the current protein pair from the protein pairs model."""
    self._protein_pair_model.remove_protein_pair(
        self.get_current_protein_pair_tree_index()
    )

  # <editor-fold desc="Menu Options">
  def show_menu_options_with_protein_pair(self) -> None:
    """Sets the model for the protein pairs tree view in the main UI to display the menu options."""
    self._main_view.ui.protein_pairs_tree_view.setModel(
        self._protein_pair_model
    )
    self._main_view.ui.protein_pairs_tree_view.setHeaderHidden(True)

  def show_menu_options_without_protein_pair(self) -> None:
    """Disables the button options related to protein pairs."""
    self._main_view.ui.btn_delete_protein_pair.setEnabled(False)
    self._main_view.ui.btn_open_protein_pair_session.setEnabled(False)
    self._main_view.ui.btn_create_protein_pair_scene.setEnabled(False)
    self._main_view.ui.btn_update_protein_pair_scene.setEnabled(False)

  # </editor-fold>

  def manage_ui_of_protein_pairs_tab(
      self,
      an_object_type: str,
      the_pymol_session_manager: "pymol_session_manager.PymolSessionManager",
  ) -> None:
    """Manages the user interface of the protein pairs tab.

    Args:
        an_object_type (str): A string representing the type of object.
        the_pymol_session_manager (pymol_session_manager.PymolSessionManager): An instance of pymol_session_manager.PymolSessionManager.

    Raises:
        exception.IllegalArgumentError: If an_object_type is either None or an empty string or if the_pymol_session_manager is None.
    """
    # <editor-fold desc="Checks">
    if an_object_type is None or an_object_type == "":
      logger.error("an_object_type is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "an_object_type is either None or an empty string."
      )
    if the_pymol_session_manager is None:
      logger.error("the_pymol_session_manager is None.")
      raise exception.IllegalArgumentError("the_pymol_session_manager is None.")

    # </editor-fold>

    self._main_view.ui.lbl_info_4.hide()

    tmp_is_protein_pair_in_session_flag: bool = (
        the_pymol_session_manager.is_the_current_protein_pair_in_session(
            self.get_current_active_protein_pair_object().name
        )
    )
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
        self._main_view.ui.lbl_info_3.setText(
            "Please load the PyMOL session of the\nselected protein pair."
        )
      self._main_view.ui.action_protein_regions.setEnabled(False)

    elif an_object_type == "protein":
      self._main_view.ui.btn_delete_protein_pair.setEnabled(False)
      self._main_view.ui.btn_open_protein_pair_session.setEnabled(False)
      self.hide_protein_pair_pymol_scene_configuration()

      if tmp_is_protein_pair_in_session_flag and tmp_current_scene_name == "":
        self._main_view.ui.lbl_info_3.setText("Please select a scene.")
      elif tmp_is_protein_pair_in_session_flag and tmp_current_scene_name != "":
        self._main_view.ui.lbl_info_3.setText("Please select a chain.")
      else:
        self._main_view.ui.lbl_info_3.setText(
            "Please load the PyMOL session of the\nselected protein pair."
        )
      self._main_view.ui.action_protein_regions.setEnabled(False)

    elif an_object_type == "scene":
      self._main_view.ui.btn_delete_protein_pair.setEnabled(False)
      self._main_view.ui.btn_open_protein_pair_session.setEnabled(False)
      self.hide_protein_pair_pymol_scene_configuration()

      if (
          tmp_is_protein_pair_in_session_flag
          and tmp_current_scene_name == "base"
      ):
        self._main_view.ui.btn_delete_protein_pair_scene.setEnabled(False)
        self._main_view.ui.lbl_info_3.setText("Please select a chain.")
      elif (
          tmp_is_protein_pair_in_session_flag
          and tmp_current_scene_name != "base"
      ):
        self._main_view.ui.btn_delete_protein_pair_scene.setEnabled(True)
        self._main_view.ui.lbl_info_3.setText("Please select a chain.")
      else:
        self._main_view.ui.lbl_info_3.setText(
            "Please load the PyMOL session of the \nselected protein pair."
        )

      if tmp_is_protein_pair_in_session_flag:
        self._main_view.ui.lbl_info_3.setText("Please select a chain.")
      else:
        self._main_view.ui.lbl_info_3.setText(
            "Please load the PyMOL session of the \nselected protein pair."
        )
      self._main_view.ui.action_protein_regions.setEnabled(False)

    elif an_object_type == "chain":
      self._main_view.ui.btn_delete_protein_pair.setEnabled(False)
      self._main_view.ui.btn_open_protein_pair_session.setEnabled(False)
      if tmp_is_protein_pair_in_session_flag:
        self.show_protein_pair_pymol_scene_configuration()
        self._main_view.ui.action_protein_regions.setEnabled(True)
      else:
        self.hide_protein_pair_pymol_scene_configuration()
        self._main_view.ui.action_protein_regions.setEnabled(False)

    elif an_object_type == "header":
      self._main_view.ui.btn_delete_protein_pair.setEnabled(False)
      self._main_view.ui.btn_open_protein_pair_session.setEnabled(False)
      self.hide_protein_pair_pymol_scene_configuration()
      if tmp_is_protein_pair_in_session_flag and tmp_current_scene_name == "":
        self._main_view.ui.lbl_info_3.setText("Please select a scene.")
      elif tmp_is_protein_pair_in_session_flag and tmp_current_scene_name != "":
        self._main_view.ui.lbl_info_3.setText("Please select a chain.")
      else:
        self._main_view.ui.lbl_info_3.setText(
            "Please load the PyMOL session of the \nselected protein pair."
        )
      self._main_view.ui.action_protein_regions.setEnabled(False)
    else:
      constants.PYSSA_LOGGER.warning(
          "Unknown object type on proteins tab selected."
      )
      self._main_view.ui.action_protein_regions.setEnabled(False)

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
    else:
      self._main_view.ui.btn_create_protein_pair_scene.setEnabled(False)
      self._main_view.ui.btn_update_protein_pair_scene.setEnabled(False)
      self._main_view.ui.btn_delete_protein_pair_scene.setEnabled(False)

  def manage_coloring_by_element_option_for_protein_pair_chain(self) -> None:
    """Manages coloring by element option for protein pair chain."""
    if self.get_protein_pair_repr_toggle_flag() == 1:
      if (
          self._main_view.tg_protein_pair_sticks.toggle_button.isChecked()
          or self._main_view.tg_protein_pair_lines.toggle_button.isChecked()
          or self._main_view.tg_protein_pair_spheres.toggle_button.isChecked()
          or self._main_view.tg_protein_pair_dots.toggle_button.isChecked()
          or self._main_view.tg_protein_pair_mesh.toggle_button.isChecked()
          or self._main_view.tg_protein_pair_surface.toggle_button.isChecked()
      ):
        # self._main_view.ui.btn_protein_pair_color_atoms.setEnabled(True)
        # self._main_view.ui.btn_protein_pair_reset_atoms.setEnabled(True)
        self._main_view.tg_protein_pair_color_atoms.setEnabled(True)
      else:
        # self._main_view.ui.btn_protein_pair_color_atoms.setEnabled(False)
        # self._main_view.ui.btn_protein_pair_reset_atoms.setEnabled(False)
        self._main_view.tg_protein_pair_color_atoms.setEnabled(False)
    else:
      if (
          self._main_view.ui.cb_protein_pair_sticks.isChecked()
          or self._main_view.ui.cb_protein_pair_lines.isChecked()
          or self._main_view.ui.cb_protein_pair_spheres.isChecked()
          or self._main_view.ui.cb_protein_pair_dots.isChecked()
          or self._main_view.ui.cb_protein_pair_mesh.isChecked()
          or self._main_view.ui.cb_protein_pair_surface.isChecked()
      ):
        self._main_view.ui.btn_protein_pair_color_atoms.setEnabled(True)
        self._main_view.ui.btn_protein_pair_reset_atoms.setEnabled(True)
      else:
        self._main_view.ui.btn_protein_pair_color_atoms.setEnabled(False)
        self._main_view.ui.btn_protein_pair_reset_atoms.setEnabled(False)

  def manage_hydrogen_representation_for_protein_pair_chain(self) -> None:
    """Manages the representation of hydrogen atoms for a protein pair chain."""
    if (
        self._main_view.tg_protein_pair_sticks.toggle_button.isChecked()
        or self._main_view.tg_protein_pair_lines.toggle_button.isChecked()
        or self._main_view.tg_protein_pair_spheres.toggle_button.isChecked()
    ):
      # self._main_view.ui.btn_protein_pair_show_hydrogens.setEnabled(True)
      # self._main_view.ui.btn_protein_pair_hide_hydrogens.setEnabled(True)
      # self._main_view.tg_protein_pair_hydrogen_atoms.setEnabled(True)
      pass
    else:
      # self._main_view.ui.btn_protein_pair_show_hydrogens.setEnabled(False)
      # self._main_view.ui.btn_protein_pair_hide_hydrogens.setEnabled(False)
      # self._main_view.tg_protein_pair_hydrogen_atoms.setEnabled(False)
      pass

  def set_current_active_chain_color_of_protein_pair(
      self, a_color: str
  ) -> None:
    """Sets the current active chain color for a protein pair.

    Args:
        a_color: The color to set as the current active chain color.

    Raises:
        exception.IllegalArgumentError: If `a_color` is None or an empty string.
        ValueError: If the type of the current index in the protein_pairs_tree_view is not one of the expected types.
    """
    # <editor-fold desc="Checks">
    if a_color is None or a_color == "":
      logger.error("a_color is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_color is either None or an empty string."
      )

    # </editor-fold>

    tmp_type = self._main_view.ui.protein_pairs_tree_view.currentIndex().data(
        enums.ModelEnum.TYPE_ROLE
    )
    if tmp_type == "protein_pair":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "protein":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "header":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "scene":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "chain":
      self._main_view.ui.protein_pairs_tree_view.model().setData(
          self._main_view.ui.protein_pairs_tree_view.currentIndex(),
          a_color,
          enums.ModelEnum.CHAIN_COLOR_ROLE,
      )
    else:
      raise ValueError("Unknown type!")

  def get_current_active_chain_color_of_protein_pair(self) -> str:
    """Returns the chain color of the current active branch.

    Returns:
        A chain color.
    """
    tmp_type = self._main_view.ui.protein_pairs_tree_view.currentIndex().data(
        enums.ModelEnum.TYPE_ROLE
    )
    if tmp_type == "protein_pair":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "protein":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "header":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "scene":
      raise ValueError(f"Cannot get a chain object if the type is: {tmp_type}!")
    elif tmp_type == "chain":
      return self._main_view.ui.protein_pairs_tree_view.currentIndex().data(
          enums.ModelEnum.CHAIN_COLOR_ROLE
      )
    else:
      raise ValueError("Unknown type!")

  def get_current_protein_pair_representation_states(
      self,
  ) -> list[tuple[enums.PyMOLRepresentation, bool]]:
    """Gets the representation toggle states of a protein chain on the Protein Pairs tab.

    Returns:
        A list of tuples representing the current representation states of the protein pair. Each tuple
        contains a PyMOLRepresentation enum value and a boolean indicating whether the representation
        is enabled or disabled.
    """
    tmp_representation_states: list[tuple[enums.PyMOLRepresentation, bool]] = []
    if self._main_view.tg_protein_pair_cartoon.toggle_button.isChecked():
      tmp_representation_states.append(
          (enums.PyMOLRepresentation.CARTOON, True)
      )
    else:
      tmp_representation_states.append(
          (enums.PyMOLRepresentation.CARTOON, False)
      )
    if self._main_view.tg_protein_pair_sticks.toggle_button.isChecked():
      tmp_representation_states.append((enums.PyMOLRepresentation.STICKS, True))
    else:
      tmp_representation_states.append(
          (enums.PyMOLRepresentation.STICKS, False)
      )
    if self._main_view.tg_protein_pair_ribbon.toggle_button.isChecked():
      tmp_representation_states.append((enums.PyMOLRepresentation.RIBBON, True))
    else:
      tmp_representation_states.append(
          (enums.PyMOLRepresentation.RIBBON, False)
      )
    if self._main_view.tg_protein_pair_lines.toggle_button.isChecked():
      tmp_representation_states.append((enums.PyMOLRepresentation.LINES, True))
    else:
      tmp_representation_states.append((enums.PyMOLRepresentation.LINES, False))
    if self._main_view.tg_protein_pair_spheres.toggle_button.isChecked():
      tmp_representation_states.append(
          (enums.PyMOLRepresentation.SPHERES, True)
      )
    else:
      tmp_representation_states.append(
          (enums.PyMOLRepresentation.SPHERES, False)
      )
    if self._main_view.tg_protein_pair_dots.toggle_button.isChecked():
      tmp_representation_states.append((enums.PyMOLRepresentation.DOTS, True))
    else:
      tmp_representation_states.append((enums.PyMOLRepresentation.DOTS, False))
    if self._main_view.tg_protein_pair_mesh.toggle_button.isChecked():
      tmp_representation_states.append((enums.PyMOLRepresentation.MESH, True))
    else:
      tmp_representation_states.append((enums.PyMOLRepresentation.MESH, False))
    if self._main_view.tg_protein_pair_surface.toggle_button.isChecked():
      tmp_representation_states.append(
          (enums.PyMOLRepresentation.SURFACE, True)
      )
    else:
      tmp_representation_states.append(
          (enums.PyMOLRepresentation.SURFACE, False)
      )
    return tmp_representation_states

  # <editor-fold desc="Scene">
  def add_scene_to_protein_pairs_model(self, a_scene_name: str) -> None:
    """Adds a scene to the protein pairs model.

    Args:
        a_scene_name (str): The name of the scene to be added.

    Raises:
        exception.IllegalArgumentError: if a_scene_name is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_scene_name is None or a_scene_name == "":
      logger.error("a_scene_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_scene_name is either None or an empty string."
      )

    # </editor-fold>

    tmp_scene_item = QtGui.QStandardItem(a_scene_name)
    tmp_scene_item.setData("scene", enums.ModelEnum.TYPE_ROLE)
    self._protein_pair_model.add_scene(
        self.get_current_protein_pair_tree_index(),
        tmp_scene_item,
    )

  def show_protein_pair_pymol_scene_configuration(self) -> None:
    """Shows the configuration of the protein pair scene in PyMOL."""
    self._main_view.ui.frame_protein_pair_color.show()
    self._main_view.ui.frame_protein_pair_repr.show()
    if (
        self._settings_manager.settings.protein_pairs_tab_use_combobox_for_colors
        == 1
    ):
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
      self._main_view.ui.verticalLayout_18.setSpacing(
          0
      )  # layout for the representation section
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

  def hide_protein_pair_pymol_scene_configuration(self) -> None:
    """Hides Protein Pair PyMOL Scene Configuration."""
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

  def remove_scene_from_protein_pairs_model(
      self, the_model_index_of_the_scene: QtCore.QModelIndex
  ) -> None:
    """Removes the scene from the protein pairs model.

    Args:
        the_model_index_of_the_scene (QtCore.QModelIndex): The index of the scene to be removed from the protein pairs model.

    Raises:
        exception.IllegalArgumentError: If `the_model_index_of_the_scene` is None.
    """
    # <editor-fold desc="Checks">
    if the_model_index_of_the_scene is None:
      logger.error("the_model_index_of_the_scene is None.")
      raise exception.IllegalArgumentError(
          "the_model_index_of_the_scene is None."
      )

    # </editor-fold>

    self._protein_pair_model.remove_scene(the_model_index_of_the_scene)

  # </editor-fold>

  # </editor-fold>

  def update_settings(self) -> None:
    """Deserializes the settings json file."""
    self._settings_manager.settings = (
        self._settings_manager.settings.deserialize_settings()
    )

  def restore_default_main_view(self) -> None:
    """Restores the default main view by clearing the sequences table and setting its row count to 0."""
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

  def block_gui(self, with_wait_cursor: bool = False) -> None:
    """Starts the wait cursor.

    Args:
        with_wait_cursor (bool): A boolean indicating whether to display a wait cursor.

    Raises:
        exception.IllegalArgumentError: If `with_wait_cursor` is None.
    """
    # <editor-fold desc="Checks">
    if with_wait_cursor is None:
      logger.error("with_wait_cursor is None.")
      raise exception.IllegalArgumentError("with_wait_cursor is None.")

    # </editor-fold>

    if with_wait_cursor is True:
      QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
    self._main_view.disable_menu_bar_without_exit_application()
    self._main_view.disable_tab_widget()
    self._main_view.disable_job_panels()

  def stop_wait_cursor(self) -> None:
    """Stops the cursor."""
    QtWidgets.QApplication.restoreOverrideCursor()
    self._main_view.enable_job_panels()

  # <editor-fold desc="Job related methods">
  def add_job_entry_to_job_overview_layout(
      self, a_job_entry_widget: "job_entry.JobEntryWidget"
  ) -> None:
    """Adds a job entry widget to the job overview layout.

    Args:
        a_job_entry_widget (job_entry.JobEntryWidget): The job entry widget to add.

    Raises:
        exception.IllegalArgumentError: If `a_job_entry_widget` is None.
    """
    # <editor-fold desc="Checks">
    if a_job_entry_widget is None:
      logger.error("a_job_entry_widget is None.")
      raise exception.IllegalArgumentError("a_job_entry_widget is None.")

    # </editor-fold>

    self.job_entry_widgets.append(a_job_entry_widget)
    self._main_view.ui.job_overview_layout.insertWidget(
        self._main_view.ui.job_overview_layout.count() - 1, a_job_entry_widget
    )
    self._main_view.lbl_job_overview.hide()
    self._main_view.btn_open_job_overview.setIcon(
        self._main_view.icon_jobs_running
    )

  def update_job_entry(self, update_job_entry_signal_values: tuple) -> None:
    """Updates a job entry with the given signal values.

    Args:
        update_job_entry_signal_values (tuple): A tuple containing information about the job entry update.
            The tuple should have the format (a_job_entry_widget, a_description, a_value), where:
                - a_job_entry_widget: A reference to the job entry widget to update.
                - a_description: The description of the job entry update.
                - a_value: The value of the job entry update.

    Notes:
        Gets signal from the classes of job.py with the signal "update_job_entry_signal"
    """
    # <editor-fold desc="Checks">
    if update_job_entry_signal_values is None:
      logger.error("update_job_entry_signal_values is None.")
      raise exception.IllegalArgumentError(
          "update_job_entry_signal_values is None."
      )

    # </editor-fold>

    _, a_description, a_value = update_job_entry_signal_values
    a_job_entry_widget: "job_entry.JobEntryWidget" = (
        update_job_entry_signal_values[0]
    )  # Unpack of 0 because of type annotation need
    a_job_entry_widget.ui.progress_bar_job.setValue(a_value)
    if a_value == 100:
      tmp_type = a_job_entry_widget.job_base_information.job_type
      tmp_progress = a_job_entry_widget.job_base_information.job_progress
      if (
          tmp_type == enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS
          and tmp_progress == enums.JobProgress.RUNNING
      ):
        return
      self._create_job_notification(a_job_entry_widget, a_description)
    elif a_value == 0 or a_value == 25:
      a_job_entry_widget.ui.btn_cancel_job.setEnabled(True)
    else:
      a_job_entry_widget.ui.btn_cancel_job.setEnabled(False)

  def _create_job_notification(
      self, a_job_entry_widget: "job_entry.JobEntryWidget", a_description: str
  ) -> None:
    """Creates job notification widget after a job finished.

    Args:
        a_job_entry_widget (job_entry.JobEntryWidget): The job entry widget that triggered the job notification creation.
        a_description (str): The description of the job notification.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None or if `a_description` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_job_entry_widget is None:
      logger.error("a_job_entry_widget is None.")
      raise exception.IllegalArgumentError("a_job_entry_widget is None.")
    if a_description is None:
      logger.error("a_description is None.")
      raise exception.IllegalArgumentError("a_description is None.")

    # </editor-fold>

    if (
        a_job_entry_widget.job_base_information.project_name
        == self._current_project.get_project_name()
    ):
      tmp_job_is_from_current_project = True
    else:
      tmp_job_is_from_current_project = False

    tmp_job_notification_widget = job_entry.JobNotificationWidget(
        a_description,
        a_job_entry_widget.job_base_information,
        tmp_job_is_from_current_project,
        self.refresh_after_job_finished_signal,
    )
    # remove job entry widget
    if a_job_entry_widget.parent():
      a_job_entry_widget.setParent(None)
    a_job_entry_widget.deleteLater()
    # Add new notification widget to notification panel
    self._main_view.ui.job_notification_layout.insertWidget(
        self._main_view.ui.job_notification_layout.count() - 1,
        tmp_job_notification_widget,
    )
    self._main_view.btn_open_job_notification.setIcon(
        self._main_view.icon_notify_unread
    )
    self._main_view.lbl_job_notification.hide()
    if self._main_view.ui.job_overview_layout.count() == 2:
      self._main_view.lbl_job_overview.show()
      self._main_view.btn_open_job_overview.setIcon(self._main_view.icon_jobs)

  def remove_job_notification_widget(
      self, a_job_notification_widget: "job_entry.JobNotificationWidget"
  ) -> None:
    """Removes a job notification after pressing open or refresh.

    Args:
        a_job_notification_widget (job_entry.JobNotificationWidget): The job notification widget to be removed.

    Raises:
        exception.IllegalArgumentError: If `a_job_notification_widget` is None.

    Notes:
        This function is called from the main_view_controller.
    """
    # <editor-fold desc="Checks">
    if a_job_notification_widget is None:
      logger.error("a_job_notification_widget is None.")
      raise exception.IllegalArgumentError("a_job_notification_widget is None.")

    # </editor-fold>

    logger.debug(a_job_notification_widget)
    self._main_view.ui.job_notification_layout.removeWidget(
        a_job_notification_widget
    )
    if self._main_view.ui.job_notification_layout.count() == 2:
      self._main_view.lbl_job_notification.show()
      self._main_view.btn_open_job_notification.setIcon(
          self._main_view.icon_notify
      )

  def close_job_notification_panel(self) -> None:
    """Closes the job notification panel after pressing open or refresh.

    Notes:
        This function is called from the main_view_controller.
    """
    self._main_view.ui.frame_job_notification.hide()

  def close_job_overview_panel(self) -> None:
    """Closes the job overview panel."""
    self._main_view.ui.frame_job_overview.hide()

  def cancel_job(self, signal_tuple: tuple) -> None:
    """Cancels a job.

    Args:
        signal_tuple (tuple): A tuple that contains information about the job.
            signal_tuple[0] (enums.JobType): The type of job.
            signal_tuple[1] (QWidget): The widget where the job entry is displayed.
            signal_tuple[2] (Job): The job object.
            signal_tuple[3] (list): A list of protein prediction info objects.
    """
    # <editor-fold desc="Checks">
    if signal_tuple is None:
      logger.error("signal_tuple is None.")
      raise exception.IllegalArgumentError("signal_tuple is None.")

    # </editor-fold>

    if signal_tuple[0] == enums.JobType.PREDICTION:
      # signal_tuple = (job_type, job_entry_widget, job object, protein_prediction_infos)
      for tmp_protein_info in signal_tuple[3]:
        self.watcher.remove_protein(tmp_protein_info.name)
      if signal_tuple[1].parent():
        signal_tuple[1].setParent(None)
      signal_tuple[1].deleteLater()
      if self.job_manager.is_prediction_job_currently_running(signal_tuple[2]):
        constants.PYSSA_LOGGER.info(
            "Structure prediction process was aborted manually."
        )
        subprocess.run(["wsl", "--shutdown"])
        constants.PYSSA_LOGGER.info("Shutdown of wsl environment.")
        filesystem_io.FilesystemCleaner.clean_prediction_scratch_folder()
        constants.PYSSA_LOGGER.info("Cleaned scratch directory.")
        self.job_manager.stop_prediction_queue()
      else:
        self.job_manager.pop_job_from_queue(signal_tuple[2])

  # </editor-fold>
