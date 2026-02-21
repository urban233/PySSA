"""
AppState: holds all mutable application-level state, including the
currently open (hot) project and any cold projects with running jobs.

The MainWindowController owns one AppState instance and reads from it
to drive refresh_ui().  AppState itself has no Qt dependency and no
knowledge of the UI — it is a plain Python object so it can be tested
without a running QApplication.

Hot project:  The project the user currently has open in the UI.
              At most one at any time.

Cold project: A project not open in the UI but with a background job
              still writing results to its database.  Multiple allowed.
"""
from __future__ import annotations

import glob
import logging
import os
import pathlib
from typing import Any, Callable, TYPE_CHECKING

from src.pyssa.gui.qt import QtGui
from src.pyssa.controller import settings_manager
from src.pyssa.controller.job_scheduler import JobScheduler
from src.pyssa.internal.data_structures.data_classes import job_descriptor
from src.pyssa.io_pyssa.db_pyssa import ProjectDatabase
from src.pyssa.io_pyssa.db_pyssa import ColdProjectHandle
from src.pyssa.internal.data_structures import workspace, settings
from src.pyssa.model import job_model
from src.pyssa.model import psa_objects_model
from src.pyssa.util import enums

if TYPE_CHECKING:
  from src.pyssa.internal.data_structures import project

logger = logging.getLogger(__name__)


class AppState:
  """All mutable application state: open project + database lifecycle.

  Args:
      on_state_changed: A zero-argument callable invoked whenever state
                        changes that should trigger a UI refresh.
                        MainWindowController passes self.refresh_ui here.
                        Defaults to a no-op so AppState can be used
                        standalone in tests.
  """

  def __init__(
          self,
          a_settings_manager: "settings_manager.SettingsManager",
          on_state_changed: Callable[[], None] = lambda: None,
  ) -> None:
    self._settings_manager = a_settings_manager
    self._on_state_changed = on_state_changed

    self._workspace: "workspace.Workspace" = workspace.Workspace(
      self._settings_manager.settings.workspace_path
    )

    self._project: "project.Project | None" = None
    self._first_pass: bool = False
    """
    A boolean flag that is used to indicate whether the refresh_ui method
    should set the PySSAObjectsModel to the panel because it would be the 
    'first pass' of the refresh_ui method since the loading of the project.
    """
    self._hot_db: ProjectDatabase | None = None
    self._cold_dbs: dict[str, ProjectDatabase] = {}

    self._pyssa_objects_model = psa_objects_model.PSAObjectsModel()

    self._job_model = job_model.JobModel()
    self._job_scheduler = JobScheduler(self._job_model)
    self._job_model.job_finished.connect(self._on_job_finished)

    self._build_workspace_model()

  @property
  def workspace(self) -> "workspace.Workspace":
    return self._workspace

  def _build_workspace_model(self) -> None:
    """Builds the workspace model.

    This method populates the workspace model with project items based on the database files found in the workspace directory.
    """
    db_pattern = os.path.join(self._workspace.path, "*.db")
    tmp_root_item = self._workspace.get_model().invisibleRootItem()
    for tmp_filename in [
      os.path.basename(file).replace(".db", "")
      for file in glob.glob(db_pattern)
    ]:
      tmp_project_item = QtGui.QStandardItem(tmp_filename)
      # tmp_filepath = pathlib.Path(f"{self._workspace.path}/{tmp_filename}.db")
      tmp_filepath = self._workspace.construct_project_db_path(tmp_filename)
      tmp_project_item.setData(tmp_filepath, enums.ModelEnum.FILEPATH_ROLE)
      tmp_root_item.appendRow(tmp_project_item)

  # ------------------------------------------------------------------
  # Read-only project access
  # ------------------------------------------------------------------

  @property
  def project(self) -> "project.Project | None":
    return self._project

  def has_open_project(self) -> bool:
    return self._project is not None

  def require_project(self) -> "project.Project":
    if self._project is None:
      raise RuntimeError("No project is currently open.")
    return self._project

  # ------------------------------------------------------------------
  # Hot project
  # ------------------------------------------------------------------

  @property
  def hot_db(self) -> ProjectDatabase | None:
    return self._hot_db

  def open_project(
          self,
          project: "project.Project",
          db: ProjectDatabase,
          on_result: Callable[[Any], None] | None = None,
  ) -> None:
    """Store a newly loaded project and its database.

    Closes any previously open hot project first.  If the project being
    opened has cold background jobs (i.e. it was closed while jobs were
    running), those jobs are transitioned back to hot so their results
    update the UI when they finish.

    Args:
        project: The project domain object.
        db: The ``ProjectDatabase`` for this project.
        on_result: Optional callback applied to every transitioned
                   cold→hot job.  Invoked on the main thread when
                   the job finishes.  If ``None``, the jobs will
                   still emit ``JobModel.job_finished`` but no
                   direct callback will fire.

    Notifies the UI via on_state_changed.
    """
    self._close_hot_db()
    self._hot_db = db
    self._project = project

    project_name = project.get_project_name()

    transitioned = self._job_scheduler.transition_to_hot(
      project_name, on_result=on_result,
    )
    if transitioned > 0:
      cold_db = self._cold_dbs.pop(project_name, None)
      if cold_db is not None:
        cold_db.close()
      logger.info(
        "Transitioned %d cold job(s) to hot for project '%s'.",
        transitioned, project_name,
      )

    logger.info("Hot project opened: '%s'.", project_name)
    self._first_pass = True
    self._on_state_changed()
    self._first_pass = False

  def close_project(self) -> None:
    """Close the hot project and notify the UI.

    If the scheduler still has running or queued jobs for this project,
    they are transitioned to cold first so that their results are
    persisted to the database instead of lost.
    """
    if self._project is not None and self._job_scheduler.has_running_jobs():
      project_name = self._project.get_project_name()
      db_path = str(
        self._workspace.construct_project_db_path(project_name),
      )
      cold_handle = self.open_cold_project(db_path, project_name)
      self._job_scheduler.transition_to_cold(project_name, cold_handle)

    self._close_hot_db()
    self._project = None
    logger.info("Hot project closed.")
    self._on_state_changed()

  def _close_hot_db(self) -> None:
    if self._hot_db is not None:
      self._hot_db.close()      # drains write queue internally
      self._hot_db = None

  # ------------------------------------------------------------------
  # Cold projects
  # ------------------------------------------------------------------

  def open_cold_project(self, db_path: str, project_id: str) -> ColdProjectHandle:
    """Open a cold project and return a write-only handle for a job.

    Idempotent: if this project_id is already open, a new handle
    wrapping the same database is returned.

    The handle's close() calls _on_cold_closed, which removes the
    entry from _cold_dbs.  No UI notification is needed for close
    (the job's on_success callback handles any UI update).
    """
    if project_id not in self._cold_dbs:
      self._cold_dbs[project_id] = ProjectDatabase(
        db_path=db_path, project_id=project_id
      )
    return ColdProjectHandle(
      project_id=project_id,
      db=self._cold_dbs[project_id],
      on_closed=self._on_cold_closed,
    )

  def _on_cold_closed(self, project_id: str) -> None:
    """Called from the worker thread when a cold job finishes.

    Must not touch Qt objects — only plain dict mutation here.
    """
    self._cold_dbs.pop(project_id, None)
    logger.info("Cold project '%s' removed.", project_id)

  # ------------------------------------------------------------------
  # Application exit
  # ------------------------------------------------------------------

  def close_all(self) -> None:
    """Drain and close every open database.  Call at application exit."""
    self._job_scheduler.shutdown()
    self._close_hot_db()
    for project_id in list(self._cold_dbs):
      db = self._cold_dbs.pop(project_id)
      db.close()

  @property
  def pyssa_objects_model(self) -> "psa_objects_model.PSAObjectsModel":
    return self._pyssa_objects_model

  @pyssa_objects_model.setter
  def pyssa_objects_model(self, value):
    self._pyssa_objects_model = value

  def clear_pyssa_objects_model(self) -> None:
    self._pyssa_objects_model.clear()

  def get_settings(self) -> "settings.Settings":
    return self._settings_manager.settings

  def is_first_pass(self) -> bool:
    return self._first_pass

  # ------------------------------------------------------------------
  # Job system
  # ------------------------------------------------------------------

  @property
  def job_model(self) -> "job_model.JobModel":
    """The application-wide job table model."""
    return self._job_model

  @property
  def job_scheduler(self) -> JobScheduler:
    """The application-wide job scheduler."""
    return self._job_scheduler

  def _on_job_finished(
      self,
      descriptor: "job_descriptor.JobDescriptor",
      result: object,
  ) -> None:
    """Handle a completed hot-project job by invoking its result callback.

    Connected to ``JobModel.job_finished``.  Only acts when the job's
    project is still the hot project.  Does **not** show any popup.

    Args:
        descriptor: The job's immutable metadata.
        result: The return value of the worker function.
    """
    if not descriptor.is_hot:
      return
    if descriptor.on_result is not None:
      try:
        descriptor.on_result(result)
      except Exception:
        logger.exception(
          "on_result callback failed for job '%s'.",
          descriptor.display_name or descriptor.job_type.value,
        )
