import datetime
import logging
import pathlib
from typing import Optional, Callable

from PyQt6 import QtGui
from PyQt6 import QtCore

from pyssa.model.data_classes import workspace_project
from pyssa.model.preference import model_definitions
from pyssa.model.util import exception
from pyssa.model.qmodel import base_tree_model
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)

__docformat__ = "google"


class ProjectModel(base_tree_model.BaseTreeModel):
  """Class for storing projects of the active workspace."""

  def __init__(self):
    """Constructor."""
    super().__init__()

  @staticmethod
  def from_workspace_path(a_workspace_path: pathlib.Path, a_callable: Callable) -> "ProjectModel":
    """Alternative constructor.

    Raises:
      exception.NoneValueError: If `a_workspace_path` is None.
      IllegalArgumentError: If `a_workspace_path` is an empty string or could not be found.
    """
    # <editor-fold desc="Checks">
    if a_workspace_path is None:
      default_logging.append_to_log_file(logger, "a_workspace_path is None.", logging.ERROR)
      raise exception.NoneValueError("a_workspace_path is None.")
    if a_workspace_path == "":
      default_logging.append_to_log_file(logger, "a_workspace_path is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_workspace_path is an empty string.")
    if not a_workspace_path.exists():
      default_logging.append_to_log_file(logger, "a_workspace_path could not be found.", logging.ERROR)
      raise exception.IllegalArgumentError("a_workspace_path could not be found.")
    # </editor-fold>
    tmp_project_model = ProjectModel()
    tmp_project_model.create_root_node()

    tmp_workspace_projects = []
    for tmp_file in a_workspace_path.iterdir():
      if tmp_file.suffix == ".db":
        tmp_database_filepath = pathlib.Path(f"{a_workspace_path}/{tmp_file.name}")
        tmp_workspace_projects.append(
          workspace_project.WorkspaceProject(
            tmp_database_filepath,
            datetime.datetime.fromtimestamp(tmp_database_filepath.stat().st_mtime),
            a_callable
          )
        )
        """
        As object the filepath is used because this will save memory.
        Another idea would be to load the entire project object, but this 
        could lead to an OOM error because all projects of a given workspace
        where then loaded into memory
        """
    tmp_workspace_projects.sort(key=lambda x: x.date_modified, reverse=True)
    for tmp_workspace_project in tmp_workspace_projects:
      tmp_project_model.add_node(
        tmp_project_model.root_node,
        tmp_workspace_project.filepath.name.replace(".db", ""),
        model_definitions.TypesEnum.PROJECT_TYPE,
        tmp_workspace_project,  # Object role
      )
    return tmp_project_model
