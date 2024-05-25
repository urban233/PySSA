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
"""Module for the watcher class."""
import logging
import queue
from typing import Optional

from pyssa.internal.data_structures import job, project
from pyssa.internal.data_structures.data_classes import job_summary
from pyssa.logging_pyssa import log_handlers
from pyssa.util import exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class Watcher:
  """Holds information about the names that cannot be used in the current opened project."""

  def __init__(self) -> None:
    """Constructor."""
    self.project = None
    self.protein_names_blacklist: list[str] = []
    self.protein_pair_names_blacklist: list[str] = []

  # <editor-fold desc="Private methods">
  def _retrieve_all_protein_names_from_project(self) -> list[str]:
    """Retrieve all protein names from project."""
    tmp_all_protein_names = []
    for tmp_protein in self.project.proteins:
      tmp_all_protein_names.append(tmp_protein.get_molecule_object())
    return tmp_all_protein_names

  def _retrieve_all_protein_pair_names_from_project(self) -> list[str]:
    """Retrieve all protein names from project."""
    tmp_all_protein_pair_names = []
    for tmp_protein_pair in self.project.protein_pairs:
      tmp_all_protein_pair_names.append(tmp_protein_pair.name)
    return tmp_all_protein_pair_names

  def _retrieve_all_prediction_names_of_the_current_project(
      self,
      a_prediction_queue: queue.Queue,
      the_current_prediction_job: Optional["job_summary.PredictionJobSummary"],
  ) -> list[str]:
    """Retrieves all prediction names from the prediction queue that belong to the current project.

    Args:
        a_prediction_queue (queue.Queue): A queue containing prediction jobs.
        the_current_prediction_job (Optional[job_summary.PredictionJobSummary]): The current prediction job summary.

    Returns:
        A list of prediction names of the current project.

    Raises:
        exception.IllegalArgumentError: If `a_prediction_queue` is None.
    """
    # <editor-fold desc="Checks">
    if a_prediction_queue is None:
      logger.error("a_prediction_queue is None.")
      raise exception.IllegalArgumentError("a_prediction_queue is None.")

    # </editor-fold>

    elements = []
    while not a_prediction_queue.empty():
      tmp_job: "job.PredictionJob" = a_prediction_queue.get()
      logger.debug(
          f"Watcher project name: {self.project.get_project_name()} vs. Job project name: {tmp_job.frozen_project.get_project_name()}"
      )
      if (
          tmp_job.frozen_project.get_project_name()
          == self.project.get_project_name()
      ):
        for tmp_protein_info in tmp_job.prediction_protein_infos:
          elements.append(tmp_protein_info.name)
    if the_current_prediction_job is not None:
      elements.extend(the_current_prediction_job.get_protein_names())
    return elements

  def _retrieve_all_distance_analysis_names_of_the_current_project(
      self,
      a_distance_analysis_queue: queue.Queue,
      the_current_distance_analysis_job: Optional[
          "job_summary.DistanceAnalysisJobSummary"
      ],
  ) -> list[str]:
    """Retrieves all analysis names from the distance analysis queue that belong to the current project.

    Args:
        a_distance_analysis_queue (queue.Queue): A queue containing distance analysis jobs.
        the_current_distance_analysis_job (Optional[job_summary.DistanceAnalysisJobSummary]): The current distance job.

    Returns:
        A list of analysis names of the current project.

    Raises:
        exception.IllegalArgumentError: If `a_distance_analysis_queue` is None.
    """
    # <editor-fold desc="Checks">
    if a_distance_analysis_queue is None:
      logger.error("a_distance_analysis_queue is None.")
      raise exception.IllegalArgumentError("a_distance_analysis_queue is None.")

    # </editor-fold>

    elements = []
    while not a_distance_analysis_queue.empty():
      tmp_job: "job.DistanceAnalysisJob" = a_distance_analysis_queue.get()
      if (
          tmp_job.frozen_project.get_project_name()
          == self.project.get_project_name()
      ):
        for tmp_analysis_name in tmp_job.list_with_analysis_names:
          elements.append(tmp_analysis_name)
    if the_current_distance_analysis_job is not None:
      elements.extend(the_current_distance_analysis_job.analysis_names)
    return elements

  # </editor-fold>

  # <editor-fold desc="Public methods">
  def setup_blacklists(
      self,
      a_project: "project.Project",
      a_prediction_queue: queue.Queue,
      a_distance_analysis_queue: queue.Queue,
      the_current_prediction_job: Optional["job_summary.PredictionJobSummary"],
      the_current_distance_analysis_job: Optional[
          "job_summary.DistanceAnalysisJobSummary"
      ],
  ) -> None:
    """Sets up the blacklists for proteins and protein pairs.

    Args:
        a_project (project.Project): The project object.
        a_prediction_queue (queue.Queue): The prediction job queue.
        a_distance_analysis_queue (queue.Queue): The distance analysis job queue.
        the_current_prediction_job (Optional[job_summary.PredictionJobSummary]): The current prediction job summary.
        the_current_distance_analysis_job (Optional[job_summary.DistanceAnalysisJobSummary]): The current distance analysis job summary.

    Raises:
        exception.IllegalArgumentError: If any of the arguments except `the_current_prediction_job` and `the_current_distance_analysis_job` are None.
    """
    # <editor-fold desc="Checks">
    if a_project is None:
      logger.error("a_project is None.")
      raise exception.IllegalArgumentError("a_project is None.")
    if a_prediction_queue is None:
      logger.error("a_prediction_queue is None.")
      raise exception.IllegalArgumentError("a_prediction_queue is None.")
    if a_distance_analysis_queue is None:
      logger.error("a_distance_analysis_queue is None.")
      raise exception.IllegalArgumentError("a_distance_analysis_queue is None.")

    # </editor-fold>

    # Reset lists
    self.protein_names_blacklist: list[str] = []
    self.protein_pair_names_blacklist: list[str] = []
    # Setup process
    self.project = a_project
    self.protein_names_blacklist: list[str] = (
        self._retrieve_all_protein_names_from_project()
    )
    self.protein_pair_names_blacklist: list[str] = (
        self._retrieve_all_protein_pair_names_from_project()
    )

    # Extend blacklist with job information
    self.protein_names_blacklist.extend(
        self._retrieve_all_prediction_names_of_the_current_project(
            a_prediction_queue,
            the_current_prediction_job,
        ),
    )
    self.protein_pair_names_blacklist.extend(
        self._retrieve_all_distance_analysis_names_of_the_current_project(
            a_distance_analysis_queue,
            the_current_distance_analysis_job,
        ),
    )
    logger.info(
        f"After setting up the blacklist, the protein blacklist is      {self.protein_names_blacklist}."
    )
    logger.info(
        f"After setting up the blacklist, the protein pair blacklist is {self.protein_pair_names_blacklist}."
    )

  def add_proteins_from_new_job(self, prediction_protein_infos: list) -> None:
    """Adds the protein names from the new job to the project's blacklist.

    Args:
        prediction_protein_infos (list): A list of protein information.

    Raises:
        exception.IllegalArgumentError: If `prediction_protein_infos` is None.
    """
    # <editor-fold desc="Checks">
    if prediction_protein_infos is None:
      logger.error("prediction_protein_infos is None.")
      raise exception.IllegalArgumentError("prediction_protein_infos is None.")

    # </editor-fold>

    logger.info(
        f"Adding the protein names from the new job {prediction_protein_infos} to the project's blacklist."
    )
    for tmp_protein_info in prediction_protein_infos:
      self.protein_names_blacklist.append(tmp_protein_info.name)
    logger.info(
        f"After adding the new names, the protein blacklist is {self.protein_names_blacklist}."
    )

  def add_protein(self, a_protein_name: str) -> None:
    """Adds a protein name to the project's blacklist.

    Args:
        a_protein_name (str): The name of the protein to be added to the blacklist.

    Raises:
        raise exception.IllegalArgumentError: If `a_protein_name` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_protein_name is None or a_protein_name == "":
      logger.error("a_protein_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_protein_name is either None or an empty string."
      )

    # </editor-fold>

    logger.info(
        f"Adding the protein name {a_protein_name} to the project's blacklist."
    )
    self.protein_names_blacklist.append(a_protein_name)
    logger.info(
        f"After adding the name, the protein blacklist is {self.protein_names_blacklist}."
    )

  def remove_protein(self, a_protein_name: str) -> None:
    """Removes a protein name from the project's blacklist.

    Args:
        a_protein_name: The name of the protein to be removed from the project's blacklist.

    Raises:
        raise exception.IllegalArgumentError: If `a_protein_name` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_protein_name is None or a_protein_name == "":
      logger.error("a_protein_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_protein_name is either None or an empty string."
      )

    # </editor-fold>

    logger.info(
        f"Removing the protein name {a_protein_name} from the project's blacklist."
    )
    self.protein_names_blacklist.pop(
        self.protein_names_blacklist.index(a_protein_name)
    )
    logger.info(
        f"After removing the name, the protein blacklist is {self.protein_names_blacklist}."
    )

  def add_protein_pairs_from_new_job(self, analysis_run_names: list) -> None:
    """Adds protein pair names from a new job to the project's blacklist.

    Args:
        analysis_run_names (list): A list of analysis run names.

    Raises:
        exception.IllegalArgumentError: If `analysis_run_names` is None.
    """
    # <editor-fold desc="Checks">
    if analysis_run_names is None:
      logger.error("analysis_run_names is None.")
      raise exception.IllegalArgumentError("analysis_run_names is None.")

    # </editor-fold>

    logger.info(
        f"Adding the protein pair names from the new job {analysis_run_names} to the project's blacklist."
    )
    for tmp_analysis_run_name in analysis_run_names:
      tmp_string = tmp_analysis_run_name.replace(";", "_")
      tmp_string_2 = tmp_string.replace(",", "_")
      self.protein_pair_names_blacklist.append(tmp_string_2)
    logger.info(
        f"After adding the new names, the protein pair blacklist is {self.protein_pair_names_blacklist}."
    )

  def add_protein_pair(self, a_protein_pair_name: str) -> None:
    """Adds a protein pair name to the project's blacklist.

    Args:
        a_protein_pair_name (str): The protein pair name to be added to the blacklist.

    Raises:
        exception.IllegalArgumentError: If `a_protein_pair_name` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_protein_pair_name is None or a_protein_pair_name == "":
      logger.error("a_protein_pair_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_protein_pair_name is either None or an empty string."
      )

    # </editor-fold>

    logger.info(
        f"Adding the protein pair name {a_protein_pair_name} to the project's blacklist."
    )
    self.protein_pair_names_blacklist.append(a_protein_pair_name)
    logger.info(
        f"After adding the name, the protein pair blacklist is {self.protein_pair_names_blacklist}."
    )

  def remove_protein_pair(self, a_protein_pair_name: str) -> None:
    """Removes a protein pair name from the project's blacklist.

    Args:
        a_protein_pair_name (str): The name of the protein pair to be removed from the blacklist.

    Raises:
        exception.IllegalArgumentError: If `a_protein_pair_name` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_protein_pair_name is None or a_protein_pair_name == "":
      logger.error("a_protein_pair_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_protein_pair_name is either None or an empty string."
      )

    # </editor-fold>

    logger.info(
        f"Removing the protein pair name {a_protein_pair_name} from the project's blacklist."
    )
    self.protein_pair_names_blacklist.pop(
        self.protein_pair_names_blacklist.index(a_protein_pair_name)
    )
    logger.info(
        f"After removing the name, the protein pair blacklist is {self.protein_pair_names_blacklist}."
    )

  def is_protein_name_on_blacklist(self, a_protein_name: str) -> bool:
    """Checks if a protein name is on the blacklist.

    Args:
        a_protein_name (str): The protein name to check.

    Returns:
        True if the protein name is on the blacklist, False otherwise.

    Raises:
        exception.IllegalArgumentError: If `a_protein_name` is None.
    """
    # <editor-fold desc="Checks">
    if a_protein_name is None:
      logger.error("a_protein_name is None.")
      raise exception.IllegalArgumentError("a_protein_name is None.")

    # </editor-fold>

    if a_protein_name in self.protein_names_blacklist:
      return True
    return False

  # </editor-fold>
