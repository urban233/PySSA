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
"""Module contains job summary dataclasses."""
import logging
from dataclasses import dataclass
from src.pyssa.internal.data_structures.data_classes import prediction_protein_info
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import enums
from src.pyssa.util import exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


@dataclass
class JobBaseInformation:
  """Contains all general information about a job."""

  # <editor-fold desc="Class attributes">
  job_type: "enums.JobType"
  """The type of job."""

  project_name: str
  """The name of the project."""

  protein_names: list[str]
  """The names of the proteins."""

  protein_pair_names: list[str]
  """The names of the protein pairs."""

  job_progress: "enums.JobProgress"
  """The progress of the job."""

  # </editor-fold>

  def add_image_filepath(self, a_filepath: str) -> None:
    """Adds the filepath where the image should be saved.

    Args:
        a_filepath: The filepath of the image to be added.

    Raises:
        exception.IllegalArgumentError: If a_filepath is None.
    """
    # <editor-fold desc="Checks">
    if a_filepath is None:
      logger.error("a_filepath is None.")
      raise exception.IllegalArgumentError("a_filepath is None.")

    # </editor-fold>

    self._filepath = a_filepath

  def get_image_filepath(self) -> str:
    """Gets the filepath where the image should be saved.

    Returns:
        The filepath where the image should be saved or an empty string if the `_filepath` is None.
    """
    if self._filepath is None:
      return ""
    return self._filepath


@dataclass
class PredictionJobSummary:
  """Contains all information about a prediction job."""

  # <editor-fold desc="Class attributes">
  prediction_protein_infos: list[prediction_protein_info.PredictionProteinInfo]
  """The list of prediction protein infos."""

  # </editor-fold>

  def get_protein_names(self) -> list[str]:
    """Returns a list of protein names.

    This method retrieves the protein names from the prediction protein infos stored in the class instance.

    Returns:
        list[str]: A list of protein names.
    """
    tmp_protein_names = []
    for tmp_prediction_info in self.prediction_protein_infos:
      tmp_protein_names.append(tmp_prediction_info.name)
    return tmp_protein_names


@dataclass
class DistanceAnalysisJobSummary:
  """Contains all information about a distance analysis job."""

  # <editor-fold desc="Class attributes">
  analysis_names: list[str]
  """A list of analysis names."""

  # </editor-fold>
