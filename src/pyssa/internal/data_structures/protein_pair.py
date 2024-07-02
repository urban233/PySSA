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
"""Module for the protein pair class."""
import copy
import logging
from typing import TYPE_CHECKING
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.io_pyssa import path_util
from src.pyssa.util import exception

if TYPE_CHECKING:
  from src.pyssa.internal.data_structures import protein
  from src.pyssa.internal.data_structures import structure_analysis

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class ProteinPair:
  """This class consists of two Protein objects. It is used to have a better workflow for the analysis."""

  # <editor-fold desc="Class attributes">
  _id: int
  """The unique identifier of the protein."""

  protein_1: "protein.Protein"
  """The first protein of the protein pair."""

  protein_2: "protein.Protein"
  """The second protein of the protein pair."""

  distance_analysis: "structure_analysis.DistanceAnalysis" = None
  """A directory where all results related to the protein will be stored."""

  pymol_session_filepath: path_util.FilePath
  """The full filepath where the session file is stored."""

  pymol_session: str
  """A base64 string of the pymol session."""

  db_project_id: int
  """The project id from the database."""

  # </editor-fold>

  def __init__(
      self, protein_1: "protein.Protein", protein_2: "protein.Protein"
  ) -> None:
    """Constructor.

    Args:
        protein_1 (protein.Protein): The first protein object.
        protein_2 (protein.Protein): The second protein object.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.

    Notes:
        The protein objects getting "deep copied"!
    """
    # <editor-fold desc="Checks">
    if protein_1 is None:
      logger.error("protein_1 is None.")
      raise exception.IllegalArgumentError("protein_1 is None.")
    if protein_2 is None:
      logger.error("protein_2 is None.")
      raise exception.IllegalArgumentError("protein_2 is None.")

    # </editor-fold>

    self.protein_1: "protein.Protein" = copy.deepcopy(protein_1)
    self.protein_2: "protein.Protein" = copy.deepcopy(protein_2)
    self.name = f"{self.protein_1.get_molecule_object()}_with_{self.protein_2.get_molecule_object()}"
    self.pymol_session = ""

  def set_distance_analysis(
      self, a_value: "distance_analysis.DistanceAnalysis"
  ) -> None:
    """Sets the distance analysis object.

    Args:
        a_value (distance_analysis.DistanceAnalysis): A distance analysis object.

    Raises:
        exception.IllegalArgumentError: If a_value is None.
    """
    # <editor-fold desc="Checks">
    if a_value is None:
      logger.error("a_value is None.")
      raise exception.IllegalArgumentError("a_value is None.")

    # </editor-fold>

    self.distance_analysis = a_value

  def get_id(self) -> int:
    """Returns the ID of the object.

    Returns:
        The ID of the object.
    """
    return self._id

  def set_id(self, a_value: int) -> None:
    """Set the value of the ID for the object.

    Args:
        a_value (int): The new value to assign to the ID.

    Raises:
        exception.IllegalArgumentError: If a_value is either None or a value less than 0.
    """
    # <editor-fold desc="Checks">
    if a_value is None or a_value < 0:
      logger.error("a_value is either None or a value less than 0.")
      raise exception.IllegalArgumentError(
          "a_value is either None or a value less than 0."
      )

    # </editor-fold>

    self._id = a_value

  def get_protein_name_without_chains(self) -> str:
    """Returns the name of the protein without chains.

    The method returns a string representing the name of the protein without the chains.
    It combines the molecule objects of protein_1 and protein_2 with a hyphen in between.

    Returns:
        The name of the protein pair without chains.

    """
    return f"{self.protein_1.get_molecule_object()}-{self.protein_2.get_molecule_object()}"
