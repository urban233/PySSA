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
"""Module contains the distance analysis class."""
import logging
import pathlib
from typing import TYPE_CHECKING
import numpy as np
from pyssa.logging_pyssa import log_handlers
from pyssa.io_pyssa import path_util
from pyssa.internal.data_structures import results
from pyssa.util import exception

if TYPE_CHECKING:
  from pyssa.internal.data_structures import protein_pair
  from pyssa.internal.data_structures import settings
  from pyssa.internal.data_structures import selection

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class DistanceAnalysis:
  """Contains information about the distance analysis."""

  # <editor-fold desc="Class attributes">
  _protein_pair_for_analysis: "protein_pair.ProteinPair"
  """A pair of proteins which get analyzed."""

  app_settings: "settings.Settings"
  """The settings of PySSA."""

  cutoff: float
  """The cutoff value for the align command from PyMOL."""

  cycles: int
  """The number of refinement cycles for the align command from PyMOL."""

  figure_size: tuple[float, float]
  """The size of the figure."""

  alignment_file_name: str = "aln"
  """The filename of the alignment file which gets created during the align command."""

  export_dirname: pathlib.Path
  """A directory where all results related to the protein will be stored."""

  pymol_session_filepath: path_util.FilePath
  """The session filepath for the PyMOL session of the analysis."""

  distance_analysis_data: dict[str, np.ndarray] = {}
  """A dictionary containing the distance analysis results."""

  analysis_results: "results.DistanceAnalysisResults" = None
  """An object that contains information about the distance analysis results."""

  rmsd_dict: dict = {"rmsd": 0.0, "aligned_aa": 0}
  """A dictionary containing the RMSD value and the number of aligned residues."""

  # </editor-fold>

  def __init__(
      self,
      protein_pair_for_analysis: "protein_pair.ProteinPair",
      a_cutoff: float,
      cycles: int,
  ) -> None:
    """Constructor.

    Args:
        protein_pair_for_analysis (protein_pair.ProteinPair): An instance of the `protein_pair.ProteinPair` class representing the pair of proteins to be analyzed.
        a_cutoff (float): A float representing the cutoff value for distance analysis.
        cycles (int): An integer representing the number of cycles to be performed in the analysis.

    Raises:
        exception.IllegalArgumentError: If protein_pair_for_analysis is None.
        exception.IllegalArgumentError: If a_cutoff is either None or a negative value.
        exception.IllegalArgumentError: If cycles is either None or a negative value.
    """
    # <editor-fold desc="Checks">
    if protein_pair_for_analysis is None:
      logger.error("protein_pair_for_analysis is None.")
      raise exception.IllegalArgumentError("protein_pair_for_analysis is None.")
    if a_cutoff is None or a_cutoff < 0:
      logger.error("a_cutoff is either None or a negative value.")
      raise exception.IllegalArgumentError(
          "a_cutoff is either None or a negative value."
      )
    if cycles is None or cycles < 0:
      logger.error("cycles is either None or a negative value.")
      raise exception.IllegalArgumentError(
          "cycles is either None or a negative value."
      )

    # </editor-fold>

    self._protein_pair_for_analysis: protein_pair.ProteinPair = (
        protein_pair_for_analysis
    )
    self.name = f"dist_analysis_{self._protein_pair_for_analysis.name}"
    self.cutoff: float = a_cutoff
    self.cycles: int = cycles
    self.figure_size = (11.0, 6.0)
    self.alignment_file_name = "aln"

  def create_align_selections(
      self,
      protein_1_selection: "selection.Selection",
      protein_2_selection: "selection.Selection",
  ) -> None:
    """Creates the selection which are needed for the align command.

    Args:
        protein_1_selection (selection.Selection): The selection of protein 1 to be aligned.
        protein_2_selection (selection.Selection): The selection of protein 2 to be aligned.

    Raises:
        exception.IllegalArgumentError: If `protein_1_selection` or `protein_2_selection` is either None or an empty string.
        ValueError: If the selection is illegal, i.e., it does not belong to the respective protein.

    """
    # <editor-fold desc="Checks">
    if protein_1_selection is None or protein_1_selection == "":
      logger.error("protein_1_selection is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "protein_1_selection is either None or an empty string."
      )
    if protein_2_selection is None or protein_2_selection == "":
      logger.error("protein_2_selection is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "protein_2_selection is either None or an empty string."
      )
    if (
        protein_1_selection.molecule_object
        != self._protein_pair_for_analysis.protein_1.get_molecule_object()
    ):  # pylint: disable=line-too-long
      logger.error("Selection is illegal.")
      raise ValueError("Selection is illegal.")
    if (
        protein_2_selection.molecule_object
        != self._protein_pair_for_analysis.protein_2.get_molecule_object()
    ):  # pylint: disable=line-too-long
      logger.error("Selection is illegal.")
      raise ValueError("Selection is illegal.")

    # </editor-fold>

    # logger.debug(
    #     f"1st argument of <create_align_selections>: "
    #     f"{protein_1_selection.selection_string} {protein_1_selection}",
    # )
    # logger.debug(
    #     f"2nd argument of <create_align_selections>: "
    #     f"{protein_2_selection.selection_string} {protein_2_selection}",
    # )

    self._protein_pair_for_analysis.protein_1.pymol_selection = (
        protein_1_selection
    )
    self._protein_pair_for_analysis.protein_2.pymol_selection = (
        protein_2_selection
    )
    logger.debug(
        f"Prot 1 sele in <create_align_selections>: "
        f"{self._protein_pair_for_analysis.protein_1.pymol_selection.selection_string} "
        f"{self._protein_pair_for_analysis.protein_1.pymol_selection}",
    )
    logger.debug(
        f"Prot 2 sele in <create_align_selections>: "
        f"{self._protein_pair_for_analysis.protein_2.pymol_selection.selection_string} "
        f"{self._protein_pair_for_analysis.protein_2.pymol_selection}",
    )
