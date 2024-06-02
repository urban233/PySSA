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
"""Module for the results class."""
import csv
import pathlib
from xml.etree import ElementTree
import logging
import numpy as np

from src.pyssa.io_pyssa import path_util
from src.pyssa.io_pyssa import binary_data
from src.pyssa.io_pyssa.xml_pyssa import element_names
from src.pyssa.io_pyssa.xml_pyssa import attribute_names
from src.pyssa.util import pyssa_keys, exception
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import constants

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class DistanceAnalysisResults:
  """Contains the results of the distance analysis done in PyMOL."""

  # <editor-fold desc="Class attributes">
  distance_data: dict[str, np.ndarray]
  """The data of the distance analysis including: position, chain, residue, residue number and distance."""

  pymol_session: str
  """The base64 string of the pymol session used for the distance analysis."""

  rmsd: float
  """The RMSD value of the protein pair."""

  aligned_aa: str
  """The aligned residues of the protein pair."""

  structure_aln_image: tuple[str, str] = ()
  """A tuple with the basename and the base64 string of the structure alignment image."""

  interesting_regions_images: list[tuple[str, str]] = []
  """A list of tuples with the basename and the base64 string of the interesting region image."""

  # </editor-fold>

  def __init__(
      self,
      distance_data: dict,
      pymol_session: str,
      rmsd: float,
      aligned_aa: str,
  ) -> None:
    """Constructor.

    Args:
        distance_data (dict): The distance data containing information about the distances between atoms.
        pymol_session (str): The PyMOL session file name or path.
        rmsd (float): The root-mean-square deviation (RMSD) value.
        aligned_aa (str): The aligned amino acid sequence.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if distance_data is None:
      logger.error("distance_data is None.")
      raise exception.IllegalArgumentError("distance_data is None.")
    if pymol_session is None:
      logger.error("pymol_session is None.")
      raise exception.IllegalArgumentError("pymol_session is None.")
    if rmsd is None:
      logger.error("rmsd is None.")
      raise exception.IllegalArgumentError("rmsd is None.")
    if aligned_aa is None:
      logger.error("aligned_aa is None.")
      raise exception.IllegalArgumentError("aligned_aa is None.")

    # </editor-fold>

    self.distance_data = distance_data
    self.pymol_session = pymol_session
    self.rmsd = rmsd
    self.aligned_aa = aligned_aa

  def export_distance_data_as_csv(self, a_filepath: str) -> None:
    """Exports distance data as a CSV file.

    Args:
        a_filepath (str): The file path to save the CSV file.

    Raises:
        exception.IllegalArgumentError: If a_filepath is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_filepath is None or a_filepath == "":
      logger.error("a_filepath is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_filepath is either None or an empty string."
      )

    # </editor-fold>

    tmp_data_to_write = np.transpose(
        [
            self.distance_data["index"],
            self.distance_data["ref_chain"],
            self.distance_data["ref_pos"],
            self.distance_data["ref_resi"],
            self.distance_data["model_chain"],
            self.distance_data["model_pos"],
            self.distance_data["model_resi"],
            self.distance_data["distance"],
        ]
    )
    # Writing to the TSV file
    with open(a_filepath, mode="w", newline="", encoding="utf-8") as tsv_file:
      csv_writer = csv.writer(tsv_file, delimiter=",")
      # Write the header
      csv_writer.writerow(
          [
              "Residue_Pair_No",
              "Protein_1_Chain",
              "Protein_1_Position",
              "Protein_1_Residue",
              "Protein_2_Chain",
              "Protein_2_Position",
              "Protein_2_Residue",
              "Distance",
          ]
      )
      # Write the data to the TSV file
      csv_writer.writerows(tmp_data_to_write)
