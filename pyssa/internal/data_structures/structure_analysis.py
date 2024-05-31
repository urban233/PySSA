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
"""Module for structure analysis class."""
import logging
import pathlib
import numpy as np
from typing import TYPE_CHECKING

import zmq

from auxiliary_pymol import auxiliary_pymol_client
from pyssa.controller import database_manager
from pyssa.io_pyssa import filesystem_helpers, bio_data
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants, pyssa_keys, enums, constant_messages
from pyssa.util import exception
from pyssa.internal.data_structures import results, job

if TYPE_CHECKING:
  from pyssa.internal.data_structures import project
  from pyssa.internal.data_structures import protein_pair
  from pyssa.internal.data_structures import settings


logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class Analysis:
  """Contains information about the type of analysis."""

  # <editor-fold desc="Class attributes">
  analysis_list: list["protein_pair.ProteinPair"] = []
  """A list of protein pairs to analyze."""

  app_project: "project.Project"
  """The project object for which analysis is being performed."""

  # </editor-fold>

  def __init__(self, app_project: "project.Project") -> None:
    """Constructor.

    Args:
        app_project(project.Project): The project that this analysis is used.

    Raises:
        exception.IllegalArgumentError: If app_project is None.
    """
    # <editor-fold desc="Checks">
    if app_project is None:
      logger.error("app_project is None.")
      raise exception.IllegalArgumentError("app_project is None.")

    # </editor-fold>

    self.app_project: "project.Project" = app_project

  def _fetch_pdb_atoms_for_all_proteins(self) -> None:
    """Fetches PDB atom data for all proteins in the analysis list from the database and sets the retrieved data to the respective proteins' pdb_data attribute."""
    with database_manager.DatabaseManager(
        str(self.app_project.get_database_filepath())
    ) as db_manager:
      for tmp_protein_pair in self.analysis_list:
        tmp_pdb_atom_db_data_1 = db_manager.get_pdb_atoms_of_protein(
            tmp_protein_pair.protein_1.get_id()
        )
        tmp_pdb_atom_db_data_2 = db_manager.get_pdb_atoms_of_protein(
            tmp_protein_pair.protein_2.get_id()
        )
        pdb_atom_dict_1 = [
            {key.value: value for key, value in zip(enums.PdbAtomEnum, t)}
            for t in tmp_pdb_atom_db_data_1
        ]
        pdb_atom_dict_2 = [
            {key.value: value for key, value in zip(enums.PdbAtomEnum, t)}
            for t in tmp_pdb_atom_db_data_2
        ]
        tmp_protein_pair.protein_1.set_pdb_data(pdb_atom_dict_1)
        tmp_protein_pair.protein_2.set_pdb_data(pdb_atom_dict_2)

  def _get_last_id_of_protein_pairs(self) -> int:
    """Get the last ID of protein pairs.

    Returns:
        The last ID of protein pairs.
    """
    if len(self.app_project.protein_pairs) == 0:
      return 1
    return max(
        self.app_project.protein_pairs, key=lambda obj: obj.get_id()
    ).get_id()

  def run_distance_analysis(
      self,
      the_image_creation_option: bool,
      the_main_socket: zmq.Socket,
      a_socket: zmq.Socket,
  ) -> None:
    """Runs the distance analysis for all protein pairs of the analysis job.

    Args:
        the_image_creation_option (bool): A flag indicating whether to create images.
        the_main_socket (zmq.Socket): The main socket for communication.
        a_socket (zmq.Socket): A secondary socket for communication.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
        UnableToSetImageError: If no image was setting.
        UnableToOpenFileError: If file could not be opened.
        exception.UnableToDoAnalysisError: If the analysis fails.
    """
    # <editor-fold desc="Checks">
    if the_image_creation_option is None:
      logger.error("the_image_creation_option is None.")
      raise exception.IllegalArgumentError("the_image_creation_option is None.")
    if the_main_socket is None:
      logger.error("the_main_socket is None.")
      raise exception.IllegalArgumentError("the_main_socket is None.")
    if a_socket is None:
      logger.error("a_socket is None.")
      raise exception.IllegalArgumentError("a_socket is None.")

    # </editor-fold>

    logger.debug(
        f"The function argument of the value for the image_creation_option is: {the_image_creation_option}",
    )
    self._fetch_pdb_atoms_for_all_proteins()
    tmp_latest_protein_pair_id = self._get_last_id_of_protein_pairs()
    # create scratch dirs
    filesystem_helpers.create_directory(constants.SCRATCH_DIR_ANALYSIS)

    for tmp_protein_pair in self.analysis_list:
      logger.info(f"The protein pair: {tmp_protein_pair.name} gets analyzed.")
      try:
        # do distance analysis in PyMOL
        if tmp_protein_pair is None:
          logger.error(constant_messages.ARGUMENT_IS_ILLEGAL)
          raise exception.IllegalArgumentError(
              constant_messages.ARGUMENT_IS_ILLEGAL
          )

        tmp_protein_1_name = tmp_protein_pair.protein_1.get_molecule_object()
        tmp_protein_1_pdb_cache_filepath = pathlib.Path(
            f"{constants.CACHE_PROTEIN_DIR}/{tmp_protein_1_name}.pdb",  # fixme: Must be a string, but is a path!
        )
        tmp_protein_2_name = tmp_protein_pair.protein_2.get_molecule_object()
        tmp_protein_2_pdb_cache_filepath = pathlib.Path(
            f"{constants.CACHE_PROTEIN_DIR}/{tmp_protein_2_name}.pdb",  # fixme: Must be a string, but is a path!
        )
        try:
          bio_data.build_pdb_file(
              tmp_protein_pair.protein_1.get_pdb_data(),
              tmp_protein_1_pdb_cache_filepath,
          )
          bio_data.build_pdb_file(
              tmp_protein_pair.protein_2.get_pdb_data(),
              tmp_protein_2_pdb_cache_filepath,
          )
        except Exception as e:
          logger.error(f"PDB file could not be built. Error: {e}")

        tmp_reply = auxiliary_pymol_client.send_request_to_auxiliary_pymol(
            the_main_socket,
            a_socket,
            job.DistanceAnalysisJobDescription(
                tmp_protein_pair.name,
                tmp_protein_1_pdb_cache_filepath,
                tmp_protein_2_pdb_cache_filepath,
                tmp_protein_pair.protein_1.pymol_selection.selection_string,
                tmp_protein_pair.protein_2.pymol_selection.selection_string,
                tmp_protein_pair.distance_analysis.cutoff,
                tmp_protein_pair.distance_analysis.cycles,
            ),
        )
        if tmp_reply["result"] == "error" and tmp_reply["data"] == "More than one state":
          raise exception.UnableToDoAnalysisError(tmp_reply["data"])
        if tmp_reply["data"] is not None:
          print(tmp_reply["data"])
          distance_analysis_results_object_values, base64_string = tmp_reply[
              "data"
          ]
          distances, base64_string, rmsd, aligned_residues = (
              distance_analysis_results_object_values
          )
          result_hashtable: dict[str, np.ndarry] = {
              pyssa_keys.ARRAY_DISTANCE_INDEX: np.array(
                  distances[pyssa_keys.ARRAY_DISTANCE_INDEX]
              ),
              pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN: np.array(
                  distances[pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN]
              ),
              pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION: np.array(
                  distances[pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION]
              ),
              pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI: np.array(
                  distances[pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI]
              ),
              pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN: np.array(
                  distances[pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN]
              ),
              pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION: np.array(
                  distances[pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION]
              ),
              pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI: np.array(
                  distances[pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI]
              ),
              pyssa_keys.ARRAY_DISTANCE_DISTANCES: np.array(
                  distances[pyssa_keys.ARRAY_DISTANCE_DISTANCES]
              ),
          }
          print(result_hashtable)

          tmp_protein_pair.distance_analysis.analysis_results = (
              results.DistanceAnalysisResults(
                  result_hashtable,
                  base64_string,
                  rmsd,
                  aligned_residues,
              )
          )
          tmp_protein_pair.pymol_session = base64_string
        else:
          logger.warning("Returning data was None!")
      except exception.IllegalArgumentError:
        logger.error("The argument filename is illegal.")
        raise exception.UnableToOpenFileError(
            f"filename: {tmp_protein_pair.name}"
        )
      except exception.UnableToDoAnalysisError:
        raise exception.UnableToDoAnalysisError("Ambiguous selection during structure analysis.")
      except Exception as e:
        logger.error(f"Unknown error: {e}")
        raise exception.UnableToSetImageError("")
      filesystem_helpers.delete_directory(constants.SCRATCH_DIR_ANALYSIS)

  def run_analysis(
      self,
      the_analysis_type: str,
      the_image_option: bool,
      the_main_socket: zmq.Socket,
      a_socket: zmq.Socket,
  ) -> None:
    """Starts the distance analysis.

    Args:
        the_analysis_type (str): The type of analysis to run. Possible values are "distance".
        the_image_option (bool): Option for image analysis.
        the_main_socket (zmq.Socket): The main socket for communication.
        a_socket (zmq.Socket): A socket for communication.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
        ValueError: If `self.analysis_list` is empty.
        ValueError: If analysis type is unknown.
    """
    # <editor-fold desc="Checks">
    if the_analysis_type is None:
      logger.error("the_analysis_type is None.")
      raise exception.IllegalArgumentError("the_analysis_type is None.")
    if the_image_option is None:
      logger.error("the_image_option is None.")
      raise exception.IllegalArgumentError("the_image_option is None.")
    if the_main_socket is None:
      logger.error("the_main_socket is None.")
      raise exception.IllegalArgumentError("the_main_socket is None.")
    if a_socket is None:
      logger.error("a_socket is None.")
      raise exception.IllegalArgumentError("a_socket is None.")

    if len(self.analysis_list) == 0:
      logger.error("Analysis list is empty.")
      logger.debug(f"self.analysis_list: {self.analysis_list}")
      raise ValueError("Analysis list is empty.")

    # </editor-fold>

    if the_analysis_type == "distance":
      self.run_distance_analysis(the_image_option, the_main_socket, a_socket)
    else:
      tmp_msg: str = f"Unknown analysis type: {the_analysis_type}"
      logger.error(tmp_msg)
      raise ValueError(tmp_msg)


class DistanceAnalysis:
  """Contains all information about the distance analysis of a certain protein pair."""

  # <editor-fold desc="Checks">
  name: str
  """The name of the distance analysis."""

  cutoff: float
  """The cutoff value for the align command from PyMOL."""

  cycles: int
  """The number of refinement cycles for the align command from PyMOL."""

  figure_size: tuple[float, float]
  """The size of the figures."""

  analysis_results: "results.DistanceAnalysisResults" = None
  """The object which contains all results from the distance analysis done in PyMOL."""

  # </editor-fold>

  def __init__(
      self,
      the_app_settings: "settings.Settings",
      protein_pair_name: str = "",
      distance_analysis_name: str = "",
  ) -> None:
    """Constructor.

    Args:
        the_app_settings (settings.Settings): The settings object that contains the cutoff and cycles values.
        protein_pair_name (str, optional): The name of the protein pair. Defaults to an empty string.
        distance_analysis_name (str, optional): The name of the distance analysis. Defaults to an empty string.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_app_settings is None:
      logger.error("the_app_settings is None.")
      raise exception.IllegalArgumentError("the_app_settings is None.")
    if protein_pair_name is None:
      logger.error("protein_pair_name is None.")
      raise exception.IllegalArgumentError("protein_pair_name is None.")
    if distance_analysis_name is None:
      logger.error("distance_analysis_name is None.")
      raise exception.IllegalArgumentError("distance_analysis_name is None.")

    # </editor-fold>

    self.name: str = f"dist_analysis_{protein_pair_name}"
    self.cutoff: float = the_app_settings.cutoff
    self.cycles: int = the_app_settings.cycles
    self.figure_size: tuple[float, float] = (11.0, 6.0)
