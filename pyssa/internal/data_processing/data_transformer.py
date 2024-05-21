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
"""Module contains helper functions for specific data transformations."""
import copy
import logging
from pyssa.internal.data_structures.data_classes import prediction_protein_info
from pyssa.internal.data_structures import sequence
from pyssa.logging_pyssa import log_handlers
from pyssa.internal.data_structures.data_classes import analysis_run_info
from pyssa.internal.data_structures import protein
from pyssa.internal.data_structures import protein_pair
from pyssa.internal.analysis_types import distance_analysis
from pyssa.util import analysis_util
from pyssa.util import exception
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
  from pyssa.internal.data_structures import project
  from pyssa.internal.data_structures import settings

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


def transform_protein_name_seq_tuple_to_sequence_obj(
    proteins_to_predict: list["prediction_protein_info.PredictionProteinInfo"],
) -> list["protein.Protein"]:
  """Transforms the list of proteins to a list of protein objects.

  Args:
      proteins_to_predict: A list of tuples which consists of the protein name and sequence.

  Returns:
      A list of protein objects.

  Raises:
      IllegalArgumentError: If `proteins_to_predict` is None or the first element in an empty string.
  """
  # <editor-fold desc="Checks">
  if proteins_to_predict is None:
    logger.error("proteins_to_predict is None.")
    raise exception.IllegalArgumentError("proteins_to_predict is None.")
  if proteins_to_predict[0].name == "":
    logger.error("proteins_to_predict[0] is an empty string.")
    raise exception.IllegalArgumentError(
        "proteins_to_predict[0] is an empty string."
    )
  if not proteins_to_predict[0].sequences:
    logger.error("proteins_to_predict[0] has no sequence.")
    raise exception.IllegalArgumentError(
        "proteins_to_predict[0] has no sequence."
    )

  # </editor-fold>

  protein_objects: list["protein.Protein"] = []

  for tmp_protein_to_predict in proteins_to_predict:
    protein_with_sequence = protein.Protein(tmp_protein_to_predict.name)
    protein_with_sequence.chains = []
    for tmp_seq in tmp_protein_to_predict.sequences:
      logger.debug(
          "tmp_seq: %s, tmp_seqs: %r", tmp_seq, tmp_protein_to_predict.sequences
      )
      protein_seq = sequence.Sequence(tmp_protein_to_predict.name, tmp_seq)
      protein_with_sequence.append_chain(
          chain_name="", chain_sequence=protein_seq, chain_type="protein_chain"
      )
    protein_with_sequence.add_chain_names_to_chains()
    protein_objects.append(protein_with_sequence)

  logger.debug(protein_objects)
  return protein_objects


class DistanceAnalysisDataTransformer:
  """Transforms data from the gui to a manageable format."""

  # <editor-fold desc="Class attributes">
  analysis_run_name: str
  """The name of a single analysis run."""

  current_project: "project.Project"
  """The current project of the main window."""

  settings: "settings.Settings"
  """The settings of PySSA."""

  analysis_run: "analysis_run_info.AnalysisRunInfo"
  """The information about the analysis run, includes the names and chains and analysis name."""

  proteins: tuple["protein.Protein", "protein.Protein"]
  """A tuple of two proteins."""

  analysis_protein_pair: "protein_pair.ProteinPair"
  """The protein pair for the distance analysis."""

  analysis_distance: "distance_analysis.DistanceAnalysis"
  """The distance analysis object."""

  # </editor-fold>

  def __init__(
      self,
      an_analysis_run_name: str,
      the_app_project: "project.Project",
      a_cutoff: float,
      cycles: int,
  ) -> None:
    """Constructor.

    Args:
        an_analysis_run_name (str): The name of the analysis run.
        the_app_project (project.Project): The project object for the app.
        a_cutoff (float): The cutoff value.
        cycles (int): The number of cycles.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if an_analysis_run_name is None:
      logger.error("an_analysis_run_name is None.")
      raise exception.IllegalArgumentError("an_analysis_run_name is None.")
    if the_app_project is None:
      logger.error("the_app_project is None.")
      raise exception.IllegalArgumentError("the_app_project is None.")
    if a_cutoff is None:
      logger.error("a_cutoff is None.")
      raise exception.IllegalArgumentError("a_cutoff is None.")
    if cycles is None:
      logger.error("cycles is None.")
      raise exception.IllegalArgumentError("cycles is None.")
    # </editor-fold>

    self.analysis_run_name: str = an_analysis_run_name
    self.current_project: "project.Project" = the_app_project
    self.cutoff = a_cutoff
    self.cycles = cycles

  def _create_analysis_run_info(self) -> None:
    """Create the analysis run info based on the analysis run name.

    Raises:
        UnableToCreateAnalysisRunInfoError: If the analysis run name is illegal and no analysis run info can be created.
    """
    try:
      tmp_analysis_run_infos: list[Union[str, list[str]]] = (
          analysis_util.split_analysis_run_name_in_protein_name_and_chain(
              self.analysis_run_name
          )
      )
      logger.debug(f"tmp_analysis_run_infos: {tmp_analysis_run_infos}")
      self.analysis_run = analysis_run_info.AnalysisRunInfo(
          tmp_analysis_run_infos[0],
          tmp_analysis_run_infos[1],
          tmp_analysis_run_infos[2],
          tmp_analysis_run_infos[3],
          self.analysis_run_name,
      )
    except exception.IllegalArgumentError:
      logger.error("Creating the analysis run info failed.")
      raise exception.UnableToCreateAnalysisRunInfoError("")

  def _create_proteins_for_analysis(self) -> None:
    """Creates protein objects for distance analysis.

    Raises:
        ProteinNotFoundInCurrentProjectError: If a protein is not found in the current project.
    """
    protein_1: "protein.Protein" = self.current_project.search_protein(
        self.analysis_run.get_protein_name_1()
    )
    if protein_1 is None:
      logger.error(
          "No protein with the given protein name: "
          f"{self.analysis_run.get_protein_name_1()} found in the current project.",
      )
      raise exception.ProteinNotFoundInCurrentProjectError("")
    else:
      logger.debug("Protein 1 object successfully created.")

    if self.analysis_run.are_protein_names_identical():
      protein_2: "protein.Protein" = copy.deepcopy(protein_1)
      protein_1: "protein.Protein" = copy.deepcopy(protein_2)
      protein_1.set_molecule_object(f"{protein_1.get_molecule_object()}_1")
      protein_2.set_molecule_object(f"{protein_2.get_molecule_object()}_2")
    else:
      protein_2: "protein.Protein" = self.current_project.search_protein(
          self.analysis_run.get_protein_name_2()
      )
      if protein_2 is None:
        logger.error(
            "No protein with the given protein name: "
            f"{self.analysis_run.get_protein_name_2()} found in the current project.",
        )
        raise exception.ProteinNotFoundInCurrentProjectError("")
      else:
        logger.debug("Protein 2 object successfully created.")
    self.proteins = (protein_1, protein_2)

  def _create_analysis_name(self) -> str:
    """Creates the name of the analysis.

    Returns:
        The analysis name.
    """
    if len(self.proteins[0].chains) != 0:
      analysis_name = f"{self.proteins[0].get_molecule_object()};{self.analysis_run.protein_chains_1}_vs_{self.proteins[1].get_molecule_object()};{self.analysis_run.protein_chains_2}"  # noqa
      analysis_name = analysis_name.replace(";", "_")
      analysis_name = analysis_name.replace(",", "_")
      analysis_name = analysis_name.replace("[", "")
      analysis_name = analysis_name.replace("]", "")
      analysis_name = analysis_name.replace("'", "")
      analysis_name = analysis_name.replace(" ", "")
    else:
      analysis_name = f"{self.proteins[0].get_molecule_object()}_vs_{self.proteins[1].get_molecule_object()}"
    return analysis_name

  def _create_protein_pair(self) -> None:
    """Creates the protein pair for the distance analysis.

    Raises:
        UnableToCreateProteinPairError: If the creation of the protein pair failed.
    """
    try:
      self.analysis_protein_pair = protein_pair.ProteinPair(
          self.proteins[0], self.proteins[1]
      )
      self.analysis_protein_pair.name = self._create_analysis_name()
    except exception.IllegalArgumentError:
      logger.error("Creating the protein pair for the analysis failed.")
      raise exception.UnableToCreateProteinPairError("")

  def _create_distance_analysis(self) -> None:
    """Creates distance analysis object and sets the object into the protein pair.

    Raises:
        UnableToCreateDistanceAnalysis: If the creation of the distance analysis or setting the object failed.
    """
    try:
      self.analysis_protein_pair.set_distance_analysis(
          distance_analysis.DistanceAnalysis(
              self.analysis_protein_pair, self.cutoff, self.cycles
          ),
      )
    except exception.IllegalArgumentError:
      logger.error(
          "Creating and setting the distance analysis into the protein pair faild."
      )
      raise exception.UnableToCreateDistanceAnalysisObjectError("")

  def _set_selection_for_analysis(self) -> None:
    """Sets the selection for the analysis into the distance analysis object."""
    prot_1_selection = (
        analysis_util.create_selection_strings_for_structure_alignment(
            self.analysis_protein_pair.protein_1,
            self.analysis_run.protein_chains_1,
        )
    )
    prot_2_selection = (
        analysis_util.create_selection_strings_for_structure_alignment(
            self.analysis_protein_pair.protein_2,
            self.analysis_run.protein_chains_2,
        )
    )
    self.analysis_protein_pair.distance_analysis.create_align_selections(
        prot_1_selection, prot_2_selection
    )

  def transform_gui_input_to_distance_analysis_object(
      self,
  ) -> "protein_pair.ProteinPair":
    """Transforms the gui data into a distance analysis object.

    Returns:
        The protein pair object that gets analyzed.

    Raises:
        UnableToTransformDataForAnalysisError: If an error occurs during the transformation process.
    """
    try:
      self._create_analysis_run_info()
      self._create_proteins_for_analysis()
      self._create_protein_pair()
      self._create_distance_analysis()
      self._set_selection_for_analysis()
    except exception.UnableToCreateAnalysisRunInfoError:
      logger.error("Unable to create analysis run info for distance analysis.")
      raise exception.UnableToTransformDataForAnalysisError("")
    except exception.ProteinNotFoundInCurrentProjectError:
      logger.error(
          "Unable to find one of the proteins for the distance analysis."
      )
      raise exception.UnableToTransformDataForAnalysisError("")
    except exception.UnableToCreateProteinPairError:
      logger.error("Unable to create protein pair for distance analysis.")
      raise exception.UnableToTransformDataForAnalysisError("")
    except exception.UnableToCreateDistanceAnalysisObjectError:
      logger.error(
          "Unable to create or set the distance analysis object for distance analysis."
      )
      raise exception.UnableToTransformDataForAnalysisError("")
    except Exception as e:
      logger.error(f"Unknown error: {e}.")
      raise exception.UnableToTransformDataForAnalysisError("")
    return self.analysis_protein_pair
