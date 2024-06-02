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
"""Module contains helper function for the analysis process."""
import copy
import logging
from typing import Union
from typing import TYPE_CHECKING
from src.pyssa.internal.data_processing import data_transformer
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.internal.data_structures import protein
from src.pyssa.internal.data_structures import selection
from src.pyssa.util import protein_util, exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"

if TYPE_CHECKING:
  from src.pyssa.internal.data_structures import project


def transform_gui_input_to_practical_data(
    a_list_with_analysis_names: list,
    a_project: "project.Project",
    a_cutoff: float,
    cycles: int,
) -> list:
  """Transforms GUI input into practical data for analysis.

  Args:
      a_list_with_analysis_names (list): A list with analysis names.
      a_project (project.Project): The project object.
      a_cutoff (float): The cutoff value.
      cycles (int): The number of cycles.

  Returns:
      list: The transformed distance analysis runs.

  Raises:
      exception.IllegalArgumentError: If an argument is illegal.
      exception.UnableToTransformDataForAnalysisError: If the transformation process fails or an unknown error occurs.
  """
  # <editor-fold desc="Checks">
  if a_list_with_analysis_names is None or len(a_list_with_analysis_names) == 0:
    logger.error("a_list_with_analysis_names is either None or an empty list.")
    raise exception.IllegalArgumentError(
        "a_list_with_analysis_names is either None or an empty list."
    )
  if a_project is None:
    logger.error("a_project is None.")
    raise exception.IllegalArgumentError("a_project is None.")
  if a_cutoff is None or a_cutoff < 0:
    logger.error("a_cutoff is either None or a negative number.")
    raise exception.IllegalArgumentError(
        "a_cutoff is either None or a negative number."
    )
  if cycles is None or cycles < 0:
    logger.error("cycles is either None or a negative number.")
    raise exception.IllegalArgumentError(
        "cycles is either None or a negative number."
    )

  # </editor-fold>

  distance_analysis_runs = []
  try:
    for tmp_analysis_name in a_list_with_analysis_names:
      input_transformer = data_transformer.DistanceAnalysisDataTransformer(
          tmp_analysis_name,
          a_project,
          a_cutoff,
          cycles,
      )
      protein_pair_for_analysis = (
          input_transformer.transform_gui_input_to_distance_analysis_object()
      )
      new_protein_pair = copy.deepcopy(protein_pair_for_analysis)
      distance_analysis_runs.append(new_protein_pair)
    logger.debug(
        f"These are the distance analysis runs, after the data transformation: {distance_analysis_runs}.",
    )
  except exception.IllegalArgumentError:
    logger.error(
        "Transformation of data failed because an argument was illegal."
    )
    raise exception.UnableToTransformDataForAnalysisError("")
  except exception.UnableToTransformDataForAnalysisError:
    logger.error(
        "Transformation of data failed because the transformation process failed."
    )
    raise exception.UnableToTransformDataForAnalysisError("")
  except Exception as e:
    logger.error(f"Unknown error: {e}")
    raise exception.UnableToTransformDataForAnalysisError("")
  return distance_analysis_runs


def split_analysis_run_name_in_protein_name_and_chain(
    an_analysis_run_name: str,
) -> list[Union[str, list[str]]]:
  """Splits the analysis run name into protein name and chain.

  Args:
      an_analysis_run_name (str): The analysis run name that will be split.

  Returns:
      A list of the analysis run name split into protein name and chain.

  Raises:
      IllegalArgumentException: If the argument is None or an empty string or does not contain the _vs_.
  """
  # <editor-fold desc="Checks">
  if an_analysis_run_name is None or an_analysis_run_name == "":
    logger.error("an_analysis_run_name is either None or an empty string.")
    raise exception.IllegalArgumentError(
        "an_analysis_run_name is either None or an empty string."
    )
  if an_analysis_run_name.find("_vs_") == -1:
    logger.error(
        "an_analysis_run_name is invalid because no '_vs_' is present."
    )
    raise exception.IllegalArgumentError(
        "an_analysis_run_name is invalid because no '_vs_' is present."
    )

  # </editor-fold>

  tmp_separator_index: int = an_analysis_run_name.find("_vs_")

  tmp_prot_1: str = an_analysis_run_name[:tmp_separator_index]
  if tmp_prot_1.find(";") != -1:
    tmp_prot_1_name: str = tmp_prot_1[: tmp_prot_1.find(";")]
    tmp_prot_1_chains: list[str] = tmp_prot_1[tmp_prot_1.find(";") + 1 :].split(
        ","
    )
  else:
    tmp_prot_1_name: str = tmp_prot_1
    tmp_prot_1_chains: list[str] = []

  tmp_prot_2 = an_analysis_run_name[tmp_separator_index + 4 :]
  if tmp_prot_2.find(";") != -1:
    tmp_prot_2_name: str = tmp_prot_2[: tmp_prot_2.find(";")]
    tmp_prot_2_chains: list[str] = tmp_prot_2[tmp_prot_2.find(";") + 1 :].split(
        ","
    )
  else:
    tmp_prot_2_name: str = tmp_prot_2
    tmp_prot_2_chains: list[str] = []
  return [
      tmp_prot_1_name,
      tmp_prot_1_chains,
      tmp_prot_2_name,
      tmp_prot_2_chains,
  ]


def create_selection_strings_for_structure_alignment(
    a_protein: "protein.Protein",
    selected_chains: list,
) -> selection.Selection:
  """Creates the selection string for the structure alignment for a single protein.

  Args:
      a_protein (protein.Protein): The protein object to create selection strings from.
      selected_chains (list): A list of chain names to select.

  Returns:
      A selection object containing the selected chains from the protein.

  Raises:
      IllegalArgumentError: If `a_protein` or `selected_chains` is None.
  """
  # <editor-fold desc="Checks">
  if a_protein is None:
    logger.error("a_protein is None")
    raise exception.IllegalArgumentError("a_protein is None")
  if selected_chains is None:
    logger.error("selected_chains is None")
    raise exception.IllegalArgumentError("selected_chains is None")

  # </editor-fold>

  tmp_protein_selection = selection.Selection(a_protein.get_molecule_object())
  if len(selected_chains) != 0:
    protein_chains = []
    for tmp_chain in selected_chains:
      protein_chains.append(tmp_chain)
    tmp_protein_selection.set_selections_from_chains_ca(
        protein_util.get_chains_from_list_of_chain_names(
            a_protein, protein_chains
        ),
    )
  else:
    all_protein_chains = []
    for tmp_chain in a_protein.chains:
      if tmp_chain.chain_type != "non_protein_chain":
        all_protein_chains.append(tmp_chain)
    tmp_protein_selection.set_selections_from_chains_ca(all_protein_chains)
  return tmp_protein_selection


def get_highest_start_index(ref_index: int, model_index: int) -> int:
  """Gets the highest start index between the reference index and the model index.

  Args:
      ref_index (int): The reference index.
      model_index (int): The model index.

  Returns:
      An integer that is the highest start index between the reference index and the model index.

  Raises:
      exception.IllegalArgumentError: If an argument is None.
      ValueError: If the reference index is not greater than, less than, or equal to the model index.
  """
  # <editor-fold desc="Checks">
  if ref_index is None:
    logger.error("ref_index is None.")
    raise exception.IllegalArgumentError("ref_index is None.")
  if model_index is None:
    logger.error("model_index is None.")
    raise exception.IllegalArgumentError("model_index is None.")

  # </editor-fold>

  if ref_index > model_index:
    return ref_index
  if ref_index < model_index:
    return model_index
  if ref_index == model_index:
    return ref_index
  raise ValueError()
