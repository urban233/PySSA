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
from pyssa.internal.data_processing import data_transformer
from pyssa.logging_pyssa import log_handlers
from pyssa.internal.data_structures import protein
from pyssa.internal.data_structures import selection
from pyssa.util import protein_util, exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)

if TYPE_CHECKING:
    from pyssa.internal.data_structures import project
    from pyssa.internal.data_structures import settings


def transform_gui_input_to_practical_data(
    a_list_with_analysis_names: list,
    a_project: "project.Project",
    a_cutoff: float, cycles: int,
) -> list:
    """Transforms the input from the gui to a practical data basis that can be used to set up analysis runs.

    Raises:
        UnableToTransformDataForAnalysis: If the transformation process failed.
    """
    distance_analysis_runs = []
    try:
        for tmp_analysis_name in a_list_with_analysis_names:
            input_transformer = data_transformer.DistanceAnalysisDataTransformer(
                tmp_analysis_name,
                a_project,
                a_cutoff,
                cycles
            )
            protein_pair_for_analysis = input_transformer.transform_gui_input_to_distance_analysis_object()
            new_protein_pair = copy.deepcopy(protein_pair_for_analysis)
            distance_analysis_runs.append(new_protein_pair)
        logger.debug(
            f"These are the distance analysis runs, after the data transformation: {distance_analysis_runs}.",
        )
    except exception.IllegalArgumentError:
        logger.error("Transformation of data failed because an argument was illegal.")
        raise exception.UnableToTransformDataForAnalysisError("")
    except exception.UnableToTransformDataForAnalysisError:
        logger.error("Transformation of data failed because the transformation process failed.")
        raise exception.UnableToTransformDataForAnalysisError("")
    except Exception as e:
        logger.error(f"Unknown error: {e}")
        raise exception.UnableToTransformDataForAnalysisError("")
    return distance_analysis_runs


def split_analysis_run_name_in_protein_name_and_chain(an_analysis_run_name: str) -> list[Union[str, list[str]]]:
    """Splits the analysis run name into protein name and chain.

    Raises:
        IllegalArgumentException: If the argument is None or an empty string or does not contain the _vs_.
    """
    # <editor-fold desc="Checks">
    if an_analysis_run_name is None or an_analysis_run_name == "":
        logger.error(f"The argument 'an_analysis_run_name' is illegal: {an_analysis_run_name}!")
        raise exception.IllegalArgumentError("An argument is illegal.")
    if an_analysis_run_name.find("_vs_") == -1:
        logger.error("The argument 'an_analysis_run_name' is invalid because no '_vs_' is present.")
        raise exception.IllegalArgumentError("An argument is illegal")

    # </editor-fold>

    tmp_separator_index: int = an_analysis_run_name.find("_vs_")

    tmp_prot_1: str = an_analysis_run_name[:tmp_separator_index]
    if tmp_prot_1.find(";") != -1:
        tmp_prot_1_name: str = tmp_prot_1[: tmp_prot_1.find(";")]
        tmp_prot_1_chains: list[str] = tmp_prot_1[tmp_prot_1.find(";") + 1 :].split(",")
    else:
        tmp_prot_1_name: str = tmp_prot_1
        tmp_prot_1_chains: list[str] = []

    tmp_prot_2 = an_analysis_run_name[tmp_separator_index + 4 :]
    if tmp_prot_2.find(";") != -1:
        tmp_prot_2_name: str = tmp_prot_2[: tmp_prot_2.find(";")]
        tmp_prot_2_chains: list[str] = tmp_prot_2[tmp_prot_2.find(";") + 1 :].split(",")
    else:
        tmp_prot_2_name: str = tmp_prot_2
        tmp_prot_2_chains: list[str] = []
    return [tmp_prot_1_name, tmp_prot_1_chains, tmp_prot_2_name, tmp_prot_2_chains]


def create_selection_strings_for_structure_alignment(
    a_protein: "protein.Protein",
    selected_chains: list,
) -> selection.Selection:
    """Creates the selection string for the structure alignment for a single protein.

    Raises:
        IllegalArgumentError: If an argument is None.
    """
    # <editor-fold desc="Checks">
    if a_protein is None:
        logger.error(f"The argument 'a_value' is illegal: {a_protein}!")
        raise exception.IllegalArgumentError("")
    if selected_chains is None:
        logger.error(f"The argument 'a_value' is illegal: {selected_chains}!")
        raise exception.IllegalArgumentError("")

    # </editor-fold>

    tmp_protein_selection = selection.Selection(a_protein.get_molecule_object())
    if len(selected_chains) != 0:
        protein_chains = []
        for tmp_chain in selected_chains:
            protein_chains.append(tmp_chain)
        tmp_protein_selection.set_selections_from_chains_ca(
            protein_util.get_chains_from_list_of_chain_names(a_protein, protein_chains),
        )
    else:
        all_protein_chains = []
        for tmp_chain in a_protein.chains:
            if tmp_chain.chain_type != "non_protein_chain":
                all_protein_chains.append(tmp_chain)
        tmp_protein_selection.set_selections_from_chains_ca(all_protein_chains)
    return tmp_protein_selection


def get_highest_start_index(ref_index: int, model_index: int) -> int:
    """Gets the highest index of two indices.

    Args:
        ref_index: the index of the reference protein.
        model_index: the index of the model protein.
    """
    if ref_index > model_index:
        return ref_index
    if ref_index < model_index:
        return model_index
    if ref_index == model_index:
        return ref_index
    raise ValueError()
