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
"""Module for all asynchronous functions used in the main presenter."""
import logging
from pyssa.internal.data_structures import project
from pyssa.internal.portal import protein_operations
from pyssa.logging_pyssa import log_handlers
from pyssa.util import exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


def check_chains_for_analysis(the_protein_1_name: str, the_protein_2_name: str, a_project: "project.Project") -> tuple:
    """Checks if the two proteins have only one chain.

    Args:
        the_protein_1_name (str): The name of the first protein from the analysis.
        the_protein_2_name (str): The name of the second protein from the analysis.
        a_project (project.Project): The current project.

    Returns:
        A tuple with the form ("result", tmp_is_only_one_chain, tmp_analysis_run_name, tmp_protein_1, tmp_protein_2).
    
    Raises:
        exception.IllegalArgumentError: If the_protein_1_name is either None or an empty string.
        exception.IllegalArgumentError: If the_protein_2_name is either None or an empty string.
        exception.IllegalArgumentError: If a_project is None.
    """
    # <editor-fold desc="Checks">
    if the_protein_1_name is None or the_protein_1_name == "":
        logger.error("the_protein_1_name is either None or an empty string.")
        raise exception.IllegalArgumentError("the_protein_1_name is either None or an empty string.")
    if the_protein_2_name is None or the_protein_2_name == "":
        logger.error("the_protein_2_name is either None or an empty string.")
        raise exception.IllegalArgumentError("the_protein_2_name is either None or an empty string.")
    if a_project is None:
        logger.error("a_project is None.")
        raise exception.IllegalArgumentError("a_project is None.")
    
    # </editor-fold>
    
    tmp_is_only_one_chain: bool = False
    tmp_analysis_run_name: str = ""
    tmp_protein_1 = a_project.search_protein(the_protein_1_name)
    tmp_protein_2 = a_project.search_protein(the_protein_2_name)
    if len(tmp_protein_1.get_protein_sequences()) == 1: # or len(tmp_protein_2.get_protein_sequences()) == 1:
        tmp_is_only_one_chain = True
        tmp_protein_1_first_protein_chain_letter = protein_operations.get_chain_letter_of_first_protein_sequence(
            tmp_protein_1.chains,
        )
        tmp_protein_2_first_protein_chain_letter = protein_operations.get_chain_letter_of_first_protein_sequence(
            tmp_protein_2.chains,
        )
        tmp_analysis_run_name = (
            f"{tmp_protein_1.get_molecule_object()};{tmp_protein_1_first_protein_chain_letter}"
            f"_vs_{tmp_protein_2.get_molecule_object()};{tmp_protein_2_first_protein_chain_letter}"
        )
    return "result", tmp_is_only_one_chain, tmp_analysis_run_name, tmp_protein_1, tmp_protein_2


def check_chains_for_subsequent_analysis(
    the_protein_1_name: str,
    the_protein_2_name: str,
    a_project: "project.Project",
    a_list_with_proteins_to_predict: list[str],
) -> tuple:
    """Checks if the two proteins have only one chain.

    Args:
        the_protein_1_name (str): The name of the first protein from the analysis.
        the_protein_2_name (str): The name of the second protein from the analysis.
        a_project (project.Project): The current project.
        a_list_with_proteins_to_predict (list[str]): A list of all proteins that should be predicted.
    
    Raises:
        exception.IllegalArgumentError: If the_protein_1_name is either None or an empty string.
        exception.IllegalArgumentError: If the_protein_2_name is either None or an empty string.
        exception.IllegalArgumentError: If a_project is None.
        exception.IllegalArgumentError: If a_list_with_proteins_to_predict is None.
    """
    # <editor-fold desc="Checks">
    if the_protein_1_name is None or the_protein_1_name == "":
        logger.error("the_protein_1_name is either None or an empty string.")
        raise exception.IllegalArgumentError("the_protein_1_name is either None or an empty string.")
    if the_protein_2_name is None or the_protein_2_name == "":
        logger.error("the_protein_2_name is either None or an empty string.")
        raise exception.IllegalArgumentError("the_protein_2_name is either None or an empty string.")
    if a_project is None:
        logger.error("a_project is None.")
        raise exception.IllegalArgumentError("a_project is None.")
    if a_list_with_proteins_to_predict is None:
        logger.error("a_list_with_proteins_to_predict is None.")
        raise exception.IllegalArgumentError("a_list_with_proteins_to_predict is None.")
    
    # </editor-fold>
    
    tmp_protein_1_only_one_chain: bool = False
    tmp_protein_2_only_one_chain: bool = False
    tmp_sub_analysis_name_1: str = ""
    tmp_sub_analysis_name_2: str = ""
    if the_protein_1_name not in a_list_with_proteins_to_predict:
        tmp_protein_1 = a_project.search_protein(the_protein_1_name)
        if len(tmp_protein_1.chains) == 1:
            tmp_protein_1_only_one_chain = True
            tmp_sub_analysis_name_1 = f"{tmp_protein_1.get_molecule_object()};{tmp_protein_1.chains[0].chain_letter}"
    else:
        tmp_protein_1_only_one_chain = True
        tmp_sub_analysis_name_1 = f"{the_protein_1_name};A"

    if the_protein_2_name not in a_list_with_proteins_to_predict:
        tmp_protein_2 = a_project.search_protein(the_protein_2_name)
        if len(tmp_protein_2.chains) == 1:
            tmp_protein_2_only_one_chain = True
            tmp_sub_analysis_name_2 = f"{tmp_protein_2.get_molecule_object()};{tmp_protein_2.chains[0].chain_letter}"
    else:
        tmp_protein_2_only_one_chain = True
        tmp_sub_analysis_name_2 = f"{the_protein_2_name};A"

    if tmp_protein_1_only_one_chain and tmp_protein_2_only_one_chain:
        tmp_analysis_run_name = f"{tmp_sub_analysis_name_1}_vs_{tmp_sub_analysis_name_2}"
    else:
        tmp_analysis_run_name = ""
    return "result", tmp_analysis_run_name
