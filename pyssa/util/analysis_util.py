#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2022
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
"""This module contains helper function for the analysis process."""
import pathlib
import logging
from typing import Union
from pyssa.logging_pyssa import log_handlers
from pyssa.internal.data_structures.data_classes import protein_analysis_info
from pyssa.internal.data_structures import protein
from pyssa.internal.data_structures import selection
from pyssa.util import protein_util, exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


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
        tmp_prot_1_name: str = tmp_prot_1[:tmp_prot_1.find(";")]
        tmp_prot_1_chains: list[str] = tmp_prot_1[tmp_prot_1.find(";") + 1:].split(",")
    else:
        tmp_prot_1_name: str = tmp_prot_1
        tmp_prot_1_chains: list[str] = []

    tmp_prot_2 = an_analysis_run_name[tmp_separator_index + 4:]
    if tmp_prot_2.find(";") != -1:
        tmp_prot_2_name: str = tmp_prot_2[:tmp_prot_2.find(";")]
        tmp_prot_2_chains: list[str] = tmp_prot_2[tmp_prot_2.find(";") + 1:].split(",")
    else:
        tmp_prot_2_name: str = tmp_prot_2
        tmp_prot_2_chains: list[str] = []
    return [tmp_prot_1_name, tmp_prot_1_chains, tmp_prot_2_name, tmp_prot_2_chains]


def create_selection_strings_for_structure_alignment(a_protein: 'protein.Protein', selected_chains: list) -> selection.Selection:
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
            protein_util.get_chains_from_list_of_chain_names(a_protein, protein_chains))
        logger.debug(
            f"This is one selection created with <create_selection_strings_for_structure_alignment>: {tmp_protein_selection.selection_string}")
    else:
        all_protein_chains = []
        for tmp_chain in a_protein.chains:
            if tmp_chain.chain_type != "non_protein_chain":
                all_protein_chains.append(tmp_chain)
        tmp_protein_selection.set_selections_from_chains_ca(all_protein_chains)
    return tmp_protein_selection


def transform_data_for_analysis(proteins_for_analysis, project) -> list[tuple['protein.Protein', 'protein.Protein', pathlib.Path, str]]:
    """This function transforms the data from a protein list to a usable format for the analysis algorithm

    Returns:
        a list which consists of the two proteins, the export path for the results and the name of the analysis run
    """
    tmp_data = []
    for protein_pair_tuple in proteins_for_analysis:
        prot_1_name = protein_pair_tuple[0].protein_name.replace(".pdb", "")
        prot_2_name = protein_pair_tuple[1].protein_name.replace(".pdb", "")
        # create new protein objects
        prot_1: protein.Protein = project.search_protein(prot_1_name)
        if prot_1_name == prot_2_name:
            prot_2: protein.Protein = prot_1.duplicate_protein()
            prot_1: protein.Protein = prot_2.duplicate_protein()
            prot_1.molecule_object = f"{prot_1.molecule_object}_1"
            prot_2.molecule_object = f"{prot_2.molecule_object}_2"
        else:
            prot_2: protein.Protein = project.search_protein(prot_2_name)

        create_selection_strings_for_structure_alignment(prot_1, protein_pair_tuple)
        create_selection_strings_for_structure_alignment(prot_2, protein_pair_tuple)
        # set selection
        prot_1_selection = selection.Selection(prot_1.molecule_object)
        prot_1_chains_selected = protein_pair_tuple[0].protein_chains
        if prot_1_chains_selected is not None:
            prot_1_chains = []
            for tmp_chain in prot_1_chains_selected:
                prot_1_chains.append(tmp_chain)
            prot_1_selection.set_selections_from_chains_ca(protein_util.get_chains_from_list_of_chain_names(prot_1, prot_1_chains))
        else:
            prot_1_selection.set_selections_without_chains_ca()
        prot_1.pymol_selection.selection_string = prot_1_selection.selection_string



        prot_2_chains_selected = protein_pair_tuple[1].protein_chains
        if prot_2_chains_selected is not None:
            prot_2_chains = []
            for tmp_chain in prot_2_chains_selected:
                prot_2_chains.append(tmp_chain)
            prot_2.set_chains(prot_2_chains)
        else:
            prot_2_chains = []

        if len(prot_1.chains) != 0:
            analysis_name = f"{prot_1.molecule_object};{prot_1_chains}_vs_{prot_2.molecule_object};{prot_2_chains}"
            analysis_name = analysis_name.replace(";", "_")
            analysis_name = analysis_name.replace(",", "_")
            analysis_name = analysis_name.replace("[", "")
            analysis_name = analysis_name.replace("]", "")
            analysis_name = analysis_name.replace("'", "")
        else:
            analysis_name = f"{prot_1.molecule_object}_vs_{prot_2.molecule_object}"
        export_dir = pathlib.Path(
            f"{project.get_protein_pairs_path()}/{analysis_name}")
        transformed_data: tuple = prot_1, prot_2, export_dir, analysis_name
        tmp_data.append(transformed_data)
    return tmp_data


def count_atoms_in_selection(pymol_obj):
    """This function counts the atoms within a pymol molecule object."""
    count = 0
    first_seq_index = 0
    for atom in pymol_obj.atom:
        if count == 0:
            first_seq_index = atom.resi
        count += 1
    return count, int(first_seq_index)


def get_highest_start_index(ref_index, model_index):
    if ref_index > model_index:
        return ref_index
    if ref_index < model_index:
        return model_index
    if ref_index == model_index:
        return ref_index
    else:
        raise ValueError


def get_ref_gap(a_ref_index, a_model_index):
    return get_highest_start_index(a_ref_index, a_model_index) - a_ref_index


def get_model_gap(a_ref_index, a_model_index):
    return get_highest_start_index(a_ref_index, a_model_index) - a_model_index


def get_lowest_count(a_ref_count, a_model_count):
    if a_ref_count > a_model_count:
        return a_model_count
    if a_ref_count < a_model_count:
        return a_ref_count
    if a_ref_count == a_model_count:
        return a_ref_count
    else:
        raise ValueError
