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
"""This module contains helper functions for specific data transformations."""
import logging
import pathlib

from pyssa.internal.data_structures import sequence
from pyssa.logging_pyssa import log_handlers
from pyssa.internal.data_structures.data_classes import analysis_run_info
from pyssa.internal.data_structures import protein
from pyssa.internal.data_structures import protein_pair
from pyssa.internal.analysis_types import distance_analysis
from pyssa.util import analysis_util
from pyssa.util import protein_util
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pyssa.internal.data_structures import project
    from pyssa.internal.data_structures import settings

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def transform_protein_name_seq_tuple_to_sequence_obj(protein_name_seq_tuples: list[tuple[str, str]]) -> list[protein.Protein]:
    """

    Args:
        protein_name_seq_tuples:
            list of tuples which consists of the protein name and sequence
    Raises:
        ValueError: raised if an argument is illegal
    Returns:

    """
    # <editor-fold desc="Checks">
    if protein_name_seq_tuples[0][0] == "":
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")
    if protein_name_seq_tuples[0][1] == "":
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")

    # </editor-fold>

    protein_objects: list[protein.Protein] = []

    # create an empty dictionary
    groups = {}

    # loop over the tuples in the list
    for t in protein_name_seq_tuples:
        key = t[0]
        value = t[1]
        # if the key is not already in the dictionary, add it with an empty list
        if key not in groups:
            groups[key] = []
        # add the tuple's second element to the list associated with the key
        groups[key].append(value)

    for tmp_prot_name, tmp_seqs in groups.items():
        logger.debug("tmp_prot_name: %s, tmp_seqs: %r", tmp_prot_name, tmp_seqs)
        protein_with_sequence = protein.Protein(tmp_prot_name)
        protein_with_sequence.chains = []
        logger.debug("Chains of protein_with_sequence: %r", protein_with_sequence.chains)
        logger.debug("Name of protein_with_sequence: %s", protein_with_sequence.get_molecule_object())
        logger.debug("Memory address of protein_with_sequence: %r", protein_with_sequence)
        for tmp_seq in tmp_seqs:
            logger.debug("tmp_seq: %s, tmp_seqs: %r", tmp_seq, tmp_seqs)
            protein_seq = sequence.Sequence(tmp_prot_name, tmp_seq)
            protein_with_sequence.append_chain(chain_name="", chain_sequence=protein_seq, chain_type="protein_chain")
        logger.debug("Chains of protein_with_sequence after seq loop: %r", protein_with_sequence.chains)
        protein_with_sequence.add_chain_names_to_chains()
        protein_objects.append(protein_with_sequence)
    logger.debug(protein_objects)

    # last_protein_name = protein_name_seq_tuples[0][0]
    # protein_sequence = sequence.Sequence("", "")
    # protein_with_sequence = protein.Protein("generic")
    # if len(protein_name_seq_tuples) == 1:
    #     # single prediction
    #     logger.debug("Create protein object for single prediction.")
    #     protein_sequence.name = protein_name_seq_tuples[0][0]
    #     protein_sequence.sequence = protein_name_seq_tuples[0][1]
    #     logger.debug("protein_sequence.sequence: %s", protein_sequence.sequence)
    #     protein_with_sequence.set_molecule_object(protein_sequence.name)
    #     protein_seq = sequence.Sequence(protein_sequence.name, protein_sequence.sequence)
    #     protein_with_sequence.append_chain(chain_name="", chain_sequence=protein_seq, chain_type="protein_chain")
    #     protein_with_sequence.add_chain_names_to_chains()
    #     protein_objects.append(protein_with_sequence)
    # else:
    #     # batch prediction
    #     logger.debug("Create protein objects for batch prediction.")
    #     for tmp_protein_name_seq_tuple in protein_name_seq_tuples:
    #         logger.debug("tmp_protein_name_seq_tuple: %r", tmp_protein_name_seq_tuple)
    #         logger.debug("protein_name_seq_tuples: %r", protein_name_seq_tuples)
    #         current_protein_name = tmp_protein_name_seq_tuple[0]
    #         if last_protein_name == current_protein_name:
    #             protein_sequence.name = tmp_protein_name_seq_tuple[0]
    #             protein_sequence.sequence = tmp_protein_name_seq_tuple[1]
    #             protein_with_sequence = protein.Protein(protein_sequence.name)
    #             #protein_with_sequence.set_molecule_object(protein_sequence.name)
    #             logger.debug(protein_with_sequence.get_molecule_object())
    #             logger.debug(protein_name_seq_tuples)
    #             protein_seq = sequence.Sequence(protein_sequence.name, protein_sequence.sequence)
    #             protein_with_sequence.append_chain("", protein_seq, "protein_chain")
    #
    #             # protein_with_sequence.set_molecule_object(protein_name_seq_tuples[0])
    #             # logger.debug(protein_with_sequence.get_molecule_object())
    #             # logger.debug(protein_name_seq_tuples)
    #             # protein_seq = sequence.Sequence(protein_name_seq_tuples[0][0], protein_name_seq_tuples[0][1])
    #             # protein_with_sequence.append_chain("", protein_seq, "protein_chain")
    #         else:
    #             protein_with_sequence.add_chain_names_to_chains()
    #             protein_objects.append(protein_with_sequence)
    #             protein_sequence = sequence.Sequence("", "")
    #             protein_with_sequence = protein.Protein("generic")
    #             last_protein_name = current_protein_name
    # logger.debug(protein_objects)
    return protein_objects


# class DataTransformer:
#
#     def __init__(self,
#                  project,
#                  proteins_for_analysis: list[tuple[types.PROTEIN_ANALYSIS_INFO, types.PROTEIN_ANALYSIS_INFO]]):
#         self.project = project
#         self.proteins_for_analysis = proteins_for_analysis
#         print(self.project.proteins)
#     # def transform_to_analysis(self, project):
#     #     prot_1_name = self.ui.lbl_analysis_prot_struct_1.text().replace(".pdb", "")
#     #     prot_2_name = self.ui.lbl_analysis_prot_struct_2.text().replace(".pdb", "")
#     #     prot_1: protein.Protein = project.search_protein(prot_1_name)
#     #     if prot_1_name == prot_2_name:
#     #         prot_2: protein.Protein = prot_1.duplicate_protein()
#     #         prot_1.molecule_object = f"{prot_1.molecule_object}_1"
#     #         prot_2.molecule_object = f"{prot_2.molecule_object}_2"
#     #     else:
#     #         prot_2: protein.Protein = project.search_protein(prot_2_name)
#     #
#     #     prot_1_chains_selected = self.ui.list_analysis_ref_chains.selectedItems()
#     #     prot_1_chains = []
#     #     for tmp_chain in prot_1_chains_selected:
#     #         prot_1_chains.append(tmp_chain.text())
#     #     prot_1.set_chains(prot_1_chains)
#     #
#     #     prot_2_chains_selected = self.ui.list_analysis_model_chains.selectedItems()
#     #     prot_2_chains = []
#     #     for tmp_chain in prot_2_chains_selected:
#     #         prot_2_chains.append(tmp_chain.text())
#     #     prot_2.set_chains(prot_2_chains)
#     #
#     #     if len(prot_1.chains) != 0:
#     #         analysis_name = f"{prot_1.molecule_object};{prot_1_chains}_vs_{prot_2.molecule_object};{prot_2_chains}"
#     #         analysis_name = analysis_name.replace(";", "_")
#     #         analysis_name = analysis_name.replace(",", "_")
#     #         analysis_name = analysis_name.replace("[", "")
#     #         analysis_name = analysis_name.replace("]", "")
#     #         analysis_name = analysis_name.replace("'", "")
#     #     else:
#     #         analysis_name = f"{prot_1.molecule_object}_vs_{prot_2.molecule_object}"
#     #     export_dir = pathlib.Path(
#     #         f"{project.get_results_path()}/{analysis_name}")
#     #     return prot_1, prot_2, export_dir, analysis_name
#


class DistanceAnalysisDataTransformer:
    """This class is used to transform data from the gui to a manageable format."""

    # <editor-fold desc="Class attributes">
    """
    the name of a single analysis run
    """
    analysis_run_name: str
    """
    the current project of the main window
    """
    current_project: 'project.Project'
    """
    the settings of pyssa
    """
    settings: 'settings.Settings'
    """
    the information about the analysis run, includes the names and chains and analysis name
    """
    analysis_run: 'analysis_run_info.AnalysisRunInfo'
    """
    a tuple of two proteins
    """
    proteins: tuple['protein.Protein', 'protein.Protein']
    """
    the protein pair for the distance analysis
    """
    analysis_protein_pair: 'protein_pair.ProteinPair'
    """
    the distance analysis object
    """
    analysis_distance: 'distance_analysis.DistanceAnalysis'

    # </editor-fold>

    def __init__(self, analysis_run_name: str, app_project: 'project.Project', app_settings: 'settings.Settings'):
        self.analysis_run_name = analysis_run_name
        self.current_project = app_project
        self.settings = app_settings

    def _create_analysis_run_info(self):
        tmp_analysis_run_infos: list = analysis_util.split_analysis_run_name_in_protein_name_and_chain(self.analysis_run_name)
        self.analysis_run = analysis_run_info.AnalysisRunInfo(tmp_analysis_run_infos[0], tmp_analysis_run_infos[1],
                                                              tmp_analysis_run_infos[2], tmp_analysis_run_infos[3],
                                                              self.analysis_run_name)

    def _create_proteins_for_analysis(self):
        protein_1: protein.Protein = self.current_project.search_protein(self.analysis_run.get_protein_name_1())

        if self.analysis_run.are_protein_names_identical():
            protein_2: protein.Protein = protein_1.duplicate_protein()
            protein_1: protein.Protein = protein_2.duplicate_protein()
            protein_1.set_molecule_object(f"{protein_1.get_molecule_object()}_1")
            protein_2.set_molecule_object(f"{protein_2.get_molecule_object()}_2")
        else:
            protein_2: protein.Protein = self.current_project.search_protein(self.analysis_run.get_protein_name_2())
        self.proteins = (protein_1, protein_2)

    def _create_analysis_name(self):
        if len(self.proteins[0].chains) != 0:
            analysis_name = f"{self.proteins[0].get_molecule_object()};{self.analysis_run.protein_chains_1}_vs_{self.proteins[1].get_molecule_object()};{self.analysis_run.protein_chains_1}"
            analysis_name = analysis_name.replace(";", "_")
            analysis_name = analysis_name.replace(",", "_")
            analysis_name = analysis_name.replace("[", "")
            analysis_name = analysis_name.replace("]", "")
            analysis_name = analysis_name.replace("'", "")
        else:
            analysis_name = f"{self.proteins[0].get_molecule_object()}_vs_{self.proteins[1].get_molecule_object()}"
        return analysis_name

    def _create_protein_pair(self):
        self.analysis_protein_pair = protein_pair.ProteinPair(self.proteins[0], self.proteins[1])

    def _create_distance_analysis(self):
        self.analysis_protein_pair.set_distance_analysis(distance_analysis.DistanceAnalysis(self.analysis_protein_pair, self.settings))

    def _set_selection_for_analysis(self):
        self.analysis_protein_pair.distance_analysis.create_align_selections(
            analysis_util.create_selection_strings_for_structure_alignment(self.analysis_protein_pair.protein_1,
                                                                           self.analysis_run.protein_chains_1),
            analysis_util.create_selection_strings_for_structure_alignment(self.analysis_protein_pair.protein_2,
                                                                           self.analysis_run.protein_chains_2)
        )

    def transform_gui_input_to_distance_analysis_object(self):
        self._create_analysis_run_info()
        self._create_proteins_for_analysis()
        self._create_analysis_name()
        self._create_protein_pair()
        self._create_distance_analysis()
        self._set_selection_for_analysis()
        return self.analysis_protein_pair
