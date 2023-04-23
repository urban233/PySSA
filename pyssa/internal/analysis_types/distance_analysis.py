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
"""This module contains the distance analysis class."""
import os
import base64
import logging
import pathlib
import json
from pymol import cmd
from typing import TYPE_CHECKING

from xml.etree import ElementTree
from pyssa.io_pyssa.xml_pyssa import element_names
from pyssa.io_pyssa.xml_pyssa import attribute_names
from io_pyssa import filesystem_io
from pyssa.logging_pyssa import log_handlers
from pyssa.util import protein_pair_util
from pyssa.internal.portal import protein_pair_operations
from pyssa.io_pyssa import path_util
from pyssa.internal.portal import pymol_io
from pyssa.util import pyssa_keys
from pyssa.util import analysis_util
from pyssa.internal.data_structures import results

if TYPE_CHECKING:
    from pyssa.internal.data_structures import protein_pair
    from pyssa.internal.data_structures import settings
    from pyssa.internal.data_structures import selection
    from pyssa.internal.data_structures import project

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class DistanceAnalysis:
    """This class contains information about the distance analysis"""

    # <editor-fold desc="Class attributes">
    """
    a pair of proteins which get analyzed 
    """
    _protein_pair_for_analysis: 'protein_pair.ProteinPair'
    """
    the settings from pyssa
    """
    app_settings: 'settings.Settings'
    """
    the cutoff value for the align command from PyMOL
    """
    cutoff: float
    """
    the number of refinement cycles for the align command from PyMOL
    """
    cycles: int
    """
    the filename of the alignement file which gets created during the align command
    """
    alignment_file_name: str
    """
    a directory where all results related to the protein will be stored
    """
    export_dirname: pathlib.Path
    """
    the session filepath for the PyMOL session of the analysis
    """
    pymol_session_filepath: path_util.FilePath
    """
    
    """
    analysis_results: 'results.DistanceAnalysisResults' = None
    """
    """
    rmsd_dict: dict = {
        'rmsd': 0.0,
        'aligned_aa': 0,
    }

    # </editor-fold>

    def __init__(self, protein_pair_for_analysis: 'protein_pair.ProteinPair', app_settings: 'settings.Settings'):
        """Constructor

        Args:
            protein_pair_for_analysis:
                a pair of proteins which get analyzed
            app_settings:
                the settings from pyssa
            export_dirname:
                a directory where all results will be stored
        """
        self.distance_analysis_data = []
        self._protein_pair_for_analysis: protein_pair.ProteinPair = protein_pair_for_analysis
        self.name = f"dist_analysis_{self._protein_pair_for_analysis.name}"
        self.app_settings: settings.Settings = app_settings
        self.cutoff: float = app_settings.cutoff
        self.cycles: int = app_settings.cycles
        self.alignment_file_name = f"{self._protein_pair_for_analysis.name}_alignment"
        # self.distance_analysis_subdirs = {
        #     pyssa_keys.DISTANCE_ANALYSIS_SUBDIR: pathlib.Path(f"{distance_analysis_dirname}"),
        #     pyssa_keys.DISTANCE_ANALYSIS_SESSION_SUBDIR: pathlib.Path(f"{distance_analysis_dirname}/session"),
        #     pyssa_keys.DISTANCE_ANALYSIS_RESULTS_SUBDIR: pathlib.Path(f"{distance_analysis_dirname}/results"),
        #     pyssa_keys.DISTANCE_ANALYSIS_OBJECTS_SUBDIR: pathlib.Path(f"{distance_analysis_dirname}/.objects"),
        # }
        # for key in self.distance_analysis_subdirs:
        #     if not os.path.exists(self.distance_analysis_subdirs.get(key)):
        #         os.mkdir(self.distance_analysis_subdirs.get(key))
        #
        # self.export_dirname = self.distance_analysis_subdirs.get(pyssa_keys.DISTANCE_ANALYSIS_RESULTS_SUBDIR)
        # self.pymol_session_filepath = path_util.FilePath(f"{self.distance_analysis_subdirs.get(pyssa_keys.DISTANCE_ANALYSIS_SESSION_SUBDIR)}/{self.protein_pair_for_analysis.name}_analysis_session.pse")

    def save_distance_analysis_session(self) -> None:
        """This function saves the pymol session of the Protein pair distance analysis.

        """
        protein_pair_operations.save_session_of_protein_pair(self._protein_pair_for_analysis.name)

    def create_align_selections(self,
                                protein_1_selection: 'selection.Selection',
                                protein_2_selection: 'selection.Selection') -> None:
        """This function creates the selection which are needed for the align command.

        """
        # <editor-fold desc="Checks">
        if protein_1_selection.molecule_object != self._protein_pair_for_analysis.protein_1.get_molecule_object():
            logger.error("Selection is illegal.")
            raise ValueError("Selection is illegal.")
        if protein_2_selection.molecule_object != self._protein_pair_for_analysis.protein_2.get_molecule_object():
            logger.error("Selection is illegal.")
            raise ValueError("Selection is illegal.")

        # </editor-fold>

        self._protein_pair_for_analysis.protein_1.pymol_selection = protein_1_selection
        self._protein_pair_for_analysis.protein_2.pymol_selection = protein_2_selection

    def align_protein_pair_for_analysis(self) -> tuple:
        """This function aligns the protein pair with the PyMOL align command.

        Returns:
            a tuple with the rmsd and aligned amino acids
        """
        results = protein_pair_operations.align_protein_pair(self._protein_pair_for_analysis.protein_1.pymol_selection.selection_string,
                                                             self._protein_pair_for_analysis.protein_2.pymol_selection.selection_string,
                                                             self.alignment_file_name)
        # save the align object from pymol as alignment file
        # if not os.path.exists(pathlib.Path(f"{self.protein_pair_for_analysis.analysis_results}/alignment_files")):
        #     os.mkdir(pathlib.Path(f"{self.protein_pair_for_analysis.analysis_results}/alignment_files"))
        # cmd.save(pathlib.Path(f"{self.protein_pair_for_analysis.analysis_results}/alignment_files/{self.alignment_file_name}.aln"))
        return results[0], results[1]

    def do_analysis_in_pymol(self, app_project: 'project.Project'):
        """This function does the distance analysis of the protein pair.

        """
        self._protein_pair_for_analysis.load_protein_pair_in_pymol() # This creates a new pymol session
        self._protein_pair_for_analysis.color_protein_pair()
        align_results = self.align_protein_pair_for_analysis()
        self.rmsd_dict = {
            "rmsd": str(round(align_results[0], 2)),
            "aligned_residues": str(align_results[1]),
        }
        # rmsd_file = open(pathlib.Path(f"{self.protein_pair_for_analysis.analysis_results}/rmsd.json"), "w", encoding="utf-8")
        # json.dump(rmsd_dict, rmsd_file, indent=4)
        # rmsd_file.close()
        # extract single chain selections into a list
        protein_1_chain_selections = self._protein_pair_for_analysis.protein_1.pymol_selection.selection_string.split(",")
        protein_2_chain_selections = self._protein_pair_for_analysis.protein_2.pymol_selection.selection_string.split(",")
        # get the pymol chempy objects of the selections
        protein_1_ca_pymol_objects = []
        for tmp_selection_string in protein_1_chain_selections:
            protein_1_ca_pymol_objects.append(cmd.get_model(tmp_selection_string))
            logger.debug("Prot1")
            logger.debug(tmp_selection_string)
        protein_2_ca_pymol_objects = []
        for tmp_selection_string in protein_2_chain_selections:
            protein_2_ca_pymol_objects.append(cmd.get_model(tmp_selection_string))
            logger.debug("Prot2")
            logger.debug(tmp_selection_string)

        # create list which consists of a tuple (prot_1_ca, prot_2_ca)
        pymol_ca_object_pairs = []
        for i in range(len(protein_1_ca_pymol_objects)):
            pymol_ca_object_pairs.append((protein_1_ca_pymol_objects[i], protein_2_ca_pymol_objects[i]))
        self.distance_analysis_data = []
        for tmp_pair in pymol_ca_object_pairs:
            prot_ref_ca = tmp_pair[0]
            prot_model_ca = tmp_pair[1]

            ref_count, ref_start_index = analysis_util.count_atoms_in_selection(prot_ref_ca)
            model_count, model_start_index = analysis_util.count_atoms_in_selection(prot_model_ca)

            count = analysis_util.get_lowest_count(ref_count, model_count)
            distances = protein_pair_util.calculate_distance_between_ca_atoms(
                count, prot_ref_ca, prot_model_ca,
                analysis_util.get_ref_gap(ref_start_index, model_start_index),
                analysis_util.get_model_gap(ref_start_index, model_start_index),
                self._protein_pair_for_analysis.protein_1.get_molecule_object(),
                self._protein_pair_for_analysis.protein_2.get_molecule_object(),
            )
            self.distance_analysis_data.append(distances)
        #self.save_distance_analysis_session()
        self.analysis_results = results.DistanceAnalysisResults(
            self.distance_analysis_data,
            pymol_io.convert_pymol_session_to_base64_string(self._protein_pair_for_analysis.name),
            self.rmsd_dict['rmsd'],
            self.rmsd_dict['aligned_residues']
        )
        #return self._protein_pair_for_analysis

    def serialize_distance_analysis(self, xml_distance_analysis_element):
        """This function serialize the protein pair object

        """
        tmp_distance_analysis = ElementTree.SubElement(xml_distance_analysis_element, element_names.DISTANCE_ANALYSIS)
        tmp_distance_analysis.set(attribute_names.DISTANCE_ANALYSIS_NAME, str(self.name))
        tmp_distance_analysis.set(attribute_names.DISTANCE_ANALYSIS_CUTOFF, str(self.cutoff))
        tmp_distance_analysis.set(attribute_names.DISTANCE_ANALYSIS_CYCLES, str(self.cycles))

        logger.debug(self.name)
        logger.debug(self.analysis_results)
        if len(self.distance_analysis_data) == 0:
            self.analysis_results.serialize_distance_analysis_results(tmp_distance_analysis)
        else:
            tmp_results_data = ElementTree.SubElement(tmp_distance_analysis, element_names.DISTANCE_ANALYSIS_RESULTS)
            tmp_results_data.set(attribute_names.DISTANCE_ANALYSIS_RMSD, str(self.rmsd_dict['rmsd']))
            tmp_results_data.set(attribute_names.DISTANCE_ANALYSIS_ALIGNED_AA, str(self.rmsd_dict['aligned_residues']))

            for distances in self.distance_analysis_data:
                tmp_distance_data = ElementTree.SubElement(tmp_results_data, element_names.DISTANCE_ANALYSIS_DISTANCE_RESULTS)
                tmp_index_data = ElementTree.SubElement(tmp_distance_data, element_names.DISTANCE_ANALYSIS_INDEX_LIST)
                tmp_index_data.text = str(distances['index'].tolist())
                tmp_prot_1_chain_data = ElementTree.SubElement(tmp_distance_data, element_names.DISTANCE_ANALYSIS_PROT_1_CHAIN_LIST)
                tmp_prot_1_chain_data.text = str(distances['ref_chain'].tolist())
                tmp_prot_1_position_data = ElementTree.SubElement(tmp_distance_data, element_names.DISTANCE_ANALYSIS_PROT_1_POSITION_LIST)
                tmp_prot_1_position_data.text = str(distances['ref_pos'].tolist())
                tmp_prot_1_residue_data = ElementTree.SubElement(tmp_distance_data, element_names.DISTANCE_ANALYSIS_PROT_1_RESIDUE_LIST)
                tmp_prot_1_residue_data.text = str(distances['ref_resi'].tolist())
                tmp_prot_2_chain_data = ElementTree.SubElement(tmp_distance_data, element_names.DISTANCE_ANALYSIS_PROT_2_CHAIN_LIST)
                tmp_prot_2_chain_data.text = str(distances['model_chain'].tolist())
                tmp_prot_2_position_data = ElementTree.SubElement(tmp_distance_data, element_names.DISTANCE_ANALYSIS_PROT_2_POSITION_LIST)
                tmp_prot_2_position_data.text = str(distances['model_pos'].tolist())
                tmp_prot_2_residue_data = ElementTree.SubElement(tmp_distance_data, element_names.DISTANCE_ANALYSIS_PROT_2_RESIDUE_LIST)
                tmp_prot_2_residue_data.text = str(distances['model_resi'].tolist())
                tmp_distances_data = ElementTree.SubElement(tmp_distance_data, element_names.DISTANCE_ANALYSIS_DISTANCES_LIST)
                tmp_distances_data.text = str(distances['distance'].tolist())

        tmp_session_data = ElementTree.SubElement(tmp_distance_analysis, element_names.DISTANCE_ANALYSIS_SESSION)
        tmp_session_data.set(attribute_names.PROTEIN_PAIR_SESSION, pymol_io.convert_pymol_session_to_base64_string(self.name))
