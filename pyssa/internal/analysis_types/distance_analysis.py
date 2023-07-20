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
import logging
import pathlib
import numpy as np
from pymol import cmd
from typing import TYPE_CHECKING
from xml.etree import ElementTree
from pyssa.io_pyssa.xml_pyssa import element_names
from pyssa.io_pyssa.xml_pyssa import attribute_names
from pyssa.logging_pyssa import log_handlers
from pyssa.util import protein_pair_util
from pyssa.internal.portal import protein_pair_operations
from pyssa.io_pyssa import path_util
from pyssa.internal.portal import pymol_io
from pyssa.internal.portal import graphic_operations
from pyssa.util import constants
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
    """
    figure_size: tuple[float, float]
    """
    the filename of the alignement file which gets created during the align command
    """
    alignment_file_name: str = "aln"
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
    distance_analysis_data: dict[str, np.ndarray] = {}
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
        self._protein_pair_for_analysis: protein_pair.ProteinPair = protein_pair_for_analysis
        self.name = f"dist_analysis_{self._protein_pair_for_analysis.name}"
        self.app_settings: settings.Settings = app_settings
        self.cutoff: float = app_settings.cutoff
        self.cycles: int = app_settings.cycles
        self.figure_size = (11.0, 6.0)
        #self.alignment_file_name = f"{self._protein_pair_for_analysis.name}_alignment"
        self.alignment_file_name = "aln"

    def get_protein_pair(self):
        return self._protein_pair_for_analysis

    def save_distance_analysis_session(self) -> None:
        """This function saves the pymol session of the Protein pair distance analysis.

        """
        protein_pair_operations.save_session_of_protein_pair(self._protein_pair_for_analysis.name)

    def create_align_selections(self,
                                protein_1_selection: 'selection.Selection',
                                protein_2_selection: 'selection.Selection') -> None:
        """This function creates the selection which are needed for the align command.

        """
        logger.debug(f"1st argument of <create_align_selections>: {protein_1_selection.selection_string} {protein_1_selection}")
        logger.debug(f"2nd argument of <create_align_selections>: {protein_2_selection.selection_string} {protein_2_selection}")
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
        logger.debug(f"Prot 1 sele in <create_align_selections>: {self._protein_pair_for_analysis.protein_1.pymol_selection.selection_string} {self._protein_pair_for_analysis.protein_1.pymol_selection}")
        logger.debug(
            f"Prot 2 sele in <create_align_selections>: {self._protein_pair_for_analysis.protein_2.pymol_selection.selection_string} {self._protein_pair_for_analysis.protein_2.pymol_selection}")

    def align_protein_pair_for_analysis(self) -> tuple:
        """This function aligns the protein pair with the PyMOL align command.

        Returns:
            a tuple with the rmsd and aligned amino acids
        """
        logger.debug(
            f"Prot 1 sele in <align_protein_pair_for_analysis>: {self._protein_pair_for_analysis.protein_1.pymol_selection.selection_string} {self._protein_pair_for_analysis.protein_1.pymol_selection}")
        logger.debug(
            f"Prot 2 sele in <align_protein_pair_for_analysis>: {self._protein_pair_for_analysis.protein_2.pymol_selection.selection_string} {self._protein_pair_for_analysis.protein_2.pymol_selection}")
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
        logger.info("Start of do_analysis_in_pymol() method.")
        self._protein_pair_for_analysis.load_protein_pair_in_pymol()
        logger.info(f"Loaded protein pair: {self._protein_pair_for_analysis.name} in pymol session.")
        self._protein_pair_for_analysis.color_protein_pair()
        logger.info(f"Colored protein pair: {self._protein_pair_for_analysis.name} in pymol session.")

        logger.debug(f"Protein 1 selection string: {self._protein_pair_for_analysis.protein_1.pymol_selection.selection_string} {self._protein_pair_for_analysis.protein_1.pymol_selection}")
        logger.debug(f"Protein 2 selection string: {self._protein_pair_for_analysis.protein_2.pymol_selection.selection_string} {self._protein_pair_for_analysis.protein_2.pymol_selection}")
        # # extract single chain selections into a list
        # protein_1_chain_selections = self._protein_pair_for_analysis.protein_1.pymol_selection.selection_string.split(",")
        # protein_2_chain_selections = self._protein_pair_for_analysis.protein_2.pymol_selection.selection_string.split(",")
        # # get the pymol chempy objects of the selections
        # protein_1_ca_pymol_objects = []
        # for tmp_selection_string in protein_1_chain_selections:
        #     protein_1_ca_pymol_objects.append(cmd.get_model(tmp_selection_string))
        #     logger.debug(f"Prot1: selection {tmp_selection_string}")
        # protein_2_ca_pymol_objects = []
        # for tmp_selection_string in protein_2_chain_selections:
        #     protein_2_ca_pymol_objects.append(cmd.get_model(tmp_selection_string))
        #     logger.debug(f"Prot2: selection {tmp_selection_string}")

        align_results = self.align_protein_pair_for_analysis()
        logger.info(f"Aligned protein pair: {self._protein_pair_for_analysis.name} in pymol session.")
        self.rmsd_dict = {
            "rmsd": str(round(align_results[0], 2)),
            "aligned_residues": str(align_results[1]),
        }

        distances = protein_pair_util.calculate_distance_between_ca_atoms(
            self._protein_pair_for_analysis.protein_1.get_molecule_object(),
            self._protein_pair_for_analysis.protein_2.get_molecule_object(),
        )
        logger.info(f"Calculated distances of protein pair: {self._protein_pair_for_analysis.name}.")
        self.distance_analysis_data = distances
        self.analysis_results = results.DistanceAnalysisResults(
            self.distance_analysis_data,
            pymol_io.convert_pymol_session_to_base64_string(self._protein_pair_for_analysis.name),
            self.rmsd_dict['rmsd'],
            self.rmsd_dict['aligned_residues']
        )

    def take_image_of_protein_pair(self,
                                   representation: str,
                                   filename: str,
                                   selection: str = "",
                                   ray_shadows: bool = False,
                                   opaque_background: int = 0,
                                   take_images=True) -> None:
        """This function takes an image of the whole Protein/Protein pair.

        Note:
            The png file will be saved under the relative path
            (if export_data_dir = "data/results"):
            ``data/results/images``

        Args:
            alignment_filename (str):
                name of the alignment file from the structure alignment
            representation (str):
                defines the type of molecular representation
                like cartoon or ribbon
            filename (str):
                name of the png image file
            selection (str, optional):
                the atoms which MUST NOT displayed in the image
            ray_shadows (bool, optional):
                false if no shadows, true if shadows should be displayed
            opaque_background (int, optional):
                0 for a transparent background and 1 for a white background

        Raises:
            ValueError: If opaque_background is not 0 or 1.
        """
        # determine the option for ray_shadows
        if not ray_shadows:
            opt_ray_shadows: str = "off"
        else:
            opt_ray_shadows: str = "on"

        REPRESENTATION: str = "cartoon"
        cmd.show(REPRESENTATION)

        if selection != "":
            cmd.hide(representation, selection)

        aln_obj_representation: str = "cgo"
        cmd.hide(aln_obj_representation, "all")
        cmd.orient()
        cmd.center()
        # set image parameters
        cmd.bg_color(constants.PYMOL_DEFAULT_BACKGROUND_COLOR)
        cmd.set("ray_trace_mode", constants.PYMOL_DEFAULT_RAY_TRACE_MODE)
        cmd.set("antialias", constants.PYMOL_DEFAULT_ANTIALIAS)
        cmd.set("ray_shadows", opt_ray_shadows)
        cmd.set('ray_opaque_background', opaque_background)
        graphic_operations.setup_default_session_graphic_settings()
        self._protein_pair_for_analysis.color_protein_pair()

        cmd.scene(key=f"{self._protein_pair_for_analysis.protein_1.get_molecule_object()}-"
                      f"{self._protein_pair_for_analysis.protein_2.get_molecule_object()}",
                  action="store")

        if take_images is True:
            # check if path exists where the data will be exported,
            # if not the directory will be created
            if not os.path.exists(f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR}"):
                os.mkdir(f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR}")
            # save image as 300 dpi png image
            cmd.ray(2400, 2400, renderer=0)
            cmd.png(f'{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR}/{filename}.png',
                    dpi=300)

    def take_image_of_interesting_regions(self,
                                          cutoff: float,
                                          filename: str,
                                          ray_shadows: bool = False,
                                          opaque_background: int = 0,
                                          take_images=True):
        """This function takes images of interesting regions of the alignment

        Args:
            cutoff (float):
                defines a border of which the specific regions begin,
                if the distance is greater than the cutoff, the amino acid is
                categorized as "interesting"
            filename (str):
                name of the png image file
            ray_shadows (bool, optional):
                false if no shadows, true if shadows should be displayed
            opaque_background (int, optional):
                0 for a transparent background and 1 for a white background

        """

        # set default parameters
        cmd.set("label_size", 14)
        cmd.set("label_font_id", 13)
        cmd.set("label_color", "hotpink")
        cmd.set("depth_cue", 0)
        # cmd.set("fog_start", 0.6)
        cmd.bg_color(constants.PYMOL_DEFAULT_BACKGROUND_COLOR)

        j: int = 0
        measurement_obj: str = f"measure{j}"
        REPRESENTATION: str = "ribbon"

        cmd.hide("cartoon", "all")
        cmd.show(REPRESENTATION, "all")
        graphic_operations.setup_default_session_graphic_settings()
        self._protein_pair_for_analysis.color_protein_pair()

        i: int = 0
        for distance_value in self.analysis_results.distance_data.get("distance"):
            if float(distance_value) > float(cutoff):
                # here algorithm for image
                ref_pos_array: np.ndarray = self.analysis_results.distance_data.get("ref_pos")
                ref_pos: int = ref_pos_array[i]
                ref_chain_array: np.ndarray = self.analysis_results.distance_data.get("ref_chain")
                ref_chain: str = ref_chain_array[i]
                model_pos_array: np.ndarray = self.analysis_results.distance_data.get("model_pos")
                model_pos: int = model_pos_array[i]
                model_chain_array: np.ndarray = self.analysis_results.distance_data.get("model_chain")
                model_chain: str = model_chain_array[i]

                measurement_obj: str = f"measure{j}"
                # create two atoms for the get_distance command
                atom1: str = f"/{self._protein_pair_for_analysis.protein_1.get_molecule_object()}//" \
                             f"{ref_chain}/{ref_pos}/CA"
                atom2: str = f"/{self._protein_pair_for_analysis.protein_2.get_molecule_object()}//" \
                             f"{model_chain}/{model_pos}/CA"
                # zoom to reference amino acid
                cmd.zoom(f"/{self._protein_pair_for_analysis.protein_1.get_molecule_object()}//"
                         f"{ref_chain}/{ref_pos}", 10)
                # create distance object with get_distance command
                cmd.distance(measurement_obj, atom1, atom2)
                cmd.label(atom1, "'%s-%s' % (resn, resi)")
                cmd.label(atom2, "'%s-%s' % (resn, resi)")
                cmd.set("label_position", (0, 0, 10))
                # determine the option for ray_shadows
                if ray_shadows == False:
                    opt_ray_shadows: str = "off"
                else:
                    opt_ray_shadows: str = "on"
                # set image parameters
                cmd.bg_color(constants.PYMOL_DEFAULT_BACKGROUND_COLOR)
                cmd.set("ray_trace_mode", 3)
                cmd.set("antialias", 2)
                cmd.set("ray_shadows", opt_ray_shadows)
                cmd.set('ray_opaque_background', opaque_background)

                cmd.scene(key=f"{ref_pos}-{model_pos}", action="store")

                if take_images is True:
                    # check if path exists where the data will be exported,
                    # if not the directory will be created
                    if not os.path.exists(
                           constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR):
                        os.mkdir(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR)
                    # save image as 300 dpi png image
                    cmd.ray(2400, 2400, renderer=0)
                    cmd.png(
                        f'{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR}/{filename}_{ref_pos}.png',
                        dpi=300)
                # hide created labels
                cmd.hide("labels", atom1)
                cmd.hide("labels", atom2)
                cmd.hide("labels", measurement_obj)
                cmd.hide("dashes", measurement_obj)
                i += 1
                j += 1
            else:
                i += 1

    def serialize_distance_analysis(self, xml_distance_analysis_element):
        """This function serialize the protein pair object

        """
        tmp_distance_analysis = ElementTree.SubElement(xml_distance_analysis_element, element_names.DISTANCE_ANALYSIS)
        tmp_distance_analysis.set(attribute_names.DISTANCE_ANALYSIS_NAME, str(self.name))
        tmp_distance_analysis.set(attribute_names.DISTANCE_ANALYSIS_CUTOFF, str(self.cutoff))
        tmp_distance_analysis.set(attribute_names.DISTANCE_ANALYSIS_CYCLES, str(self.cycles))

        logger.debug(f"Serialization of: {self.name}")
        logger.debug(self.analysis_results)
        # if len(self.distance_analysis_data) == 0:
        #     self.analysis_results.serialize_distance_analysis_results(tmp_distance_analysis)
        # else:
        # tmp_results_data = ElementTree.SubElement(tmp_distance_analysis, element_names.DISTANCE_ANALYSIS_RESULTS)

        self.analysis_results.serialize_distance_analysis_results(tmp_distance_analysis)
        #
        # # <editor-fold desc="RMSD and aligned AA">
        # tmp_results_data.set(attribute_names.DISTANCE_ANALYSIS_RMSD, str(self.analysis_results.rmsd))
        # tmp_results_data.set(attribute_names.DISTANCE_ANALYSIS_ALIGNED_AA, str(self.analysis_results.aligned_aa))
        #
        # # </editor-fold>
        #
        # # <editor-fold desc="Distance hashtable">
        # tmp_distance_data = ElementTree.SubElement(tmp_results_data, element_names.DISTANCE_ANALYSIS_DISTANCE_RESULTS)
        # tmp_index_data = ElementTree.SubElement(tmp_distance_data, element_names.DISTANCE_ANALYSIS_INDEX_LIST)
        # tmp_index_data.text = str(self.distance_analysis_data['index'].tolist())
        # tmp_prot_1_chain_data = ElementTree.SubElement(tmp_distance_data, element_names.DISTANCE_ANALYSIS_PROT_1_CHAIN_LIST)
        # tmp_prot_1_chain_data.text = str(self.distance_analysis_data['ref_chain'].tolist())
        # tmp_prot_1_position_data = ElementTree.SubElement(tmp_distance_data, element_names.DISTANCE_ANALYSIS_PROT_1_POSITION_LIST)
        # tmp_prot_1_position_data.text = str(self.distance_analysis_data['ref_pos'].tolist())
        # tmp_prot_1_residue_data = ElementTree.SubElement(tmp_distance_data, element_names.DISTANCE_ANALYSIS_PROT_1_RESIDUE_LIST)
        # tmp_prot_1_residue_data.text = str(self.distance_analysis_data['ref_resi'].tolist())
        # tmp_prot_2_chain_data = ElementTree.SubElement(tmp_distance_data, element_names.DISTANCE_ANALYSIS_PROT_2_CHAIN_LIST)
        # tmp_prot_2_chain_data.text = str(self.distance_analysis_data['model_chain'].tolist())
        # tmp_prot_2_position_data = ElementTree.SubElement(tmp_distance_data, element_names.DISTANCE_ANALYSIS_PROT_2_POSITION_LIST)
        # tmp_prot_2_position_data.text = str(self.distance_analysis_data['model_pos'].tolist())
        # tmp_prot_2_residue_data = ElementTree.SubElement(tmp_distance_data, element_names.DISTANCE_ANALYSIS_PROT_2_RESIDUE_LIST)
        # tmp_prot_2_residue_data.text = str(self.distance_analysis_data['model_resi'].tolist())
        # tmp_distances_data = ElementTree.SubElement(tmp_distance_data, element_names.DISTANCE_ANALYSIS_DISTANCES_LIST)
        # tmp_distances_data.text = str(self.distance_analysis_data['distance'].tolist())
        #
        # # </editor-fold>
        #
        # # <editor-fold desc="Images">
        # tmp_image_data = ElementTree.SubElement(tmp_results_data, element_names.DISTANCE_ANALYSIS_IMAGES)
        # tmp_image_structure_aln = ElementTree.SubElement(tmp_image_data, element_names.DISTANCE_ANALYSIS_STRUCTURE_ALN_IMAGE)
        # #tmp_image_filename = os.listdir(constants.CACHE_STRUCTURE_ALN_IMAGES_DIR)
        # #tmp_image_structure_aln.text = binary_data.create_base64_string_from_file(path_util.FilePath(pathlib.Path(f"{constants.CACHE_STRUCTURE_ALN_IMAGES_DIR}/{tmp_image_filename[1]}")))
        # tmp_image_structure_aln.text = self.analysis_results.structure_aln_image
        # tmp_image_interesting_reg = ElementTree.SubElement(tmp_image_data, element_names.DISTANCE_ANALYSIS_ALN_IMAGES_INTERESTING_REGIONS)
        # for tmp_base64_image in self.analysis_results.interesting_regions_images:
        #     tmp_image_interesting_reg.text = tmp_base64_image
        # #for tmp_filename in os.listdir(constants.CACHE_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR):
        # #    tmp_image_interesting_reg.text = binary_data.create_base64_string_from_file(path_util.FilePath(pathlib.Path(f"{constants.CACHE_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR}/{tmp_filename}")))
        #
        # # </editor-fold>

        tmp_session_data = ElementTree.SubElement(tmp_distance_analysis, element_names.DISTANCE_ANALYSIS_SESSION)
        tmp_session_data.set(attribute_names.PROTEIN_PAIR_SESSION, pymol_io.convert_pymol_session_to_base64_string(self.name))
