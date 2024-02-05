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
"""Module for structure analysis class."""
import logging
import os.path
import pathlib
from xml.etree import ElementTree

import pymol
from pymol import cmd
from typing import TYPE_CHECKING
from PyQt5 import QtCore

from pyssa.controller import database_manager
from pyssa.internal.portal import pymol_io
from pyssa.io_pyssa import path_util, filesystem_helpers
from pyssa.io_pyssa.xml_pyssa import element_names, attribute_names
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants, distance_analysis_util, pyssa_keys, enums
from pyssa.util import exception

if TYPE_CHECKING:
    from pyssa.internal.data_structures import project
    from pyssa.internal.data_structures import protein_pair
    from pyssa.internal.data_structures import settings
    from pyssa.internal.data_structures import results

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class Analysis:
    """Contains information about the type of analysis."""

    analysis_list: list["protein_pair.ProteinPair"] = []
    app_project: "project.Project"

    def __init__(self, app_project: "project.Project") -> None:
        """Initialize the app project.

        Args:
            app_project(project.Project): The project that this analysis is used.
        """
        self.app_project: "project.Project" = app_project

    def _fetch_pdb_atoms_for_all_proteins(self):
        with database_manager.DatabaseManager(str(self.app_project.get_database_filepath())) as db_manager:
            db_manager.open_project_database()
            for tmp_protein_pair in self.analysis_list:
                tmp_pdb_atom_db_data_1 = db_manager.get_pdb_atoms_of_protein(tmp_protein_pair.protein_1.get_id())
                tmp_pdb_atom_db_data_2 = db_manager.get_pdb_atoms_of_protein(tmp_protein_pair.protein_2.get_id())
                pdb_atom_dict_1 = [{key.value: value for key, value in zip(enums.PdbAtomEnum, t)} for t in
                                   tmp_pdb_atom_db_data_1]
                pdb_atom_dict_2 = [{key.value: value for key, value in zip(enums.PdbAtomEnum, t)} for t in
                                   tmp_pdb_atom_db_data_2]
                tmp_protein_pair.protein_1.set_pdb_data(pdb_atom_dict_1)
                tmp_protein_pair.protein_2.set_pdb_data(pdb_atom_dict_2)
            db_manager.close_project_database()

    def _get_last_id_of_protein_pairs(self) -> int:
        return max(self.app_project.protein_pairs, key=lambda obj: obj.get_id()).get_id()

    def run_distance_analysis(self, the_image_creation_option: bool) -> None:
        """Runs the distance analysis for all protein pairs of the analysis job.

        Args:
            the_image_creation_option: Is a boolean, indicating to take images or not.

        Raises:
            UnableToReinitializePymolSessionError: If pymol session could not be reinitialized.
            UnableToTakeImageError: If no image was taken.
            UnableToSetImageError: If no image was setting.
            UnableToSafeSessionError: If session could not be saved.
            UnableToOpenFileError: If file could not be opened.
            FileNotFoundError: If file could not be founded.
            DirectoryNotFoundError: If directory could not be founded.
            IllegalArgumentError: If an invalid argument was used.
        """
        logger.debug(
            f"The function argument of the value for the image_creation_option is: {the_image_creation_option}",
        )
        self._fetch_pdb_atoms_for_all_proteins()
        tmp_latest_protein_pair_id = self._get_last_id_of_protein_pairs()
        # create scratch dirs
        filesystem_helpers.create_directory(constants.SCRATCH_DIR_ANALYSIS)
        filesystem_helpers.create_directory(constants.SCRATCH_DIR_IMAGES)
        filesystem_helpers.create_directory(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR)
        filesystem_helpers.create_directory(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR)

        for tmp_protein_pair in self.analysis_list:
            logger.info(f"The protein pair: {tmp_protein_pair.name} gets analyzed.")
            try:
                cmd.reinitialize()
            except pymol.CmdException:
                tmp_msg: str = "Unable to reinitialize the pymol session."
                logger.error(tmp_msg)
                raise exception.UnableToReinitializePymolSessionError(tmp_msg)
            try:
                # do distance analysis in PyMOL
                tmp_protein_pair.distance_analysis.analysis_results = (
                    distance_analysis_util.do_distance_analysis_in_pymol(
                        tmp_protein_pair,
                    )
                )
                # create scene for structure alignment // take images of structure alignment if necessary
                distance_analysis_util.create_scene_of_protein_pair(
                    a_protein_pair=tmp_protein_pair,
                    filename=f"structure_aln_{tmp_protein_pair.name}",
                    take_images=the_image_creation_option,
                )
                # create scenes for interesting regions // take image of interesting regions if necessary
                distance_analysis_util.create_scenes_of_interesting_regions(
                    tmp_protein_pair.distance_analysis.analysis_results.distance_data,
                    tmp_protein_pair.protein_1.get_molecule_object(),
                    tmp_protein_pair.protein_2.get_molecule_object(),
                    tmp_protein_pair.distance_analysis.cutoff,
                    take_images=the_image_creation_option,
                    filename=f"interesting_reg_{tmp_protein_pair.name}",
                )
                if the_image_creation_option is True:
                    logger.info("Setting the structure alignment image into the results object ...")
                    tmp_protein_pair.distance_analysis.analysis_results.set_structure_aln_image(
                        path_util.FilePath(
                            pathlib.Path(
                                f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR}/"
                                f"structure_aln_{tmp_protein_pair.name}.png",
                            ),
                        ),
                    )
                    logger.info("Setting all image of the interesting regions into the results object ...")
                    interesting_region_filepaths = []
                    for tmp_filename in os.listdir(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR):
                        interesting_region_filepaths.append(
                            path_util.FilePath(
                                pathlib.Path(
                                    f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR}/"
                                    f"{tmp_filename}",
                                ),
                            ),
                        )
                    (
                        tmp_protein_pair.distance_analysis.analysis_results.set_interesting_region_images(
                            interesting_region_filepaths,
                        )
                    )
            except FileNotFoundError:
                tmp_path: str = str(
                    pathlib.Path(
                        f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR}/structure_aln_{tmp_protein_pair.name}.png",
                    ),
                )
                logger.error(f"Image file could not be found! {tmp_path}")
                raise exception.UnableToOpenFileError(f"Image file: {tmp_path}")
            except exception.IllegalArgumentError:
                logger.error("The argument filename is illegal.")
                raise exception.UnableToOpenFileError(f"filename: {tmp_protein_pair.name}")
            except Exception as e:
                logger.error(f"Unknown error: {e}")
                raise exception.UnableToSetImageError("")
            filesystem_helpers.delete_directory(constants.SCRATCH_DIR_IMAGES)
            cmd.scene(
                f"{tmp_protein_pair.protein_1.get_molecule_object()}-"
                f"{tmp_protein_pair.protein_2.get_molecule_object()}",
                action="recall",
            )
            # save pymol session of distance analysis
            tmp_protein_pair.save_session_of_protein_pair()

    def run_analysis(self, the_analysis_type: str, the_image_option: bool) -> None:
        """This function is used to run the analysis.

        Args:
            the_analysis_type: Defines which type of analysis should be performed.
            the_image_option: Is a boolean, indicating to take images or not.

        Raises:
            ValueError: If analysis list is empty.
        """
        # <editor-fold desc="Checks">
        if len(self.analysis_list) == 0:
            logger.error("Analysis list is empty.")
            logger.debug(f"self.analysis_list: {self.analysis_list}")
            raise ValueError("Analysis list is empty.")

        # </editor-fold>

        if the_analysis_type == "distance":
            self.run_distance_analysis(the_image_option)
        else:
            tmp_msg: str = f"Unknown analysis type: {the_analysis_type}"
            logger.error(tmp_msg)
            raise ValueError(tmp_msg)


class DistanceAnalysis:
    """Contains all information about the distance analysis of a certain protein pair."""

    """
    The name of the distance analysis.
    """
    name: str

    """
    The cutoff value for the align command from PyMOL
    """
    cutoff: float

    """
    The number of refinement cycles for the align command from PyMOL
    """
    cycles: int

    """
    The size of the figures.
    """
    figure_size: tuple[float, float]

    """
    The object which contains all results from the distance analysis done in PyMOL.
    """
    analysis_results: "results.DistanceAnalysisResults" = None

    def __init__(self, the_app_settings: "settings.Settings", protein_pair_name: str = "",
                 distance_analysis_name: str = "") -> None:
        """Constructor."""
        self.name: str = f"dist_analysis_{protein_pair_name}"
        self.cutoff: float = the_app_settings.cutoff
        self.cycles: int = the_app_settings.cycles
        self.figure_size: tuple[float, float] = (11.0, 6.0)

    def serialize_to_xml_structure(self, an_xml_writer: QtCore.QXmlStreamWriter):
        an_xml_writer.writeStartElement("distance_analysis")
        an_xml_writer.writeAttribute("name", self.name)
        an_xml_writer.writeAttribute("cutoff", str(self.cutoff))
        an_xml_writer.writeAttribute("cycles", str(self.cycles))
        # Result
        an_xml_writer.writeStartElement("result")
        an_xml_writer.writeAttribute("rmsd", str(self.analysis_results.rmsd))
        an_xml_writer.writeAttribute("aligned_amino_acids", str(self.analysis_results.aligned_aa))
        # Create index_data element
        an_xml_writer.writeStartElement('index')
        an_xml_writer.writeCharacters(str(self.analysis_results.distance_data[pyssa_keys.ARRAY_DISTANCE_INDEX]))
        an_xml_writer.writeEndElement()  # End index_data element
        # Create protein_1_chain_data element
        an_xml_writer.writeStartElement('protein_1_chain')
        an_xml_writer.writeCharacters(str(self.analysis_results.distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN]))
        an_xml_writer.writeEndElement()  # End protein_1_chain_data element
        # Create protein_1_position_data element
        an_xml_writer.writeStartElement('protein_1_position')
        an_xml_writer.writeCharacters(str(self.analysis_results.distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION]))
        an_xml_writer.writeEndElement()  # End protein_1_position_data element
        # Create protein_1_residue_data element
        an_xml_writer.writeStartElement('protein_1_residue')
        an_xml_writer.writeCharacters(str(self.analysis_results.distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI]))
        an_xml_writer.writeEndElement()  # End protein_1_residue_data element
        # Create protein_2_chain_data element
        an_xml_writer.writeStartElement('protein_2_chain')
        an_xml_writer.writeCharacters(str(self.analysis_results.distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN]))
        an_xml_writer.writeEndElement()  # End protein_2_chain_data element
        # Create protein_2_position_data element
        an_xml_writer.writeStartElement('protein_2_position')
        an_xml_writer.writeCharacters(str(self.analysis_results.distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION]))
        an_xml_writer.writeEndElement()  # End protein_2_position_data element
        # Create protein_2_residue_data element
        an_xml_writer.writeStartElement('protein_2_residue')
        an_xml_writer.writeCharacters(str(self.analysis_results.distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI]))
        an_xml_writer.writeEndElement()  # End protein_2_residue_data element
        # Create distances_data element
        an_xml_writer.writeStartElement('distances')
        an_xml_writer.writeCharacters(str(self.analysis_results.distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES]))
        an_xml_writer.writeEndElement()  # End distances_data element
        # Session data
        an_xml_writer.writeStartElement("session_data")
        an_xml_writer.writeAttribute("session", self.analysis_results.pymol_session)
        an_xml_writer.writeEndElement()
        an_xml_writer.writeEndElement()  # End result

    def serialize_distance_analysis(self, xml_distance_analysis_element) -> None:  # noqa: ANN001
        """This function serialize the protein pair object."""
        tmp_distance_analysis = ElementTree.SubElement(
            xml_distance_analysis_element,
            element_names.DISTANCE_ANALYSIS,
        )
        tmp_distance_analysis.set(attribute_names.DISTANCE_ANALYSIS_NAME, str(self.name))
        tmp_distance_analysis.set(attribute_names.DISTANCE_ANALYSIS_CUTOFF, str(self.cutoff))
        tmp_distance_analysis.set(attribute_names.DISTANCE_ANALYSIS_CYCLES, str(self.cycles))

        self.analysis_results.serialize_distance_analysis_results(tmp_distance_analysis)
        tmp_session_data = ElementTree.SubElement(
            tmp_distance_analysis,
            element_names.DISTANCE_ANALYSIS_SESSION,
        )
        tmp_session_data.set(
            attribute_names.PROTEIN_PAIR_SESSION,
            pymol_io.convert_pymol_session_to_base64_string(self.name),
        )
