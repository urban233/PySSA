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
import pathlib
import numpy as np
from typing import TYPE_CHECKING

from auxiliary_pymol import auxiliary_pymol_client
from pyssa.controller import database_manager
from pyssa.io_pyssa import path_util, filesystem_helpers, bio_data
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants, distance_analysis_util, pyssa_keys, enums, constant_messages
from pyssa.util import exception
from pyssa.internal.data_structures import results, job

if TYPE_CHECKING:
    from pyssa.internal.data_structures import project
    from pyssa.internal.data_structures import protein_pair
    from pyssa.internal.data_structures import settings


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
            for tmp_protein_pair in self.analysis_list:
                tmp_pdb_atom_db_data_1 = db_manager.get_pdb_atoms_of_protein(tmp_protein_pair.protein_1.get_id())
                tmp_pdb_atom_db_data_2 = db_manager.get_pdb_atoms_of_protein(tmp_protein_pair.protein_2.get_id())
                pdb_atom_dict_1 = [{key.value: value for key, value in zip(enums.PdbAtomEnum, t)} for t in
                                   tmp_pdb_atom_db_data_1]
                pdb_atom_dict_2 = [{key.value: value for key, value in zip(enums.PdbAtomEnum, t)} for t in
                                   tmp_pdb_atom_db_data_2]
                tmp_protein_pair.protein_1.set_pdb_data(pdb_atom_dict_1)
                tmp_protein_pair.protein_2.set_pdb_data(pdb_atom_dict_2)

    def _get_last_id_of_protein_pairs(self) -> int:
        if len(self.app_project.protein_pairs) == 0:
            return 1
        else:
            return max(self.app_project.protein_pairs, key=lambda obj: obj.get_id()).get_id()

    def run_distance_analysis(self, the_image_creation_option: bool, the_main_socket, a_socket) -> None:
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
        # filesystem_helpers.create_directory(constants.SCRATCH_DIR_IMAGES)
        # filesystem_helpers.create_directory(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR)
        # filesystem_helpers.create_directory(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR)

        for tmp_protein_pair in self.analysis_list:
            logger.info(f"The protein pair: {tmp_protein_pair.name} gets analyzed.")
            # try:
            #     cmd.reinitialize()
            # except pymol.CmdException:
            #     tmp_msg: str = "Unable to reinitialize the pymol session."
            #     logger.error(tmp_msg)
            #     raise exception.UnableToReinitializePymolSessionError(tmp_msg)
            try:
                # do distance analysis in PyMOL
                if tmp_protein_pair is None:
                    logger.error(constant_messages.ARGUMENT_IS_ILLEGAL)
                    raise exception.IllegalArgumentError(constant_messages.ARGUMENT_IS_ILLEGAL)

                tmp_protein_1_name = tmp_protein_pair.protein_1.get_molecule_object()
                tmp_protein_1_pdb_cache_filepath = pathlib.Path(
                    f"{constants.CACHE_PROTEIN_DIR}/{tmp_protein_1_name}.pdb"
                )
                tmp_protein_2_name = tmp_protein_pair.protein_2.get_molecule_object()
                tmp_protein_2_pdb_cache_filepath = pathlib.Path(
                    f"{constants.CACHE_PROTEIN_DIR}/{tmp_protein_2_name}.pdb"
                )
                try:
                    bio_data.build_pdb_file(
                        tmp_protein_pair.protein_1.get_pdb_data(),
                        tmp_protein_1_pdb_cache_filepath
                    )
                    bio_data.build_pdb_file(
                        tmp_protein_pair.protein_2.get_pdb_data(),
                        tmp_protein_2_pdb_cache_filepath
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
                    )
                )

                # the_main_socket.send_string("Distance Analysis")
                # response = the_main_socket.recv_string()
                # print(f"Received response: {response}")
                # message = {
                #     "job_type": "Distance Analysis",
                #     "the_protein_pair_name": tmp_protein_pair.name,
                #     "a_protein_1_pdb_cache_filepath": str(tmp_protein_1_pdb_cache_filepath),
                #     "a_protein_2_pdb_cache_filepath":str(tmp_protein_2_pdb_cache_filepath),
                #     "a_protein_1_pymol_selection_string": tmp_protein_pair.protein_1.pymol_selection.selection_string,
                #     "a_protein_2_pymol_selection_string": tmp_protein_pair.protein_2.pymol_selection.selection_string,
                #     "a_cutoff": tmp_protein_pair.distance_analysis.cutoff,
                #     "the_cycles": tmp_protein_pair.distance_analysis.cycles,
                # }
                # the_main_socket.send_json(message)
                # response = the_main_socket.recv_string()
                # print(f"Received response: {response}")
                # # Wait for the response from the server
                # a_socket.send_json({"job_type": "Distance Analysis"})

                result = tmp_reply["result"]
                if tmp_reply["data"] is not None:
                    distance_analysis_results_object_values, base64_string = tmp_reply["data"]
                    distances, base64_string, rmsd, aligned_residues = distance_analysis_results_object_values
                    result_hashtable: dict[str, np.ndarry] = {
                        pyssa_keys.ARRAY_DISTANCE_INDEX: np.array(distances[pyssa_keys.ARRAY_DISTANCE_INDEX]),
                        pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN: np.array(distances[pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN]),
                        pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION: np.array(distances[pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION]),
                        pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI: np.array(distances[pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI]),
                        pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN: np.array(distances[pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN]),
                        pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION: np.array(distances[pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION]),
                        pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI: np.array(distances[pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI]),
                        pyssa_keys.ARRAY_DISTANCE_DISTANCES: np.array(distances[pyssa_keys.ARRAY_DISTANCE_DISTANCES]),
                    }
                    print(result_hashtable)

                    tmp_protein_pair.distance_analysis.analysis_results = results.DistanceAnalysisResults(
                        result_hashtable, base64_string, rmsd, aligned_residues
                    )
                    tmp_protein_pair.pymol_session = base64_string
                else:
                    logger.warning("Returning data was None!")

                # tmp_protein_pair.distance_analysis.analysis_results, tmp_protein_pair.pymol_session = auxiliary_pymol.AuxiliaryPyMOL.do_distance_analysis(
                #     tmp_protein_pair
                # )
                # tmp_protein_pair.distance_analysis.analysis_results = (
                #     distance_analysis_util.do_distance_analysis_in_pymol(
                #         tmp_protein_pair,
                #     )
                # )

                # create scene for structure alignment // take images of structure alignment if necessary
                # distance_analysis_util.create_scene_of_protein_pair(
                #     a_protein_pair=tmp_protein_pair,
                #     filename=f"structure_aln_{tmp_protein_pair.name}",
                #     take_images=the_image_creation_option,
                # )
                # create scenes for interesting regions // take image of interesting regions if necessary
                # distance_analysis_util.create_scenes_of_interesting_regions(
                #     tmp_protein_pair.distance_analysis.analysis_results.distance_data,
                #     tmp_protein_pair.protein_1.get_molecule_object(),
                #     tmp_protein_pair.protein_2.get_molecule_object(),
                #     tmp_protein_pair.distance_analysis.cutoff,
                #     take_images=the_image_creation_option,
                #     filename=f"interesting_reg_{tmp_protein_pair.name}",
                # )
                # if the_image_creation_option is True:
                #     logger.info("Setting the structure alignment image into the results object ...")
                #     tmp_protein_pair.distance_analysis.analysis_results.set_structure_aln_image(
                #         path_util.FilePath(
                #             pathlib.Path(
                #                 f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR}/"
                #                 f"structure_aln_{tmp_protein_pair.name}.png",
                #             ),
                #         ),
                #     )
                #     logger.info("Setting all image of the interesting regions into the results object ...")
                #     interesting_region_filepaths = []
                #     for tmp_filename in os.listdir(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR):
                #         interesting_region_filepaths.append(
                #             path_util.FilePath(
                #                 pathlib.Path(
                #                     f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR}/"
                #                     f"{tmp_filename}",
                #                 ),
                #             ),
                #         )
                #     (
                #         tmp_protein_pair.distance_analysis.analysis_results.set_interesting_region_images(
                #             interesting_region_filepaths,
                #         )
                #     )
            # except FileNotFoundError:
            #     tmp_path: str = str(
            #         pathlib.Path(
            #             f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR}/structure_aln_{tmp_protein_pair.name}.png",
            #         ),
            #     )
            #     logger.error(f"Image file could not be found! {tmp_path}")
            #     raise exception.UnableToOpenFileError(f"Image file: {tmp_path}")
            except exception.IllegalArgumentError:
                logger.error("The argument filename is illegal.")
                raise exception.UnableToOpenFileError(f"filename: {tmp_protein_pair.name}")
            except Exception as e:
                logger.error(f"Unknown error: {e}")
                raise exception.UnableToSetImageError("")
            filesystem_helpers.delete_directory(constants.SCRATCH_DIR_ANALYSIS)
            # cmd.scene(
            #     f"{tmp_protein_pair.protein_1.get_molecule_object()}-"
            #     f"{tmp_protein_pair.protein_2.get_molecule_object()}",
            #     action="recall",
            # )
            # save pymol session of distance analysis
            # tmp_protein_pair.save_session_of_protein_pair()

    def run_analysis(self, the_analysis_type: str, the_image_option: bool, the_main_socket, a_socket) -> None:
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
            self.run_distance_analysis(the_image_option, the_main_socket, a_socket)
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
