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
import shutil
from pymol import cmd
from typing import TYPE_CHECKING
from pyssa.io_pyssa import path_util
from pyssa.logging_pyssa import log_handlers
from pyssa.internal.data_structures import protein_pair
from pyssa.util import constants
from pyssa.util import exception
from PyQt5.QtWidgets import QCheckBox

if TYPE_CHECKING:
    from pyssa.internal.data_structures import protein_pair
    from pyssa.internal.data_structures import project

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class Analysis:
    """This class contains information about the type of analysis."""

    analysis_list: list[protein_pair.ProteinPair] = []
    app_project: 'project.Project'

    def __init__(self, app_project: project.Project) -> None:
        self.app_project = app_project

    def run_analysis(self, cb_analysis_images: QCheckBox, tmp_protein_pair: protein_pair) -> None:
        """This function is used to run the analysis.

        Args:
            cb_analysis_images: Is a boolean, indicating to take images or not.
            tmp_protein_pair: Is the protein pair, used for the analysis.

        Raises:
            ValueError: If analysis list is empty.
            UnableToTakeImageError: If no image was taken.
            UnableToSetImageError: If no image was setting.
            UnableToSafeSessionError: If session could not be saved.
        """
        logger.debug(f"self.analysis_list: {self.analysis_list}")
        # <editor-fold desc="Checks">
        if len(self.analysis_list) == 0:
            logger.error("Analysis list is empty.")
            raise ValueError("Analysis list is empty.")
        # </editor-fold>
        logger.debug(self.analysis_list)
        try:
            for tmp_protein_pair in self.analysis_list:
                logger.info(f"The protein pair: {tmp_protein_pair.name} gets analyzed.")
                cmd.reinitialize()
                tmp_protein_pair.distance_analysis.do_analysis_in_pymol()
                take_images = False
        except exception.UnableToDoAnalysisError:
            logger.error("The analysis in PyMOL failed!")
            raise exception.UnableToDoAnalysisError("")

        try:
            # make images if checked
            for tmp_protein_pair in self.analysis_list:
                if cb_analysis_images.isChecked():
                    take_images = True
                if not os.path.exists(constants.SCRATCH_DIR_IMAGES):
                    os.mkdir(constants.SCRATCH_DIR_IMAGES)
                if not os.path.exists(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR):
                    os.mkdir(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR)
                if not os.path.exists(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR):
                    os.mkdir(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR)
                tmp_protein_pair.distance_analysis.take_image_of_protein_pair(filename=f"structure_aln_{tmp_protein_pair.name}",
                                                                              representation="cartoon", take_images=take_images)
        except exception.UnableToTakeImageError:
            logger.error("Could not take image of the protein pair!")
            raise exception.UnableToTakeImageError("")

        try:
            if take_images is True:
                tmp_protein_pair.distance_analysis.analysis_results.set_structure_aln_image(
                    path_util.FilePath(
                        pathlib.Path(
                            f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR}/structure_aln_{tmp_protein_pair.name}.png"),
                    ),
                )
                logger.debug(tmp_protein_pair.distance_analysis.analysis_results.structure_aln_image[0])
        except exception.UnableToTakeImageError:
            logger.error("Could not take image of the structure alignment!")
            raise exception.UnableToTakeImageError("")

        try:
            tmp_protein_pair.distance_analysis.take_image_of_interesting_regions(tmp_protein_pair.distance_analysis.cutoff,
                                                                                 f"interesting_reg_{tmp_protein_pair.name}", take_images=take_images)
        except exception.UnableToTakeImageError:
            logger.error("Could not take images of interesting regions!")
            raise exception.UnableToTakeImageError("")

        try:
            if take_images is True:
                interesting_region_filepaths = []
                for tmp_filename in os.listdir(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR):
                    interesting_region_filepaths.append(
                        path_util.FilePath(
                            pathlib.Path(f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR}/{tmp_filename}"),
                        ),
                    )
                tmp_protein_pair.distance_analysis.analysis_results.set_interesting_region_images(interesting_region_filepaths)
            shutil.rmtree(constants.SCRATCH_DIR_IMAGES)
            cmd.scene(f"{tmp_protein_pair.protein_1.get_molecule_object()}-{tmp_protein_pair.protein_2.get_molecule_object()}", 
                      action="recall")
        except exception.UnableToSetImageError:
            logger.error("Could not set images of interesting regions!")
            raise exception.UnableToSetImageError("")

        try:
            tmp_protein_pair.save_session_of_protein_pair()
        except exception.UnableToSafeSessionError:
            logger.error("Could not save session of the protein pair!")
            raise exception.UnableToSafeSessionError

        try:
            self.app_project.add_protein_pair(tmp_protein_pair)
        except exception.UnableToAddProteinPairError:
            logger.error("Could not add protein pair!")

        self.analysis_list.clear()
