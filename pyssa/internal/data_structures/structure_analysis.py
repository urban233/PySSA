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

import pymol
from pymol import cmd
from typing import TYPE_CHECKING
from pyssa.io_pyssa import path_util, filesystem_helpers
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants
from pyssa.util import exception
from PyQt5.QtWidgets import QCheckBox

if TYPE_CHECKING:
    from pyssa.internal.data_structures import project
    from pyssa.internal.data_structures import protein_pair

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class Analysis:
    """This class contains information about the type of analysis."""

    analysis_list: list['protein_pair.ProteinPair'] = []
    app_project: 'project.Project'
    directory_path: 'pathlib.Path'

    def __init__(self, app_project: 'project.Project') -> None:
        """Initialize the app project.

        Args:
            app_project(project.Project): The project that this analysis is used.
        """
        self.app_project = app_project

    def run_analysis(self, cb_analysis_images: QCheckBox) -> None:
        """This function is used to run the analysis.

        Args:
            cb_analysis_images: Is a boolean, indicating to take images or not.

        Raises:
            ValueError: If analysis list is empty.
            UnableToTakeImageError: If no image was taken.
            UnableToSetImageError: If no image was setting.
            UnableToSafeSessionError: If session could not be saved.
            UnableToOpenFileError: If file could not be opened.
            FileNotFoundError: If file could not be founded.
            DirectoryNotFoundError: If directory could not be founded.
            IllegalArgumentError: If an invalid argument was used.
        """
        logger.debug(f"self.analysis_list: {self.analysis_list}")
        # <editor-fold desc="Checks">
        if len(self.analysis_list) == 0:
            logger.error("Analysis list is empty.")
            raise ValueError("Analysis list is empty.")

        # </editor-fold>

        if cb_analysis_images.isChecked():
            # make images if checked
            take_images: bool = True
        else:
            take_images: bool = False

        # create scene/ image of structure alignment
        filesystem_helpers.create_directory(constants.SCRATCH_DIR_IMAGES)
        filesystem_helpers.create_directory(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR)
        filesystem_helpers.create_directory(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR)

        for tmp_protein_pair in self.analysis_list:
            logger.info(f"The protein pair: {tmp_protein_pair.name} gets analyzed.")
            try:
                cmd.reinitialize()
            except pymol.CmdException:
                logger.error("Unable to reinitialize the pymol session.")
                raise RuntimeError()  # Todo: needs better exception
            try:
                # run analysis in pymol
                tmp_protein_pair.distance_analysis.do_analysis_in_pymol()
                tmp_protein_pair.distance_analysis.take_image_of_protein_pair(
                    filename=f"structure_aln_{tmp_protein_pair.name}",
                    representation="cartoon", take_images=take_images)
                if take_images is True:
                    tmp_protein_pair.distance_analysis.analysis_results.set_structure_aln_image(
                        path_util.FilePath(
                            pathlib.Path(
                                f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR}/"
                                f"structure_aln_{tmp_protein_pair.name}.png"),
                        ),
                    )
                    logger.debug(tmp_protein_pair.distance_analysis.analysis_results.structure_aln_image[0])
                tmp_protein_pair.distance_analysis.take_image_of_interesting_regions(
                    tmp_protein_pair.distance_analysis.cutoff,
                    f"interesting_reg_{tmp_protein_pair.name}",
                    take_images=take_images)
                if take_images is True:
                    interesting_region_filepaths = []
                    for tmp_filename in os.listdir(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR):
                        interesting_region_filepaths.append(
                            path_util.FilePath(
                                pathlib.Path(f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR}/"
                                             f"{tmp_filename}"),
                            ),
                        )
                    (tmp_protein_pair.distance_analysis.analysis_results.
                     set_interesting_region_images(interesting_region_filepaths))
                shutil.rmtree(constants.SCRATCH_DIR_IMAGES)
                cmd.scene(f"{tmp_protein_pair.protein_1.get_molecule_object()}-{tmp_protein_pair.protein_2.get_molecule_object()}",
                          action="recall")
                tmp_protein_pair.save_session_of_protein_pair()
                self.app_project.add_protein_pair(tmp_protein_pair)
            except pymol.CmdException:
                logger.error("Unable to recall scene in pymol session.")
                raise RuntimeError()  # Todo: needs better exception
            except exception.UnableToDoAnalysisError:
                logger.error("The analysis in PyMOL failed!")
                raise exception.UnableToDoAnalysisError("")
            except exception.UnableToColorProteinPairError:
                logger.error("Could not color the protein pair!")
                raise exception.UnableToTakeImageError("")
            except exception.UnableToTakeImageError:
                logger.error("Could not take image of the structure alignment!")
                raise exception.UnableToTakeImageError("")
            except FileNotFoundError:
                logger.error("Image file could not be found!")
                raise exception.UnableToOpenFileError(f"Image file: {constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR}/"
                                                      f"structure_aln_{tmp_protein_pair.name}.png")
            except exception.IllegalArgumentError:
                logger.error("The argument filename is illegal.")
                raise exception.UnableToOpenFileError(f"filename: {tmp_protein_pair.name}")
            except exception.UnableToSetImageError:
                logger.error("Could not set images of interesting regions!")
                raise exception.UnableToSetImageError("")
            except exception.DirectoryNotFoundError:
                logger.error(f"Could not find directory: "
                             f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR}.")
                raise exception.DirectoryNotFoundError("")
            except exception.UnableToSafeSessionError:
                logger.error("Could not save session of the protein pair!")
                raise exception.UnableToSafeSessionError("")
            except exception.UnableToAddProteinPairError:
                logger.error("Could not add protein pair!")
                raise exception.UnableToAddProteinPairError("")
            except Exception as e:
                logger.error(f"Unknown error: {e}")
                raise exception.UnableToSetImageError("")
        self.analysis_list.clear()
