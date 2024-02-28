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
"""Module for the results class."""
import csv
import pathlib
from xml.etree import ElementTree
import logging
import numpy as np

from pyssa.io_pyssa import path_util
from pyssa.io_pyssa import binary_data
from pyssa.io_pyssa.xml_pyssa import element_names
from pyssa.io_pyssa.xml_pyssa import attribute_names
from pyssa.util import pyssa_keys
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class DistanceAnalysisResults:
    """Contains the results of the distance analysis done in PyMOL."""

    """
    The data of the distance analysis including: position, chain, residue, residue number and distance.
    """
    distance_data: dict[str, np.ndarray]

    """
    The base64 string of the pymol session used for the distance analysis
    """
    pymol_session: str

    """
    The RMSD value of the protein pair.
    """
    rmsd: float

    """
    The aligned residues of the protein pair.
    """
    aligned_aa: str

    """
    A tuple with the basename and the base64 string of the structure alignment image
    """
    structure_aln_image: tuple[str, str] = ()

    """
    A list of tuples with the basename and the base64 string of the interesting region image
    """
    interesting_regions_images: list[tuple[str, str]] = []

    def __init__(self, distance_data: dict, pymol_session: str, rmsd: float, aligned_aa: str) -> None:
        """Constructor."""
        self.distance_data = distance_data
        self.pymol_session = pymol_session
        self.rmsd = rmsd
        self.aligned_aa = aligned_aa

    def set_structure_aln_image(self, filepath: path_util.FilePath) -> None:
        """This function sets the structure alignment image."""
        self.structure_aln_image = (filepath.get_basename(), binary_data.create_base64_string_from_file(filepath))

    def set_interesting_region_images(self, filepaths: list[path_util.FilePath]) -> None:
        """This function sets the interesting region images.

        Args:
            filepaths (list[path_util.FilePath]): List of file paths.
        """
        for tmp_filepath in filepaths:
            self.interesting_regions_images.append(
                (tmp_filepath.get_basename(), binary_data.create_base64_string_from_file(tmp_filepath)),
            )

    def create_image_png_files_from_base64(self) -> None:
        """Creates png files from the base64 data."""
        binary_data.write_binary_file_from_base64_string(
            pathlib.Path(f"{constants.CACHE_STRUCTURE_ALN_IMAGES_DIR}/{self.structure_aln_image[0]}"),
            self.structure_aln_image[1],
        )
        for tmp_interesting_reg in self.interesting_regions_images:
            binary_data.write_binary_file_from_base64_string(
                pathlib.Path(
                    f"{constants.CACHE_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR}/{tmp_interesting_reg[0]}",
                ),
                tmp_interesting_reg[1],
            )

    def export_distance_data_as_csv(self, a_filepath):
        tmp_data_to_write = np.transpose([
            self.distance_data["index"],
            self.distance_data["ref_chain"],
            self.distance_data["ref_pos"],
            self.distance_data["ref_resi"],
            self.distance_data["model_chain"],
            self.distance_data["model_pos"],
            self.distance_data["model_resi"],
            self.distance_data["distance"],
        ])
        # Writing to the TSV file
        with open(a_filepath, mode='w', newline='', encoding='utf-8') as tsv_file:
            csv_writer = csv.writer(tsv_file, delimiter=',')
            # Write the header
            csv_writer.writerow([
                    "Residue_Pair_No",
                    "Protein_1_Chain",
                    "Protein_1_Position",
                    "Protein_1_Residue",
                    "Protein_2_Chain",
                    "Protein_2_Position",
                    "Protein_2_Residue",
                    "Distance"
            ])
            # Write the data to the TSV file
            csv_writer.writerows(tmp_data_to_write)


    def serialize_distance_analysis_results(self, tmp_distance_analysis) -> None:  # noqa: ANN001
        """Serializes the distance analysis results object."""
        tmp_results_data = ElementTree.SubElement(tmp_distance_analysis, element_names.DISTANCE_ANALYSIS_RESULTS)
        tmp_results_data.set(attribute_names.DISTANCE_ANALYSIS_RMSD, str(self.rmsd))
        tmp_results_data.set(attribute_names.DISTANCE_ANALYSIS_ALIGNED_AA, str(self.aligned_aa))

        # logger.debug(f"These are the keys of the self.distance_data hashtable: {self.distance_data.keys()}")
        # logger.debug(f"These are the values of the self.distance_data hashtable: {self.distance_data.values()}")
        # <editor-fold desc="Distance data">
        tmp_distance_data = ElementTree.SubElement(tmp_results_data, element_names.DISTANCE_ANALYSIS_DISTANCE_RESULTS)
        tmp_index_data = ElementTree.SubElement(tmp_distance_data, element_names.DISTANCE_ANALYSIS_INDEX_LIST)
        tmp_index_data.text = str(self.distance_data[pyssa_keys.ARRAY_DISTANCE_INDEX].tolist())
        tmp_prot_1_chain_data = ElementTree.SubElement(
            tmp_distance_data,
            element_names.DISTANCE_ANALYSIS_PROT_1_CHAIN_LIST,
        )
        tmp_prot_1_chain_data.text = str(self.distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN].tolist())
        tmp_prot_1_position_data = ElementTree.SubElement(
            tmp_distance_data,
            element_names.DISTANCE_ANALYSIS_PROT_1_POSITION_LIST,
        )
        tmp_prot_1_position_data.text = str(self.distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION].tolist())
        tmp_prot_1_residue_data = ElementTree.SubElement(
            tmp_distance_data,
            element_names.DISTANCE_ANALYSIS_PROT_1_RESIDUE_LIST,
        )
        tmp_prot_1_residue_data.text = str(self.distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI].tolist())
        tmp_prot_2_chain_data = ElementTree.SubElement(
            tmp_distance_data,
            element_names.DISTANCE_ANALYSIS_PROT_2_CHAIN_LIST,
        )
        tmp_prot_2_chain_data.text = str(self.distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN].tolist())
        tmp_prot_2_position_data = ElementTree.SubElement(
            tmp_distance_data,
            element_names.DISTANCE_ANALYSIS_PROT_2_POSITION_LIST,
        )
        tmp_prot_2_position_data.text = str(self.distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION].tolist())
        tmp_prot_2_residue_data = ElementTree.SubElement(
            tmp_distance_data,
            element_names.DISTANCE_ANALYSIS_PROT_2_RESIDUE_LIST,
        )
        tmp_prot_2_residue_data.text = str(self.distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI].tolist())
        tmp_distances_data = ElementTree.SubElement(tmp_distance_data, element_names.DISTANCE_ANALYSIS_DISTANCES_LIST)
        tmp_distances_data.text = str(self.distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES].tolist())

        # </editor-fold>

        # <editor-fold desc="Images">
        if len(self.structure_aln_image) != 0:
            tmp_image_data = ElementTree.SubElement(tmp_results_data, element_names.DISTANCE_ANALYSIS_IMAGES)
            tmp_image_structure_aln = ElementTree.SubElement(
                tmp_image_data,
                element_names.DISTANCE_ANALYSIS_STRUCTURE_ALN_IMAGE,
            )
            tmp_image_structure_aln.set(
                attribute_names.DISTANCE_ANALYSIS_STRUCTURE_ALN_IMAGE_BASENAME,
                self.structure_aln_image[0],
            )
            tmp_image_structure_aln.text = self.structure_aln_image[1]

            if len(self.interesting_regions_images) != 0:
                for tmp_base64_image in self.interesting_regions_images:
                    tmp_image_interesting_reg = ElementTree.SubElement(
                        tmp_image_data,
                        element_names.DISTANCE_ANALYSIS_ALN_IMAGES_INTERESTING_REGIONS,
                    )
                    tmp_image_interesting_reg.set(
                        attribute_names.DISTANCE_ANALYSIS_ALN_IMAGES_INTERESTING_REGIONS_BASENAME,
                        tmp_base64_image[0],
                    )
                    tmp_image_interesting_reg.text = tmp_base64_image[1]

        # </editor-fold>
