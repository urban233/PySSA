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
"""Module for the results class"""
from xml.etree import ElementTree

import numpy as np

from pyssa.io_pyssa.xml_pyssa import element_names
from pyssa.io_pyssa.xml_pyssa import attribute_names


class DistanceAnalysisResults:

    distance_data: dict[str, np.ndarray]
    pymol_session: str
    rmsd: float
    aligned_aa: int

    def __init__(self, distance_data: dict, pymol_session: str, rmsd: float, aligned_aa: int):
        self.distance_data = distance_data
        self.pymol_session = pymol_session
        self.rmsd = rmsd
        self.aligned_aa = aligned_aa

    def serialize_distance_analysis_results(self, tmp_distance_analysis):
        tmp_results_data = ElementTree.SubElement(tmp_distance_analysis, element_names.DISTANCE_ANALYSIS_RESULTS)
        tmp_results_data.set(attribute_names.DISTANCE_ANALYSIS_RMSD, str(self.rmsd))
        tmp_results_data.set(attribute_names.DISTANCE_ANALYSIS_ALIGNED_AA, str(self.aligned_aa))

        for distances in self.distance_data.values():
            tmp_distance_data = ElementTree.SubElement(tmp_results_data,
                                                       element_names.DISTANCE_ANALYSIS_DISTANCE_RESULTS)
            tmp_index_data = ElementTree.SubElement(tmp_distance_data, element_names.DISTANCE_ANALYSIS_INDEX_LIST)
            tmp_index_data.text = str(distances['index'].tolist())
            tmp_prot_1_chain_data = ElementTree.SubElement(tmp_distance_data,
                                                           element_names.DISTANCE_ANALYSIS_PROT_1_CHAIN_LIST)
            tmp_prot_1_chain_data.text = str(distances['ref_chain'].tolist())
            tmp_prot_1_position_data = ElementTree.SubElement(tmp_distance_data,
                                                              element_names.DISTANCE_ANALYSIS_PROT_1_POSITION_LIST)
            tmp_prot_1_position_data.text = str(distances['ref_pos'].tolist())
            tmp_prot_1_residue_data = ElementTree.SubElement(tmp_distance_data,
                                                             element_names.DISTANCE_ANALYSIS_PROT_1_RESIDUE_LIST)
            tmp_prot_1_residue_data.text = str(distances['ref_resi'].tolist())
            tmp_prot_2_chain_data = ElementTree.SubElement(tmp_distance_data,
                                                           element_names.DISTANCE_ANALYSIS_PROT_2_CHAIN_LIST)
            tmp_prot_2_chain_data.text = str(distances['model_chain'].tolist())
            tmp_prot_2_position_data = ElementTree.SubElement(tmp_distance_data,
                                                              element_names.DISTANCE_ANALYSIS_PROT_2_POSITION_LIST)
            tmp_prot_2_position_data.text = str(distances['model_pos'].tolist())
            tmp_prot_2_residue_data = ElementTree.SubElement(tmp_distance_data,
                                                             element_names.DISTANCE_ANALYSIS_PROT_2_RESIDUE_LIST)
            tmp_prot_2_residue_data.text = str(distances['model_resi'].tolist())
            tmp_distances_data = ElementTree.SubElement(tmp_distance_data,
                                                        element_names.DISTANCE_ANALYSIS_DISTANCES_LIST)
            tmp_distances_data.text = str(distances['distance'].tolist())
