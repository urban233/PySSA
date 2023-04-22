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

import os
import json
import pathlib


if __name__ == '__main__':
    import xml.etree.ElementTree as ET


    # Load XML file
    tree = ET.parse('tmp.xml')
    root = tree.getroot()

    for protein_pair in root.findall(".//protein_pair"):
        session_data = protein_pair.find("session_data")
        print(session_data.attrib)
        distance_analysis = protein_pair.find("distance_analysis")
        print(distance_analysis.attrib)
        results = distance_analysis.find("results")
        print(results.attrib)
        results_distance = results.findall("results_distance")
        for tmp_results_distance in results_distance:
            index = tmp_results_distance.find("index").text
            print(tmp_results_distance.find("protein_1_residues").text)
            print(index)


    # # iterate over all distance_analysis elements
    # for distance_analysis in root.findall(".//distance_analysis"):
    #     # get the attributes of the distance_analysis element
    #     name = distance_analysis.get("name")
    #     cutoff = distance_analysis.get("cutoff")
    #     cycles = distance_analysis.get("cycles")
    #
    #     # get the results of the distance_analysis element
    #     results = distance_analysis.find("results")
    #
    #     # get the attributes of the results element
    #     rmsd = results.get("rmsd")
    #     aligned_amino_acids = results.get("aligned_amino_acids")
    #
    #     # get the results_distance element
    #     results_distance = results.find("results_distance")
    #
    #     # get the elements of the results_distance element
    #     index = results_distance.find("index").text
    #
    #     # process the data as needed
    #     print("name:", name)
    #     print("cutoff:", cutoff)
    #     print("cycles:", cycles)
    #     print("rmsd:", rmsd)
    #     print("aligned_amino_acids:", aligned_amino_acids)
    #     print("index:", index)




