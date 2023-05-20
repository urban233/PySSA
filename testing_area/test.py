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
from xml.etree import ElementTree


class XMLDeserialzer:

    def __init__(self, filepath):
        # Read the XML file
        xml_file = open(filepath, "r")
        xml_contents = xml_file.read()
        self.xml_root = ElementTree.fromstring(xml_contents)

    def deserialize_analysis_images(self, analysis_results: 'results.DistanceAnalysisResults'):
        for tmp_protein_pair in self.xml_root.findall(".//protein_pair"):
            if tmp_protein_pair.attrib["name"] == "6OMN_1_with_6OMN_2":
                structure_aln_images = tmp_protein_pair.findall(".//auto_images/structure_aln_image")
                print(structure_aln_images[0].tag)
                interesting_reg_images = tmp_protein_pair.findall(".//auto_images/interesting_reg_image")
                for tmp_image in interesting_reg_images:
                    print(tmp_image.tag)


if __name__ == '__main__':

