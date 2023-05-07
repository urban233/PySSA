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
"""This module is used to handle the project xml file"""
from xml.etree import ElementTree
from xml.dom import minidom
from pyssa.io_pyssa.xml_pyssa import element_names
from pyssa.io_pyssa.xml_pyssa import attribute_names


class XmlHandler:

    root: ElementTree.ElementTree
    project_root: ElementTree.Element
    project_info: ElementTree.SubElement
    proteins: ElementTree.SubElement
    protein_pairs: ElementTree.SubElement

    def __init__(self):
        pass

    def write_basic_xml_file(self, project_name, workspace_path, creation_date, os, filepath):
        """This function builds the initial xml file which contains only the project information

        """
        self.project_root = ElementTree.Element(element_names.PROJECT)
        # setup project information tree
        self.project_info = ElementTree.SubElement(self.project_root, element_names.PROJECT_INFO)
        self.project_info.set(attribute_names.PROJECT_NAME, project_name)
        self.project_info.set(attribute_names.PROJECT_WORKSPACE_PATH, workspace_path)
        self.project_info.set(attribute_names.PROJECT_CREATION_DATE, creation_date)
        self.project_info.set(attribute_names.PROJECT_OS, os)
        # setup protein information tree
        self.proteins = ElementTree.SubElement(self.project_root, element_names.PROTEINS)
        self.proteins.set(attribute_names.PROTEIN_MOLECULE_OBJECT, pymol_molecule_object)
        self.proteins.set(attribute_names.PROTEIN_SELECTION, pymol_selection)
        self.proteins.set(attribute_names.PROTEIN_PDB_FILEPATH, pdb_filepath)
        self.proteins.set(attribute_names.PROTEIN_FASTA_FILEPATH, fasta_filepath)
        self.proteins.set(attribute_names.PROTEIN_RESULTS_PATH, export_dirname)
        self.proteins.set(attribute_names.PROTEIN_SESSION_FILEPATH, pymol_session_filepath)
        
        # setup protein_pair information tree
        self.protein_pairs = ElementTree.SubElement(self.project_root, element_names.PROTEIN_PAIRS)
        # write xml file to filesystem
        self.pretty_file(self.project_root, filepath)

    def build_simple_xml_file(self):
        # TODO: What is with the path and the name of the project? Is it a test or what ist it? -> Make it by default!
        self.project_root = ElementTree.Element("project")
        # --- Start of project information
        project_info = ElementTree.SubElement(self.project_root, "project_info")
        project_info.set("name", "Gelis-protein")
        project_info.set("workspace_path", "C:\\Users\\martin\\.pyssa\\default_workspace")
        project_info.set("creation_date", "2023-02-18")
        project_info.set("os", "Windows")
        # --- End of project information

        # --- Start of protein information
        self.proteins = ElementTree.SubElement(self.project_root, "proteins")
        protein = ElementTree.SubElement(self.proteins, "protein_1")
        protein.set("molecule_object", "3bmp")


        return self.project_root

    def prettify(self, elem):
        """Return a pretty-printed XML string for the Element.
        """
        rough_string = ElementTree.tostring(elem, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        return reparsed.toprettyxml(indent="  ")

    def pretty_file(self, root, filename):
        xmlstr = minidom.parseString(ElementTree.tostring(root)).toprettyxml(indent="   ")
        with open(filename, "w") as f:
            f.write(xmlstr)

    def parse_xml_file(self, src):
        self.root = ElementTree.parse(src)

    def get_project_infomation(self):

        print("-- Start of information")
        for info in self.root.iter("project_info"):
            print(info.attrib)
            return info.attrib

        # project_info_elements = self.root.findall("./project_info")
        # for info in self.root.findall("./project_info/project_name"):
        #     print(info.text)

    def add_protein(self):
        protein = ElementTree.SubElement(self.proteins, "protein_2")
        protein.set("molecule_object", "6omn")
        self.proteins.append(protein)
