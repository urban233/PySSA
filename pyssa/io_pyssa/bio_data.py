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
import typing
import pathlib
from xml.etree import ElementTree

if typing.TYPE_CHECKING:
    from pyssa.io_pyssa import path_util


def convert_pdb_file_into_xml_element(filepath: 'path_util.FilePath') -> ElementTree.Element:
    """This function creates a xml string of an existing pdb file.

    Args:
        filepath:
            path where the pdb file is stored
    Returns:
        a xml string
    """
    # TODO: write checks for function arguments
    # Read the PDB file
    pdb_file = open(filepath.get_filepath(), "r")
    pdb_contents = pdb_file.read()
    # Convert the PDB to XML format
    root = ElementTree.Element("pdb_data")
    for line in pdb_contents.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom = ElementTree.SubElement(root, "atom")
            atom.text = line
    # Write out the XML file
    return root


def convert_xml_string_to_pdb_file(xml_contents: 'ElementTree.ElementTree', path_pdb: 'pathlib.Path') -> None:
    """This function writes a pdb file based on a xml string which contains the pdb information

    Args:
        xml_contents:
            a xml element which holds the pdb_data
        path_pdb:
            full filepath to the pdb file which should get written, with filename and extension!
    """
    # TODO: write checks for function arguments
    # Write out the PDB file
    pdb_file = open(path_pdb, "w")
    for tmp_line in convert_pdb_xml_string_to_list(xml_contents):
        pdb_file.write(f"{tmp_line} \n")
    pdb_file.close()


def convert_pdb_xml_string_to_list(root):
    # TODO: write checks for function arguments
    pdb_lines = []
    for tmp_atom in root.iter("atom"):
        pdb_lines.append(tmp_atom.text)
    return pdb_lines


def convert_pdb_data_list_to_pdb_file(path_pdb, pdb_data):
    pdb_file = open(path_pdb, "w")
    for tmp_line in pdb_data:
        pdb_file.write(f"{tmp_line} \n")
    pdb_file.close()


def convert_pdb_data_list_to_xml_string(pdb_data):
    # Convert the PDB to XML format
    root = ElementTree.Element("pdb_data")
    for line in pdb_data:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom = ElementTree.SubElement(root, "atom")
            atom.text = line
    return root
