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
import logging
import os.path
import typing
import pathlib
from xml.etree import ElementTree

from pyssa.logging_pyssa import log_handlers
from pyssa.util import exception

if typing.TYPE_CHECKING:
    from pyssa.io_pyssa import path_util

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def convert_pdb_file_into_xml_element(filepath: "path_util.FilePath") -> ElementTree.Element:
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


def convert_xml_string_to_pdb_file(xml_contents: "ElementTree.ElementTree", path_pdb: "pathlib.Path") -> None:
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


def convert_pdb_data_list_to_pdb_file(a_pdb_filepath: pathlib.Path, a_pdb_data: list) -> None:
    """Convert pdb data into a pdb file.

    Args:
        a_pdb_filepath: A filepath of a pdb file.
        a_pdb_data: A list of pdb file content.

    Raises:
        IllegalArgumentError: If the argument is not usable.
        DirectoryNotFoundError: If the parent directory of the pdb file path is not found.
        PermissionError: If the parent directory of the pdb file path has no writing permission.
        UnableToOpenFileError: If the pdb file could not be open.
    """
    # <editor-fold desc="Checks">
    if a_pdb_filepath is None:
        logger.error(f"The argument a_pdb_filepath is illegal: {a_pdb_filepath}!")
        raise exception.IllegalArgumentError("An argument is illegal.")
    if not os.path.exists(a_pdb_filepath.parent):
        logger.error(f"The argument a_pdb_filepath is illegal: {a_pdb_filepath}!")
        raise exception.DirectoryNotFoundError("")
    if not os.access(a_pdb_filepath.parent, os.W_OK):
        logger.error(f"The argument a_pdb_filepath is illegal: {a_pdb_filepath}!")
        raise PermissionError()
    if a_pdb_data is None or len(a_pdb_data) == 0:
        logger.error(f"The argument a_pdb_data is illegal: {a_pdb_data}!")
        raise exception.IllegalArgumentError("An argument is illegal.")

    # </editor-fold>

    try:
        pdb_file = open(a_pdb_filepath, "w")
    except OSError:
        logger.error("pdb file could not be opened for writing.")
        raise exception.UnableToOpenFileError("")

    for tmp_line in a_pdb_data:
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
