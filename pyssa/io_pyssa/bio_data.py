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
"""Module for handling biological data."""
import logging
import os.path
import typing
import pathlib
import requests
from xml.etree import ElementTree

from pyssa.logging_pyssa import log_handlers
from pyssa.util import exception

if typing.TYPE_CHECKING:
    from pyssa.io_pyssa import path_util

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def download_pdb_file(pdb_id, save_path):
    pdb_url = f'https://files.rcsb.org/download/{pdb_id.upper()}.pdb'

    try:
        response = requests.get(pdb_url)
        response.raise_for_status()  # Check for errors
        with open(save_path, 'w') as pdb_file:
            pdb_file.write(response.text)
        print(f"PDB file {pdb_id} downloaded successfully to {save_path}")
    except requests.exceptions.HTTPError as errh:
        print(f"HTTP Error: {errh}")
    except requests.exceptions.ConnectionError as errc:
        print(f"Error Connecting: {errc}")
    except requests.exceptions.Timeout as errt:
        print(f"Timeout Error: {errt}")
    except requests.exceptions.RequestException as err:
        print(f"Error: {err}")


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
    """This function writes a pdb file based on a xml string which contains the pdb information.

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


def convert_pdb_xml_string_to_list(root) -> list:  # noqa: ANN001
    """Converts the xml string of pdb information to a list.

    Args:
        root: the xml root element.
    """
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


def convert_pdb_data_list_to_xml_string(pdb_data: list) -> ElementTree.Element:
    """Converts the list of pdb information into a xml string.

    Args:
        pdb_data: a list of pdb information.
    """
    # Convert the PDB to XML format
    root = ElementTree.Element("pdb_data")
    for line in pdb_data:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom = ElementTree.SubElement(root, "atom")
            atom.text = line
    return root


def build_pdb_file(records, a_filepath):
    with open(a_filepath, 'w') as pdb_file:
        for record in records:
            record_type = record['record_type']
            if record_type == 'ATOM':
                pdb_line = f"{record_type:6s}{record['atom_number']:5d}  {record['atom_name']:4s}" \
                           f"{record['residue_name']:3s} {record['chain_identifier']:1s}" \
                           f"{record['residue_sequence_number']:4d}{record['code_for_insertions_of_residues']:1s}" \
                           f"   {record['x_coord']:8.3f}{record['y_coord']:8.3f}{record['z_coord']:8.3f}" \
                           f"{record['occupancy']:6.2f}{record['temperature_factor']:6.2f}           " \
                           f"{record['element_symbol']:2s}{record['charge']:2s}\n"
            elif record_type == 'TER':
                try:
                    pdb_line = f"{record_type:3s}   {record['atom_number']:5d}      {record['residue_name']:3s}" \
                               f" {record['chain_identifier']:1s}{record['residue_sequence_number']:4d}\n"
                except KeyError:
                    pdb_line = ""
            elif record_type == 'HETATM':
                pdb_line = f"{record_type:6s}{record['atom_number']:5d}  {record['atom_name']:4s}" \
                           f"{record['residue_name']:3s} {record['chain_identifier']:1s}" \
                           f"{record['residue_sequence_number']:4d}    {record['x_coord']:8.3f}{record['y_coord']:8.3f}" \
                           f"{record['z_coord']:8.3f}{record['occupancy']:6.2f}{record['temperature_factor']:6.2f}           " \
                           f"{record['element_symbol']:2s}{record['charge']:2s}\n"
            else:
                raise ValueError(f"Unsupported record type: {record_type}")
            pdb_file.write(pdb_line)
    pdb_file.close()


def parse_pdb_file(a_filepath):
    records = []

    with open(a_filepath, 'r') as pdb_file:
        for line in pdb_file:
            record_type = line[0:6].strip()

            if record_type == "ATOM":
                atom_number = int(line[6:11])
                atom_name = line[12:16].strip()
                alternate_location_indicator = line[16]
                residue_name = line[17:20].strip()
                chain_identifier = line[21]
                residue_sequence_number = int(line[22:26])
                code_for_insertions_of_residues = line[26]
                x_coord = float(line[30:38])
                y_coord = float(line[38:46])
                z_coord = float(line[46:54])
                occupancy = float(line[54:60])
                temperature_factor = float(line[60:66])
                segment_identifier = line[72:76].strip()
                element_symbol = line[76:78].strip()
                charge = line[78:80].strip()

                records.append({
                    'record_type': record_type,
                    'atom_number': atom_number,
                    'atom_name': atom_name,
                    'alternate_location_indicator': alternate_location_indicator,
                    'residue_name': residue_name,
                    'chain_identifier': chain_identifier,
                    'residue_sequence_number': residue_sequence_number,
                    'code_for_insertions_of_residues': code_for_insertions_of_residues,
                    'x_coord': x_coord,
                    'y_coord': y_coord,
                    'z_coord': z_coord,
                    'occupancy': occupancy,
                    'temperature_factor': temperature_factor,
                    'segment_identifier': segment_identifier,
                    'element_symbol': element_symbol,
                    'charge': charge
                })
            elif record_type == "TER":
                try:
                    atom_number_ter = int(line[6:11])
                    residue_name_ter = line[17:20].strip()
                    chain_identifier_ter = line[21]
                    residue_sequence_number_ter = int(line[22:26])
                    records.append({
                        'record_type': record_type,
                        'atom_number': atom_number_ter,
                        'atom_name': "",
                        'alternate_location_indicator': "",
                        'residue_name': residue_name_ter,
                        'chain_identifier': chain_identifier_ter,
                        'residue_sequence_number': residue_sequence_number_ter,
                        'code_for_insertions_of_residues': "",
                        'x_coord': 0.0,
                        'y_coord': 0.0,
                        'z_coord': 0.0,
                        'occupancy': 0.0,
                        'temperature_factor': 0.0,
                        'segment_identifier': "",
                        'element_symbol': "",
                        'charge': ""
                    })
                except ValueError:
                    print("PDB file has only the record type for TER.")
                    records.append({
                        'record_type': record_type,
                    })
            elif record_type == "HETATM":
                atom_number_het = int(line[6:11])
                atom_name_het = line[12:16].strip()
                residue_name_het = line[17:20].strip()
                chain_identifier_het = line[21]
                residue_sequence_number_het = int(line[22:26])
                x_coord_het = float(line[30:38])
                y_coord_het = float(line[38:46])
                z_coord_het = float(line[46:54])
                occupancy_het = float(line[54:60])
                temperature_factor_het = float(line[60:66])
                segment_identifier_het = line[72:76].strip()
                element_symbol_het = line[76:78].strip()
                charge_het = line[78:80].strip()

                records.append({
                    'record_type': record_type,
                    'atom_number': atom_number_het,
                    'atom_name': atom_name_het,
                    'alternate_location_indicator': "",
                    'residue_name': residue_name_het,
                    'chain_identifier': chain_identifier_het,
                    'residue_sequence_number': residue_sequence_number_het,
                    'code_for_insertions_of_residues': "",
                    'x_coord': x_coord_het,
                    'y_coord': y_coord_het,
                    'z_coord': z_coord_het,
                    'occupancy': occupancy_het,
                    'temperature_factor': temperature_factor_het,
                    'segment_identifier': segment_identifier_het,
                    'element_symbol': element_symbol_het,
                    'charge': charge_het
                })
    return records
