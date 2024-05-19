#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/zielesny/PySSA>
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
import pathlib
import requests

from pyssa.logging_pyssa import log_handlers
from pyssa.util import exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


def download_pdb_file(pdb_id: str, save_path: str) -> None:
    """Downloads a pdb file with a given pdb_id and saves it to a given path.

    Args:
        pdb_id (str): A string representing the PDB ID of the file to be downloaded.
        save_path (str): A string representing the path to save the downloaded file.  # TODO: find out if the path is a filepath
    
    Raises:
        exception.IllegalArgumentError: If either `pdb_id` or `save_path` is None.")
        exception.DirectoryNotFoundError: If save_path does not exist.
    """
    # <editor-fold desc="Checks">
    if pdb_id is None:
        logger.error("pdb_id is None.")
        raise exception.IllegalArgumentError("pdb_id is None.")
    if save_path is None:
        logger.error("save_path is None.")
        raise exception.IllegalArgumentError("save_path is None.")
    if not os.path.exists(save_path):
        logger.error("save_path does not exist.")
        raise exception.DirectoryNotFoundError("save_path does not exist.")
    
    # </editor-fold>
    
    pdb_url = f'https://files.rcsb.org/download/{pdb_id.upper()}.pdb'

    try:
        response = requests.get(pdb_url)
        response.raise_for_status()  # Check for errors
        with open(save_path, 'w') as pdb_file:
            pdb_file.write(response.text)
        logger.error(f"PDB file {pdb_id} downloaded successfully to {save_path}")
    except requests.exceptions.HTTPError as errh:
        logger.error(f"HTTP Error: {errh}")
    except requests.exceptions.ConnectionError as errc:
        logger.error(f"Error Connecting: {errc}")
    except requests.exceptions.Timeout as errt:
        logger.error(f"Timeout Error: {errt}")
    except requests.exceptions.RequestException as err:
        logger.error(f"Error: {err}")


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


def build_pdb_file(records: list, a_filepath: str) -> None:
    """Builds a pdb file for the given records.

    Args:
        records (list): A list of dictionaries containing the information of each record.
        a_filepath (str): The file path of the PDB file to be created.
    
    Raises:
        exception.IllegalArgumentError: If records is either None or an empty list or a_filepath is None or an empty string.
    """
    # <editor-fold desc="Checks">
    if records is None or len(records) == 0:
        logger.error("records is either None or an empty list.")
        raise exception.IllegalArgumentError("records is either None or an empty list.")
    if a_filepath is None or a_filepath == "":
        logger.error("a_filepath is either None or an empty string.")
        raise exception.IllegalArgumentError("a_filepath is either None or an empty string.")
    if os.path.exists(a_filepath):
        logger.warning("a_filepath already exists.")
    
    # </editor-fold>
    
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
                except ValueError:  # fixme: this could be the fix for the protein import problem
                    logger.warning("PDB file has only the record type for TER.")
                    pdb_line = f"{record_type:3s}\n"
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


def parse_pdb_file(a_filepath: str) -> list[dict]:
    """Parses an existing pdb file.

    Args:
        a_filepath: The filepath of the PDB file to be parsed.

    Returns:
        A list of dictionaries containing the parsed records from the PDB file.

    Raises:
        exception.IllegalArgumentError: If `a_filepath` is either None or an empty string.
        FileNotFoundError: If `a_filepath` does not exist.
    """
    # <editor-fold desc="Checks">
    if a_filepath is None or a_filepath == "":
        logger.error("a_filepath is either None or an empty string.")
        raise exception.IllegalArgumentError("a_filepath is either None or an empty string.")
    if not os.path.exists(a_filepath):
        logger.error("a_filepath does not exist.")
        raise FileNotFoundError("a_filepath does not exist.")
    
    # </editor-fold>
    
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
                    'charge': charge,
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
                        'charge': "",
                    })
                except ValueError:
                    logger.warning("PDB file has only the record type for TER.")
                    records.append({
                        'record_type': record_type,
                        'atom_number': "",
                        'atom_name': "",
                        'alternate_location_indicator': "",
                        'residue_name': "",
                        'chain_identifier': "",
                        'residue_sequence_number': "",
                        'code_for_insertions_of_residues': "",
                        'x_coord': 0.0,
                        'y_coord': 0.0,
                        'z_coord': 0.0,
                        'occupancy': 0.0,
                        'temperature_factor': 0.0,
                        'segment_identifier': "",
                        'element_symbol': "",
                        'charge': "",
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
                    'charge': charge_het,
                })
    return records
