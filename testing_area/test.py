import re
from datetime import datetime
import sqlite3
import requests

# def parse_pdb_file(file_path):
#     pdb_entry = {}
#     chains = []
#     residues = []
#     atoms = []
#
#     with open(file_path, 'r') as pdb_file:
#         i = 0
#         for line in pdb_file:
#             record_type = line[:6].strip()
#
#             if record_type == 'HEADER':
#                 # Parse HEADER record for PDB ID and title
#                 pdb_id = line[62:66].strip()
#                 title = line[10:50].strip()
#                 pdb_entry = {
#                     'pdb_id': pdb_id,
#                     'title': title,
#                 }
#
#             elif record_type == 'ATOM':
#                 # Parse ATOM record for chain, residue, and atomic coordinates
#                 chain_id = line[21]
#                 residue_number = int(line[22:26].strip())
#                 residue_name = line[17:20].strip()
#                 atom_name = line[12:16].strip()
#                 x_coord = float(line[30:38].strip())
#                 y_coord = float(line[38:46].strip())
#                 z_coord = float(line[46:54].strip())
#                 occupancy = float(line[54:60].strip())
#                 temperature_factor = float(line[60:66].strip())
#
#                 chain = {
#                     'chain_id': i,
#                     'pdb_id': pdb_entry['pdb_id'],
#                     'chain_identifier': chain_id,
#                 }
#
#                 residue = {
#                     'residue_id': i,
#                     'pdb_id': pdb_entry['pdb_id'],
#                     'chain_identifier': chain_id,
#                     'residue_number': residue_number,
#                     'residue_name': residue_name,
#                 }
#
#                 atom = {
#                     'atom_id': i,
#                     'pdb_id': pdb_entry['pdb_id'],
#                     'chain_identifier': chain_id,
#                     'residue_number': residue_number,
#                     'residue_name': residue_name,
#                     'atom_name': atom_name,
#                     'x_coord': x_coord,
#                     'y_coord': y_coord,
#                     'z_coord': z_coord,
#                     'occupancy': occupancy,
#                     'temperature_factor': temperature_factor,
#                 }
#
#                 chains.append(chain)
#                 residues.append(residue)
#                 atoms.append(atom)
#                 i += 1
#
#     return pdb_entry, chains, residues, atoms


# Now, you can insert the parsed data into the database using the previously provided insert_data() function or a similar approach.

def parse_pdb_file(file_path):
    records = []

    with open(file_path, 'r') as pdb_file:
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
                atom_number_ter = int(line[6:11])
                residue_name_ter = line[17:20].strip()
                chain_identifier_ter = line[21]
                residue_sequence_number_ter = int(line[22:26])
                records.append({
                    'record_type': record_type,
                    'atom_number': atom_number_ter,
                    'residue_name': residue_name_ter,
                    'chain_identifier': chain_identifier_ter,
                    'residue_sequence_number': residue_sequence_number_ter
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
                    'residue_name': residue_name_het,
                    'chain_identifier': chain_identifier_het,
                    'residue_sequence_number': residue_sequence_number_het,
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

def build_pdb_file(records, output_file_path):
    with open(output_file_path, 'w') as pdb_file:
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
                pdb_line = f"{record_type:3s}   {record['atom_number']:5d}      {record['residue_name']:3s}" \
                           f" {record['chain_identifier']:1s}{record['residue_sequence_number']:4d}\n"
            elif record_type == 'HETATM':
                pdb_line = f"{record_type:6s}{record['atom_number']:5d}  {record['atom_name']:4s}" \
                           f"{record['residue_name']:3s} {record['chain_identifier']:1s}" \
                           f"{record['residue_sequence_number']:4d}    {record['x_coord']:8.3f}{record['y_coord']:8.3f}" \
                           f"{record['z_coord']:8.3f}{record['occupancy']:6.2f}{record['temperature_factor']:6.2f}           " \
                           f"{record['element_symbol']:2s}{record['charge']:2s}\n"
            else:
                raise ValueError(f"Unsupported record type: {record_type}")

            pdb_file.write(pdb_line)



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



if __name__ == '__main__':
    pdb_file_path = '3bmp_test.pdb'
    download_pdb_file("3BMP", pdb_file_path)

    # Example usage
    pdb_file_path = '3bmp.pdb'
    #parsed_pdb_entry, parsed_chains, parsed_residues, parsed_atoms = parse_pdb_file(pdb_file_path)
    #print(parse_pdb_file(pdb_file_path))

    parsed_data = parse_pdb_file(pdb_file_path)

    # Accessing the parsed data
    for record in parsed_data:
        print(record)

    output_file_path = 'output.pdb'
    build_pdb_file(parsed_data, output_file_path)