import re
from datetime import datetime
import sqlite3


def parse_pdb_file(file_path):
    pdb_entry = {}
    chains = []
    residues = []
    atoms = []

    with open(file_path, 'r') as pdb_file:
        i = 0
        for line in pdb_file:
            record_type = line[:6].strip()

            if record_type == 'HEADER':
                # Parse HEADER record for PDB ID and title
                pdb_id = line[62:66].strip()
                title = line[10:50].strip()
                pdb_entry = {
                    'pdb_id': pdb_id,
                    'title': title,
                }

            elif record_type == 'ATOM':
                # Parse ATOM record for chain, residue, and atomic coordinates
                chain_id = line[21]
                residue_number = int(line[22:26].strip())
                residue_name = line[17:20].strip()
                atom_name = line[12:16].strip()
                x_coord = float(line[30:38].strip())
                y_coord = float(line[38:46].strip())
                z_coord = float(line[46:54].strip())
                occupancy = float(line[54:60].strip())
                temperature_factor = float(line[60:66].strip())

                chain = {
                    'chain_id': i,
                    'pdb_id': pdb_entry['pdb_id'],
                    'chain_identifier': chain_id,
                }

                residue = {
                    'residue_id': i,
                    'pdb_id': pdb_entry['pdb_id'],
                    'chain_identifier': chain_id,
                    'residue_number': residue_number,
                    'residue_name': residue_name,
                }

                atom = {
                    'atom_id': i,
                    'pdb_id': pdb_entry['pdb_id'],
                    'chain_identifier': chain_id,
                    'residue_number': residue_number,
                    'residue_name': residue_name,
                    'atom_name': atom_name,
                    'x_coord': x_coord,
                    'y_coord': y_coord,
                    'z_coord': z_coord,
                    'occupancy': occupancy,
                    'temperature_factor': temperature_factor,
                }

                chains.append(chain)
                residues.append(residue)
                atoms.append(atom)
                i += 1

    return pdb_entry, chains, residues, atoms

# Now, you can insert the parsed data into the database using the previously provided insert_data() function or a similar approach.


if __name__ == '__main__':
    conn = sqlite3.connect(r"C:\Users\martin\.pyssa\db_workspace\project_1")
    cursor = conn.cursor()

    sql = ''' INSERT INTO Project(name,os)
                  VALUES(:name, :os) '''

    #cursor.execute(sql, {'name': "my_project", 'os': "win32"})
    sql = '''   UPDATE Project SET name=:name, os=:os 
                WHERE id=1
    '''
    cursor.execute(sql, {'name': "mybankkonto", 'os': "win64"})
    conn.commit()

    # Example usage
    pdb_file_path = '3bmp.pdb'
    #parsed_pdb_entry, parsed_chains, parsed_residues, parsed_atoms = parse_pdb_file(pdb_file_path)


    # cursor.execute('''
    #         SELECT *
    #         FROM atom
    #     ''')
    # # Fetch all rows
    # rows = cursor.fetchall()
    #
    # # Display the fetched data (you can customize this based on your needs)
    # for row in rows:
    #     print(row)

    conn.close()
