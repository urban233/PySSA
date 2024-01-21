import pathlib
import sqlite3

from pyssa.internal.data_structures import protein, project


class DatabaseManager:
    """Manages the project database."""
    _database_filepath: str
    _connection: sqlite3.Connection
    _cursor: sqlite3.Cursor

    def __init__(self, the_database_filepath: str) -> None:
        self._database_filepath = the_database_filepath

    def build_new_database(self) -> None:
        """Builds a new database for a new project."""
        tmp_connection = sqlite3.connect(self._database_filepath)
        tmp_cursor = tmp_connection.cursor()
        sql = """-- Project definition
            CREATE TABLE Project (
                id INTEGER NOT NULL,
                name TEXT NOT NULL,
                os TEXT,
                CONSTRAINT Project_PK PRIMARY KEY (id)
            );
        """
        tmp_cursor.execute(sql)
        sql = """-- SeqRecord definition
            CREATE TABLE SeqRecord (
                id INTEGER NOT NULL,
                seq_id TEXT,
                seq TEXT,
                name TEXT,
                project_id INTEGER,
                CONSTRAINT SeqRecord_PK PRIMARY KEY (id),
                CONSTRAINT SeqRecord_Project_FK FOREIGN KEY (project_id) REFERENCES Project(id)
            );
        """
        tmp_cursor.execute(sql)
        sql = """-- Protein definition
            CREATE TABLE Protein (
                id INTEGER NOT NULL,
                pymol_molecule_object TEXT,
                pdb_filepath TEXT,
                fasta_filepath TEXT,
                export_dirname TEXT,
                pymol_session_filepath TEXT,
                pymol_session TEXT,
                project_id INTEGER, pdb_id INTEGER,
                CONSTRAINT Protein_PK PRIMARY KEY (id),
                CONSTRAINT Protein_Project_FK FOREIGN KEY (project_id) REFERENCES Project(id)
            );
        """
        tmp_cursor.execute(sql)
        sql = """-- "Chain" definition
            CREATE TABLE "Chain" (
            id INTEGER NOT NULL,
            protein_id INTEGER,
            chain_identifier TEXT,
            chain_type TEXT, chain_sequence TEXT,
            CONSTRAINT Chain_PK PRIMARY KEY (id),
            CONSTRAINT Chain_Protein_FK FOREIGN KEY (protein_id) REFERENCES Protein(id));
        """
        tmp_cursor.execute(sql)
        sql = """-- PyMOLParameter definition
            CREATE TABLE PyMOLParameter (
                id INTEGER NOT NULL,
                color TEXT,
                representation TEXT,
                chain_id INTEGER,
                CONSTRAINT PyMOLParameter_PK PRIMARY KEY (id),
                CONSTRAINT PyMOLParameter_Chain_FK FOREIGN KEY (chain_id) REFERENCES "Chain"(id)
            );
        """
        tmp_cursor.execute(sql)
        sql = """-- PyMOLSelection definition
            CREATE TABLE PyMOLSelection (
                id INTEGER NOT NULL,
                selection_string TEXT,
                protein_id INTEGER,
                CONSTRAINT PyMOLSelection_PK PRIMARY KEY (id),
                CONSTRAINT PyMOLSelection_Protein_FK FOREIGN KEY (protein_id) REFERENCES Protein(id)
            );
        """
        tmp_cursor.execute(sql)
        sql = """-- PdbAtom definition
            CREATE TABLE PdbAtom (
                id INTEGER NOT NULL,
                record_type TEXT(6),
                atom_number INTEGER,
                atom_name TEXT(4),
                alternate_location_indicator TEXT(1),
                residue_name TEXT(3),
                chain_identifier TEXT(1),
                residue_sequence_number INTEGER,
                code_for_insertions_of_residues TEXT(1),
                x_coord REAL,
                y_coord REAL,
                z_coord REAL,
                occupancy REAL,
                temperature_factor REAL,
                segment_identifier TEXT(4),
                element_symbol TEXT(2),
                charge TEXT(2),
                protein_id INTEGER,
                CONSTRAINT PdbAtom_PK PRIMARY KEY (id),
                CONSTRAINT PdbAtom_Protein_FK FOREIGN KEY (protein_id) REFERENCES Protein(id)
            );
        """
        tmp_cursor.execute(sql)
        sql = """-- ProteinPair definition
            CREATE TABLE ProteinPair (
                id INTEGER NOT NULL,
                protein_1_id INTEGER,
                protein_2_id INTEGER,
                pymol_session_filepath TEXT,
                pymol_session TEXT,
                project_id INTEGER,
                CONSTRAINT ProteinPair_PK PRIMARY KEY (id),
                CONSTRAINT ProteinPair_Protein_FK FOREIGN KEY (protein_1_id) REFERENCES Protein(id),
                CONSTRAINT ProteinPair_Protein_FK_1 FOREIGN KEY (protein_2_id) REFERENCES Protein(id),
                CONSTRAINT ProteinPair_Project_FK FOREIGN KEY (project_id) REFERENCES Project(id)
            );
        """
        tmp_cursor.execute(sql)
        sql = """-- PyMOLParameterProteinPair definition
            CREATE TABLE PyMOLParameterProteinPair (
                id INTEGER NOT NULL,
                protein_id INTEGER,
                chain_letter TEXT,
                parameter_name TEXT,
                parameter_value TEXT,
                protein_pair_id INTEGER,
                CONSTRAINT PyMOLParameterProteinPair_PK PRIMARY KEY (id),
                CONSTRAINT PyMOLParameterProteinPair_ProteinPair_FK FOREIGN KEY (protein_pair_id) REFERENCES ProteinPair(id)
            );
        """
        tmp_cursor.execute(sql)
        sql = """-- DistanceAnalysis definition
            CREATE TABLE DistanceAnalysis (
                id INTEGER NOT NULL,
                name TEXT,
                cutoff REAL,
                cycles INTEGER,
                protein_pair_id INTEGER, figure_size_x REAL, figure_size_y REAL,
                CONSTRAINT DistanceAnalysis_PK PRIMARY KEY (id),
                CONSTRAINT DistanceAnalysis_ProteinPair_FK FOREIGN KEY (protein_pair_id) REFERENCES ProteinPair(id)
            );
        """
        tmp_cursor.execute(sql)
        sql = """-- DistanceAnalysisResults definition
            CREATE TABLE DistanceAnalysisResults (
                id INTEGER NOT NULL,
                pymol_session TEXT,
                rmsd REAL,
                aligned_aa TEXT,
                distance_analysis_id INTEGER,
                CONSTRAINT DistanceAnalysisResults_PK PRIMARY KEY (id),
                CONSTRAINT DistanceAnalysisResults_DistanceAnalysis_FK FOREIGN KEY (distance_analysis_id) REFERENCES DistanceAnalysis(id)
            );
        """
        tmp_cursor.execute(sql)
        sql = """-- DistanceAnalysisResultData definition
            CREATE TABLE DistanceAnalysisResultData (
                id INTEGER NOT NULL,
                "index" INTEGER,
                protein_1_chain TEXT(2),
                protein_1_position INTEGER,
                protein_1_residue TEXT(3),
                protein_2_chain TEXT(2),
                protein_2_position INTEGER,
                protein_2_residue TEXT(3),
                distances REAL,
                distance_analysis_results_id INTEGER,
                CONSTRAINT DistanceAnalysisResultData_PK PRIMARY KEY (id),
                CONSTRAINT DistanceAnalysisResultData_DistanceAnalysisResults_FK FOREIGN KEY (id) REFERENCES DistanceAnalysisResults(id)
            );
        """
        tmp_cursor.execute(sql)
        tmp_connection.close()

    def set_database_filepath(self, a_new_database_filepath: str) -> None:
        """Sets a new database filepath."""
        self._database_filepath = a_new_database_filepath

    def open_project_database(self):
        """Opens a project database"""
        self._connection = sqlite3.connect(self._database_filepath)
        self._cursor = self._connection.cursor()

    def close_project_database(self):
        """Closes a project database."""
        self._cursor.close()
        self._connection.close()

    def write_new_empty_project(self, a_project_name, an_os):
        """Writes a new empty project to the database."""
        sql = """   INSERT INTO Project(name, os)
                    VALUES(:name, :os)
        """
        self._cursor.execute(sql, {"name": a_project_name, "os": an_os})
        self._connection.commit()

    def write_new_protein(self, a_protein: protein.Protein):
        """Writes a new protein to the database."""
        sql = """   INSERT INTO Protein(pymol_molecule_object, pymol_session)
                    VALUES (:pymol_molecule_object, :pymol_session)
        """
        self._cursor.execute(sql, a_protein.get_object_as_dict_for_database())
        self._connection.commit()

    def get_number_of_sequences(self):
        sql = """SELECT id FROM SeqRecord"""
        self._cursor.execute(sql)
        return len(self._cursor.fetchall())

    def get_number_of_proteins(self):
        """Returns all proteins in the database."""
        sql = """SELECT id FROM Protein"""
        self._cursor.execute(sql)
        return len(self._cursor.fetchall())

    def get_number_of_protein_pairs(self):
        sql = """SELECT id FROM ProteinPair"""
        self._cursor.execute(sql)
        return len(self._cursor.fetchall())

    def get_all_protein_table_data(self):
        """Gets all records from the Protein table.

        Returns:
            a list of protein records (id, pymol_molecule_object, pdb_filepath, fasta_filepath,
            export_dirname, pymol_session_filepath, pymol_session)
        """
        sql = """
            SELECT id, pymol_molecule_object, pdb_filepath, fasta_filepath, 
            export_dirname, pymol_session_filepath, pymol_session FROM Protein"""
        self._cursor.execute(sql)
        return self._cursor.fetchall()

    def get_all_chain_table_data(self):
        """Gets all records from the Chain table.

        Returns:
            a list of chain records (id, protein_id, chain_identifier, chain_type, chain_sequence)
        """
        sql = """
            SELECT id, protein_id, chain_identifier, chain_type, chain_sequence 
            FROM Chain"""
        self._cursor.execute(sql)
        return self._cursor.fetchall()

    def get_all_pymol_parameter_table_data(self):
        """Gets all records from the PyMOLParameter table.

        Returns:
            a list of pymol parameter records (id, color, representation, chain_id)
        """
        sql = """
            SELECT id, color, representation, chain_id
            FROM PyMOLParameter"""
        self._cursor.execute(sql)
        return self._cursor.fetchall()

    def get_all_pymol_selection_table_data(self):
        """Gets all records from the PyMOLSelection table.

        Returns:
            a list of pymol selection records (id, selection_string, protein_id)
        """
        sql = """
            SELECT id, selection_string, protein_id
            FROM PyMOLSelection"""
        self._cursor.execute(sql)
        return self._cursor.fetchall()

    def get_all_pdb_atom_table_data(self):
        """Gets all records from the PdbAtom table.

        Returns:
            a list of pdb atom records (id, record_type, atom_number, atom_name, alternate_location_indicator, residue_name,
            chain_identifier, residue_sequence_number, code_for_insertions_of_residues,
            x_coord, y_coord, z_coord, occupancy, temperature_factor, segment_identifier, element_symbol,
            charge, protein_id)
        """
        sql = """
            SELECT id, record_type, atom_number, atom_name, alternate_location_indicator, residue_name,
            chain_identifier, residue_sequence_number, code_for_insertions_of_residues, 
            x_coord, y_coord, z_coord, occupancy, temperature_factor, segment_identifier, element_symbol, 
            charge, protein_id
            FROM PdbAtom"""
        self._cursor.execute(sql)
        return self._cursor.fetchall()

    def get_project_as_object(self, a_project_name: str, a_workspace_path: pathlib.Path) -> "project.Project":
        tmp_project = project.Project(a_project_name, a_workspace_path)
        if self.get_number_of_sequences() > 0:
            # There are one or more sequences in the db
            sql = """SELECT seq_id, seq, name FROM SeqRecord"""
            self._cursor.execute(sql)
            tmp_seqs_table_information = self._cursor.fetchall()
        if self.get_number_of_proteins() > 0:
            # There are one or more proteins in the db
            tmp_protein_table_data = self.get_all_protein_table_data()
            tmp_chain_table_data = self.get_all_chain_table_data()

            for tmp_protein_data in tmp_protein_table_data:
                protein_id, pymol_molecule_object, _, _, _, _, pymol_session = tmp_protein_data
                self._cursor.execute('''
                        SELECT 
                            Protein.*,
                            Chain.*,
                            PyMOLParameter.*,
                            PyMOLSelection.*,
                            PdbAtom.*
                        FROM Protein
                        JOIN Chain ON Protein.id = Chain.protein_id
                        JOIN PyMOLParameter ON Chain.id = PyMOLParameter.chain_id
                        LEFT JOIN PyMOLSelection ON Protein.id = PyMOLSelection.protein_id
                        LEFT JOIN PdbAtom ON Protein.id = PdbAtom.protein_id
                        WHERE Protein.id = ?
                    ''', (protein_id,))
                rows = self._cursor.fetchall()
                # Display the fetched data (you can customize this based on your needs)
                for row in rows:
                    print(row)
                tmp_protein = protein.Protein(pymol_molecule_object)
                tmp_protein.pymol_session = pymol_session
        if self.get_number_of_protein_pairs() > 0:
            # There are one or more protein_pairs in the db
            pass
        
        return tmp_project
