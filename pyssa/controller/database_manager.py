import logging
import pathlib
import sqlite3

from pyssa.internal.data_structures import protein, project, protein_pair, structure_analysis, results
from pyssa.logging_pyssa import log_handlers
from pyssa.util import enums, pyssa_keys

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class DatabaseManager:
    """Manages the project database."""
    _database_filepath: str
    _connection: sqlite3.Connection
    _cursor: sqlite3.Cursor

    def __init__(self, the_database_filepath: str) -> None:
        self._database_filepath = the_database_filepath

    # <editor-fold desc="Base methods">
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
                pymol_session TEXT,
                project_id INTEGER,
                pdb_id INTEGER,
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
                            "my_index" INTEGER,
                            protein_1_chain TEXT,
                            protein_1_position INTEGER,
                            protein_1_residue TEXT,
                            protein_2_chain TEXT,
                            protein_2_position INTEGER,
                            protein_2_residue TEXT,
                            distances REAL,
                            distance_analysis_results_id INTEGER,
                            CONSTRAINT DistanceAnalysisResultData_PK PRIMARY KEY (id),
                            CONSTRAINT DistanceAnalysisResultData_DistanceAnalysisResults_FK FOREIGN KEY (distance_analysis_results_id) REFERENCES DistanceAnalysisResults(id)
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

    # </editor-fold>

    def insert_new_project(self, a_project_name, an_os) -> int:
        """Writes a new empty project to the database."""
        sql = """   INSERT INTO Project(name, os)
                    VALUES(:name, :os)
        """
        self._cursor.execute(sql, {"name": a_project_name, "os": an_os})
        self._connection.commit()
        return self.get_last_id()

    # <editor-fold desc="Protein object inserts">
    def insert_new_protein(self, a_protein: "protein.Protein") -> int:
        """Writes a new protein to the database."""
        tmp_protein_id = self._insert_protein(a_protein.get_object_as_dict_for_database())
        a_protein.set_id(tmp_protein_id)
        if len(a_protein.chains) == 0:
            logger.warning("There are no chains to be inserted into the database.")
        else:
            # chains will be inserted
            for tmp_chain in a_protein.chains:
                tmp_chain_id = self._insert_chain(tmp_protein_id, tmp_chain)
                self._insert_pymol_parameter(tmp_chain_id, tmp_chain.pymol_parameters)
        self._insert_pymol_selection(tmp_protein_id, a_protein.pymol_selection.selection_string)
        for tmp_pdb_atom_dict in a_protein.get_pdb_data():
            self._insert_pdb_atom(tmp_protein_id, tmp_pdb_atom_dict)
        return tmp_protein_id

    def _insert_protein(self, a_protein_obj_dict: dict) -> int:
        sql = """   INSERT INTO Protein(pymol_molecule_object, pymol_session, project_id)
                    VALUES (:pymol_molecule_object, :pymol_session, :project_id)
        """
        self._cursor.execute(sql, a_protein_obj_dict)
        self._connection.commit()
        return self.get_last_id()

    def _insert_chain(self, the_protein_id: int, a_chain) -> int:
        sql = """   INSERT INTO Chain(protein_id, chain_identifier, chain_type, chain_sequence)
                    VALUES (?, ?, ?, ?)
        """
        tmp_params = (the_protein_id, a_chain.chain_letter, a_chain.chain_type, a_chain.chain_sequence.sequence)
        self._cursor.execute(sql, tmp_params)
        self._connection.commit()
        return self.get_last_id()

    def _insert_pymol_parameter(self, the_chain_id, a_pymol_parameter_dict: dict):
        sql = """   INSERT INTO PyMOLParameter(color, representation, chain_id)
                    VALUES (?, ?, ?)
        """
        tmp_params = (
            a_pymol_parameter_dict[enums.PymolParameterEnum.COLOR.value],
            a_pymol_parameter_dict[enums.PymolParameterEnum.REPRESENTATION.value],
            the_chain_id
        )
        self._cursor.execute(sql, tmp_params)
        self._connection.commit()

    def _insert_pymol_selection(self, the_protein_id: int, a_selection_string: str):
        sql = """   INSERT INTO PyMOLSelection(selection_string, protein_id)
                    VALUES (?, ?)
        """
        tmp_params = (
            a_selection_string,
            the_protein_id,
        )
        self._cursor.execute(sql, tmp_params)
        self._connection.commit()

    def _insert_pdb_atom(self, the_protein_id: int, a_pdb_atom_dict: dict):
        sql = """   INSERT INTO PdbAtom(record_type, atom_number, atom_name, alternate_location_indicator, residue_name,
            chain_identifier, residue_sequence_number, code_for_insertions_of_residues, 
            x_coord, y_coord, z_coord, occupancy, temperature_factor, segment_identifier, element_symbol, 
            charge, protein_id)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
        tmp_params = (
            a_pdb_atom_dict["record_type"],
            a_pdb_atom_dict["atom_number"],
            a_pdb_atom_dict["atom_name"],
            a_pdb_atom_dict["alternate_location_indicator"],
            a_pdb_atom_dict["residue_name"],
            a_pdb_atom_dict["chain_identifier"],
            a_pdb_atom_dict["residue_sequence_number"],
            a_pdb_atom_dict["code_for_insertions_of_residues"],
            a_pdb_atom_dict["x_coord"],
            a_pdb_atom_dict["y_coord"],
            a_pdb_atom_dict["z_coord"],
            a_pdb_atom_dict["occupancy"],
            a_pdb_atom_dict["temperature_factor"],
            a_pdb_atom_dict["segment_identifier"],
            a_pdb_atom_dict["element_symbol"],
            a_pdb_atom_dict["charge"],
            the_protein_id
        )
        self._cursor.execute(sql, tmp_params)
        self._connection.commit()

    # </editor-fold>

    # <editor-fold desc="Protein pair object inserts">
    def insert_new_protein_pair(self, a_protein_pair: "protein_pair.ProteinPair"):
        """Inserts a new protein pair in the project database.

        Note:
            Run this method after distance analysis finished!!!
        """
        tmp_protein_pair_id = self._insert_protein_pair(a_protein_pair)
        a_protein_pair.set_id(tmp_protein_pair_id)
        for tmp_chain in a_protein_pair.protein_1.chains:
            self._insert_pymol_parameter_protein_pair(
                a_protein_pair.protein_1.get_id(),
                tmp_chain.chain_letter,
                enums.PymolParameterEnum.COLOR.value,
                tmp_chain.pymol_parameters[enums.PymolParameterEnum.COLOR.value],
                a_protein_pair.get_id(),
            )
            self._insert_pymol_parameter_protein_pair(
                a_protein_pair.protein_1.get_id(),
                tmp_chain.chain_letter,
                enums.PymolParameterEnum.REPRESENTATION.value,
                tmp_chain.pymol_parameters[enums.PymolParameterEnum.REPRESENTATION.value],
                a_protein_pair.get_id(),
            )
        for tmp_chain in a_protein_pair.protein_2.chains:
            self._insert_pymol_parameter_protein_pair(
                a_protein_pair.protein_2.get_id(),
                tmp_chain.chain_letter,
                enums.PymolParameterEnum.COLOR.value,
                tmp_chain.pymol_parameters[enums.PymolParameterEnum.COLOR.value],
                a_protein_pair.get_id(),
            )
            self._insert_pymol_parameter_protein_pair(
                a_protein_pair.protein_2.get_id(),
                tmp_chain.chain_letter,
                enums.PymolParameterEnum.REPRESENTATION.value,
                tmp_chain.pymol_parameters[enums.PymolParameterEnum.REPRESENTATION.value],
                a_protein_pair.get_id(),
            )
        tmp_distance_analysis_id = self._insert_distance_analysis(a_protein_pair.distance_analysis, tmp_protein_pair_id)
        tmp_distance_analysis_results_id = self._insert_distance_analysis_results(
            a_protein_pair.distance_analysis.analysis_results, tmp_distance_analysis_id
        )
        self._insert_distance_data_records(tmp_distance_analysis_results_id,
                                           a_protein_pair.distance_analysis.analysis_results.distance_data)

    def _insert_protein_pair(self, a_protein_pair: "protein_pair.ProteinPair") -> int:
        sql = """   INSERT INTO ProteinPair(protein_1_id, protein_2_id, pymol_session, project_id)
                    VALUES (?, ?, ?, ?)
        """
        tmp_params = (a_protein_pair.protein_1.get_id(), a_protein_pair.protein_2.get_id(),
                      a_protein_pair.pymol_session, a_protein_pair.db_project_id)
        self._cursor.execute(sql, tmp_params)
        self._connection.commit()
        return self.get_last_id()

    def _insert_pymol_parameter_protein_pair(self, a_protein_id,
                                             a_chain_letter,
                                             a_parameter_name,
                                             a_parameter_value,
                                             the_protein_pair_id):
        sql = """   INSERT INTO PyMOLParameterProteinPair(protein_id, chain_letter, parameter_name, parameter_value, protein_pair_id)
                    VALUES (?, ?, ?, ?, ?)
        """
        tmp_params = (a_protein_id, a_chain_letter, a_parameter_name, a_parameter_value, the_protein_pair_id)
        self._cursor.execute(sql, tmp_params)
        self._connection.commit()

    def _insert_distance_analysis(self,
                                  a_distance_analysis: structure_analysis.DistanceAnalysis,
                                  the_protein_pair_id) -> int:
        sql = """   INSERT INTO DistanceAnalysis(name, cutoff, cycles, protein_pair_id, figure_size_x, figure_size_y)
                    VALUES (?, ?, ?, ?, ?, ?)
        """
        tmp_params = (
            a_distance_analysis.name,
            a_distance_analysis.cutoff,
            a_distance_analysis.cycles,
            the_protein_pair_id,
            a_distance_analysis.figure_size[0],
            a_distance_analysis.figure_size[1]
        )
        self._cursor.execute(sql, tmp_params)
        self._connection.commit()
        return self.get_last_id()

    def _insert_distance_analysis_results(self, a_distance_analysis_result: results.DistanceAnalysisResults,
                                          the_distance_analysis_id: int) -> int:
        sql = """   INSERT INTO DistanceAnalysisResults(pymol_session, rmsd, aligned_aa, distance_analysis_id)
                    VALUES (?, ?, ?, ?)
        """
        tmp_params = (
            a_distance_analysis_result.pymol_session,
            a_distance_analysis_result.rmsd,
            a_distance_analysis_result.aligned_aa,
            the_distance_analysis_id
        )
        self._cursor.execute(sql, tmp_params)
        self._connection.commit()
        return self.get_last_id()

    def _insert_distance_data_records(self, the_distance_analysis_results_id: int, distance_data: dict):
        index = list(distance_data[pyssa_keys.ARRAY_DISTANCE_INDEX])
        prot_1_chains = list(distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN])
        prot_1_position = list(distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION])
        prot_1_residue = list(distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI])
        prot_2_chains = list(distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN])
        prot_2_position = list(distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION])
        prot_2_residue = list(distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI])
        distances = list(distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES])

        sql = """   INSERT INTO DistanceAnalysisResultData(my_index, protein_1_chain, protein_1_position, protein_1_residue,
                                                            protein_2_chain, protein_2_position, protein_2_residue, 
                                                            distances, distance_analysis_results_id)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """
        for i in range(len(index)):
            tmp_params = (
                int(index[i]),
                prot_1_chains[i],
                prot_1_position[i],
                prot_1_residue[i],
                prot_2_chains[i],
                prot_2_position[i],
                prot_2_residue[i],
                distances[i],
                the_distance_analysis_results_id,
            )
            logger.info(tmp_params)
            self._cursor.execute(sql, tmp_params)
            i += 1
        self._connection.commit()

    # </editor-fold>

    def get_last_id(self):
        sql = """SELECT last_insert_rowid();"""
        self._cursor.execute(sql)
        tmp_id, = self._cursor.fetchall()[0]
        return tmp_id

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

    # <editor-fold desc="Get all data of a specific table">
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

    # </editor-fold>

    def get_project_as_object(self, a_project_name: str, a_workspace_path: pathlib.Path) -> "project.Project":
        """Creates a project object based on the data from the project database."""
        tmp_project = project.Project(a_project_name, a_workspace_path)
        if self.get_number_of_sequences() > 0:
            # There are one or more sequences in the db
            sql = """SELECT seq_id, seq, name FROM SeqRecord"""
            self._cursor.execute(sql)
            tmp_seqs_table_information = self._cursor.fetchall()
        if self.get_number_of_proteins() > 0:
            # There are one or more proteins in the db
            tmp_protein_table_data = self.get_all_protein_table_data()
            for tmp_protein_data in tmp_protein_table_data:
                protein_id, pymol_molecule_object, _, _, _, _, pymol_session = tmp_protein_data
                sql = """
                    SELECT * 
                    FROM Protein
                    JOIN Chain ON Protein.id = Chain.protein_id
                    WHERE Protein.id = ?
                """
                self._cursor.execute(sql, (protein_id,))
                # self._cursor.execute('''
                #         SELECT
                #             Protein.*,
                #             Chain.*,
                #             PyMOLParameter.*,
                #             PyMOLSelection.*,
                #             PdbAtom.*
                #         FROM Protein
                #         JOIN Chain ON Protein.id = Chain.protein_id
                #         JOIN PyMOLParameter ON Chain.id = PyMOLParameter.chain_id
                #         LEFT JOIN PyMOLSelection ON Protein.id = PyMOLSelection.protein_id
                #         LEFT JOIN PdbAtom ON Protein.id = PdbAtom.protein_id
                #         WHERE Protein.id = ?
                #     ''', (protein_id,))
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
