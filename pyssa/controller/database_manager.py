import logging
import pathlib
import sqlite3
import numpy as np
from PyQt5 import QtCore
from Bio import SeqRecord

from pyssa.internal.data_structures import protein, project, protein_pair, structure_analysis, results, chain, sequence
from pyssa.internal.thread.async_pyssa import custom_signals
from pyssa.logging_pyssa import log_handlers
from pyssa.util import enums, pyssa_keys

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class DatabaseManager:
    """Manages the project database."""
    _database_filepath: str
    _connection: sqlite3.Connection
    _cursor: sqlite3.Cursor
    _application_settings: "settings.Settings"

    def __init__(self, the_database_filepath: str) -> None:
        self._database_filepath = the_database_filepath
        self._connection = None

    # <editor-fold desc="Util methods">
    def get_last_id(self):
        sql = """SELECT last_insert_rowid();"""
        self._cursor.execute(sql)
        tmp_id, = self._cursor.fetchall()[0]
        return tmp_id

    def get_database_filepath(self):
        return self._database_filepath

    def set_application_settings(self, the_app_settings):
        self._application_settings = the_app_settings

    def get_latest_id_of_protein_table(self):
        self.open_project_database()
        # Assuming your_table_name has an auto-incrementing primary key column named 'id'
        self._cursor.execute("SELECT MAX(id) FROM Protein")
        latest_id_before_insert = self._cursor.fetchone()[0]
        self.close_project_database()
        return latest_id_before_insert if latest_id_before_insert is not None else 0

    def get_latest_id_of_a_specific_table(self, a_table_name):
        self.open_project_database()
        # Assuming your_table_name has an auto-incrementing primary key column named 'id'
        self._cursor.execute(f"SELECT MAX(id) FROM {a_table_name}")
        latest_id_before_insert = self._cursor.fetchone()[0]
        return latest_id_before_insert if latest_id_before_insert is not None else 0

    # </editor-fold>

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
                project_id INTEGER, name TEXT,
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
        logger.debug(f"Filepath of the database_manager: {self._database_filepath}")
        self._connection = sqlite3.connect(self._database_filepath)
        self._cursor = self._connection.cursor()

    def close_project_database(self):
        """Closes a project database."""
        self._cursor.close()
        self._connection.close()

    def __enter__(self):
        self._connection = sqlite3.connect(self._database_filepath)
        self._cursor = self._connection.cursor()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self._connection:
            self._connection.close()

    # </editor-fold>

    # <editor-fold desc="SeqRecord objects inserts">
    def insert_new_sequence(self, a_seq_record: SeqRecord.SeqRecord):
        self._insert_sequence(a_seq_record.id, a_seq_record.seq, a_seq_record.name)

    def _insert_sequence(self, a_seq_id, a_seq, a_name):
        sql = """   INSERT INTO SeqRecord(seq_id, seq, name, project_id)
                    VALUES (?, ?, ?, ?)
        """
        self._cursor.execute(sql, (str(a_seq_id), str(a_seq), str(a_name), 1))  # fixme: not the best solution with the id=1
        self._connection.commit()
        return self.get_last_id()
    # </editor-fold>

    # <editor-fold desc="Delete statements for SeqRecord object">
    def delete_existing_sequence(self, a_seq_record_name):
        sql = """   
            DELETE FROM SeqRecord
            WHERE name = ?
        """
        self._cursor.execute(sql, (a_seq_record_name,))
        self._connection.commit()

    # </editor-fold>

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

    # <editor-fold desc="Delete statements for protein object">
    def delete_existing_protein(self, a_protein_id: int):
        self._delete_pdb_atom(a_protein_id)
        self._delete_pymol_selection(a_protein_id)
        for tmp_chain_info in self._get_chain(a_protein_id):
            self._delete_pymol_parameter(tmp_chain_info[0])
        self._delete_chains(a_protein_id)
        self._delete_protein(a_protein_id)

    def _delete_protein(self, the_protein_id: int):
        sql = """   
            DELETE FROM Protein
            WHERE id = ?
        """
        self._cursor.execute(sql, (the_protein_id,))
        self._connection.commit()

    def _delete_chains(self, the_protein_id: int):
        sql = """   
            DELETE FROM Chain
            WHERE protein_id = ?
        """
        self._cursor.execute(sql, (the_protein_id,))
        self._connection.commit()

    def _delete_pymol_parameter(self, the_chain_id):
        sql = """   
            DELETE FROM PyMOLParameter
            WHERE chain_id = ?
        """
        self._cursor.execute(sql, (the_chain_id,))
        self._connection.commit()

    def _delete_pymol_selection(self, the_protein_id: int):
        sql = """   
            DELETE FROM PyMOLSelection
            WHERE protein_id = ?
        """
        self._cursor.execute(sql, (the_protein_id,))
        self._connection.commit()

    def _delete_pdb_atom(self, the_protein_id: int):
        sql = """   DELETE FROM PdbAtom
                    WHERE protein_id = ?
        """
        self._cursor.execute(sql, (the_protein_id,))
        self._connection.commit()

    # </editor-fold>

    # <editor-fold desc="Protein pair object inserts">
    def insert_new_protein_pair(self, a_protein_pair: "protein_pair.ProteinPair") -> int:
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
        return tmp_protein_pair_id

    def _insert_protein_pair(self, a_protein_pair: "protein_pair.ProteinPair") -> int:
        sql = """   INSERT INTO ProteinPair(protein_1_id, protein_2_id, pymol_session, project_id, name)
                    VALUES (?, ?, ?, ?, ?)
        """
        tmp_params = (a_protein_pair.protein_1.get_id(), a_protein_pair.protein_2.get_id(),
                      a_protein_pair.pymol_session, a_protein_pair.db_project_id, a_protein_pair.name)
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

    # <editor-fold desc="Protein pair delete statements">
    def delete_existing_protein_pair(self, a_protein_pair_id: int) -> None:
        # Get ids of distance analysis related objects
        tmp_distance_analysis_info = self._get_distance_analysis(a_protein_pair_id)
        tmp_distance_analysis_id, _, _, _, _, _ = tmp_distance_analysis_info
        tmp_distance_analysis_results_info = self._get_distance_analysis_results(tmp_distance_analysis_id)
        tmp_dist_analysis_results_id, _, _, _ = tmp_distance_analysis_results_info

        # delete statements
        self._delete_distance_analysis_result_data(tmp_dist_analysis_results_id)
        self._delete_distance_analysis_results(tmp_distance_analysis_id)
        self._delete_distance_analysis(a_protein_pair_id)
        self._delete_pymol_parameters_protein_pair(a_protein_pair_id)
        self._delete_protein_pair(a_protein_pair_id)

    def _delete_protein_pair(self, a_protein_pair_id: int) -> None:
        sql = """   
            DELETE FROM ProteinPair
            WHERE id = ?
        """
        self._cursor.execute(sql, (a_protein_pair_id,))
        self._connection.commit()

    def _delete_pymol_parameters_protein_pair(self, a_protein_pair_id: int) -> None:
        sql = """   
            DELETE FROM PyMOLParameterProteinPair
            WHERE protein_pair_id = ?
        """
        self._cursor.execute(sql, (a_protein_pair_id,))
        self._connection.commit()

    def _delete_distance_analysis(self, a_protein_pair_id: int) -> None:
        sql = """   
            DELETE FROM DistanceAnalysis
            WHERE protein_pair_id = ?
        """
        self._cursor.execute(sql, (a_protein_pair_id,))
        self._connection.commit()

    def _delete_distance_analysis_results(self, a_distance_analysis_id: int) -> None:
        sql = """   
            DELETE FROM DistanceAnalysisResults
            WHERE distance_analysis_id = ?
        """
        self._cursor.execute(sql, (a_distance_analysis_id,))
        self._connection.commit()

    def _delete_distance_analysis_result_data(self, a_distance_analysis_results_id: int) -> None:
        sql = """   
            DELETE FROM DistanceAnalysisResultData
            WHERE distance_analysis_results_id = ?
        """
        self._cursor.execute(sql, (a_distance_analysis_results_id,))
        self._connection.commit()
    # </editor-fold>

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

    # <editor-fold desc="Project related statements">
    def insert_new_project(self, a_project_name, an_os) -> int:
        """Writes a new empty project to the database."""
        sql = """   INSERT INTO Project(name, os)
                    VALUES(:name, :os)
        """
        self._cursor.execute(sql, {"name": a_project_name, "os": an_os})
        self._connection.commit()
        return self.get_last_id()

    def get_project_as_object(self,
                              a_project_name: str,
                              a_workspace_path: pathlib.Path,
                              the_app_settings,
                              the_progress_signal: "custom_signals.ProgressSignal" = custom_signals.ProgressSignal()) -> "project.Project":
        """Creates a project object based on the data from the project database."""
        tmp_project = project.Project(a_project_name, a_workspace_path)
        tmp_project_id, = self._get_project(a_project_name)
        tmp_project.set_id(tmp_project_id)

        tmp_progress = 0
        i = 1
        for tmp_seq_info in self._get_sequence():
            tmp_progress = min(tmp_progress + i, 40)
            the_progress_signal.emit_signal(f"Loading sequence ({i}/{len(self._get_sequence())}) ...", tmp_progress)
            tmp_seq_id, tmp_seq, tmp_seq_name = tmp_seq_info
            tmp_sequence = SeqRecord.SeqRecord(tmp_seq, id=tmp_seq_id, name=tmp_seq_name)
            tmp_project.sequences.append(tmp_sequence)
            i += 1

        i = 1
        # create protein objects
        for tmp_protein_info in self._get_protein(tmp_project_id):
            tmp_progress = 20
            tmp_message = f"Loading protein ({i}/{len(self._get_protein(tmp_project_id))}) ..."
            the_progress_signal.emit_signal(tmp_message, tmp_progress)

            tmp_protein_id, tmp_protein_name, tmp_pymol_session = tmp_protein_info
            tmp_protein = protein.Protein(tmp_protein_name)
            tmp_protein.set_id(tmp_protein_id)
            tmp_protein.pymol_session = tmp_pymol_session

            # Chains
            j = 5
            for tmp_chain_info in self._get_chain(tmp_protein_id):
                tmp_progress += j
                the_progress_signal.emit_signal(tmp_message, min(tmp_progress, 50))

                tmp_chain_id, tmp_chain_identifier, tmp_chain_type, tmp_chain_sequence = tmp_chain_info
                tmp_seq = sequence.Sequence(tmp_protein_name, tmp_chain_sequence)
                tmp_chain = chain.Chain(tmp_chain_identifier, tmp_seq, tmp_chain_type)
                tmp_chain.set_id(tmp_chain_id)
                tmp_chain.db_protein_id = tmp_protein_id
                tmp_color, tmp_representation = self._get_pymol_parameter(tmp_chain_id)
                tmp_chain.pymol_parameters = {
                    enums.PymolParameterEnum.COLOR.value: tmp_color,
                    enums.PymolParameterEnum.REPRESENTATION.value: tmp_representation,
                }
                tmp_protein.chains.append(tmp_chain)
                j += 3
            # pymol selection
            tmp_selection, = self._get_pymol_selection(tmp_protein_id)
            tmp_protein.pymol_selection.selection_string = tmp_selection
            # add protein to project
            tmp_project.proteins.append(tmp_protein)
            i += 1

        i = 1
        # create protein pair objects
        for tmp_protein_pair_info in self._get_protein_pair(tmp_project_id):
            tmp_progress = 20
            tmp_message = f"Loading protein pair ({i}/{len(self._get_protein_pair(tmp_project_id))}) ..."
            the_progress_signal.emit_signal(tmp_message, tmp_progress)
            # create protein pair
            tmp_protein_pair_id, tmp_protein_1_id, tmp_protein_2_id, tmp_pymol_session, tmp_pp_name = tmp_protein_pair_info
            tmp_protein_1_name, = self._get_protein_name_by_id(tmp_protein_1_id)
            tmp_protein_2_name, = self._get_protein_name_by_id(tmp_protein_2_id)
            tmp_protein_pair = protein_pair.ProteinPair(tmp_project.search_protein(tmp_protein_1_name),
                                                        tmp_project.search_protein(tmp_protein_2_name))
            tmp_protein_pair.set_id(tmp_protein_pair_id)
            tmp_protein_pair.db_project_id = tmp_project_id
            tmp_protein_pair.name = tmp_pp_name
            tmp_protein_pair.pymol_session = tmp_pymol_session

            tmp_progress += 20
            the_progress_signal.emit_signal(tmp_message, tmp_progress)
            # create distance analysis object
            tmp_distance_analysis_info = self._get_distance_analysis(tmp_protein_pair_id)
            tmp_distance_analysis_id, tmp_name, tmp_cutoff, tmp_cycles, tmp_figure_size_x, tmp_figure_size_y = tmp_distance_analysis_info
            tmp_distance_analysis = structure_analysis.DistanceAnalysis(the_app_settings)
            tmp_distance_analysis.name = tmp_name
            tmp_distance_analysis.cutoff = tmp_cutoff
            tmp_distance_analysis.cycles = tmp_cycles
            tmp_distance_analysis.figure_size = (tmp_figure_size_x, tmp_figure_size_y)

            tmp_protein_pair.distance_analysis = tmp_distance_analysis

            tmp_progress += 10
            the_progress_signal.emit_signal(tmp_message, tmp_progress)
            # create distance analysis results object
            tmp_distance_analysis_results_info = self._get_distance_analysis_results(tmp_distance_analysis_id)
            tmp_dist_analysis_results_id, tmp_pymol_session, tmp_rmsd, tmp_aligned_aa = tmp_distance_analysis_results_info

            index = []
            prot_1_chains = []
            prot_1_position = []
            prot_1_residue = []
            prot_2_chains = []
            prot_2_position = []
            prot_2_residue = []
            distances = []
            for tmp_distance_data in self._get_distance_analysis_result_data(tmp_dist_analysis_results_id):
                tmp_progress += 0.5
                the_progress_signal.emit_signal(tmp_message, min(round(tmp_progress, 2), 90))
                index.append(tmp_distance_data[0])
                prot_1_chains.append(tmp_distance_data[1])
                prot_1_position.append(tmp_distance_data[2])
                prot_1_residue.append(tmp_distance_data[3])
                prot_2_chains.append(tmp_distance_data[4])
                prot_2_position.append(tmp_distance_data[5])
                prot_2_residue.append(tmp_distance_data[6])
                distances.append(tmp_distance_data[7])

            tmp_distance_data_records = {
                pyssa_keys.ARRAY_DISTANCE_INDEX: np.array(index),
                pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN: np.array(prot_1_chains),
                pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION: np.array(prot_1_position),
                pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI: np.array(prot_1_residue),
                pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN: np.array(prot_2_chains),
                pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION: np.array(prot_2_position),
                pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI: np.array(prot_2_residue),
                pyssa_keys.ARRAY_DISTANCE_DISTANCES: np.array(distances),
            }
            tmp_dist_analysis_results = results.DistanceAnalysisResults(tmp_distance_data_records,
                                                                        tmp_pymol_session,
                                                                        tmp_rmsd,
                                                                        tmp_aligned_aa)
            tmp_protein_pair.distance_analysis.analysis_results = tmp_dist_analysis_results
            tmp_project.protein_pairs.append(tmp_protein_pair)
            i += 1

        return tmp_project

    def _get_project(self, a_project_name: str) -> tuple:
        sql = """SELECT id FROM Project WHERE name = ?"""
        self._cursor.execute(sql, (a_project_name,))
        return self._cursor.fetchall()[0]

    # </editor-fold>

    # <editor-fold desc="Select statements for protein objects">
    def _get_protein(self, the_project_id: int) -> list[tuple]:
        sql = """SELECT id, pymol_molecule_object, pymol_session FROM Protein WHERE project_id = ?"""
        self._cursor.execute(sql, (the_project_id,))
        return self._cursor.fetchall()

    def _get_pymol_selection(self, the_protein_id: int) -> tuple:
        sql = """SELECT selection_string FROM PyMOLSelection WHERE protein_id = ?"""
        self._cursor.execute(sql, (the_protein_id,))
        return self._cursor.fetchall()[0]

    def _get_chain(self, the_protein_id: int) -> list[tuple]:
        sql = """SELECT id, chain_identifier, chain_type, chain_sequence FROM Chain WHERE protein_id = ?"""
        self._cursor.execute(sql, (the_protein_id,))
        return self._cursor.fetchall()

    def _get_pymol_parameter(self, the_chain_id: int) -> tuple:
        sql = """SELECT color, representation FROM PyMOLParameter WHERE chain_id = ?"""
        self._cursor.execute(sql, (the_chain_id,))
        return self._cursor.fetchall()[0]

    def _get_protein_name_by_id(self, a_protein_id) -> tuple:
        sql = """SELECT pymol_molecule_object FROM Protein WHERE id = ?"""
        self._cursor.execute(sql, (a_protein_id,))
        return self._cursor.fetchall()[0]

    def get_pdb_atoms_of_protein(self, the_protein_id: int) -> list[tuple]:
        sql = """   
            SELECT record_type, atom_number, atom_name, alternate_location_indicator, residue_name,
            chain_identifier, residue_sequence_number, code_for_insertions_of_residues, 
            x_coord, y_coord, z_coord, occupancy, temperature_factor, segment_identifier, element_symbol, 
            charge, protein_id 
            FROM PdbAtom 
            WHERE protein_id = ?"""
        # tmp_params = (
        #     a_pdb_atom_dict["record_type"],
        #     a_pdb_atom_dict["atom_number"],
        #     a_pdb_atom_dict["atom_name"],
        #     a_pdb_atom_dict["alternate_location_indicator"],
        #     a_pdb_atom_dict["residue_name"],
        #     a_pdb_atom_dict["chain_identifier"],
        #     a_pdb_atom_dict["residue_sequence_number"],
        #     a_pdb_atom_dict["code_for_insertions_of_residues"],
        #     a_pdb_atom_dict["x_coord"],
        #     a_pdb_atom_dict["y_coord"],
        #     a_pdb_atom_dict["z_coord"],
        #     a_pdb_atom_dict["occupancy"],
        #     a_pdb_atom_dict["temperature_factor"],
        #     a_pdb_atom_dict["segment_identifier"],
        #     a_pdb_atom_dict["element_symbol"],
        #     a_pdb_atom_dict["charge"],
        #     the_protein_id
        # )
        self._cursor.execute(sql, (the_protein_id,))
        return self._cursor.fetchall()

    # </editor-fold>

    # <editor-fold desc="Select statements for protein pair objects">
    def _get_protein_pair(self, the_project_id) -> list[tuple]:
        sql = """SELECT id, protein_1_id, protein_2_id, pymol_session, name FROM ProteinPair WHERE project_id = ?"""
        self._cursor.execute(sql, (the_project_id,))
        return self._cursor.fetchall()

    def _get_distance_analysis(self, the_protein_pair_id: int) -> tuple:
        sql = """SELECT id, name, cutoff, cycles, figure_size_x, figure_size_y FROM DistanceAnalysis WHERE protein_pair_id = ?"""
        self._cursor.execute(sql, (the_protein_pair_id,))
        return self._cursor.fetchall()[0]

    def _get_distance_analysis_results(self, the_distance_analysis_id: int) -> tuple:
        sql = """SELECT id, pymol_session, rmsd, aligned_aa FROM DistanceAnalysisResults WHERE distance_analysis_id = ?"""
        self._cursor.execute(sql, (the_distance_analysis_id,))
        return self._cursor.fetchall()[0]

    def _get_distance_analysis_result_data(self, the_distance_analysis_results_id) -> list[tuple]:
        sql = """   SELECT my_index, protein_1_chain, protein_1_position, protein_1_residue,
                    protein_2_chain, protein_2_position, protein_2_residue, distances
                    FROM DistanceAnalysisResultData WHERE distance_analysis_results_id = ?"""
        self._cursor.execute(sql, (the_distance_analysis_results_id,))
        return self._cursor.fetchall()

    def get_pymol_parameter_for_certain_protein_chain_in_protein_pair(self,
                                                                      a_protein_pair_id: int,
                                                                      a_protein_id: int,
                                                                      a_chain_letter: str,
                                                                      a_parameter_name: str) -> tuple:
        sql = """   
            SELECT parameter_value
            FROM PyMOLParameterProteinPair WHERE protein_id = ? and chain_letter = ? and protein_pair_id = ? and parameter_name = ?"""
        self._cursor.execute(sql, (a_protein_id, a_chain_letter, a_protein_pair_id, a_parameter_name))
        return self._cursor.fetchall()[0]

    # </editor-fold>

    # <editor-fold desc="Select statements for sequences">
    def _get_sequence(self) -> list[tuple]:
        sql = """SELECT seq_id, seq, name FROM SeqRecord WHERE project_id = ?"""
        self._cursor.execute(sql, (1,))  # fixme: this is not the best solution
        return self._cursor.fetchall()

    # </editor-fold>

    # <editor-fold desc="Update statements for protein objects">
    def update_protein_chain_color(self, a_chain_id: int, a_color: str):
        sql = """
            UPDATE PyMOLParameter 
            SET color = ?
            WHERE chain_id = ?
        """
        self._cursor.execute(sql, (a_color, a_chain_id))
        self._connection.commit()

    def update_protein_chain_representation(self, a_chain_id: int, a_representation: str):
        sql = """
            UPDATE PyMOLParameter 
            SET representation = ?
            WHERE chain_id = ?
        """
        self._cursor.execute(sql, (a_representation, a_chain_id))
        self._connection.commit()

    def update_protein_name(self, the_new_protein_name: str, the_old_protein_name: str, the_protein_id: int):
        sql = """
            UPDATE Protein 
            SET pymol_molecule_object = ?
            WHERE pymol_molecule_object = ? and id = ?
        """
        self._cursor.execute(sql, (the_new_protein_name, the_old_protein_name, the_protein_id))
        self._connection.commit()

    def update_protein_pdb_atom_data(self, the_protein_id: int, a_pdb_atom_dict_list: list[dict]):
        self._delete_pdb_atom(the_protein_id)

        for tmp_pdb_atom_dict in a_pdb_atom_dict_list:
            self._insert_pdb_atom(the_protein_id, tmp_pdb_atom_dict)

        # sql = """
        #     UPDATE PdbAtom
        #     SET
        #         record_type = ?,
        #         atom_number = ?,
        #         atom_name = ?,
        #         alternate_location_indicator = ?,
        #         residue_name = ?,
        #         chain_identifier = ?,
        #         residue_sequence_number = ?,
        #         code_for_insertions_of_residues = ?,
        #         x_coord = ?,
        #         y_coord = ?,
        #         z_coord = ?,
        #         occupancy = ?,
        #         temperature_factor = ?,
        #         segment_identifier = ?,
        #         element_symbol = ?,
        #         charge = ?
        #     WHERE
        #         protein_id = ?;
        #     """
        # tmp_params = (
        #     a_pdb_atom_dict["record_type"],
        #     a_pdb_atom_dict["atom_number"],
        #     a_pdb_atom_dict["atom_name"],
        #     a_pdb_atom_dict["alternate_location_indicator"],
        #     a_pdb_atom_dict["residue_name"],
        #     a_pdb_atom_dict["chain_identifier"],
        #     a_pdb_atom_dict["residue_sequence_number"],
        #     a_pdb_atom_dict["code_for_insertions_of_residues"],
        #     a_pdb_atom_dict["x_coord"],
        #     a_pdb_atom_dict["y_coord"],
        #     a_pdb_atom_dict["z_coord"],
        #     a_pdb_atom_dict["occupancy"],
        #     a_pdb_atom_dict["temperature_factor"],
        #     a_pdb_atom_dict["segment_identifier"],
        #     a_pdb_atom_dict["element_symbol"],
        #     a_pdb_atom_dict["charge"],
        #     the_protein_id
        # )
        # self._cursor.execute(sql, tmp_params)
        # self._connection.commit()

    # </editor-fold>

    def update_pymol_parameter_for_certain_protein_chain_in_protein_pair(self,
                                                                         a_protein_pair_id: int,
                                                                         a_protein_id: int,
                                                                         a_chain_letter: str,
                                                                         the_parameter_name: str,
                                                                         the_new_parameter_value: str):
        sql = """
                    UPDATE PyMOLParameterProteinPair 
                    SET parameter_value = ?
                    WHERE protein_id = ? and chain_letter = ? and protein_pair_id = ? and parameter_name = ?
                """
        self._cursor.execute(sql, (the_new_parameter_value, a_protein_id, a_chain_letter, a_protein_pair_id, the_parameter_name))
        self._connection.commit()
        # sql = """
        #     UPDATE PyMOLParameterProteinPair
        #     FROM PyMOLParameterProteinPair WHERE protein_id = ? and chain_letter = ? and protein_pair_id = ? and parameter_name = ?"""
        # self._cursor.execute(sql, (a_protein_id, a_chain_letter, a_protein_pair_id, a_parameter_name))
        # return self._cursor.fetchall()[0]

    def update_pymol_session_of_protein_pair(self, the_protein_pair_id: int, the_new_pymol_session: str):
        sql = """
                            UPDATE ProteinPair 
                            SET pymol_session = ?
                            WHERE id = ?
                        """
        self._cursor.execute(sql, (str(the_new_pymol_session), int(the_protein_pair_id)))
        self._connection.commit()

    def update_pymol_session_of_protein(self, the_protein_id: int, the_new_pymol_session: str):
        sql = """
                            UPDATE Protein 
                            SET pymol_session = ?
                            WHERE id = ?
                        """
        self._cursor.execute(sql, (str(the_new_pymol_session), int(the_protein_id)))
        self._connection.commit()

    def update_sequence_name(self, the_new_seq_name: str, the_old_seq_name: str, the_sequence: str):
        sql = """
            UPDATE SeqRecord 
            SET name = ?
            WHERE name = ? and seq = ?
        """
        self._cursor.execute(sql, (the_new_seq_name, the_old_seq_name, the_sequence))
        self._connection.commit()

    def update_project_name(self, the_new_project_name: str):
        sql = """
            UPDATE Project 
            SET name = ?
            WHERE id = ?
        """
        self._cursor.execute(sql, (the_new_project_name, 1))
        self._connection.commit()
