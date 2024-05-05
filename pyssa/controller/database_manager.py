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
"""Module for the database interface and manager."""
import logging
import pathlib
from typing import Any

import numpy as np
from PyQt5 import QtCore
from PyQt5 import QtSql
from Bio import SeqRecord

from pyssa.internal.data_structures import protein, project, protein_pair, structure_analysis, results, chain, sequence
from pyssa.internal.thread.async_pyssa import custom_signals
from pyssa.logging_pyssa import log_handlers
from pyssa.util import enums, pyssa_keys, pyssa_exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class PyssaDatabaseInterface:
    """Custom interface layer to the project database.

    Notes:
        Uses the QSqlDatabase class to connect to the database.
        This should ensure thread-safety in the multithreaded environment (using QThreads).
    """
    def __init__(self, database_path, a_connection_name):
        self.mutex = QtCore.QMutex()
        self.db = QtSql.QSqlDatabase.addDatabase("QSQLITE", a_connection_name)
        self.db.setDatabaseName(database_path)

    def connect(self) -> bool:
        self.mutex.lock()
        try:
            if not self.db.open():
                # Handle connection error
                logger.error("Error connecting to database!")
                return False
            return True
        finally:
            self.mutex.unlock()

    def disconnect(self):
        self.mutex.lock()
        try:
            self.db.close()
        finally:
            self.mutex.unlock()

    def execute_query(self, a_qsql_query: QtSql.QSqlQuery, params: tuple):
        """Executes a given QSqlQuery object.

        Notes:
            Be aware this function has side effects!
            The given QSqlQuery will carry the results of the query after this function finished!
        """
        self.mutex.lock()
        try:
            i = 0
            for tmp_param in params:
                a_qsql_query.bindValue(i, tmp_param)
                i += 1
            if not a_qsql_query.exec():
                # Handle query execution error
                logger.error("Error executing query:", a_qsql_query.lastError().text())
        except Exception as e:
            logger.error(e)
        finally:
            self.mutex.unlock()
        return a_qsql_query


class PyssaSqlQuery:
    """Contains all SQL queries needed for the PySSA application in form of static functions."""
    def __init__(self):
        pass

    @staticmethod
    def create_sql_query(the_db: QtSql.QSqlDatabase, a_sql_statement: enums.SQLQueryStatement) -> QtSql.QSqlQuery:
        """Creates the SQL query for a given SQL statement.

        Raises:
            ConnectionError: if the database is closed.
            ValueError: if the SQL statement preparation failed.
        """
        if not the_db.isOpen():
            raise ConnectionError("The database must be open to perform this operation!")
        tmp_sql_query = QtSql.QSqlQuery(the_db)
        if tmp_sql_query.prepare(a_sql_statement.value):
            return tmp_sql_query
        else:
            raise ValueError("The SQL statement could not be prepared!")

    @staticmethod
    def create_sql_query_from_raw_statement_string(the_db: QtSql.QSqlDatabase, a_sql_statement: str) -> QtSql.QSqlQuery:
        """Creates the SQL query for a given raw SQL statement string.

        Raises:
            ConnectionError: if the database is closed.
            ValueError: if the SQL statement preparation failed.
        """
        if not the_db.isOpen():
            raise ConnectionError("The database must be open to perform this operation!")
        tmp_sql_query = QtSql.QSqlQuery(the_db)
        if tmp_sql_query.prepare(a_sql_statement):
            return tmp_sql_query
        else:
            raise ValueError("The SQL statement could not be prepared!")


class DatabaseManager:
    """Provides high-level access to the project's database by using the PyssaDatabaseInterface.

    Notes:
        This class should always use the with statement!
    """
    _database_filepath: str
    _pyssa_database_interface: "PyssaDatabaseInterface"

    def __init__(self, a_database_filepath: str, a_connection_name: str = "default") -> None:
        self._database_filepath = str(a_database_filepath)
        self._connection_name = a_connection_name

    # <editor-fold desc="Additional magic methods">
    def __enter__(self):
        self._pyssa_database_interface = PyssaDatabaseInterface(self._database_filepath, self._connection_name)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        logger.info("Tearing down the database manager object.")
    # </editor-fold>

    # <editor-fold desc="Private methods">
    def _get_next_id_of_protein_table(self):
        # <editor-fold desc="Checks">
        if not self._pyssa_database_interface.db.isOpen():
            raise pyssa_exception.DatabaseIsClosedError()
        # </editor-fold>

        tmp_sql_query = PyssaSqlQuery.create_sql_query(
            self._pyssa_database_interface.db, enums.SQLQueryStatement.GET_LATEST_ID_OF_PROTEIN_TABLE
        )
        self._pyssa_database_interface.execute_query(
            tmp_sql_query,
            params=()
        )
        tmp_sql_query.next()
        latest_id_before_insert = tmp_sql_query.value(0)
        if latest_id_before_insert is None:
            raise pyssa_exception.IllegalReturnValueError()
        if latest_id_before_insert == "":
            return 0 + 1
        return latest_id_before_insert + 1

    # <editor-fold desc="Sequences">
    def _get_data_of_all_sequences(self, a_project_id) -> list[dict[str, Any]]:
        # <editor-fold desc="Checks">
        if not self._pyssa_database_interface.db.isOpen():
            raise pyssa_exception.DatabaseIsClosedError()
        if a_project_id is None:
            raise pyssa_exception.IllegalArgumentError("a_project_id", a_project_id)
        # </editor-fold>

        tmp_sql_query_sequences = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                                 enums.SQLQueryStatement.GET_SEQUENCES)
        self._pyssa_database_interface.execute_query(
            tmp_sql_query_sequences, params=(a_project_id,)
        )
        tmp_seq_record_values = []
        while tmp_sql_query_sequences.next():
            tmp_seq_record_values.append(
                {
                    "id": tmp_sql_query_sequences.value(0),
                    "seq": tmp_sql_query_sequences.value(1),
                    "name": tmp_sql_query_sequences.value(2)
                }
            )
        return tmp_seq_record_values

    def _get_all_sequences_as_objects(self, the_seq_record_values: list[dict[str, Any]]) -> list[SeqRecord.SeqRecord]:
        # <editor-fold desc="Checks">
        if the_seq_record_values is None:
            raise pyssa_exception.IllegalArgumentError("the_seq_record_values", the_seq_record_values)
        # </editor-fold>

        tmp_sequences: list = []
        for tmp_seq_record_key_value_pairs in the_seq_record_values:
            tmp_sequence = SeqRecord.SeqRecord(
                id=tmp_seq_record_key_value_pairs["id"],
                seq=tmp_seq_record_key_value_pairs["seq"],
                name=tmp_seq_record_key_value_pairs["name"]
            )
            tmp_sequences.append(tmp_sequence)
        return tmp_sequences
    # </editor-fold>

    # <editor-fold desc="Proteins">
    def _get_protein_name_by_id(self, a_protein_id):
        # <editor-fold desc="Checks">
        if not self._pyssa_database_interface.db.isOpen():
            raise pyssa_exception.DatabaseIsClosedError()
        if a_protein_id is None:
            raise pyssa_exception.IllegalArgumentError("a_protein_id", a_protein_id)
        # </editor-fold>

        tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                       enums.SQLQueryStatement.GET_PROTEIN_NAME_BY_ID)
        self._pyssa_database_interface.execute_query(
            tmp_sql_query, params=(a_protein_id,)
        )
        tmp_sql_query.next()
        if tmp_sql_query.value(0) is None:
            raise pyssa_exception.IllegalReturnValueError()
        return tmp_sql_query.value(0)

    def _get_protein_data_by_name(self, a_protein_name: str):
        # <editor-fold desc="Checks">
        if not self._pyssa_database_interface.db.isOpen():
            raise pyssa_exception.DatabaseIsClosedError()
        if a_protein_name is None:
            raise pyssa_exception.IllegalArgumentError("a_protein_name", a_protein_name)
        # </editor-fold>

        tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                       enums.SQLQueryStatement.GET_PROTEIN_BY_NAME)
        self._pyssa_database_interface.execute_query(
            tmp_sql_query, params=(a_protein_name,)
        )
        tmp_sql_query.next()
        return {
                "id": tmp_sql_query.value(0),
                "name": tmp_sql_query.value(1),
                "pymol_session": tmp_sql_query.value(2),
        }

    def _get_data_of_all_proteins(self, a_project_id) -> list[dict[str, Any]]:
        # <editor-fold desc="Checks">
        if not self._pyssa_database_interface.db.isOpen():
            raise pyssa_exception.DatabaseIsClosedError()
        if a_project_id is None:
            raise pyssa_exception.IllegalArgumentError("a_project_id", a_project_id)
        # </editor-fold>

        tmp_sql_query_protein = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                               enums.SQLQueryStatement.GET_PROTEINS)
        self._pyssa_database_interface.execute_query(
            tmp_sql_query_protein, params=(a_project_id,)
        )
        tmp_protein_values = []
        while tmp_sql_query_protein.next():
            tmp_protein_values.append(
                {
                    "id": tmp_sql_query_protein.value(0),
                    "name": tmp_sql_query_protein.value(1),
                    "pymol_session": tmp_sql_query_protein.value(2),
                }
            )
        return tmp_protein_values

    def _get_all_proteins_as_objects(self, the_protein_values: list[dict[str, Any]]) -> list["protein.Protein"]:
        # <editor-fold desc="Checks">
        if the_protein_values is None:
            raise pyssa_exception.IllegalArgumentError("the_protein_values", the_protein_values)
        # </editor-fold>

        tmp_proteins: list["protein.Protein"] = []
        for tmp_protein_key_value_pairs in the_protein_values:
            tmp_protein = protein.Protein(tmp_protein_key_value_pairs["name"])
            tmp_protein.set_id(tmp_protein_key_value_pairs["id"])
            tmp_protein.pymol_session = tmp_protein_key_value_pairs["pymol_session"]
            tmp_protein.chains = self._get_all_chains_of_a_protein_as_objects(
                self._get_data_of_all_chains_of_a_protein(tmp_protein_key_value_pairs["id"]),
                tmp_protein.get_molecule_object(),
                tmp_protein_key_value_pairs["id"]
            )
            tmp_protein.pymol_selection.selection_string = self._get_data_of_pymol_selection_of_a_protein(
                tmp_protein_key_value_pairs["id"]
            )
            tmp_proteins.append(tmp_protein)
        return tmp_proteins

    def _get_data_of_all_chains_of_a_protein(self, a_protein_id) -> list[dict[str, Any]]:
        # <editor-fold desc="Checks">
        if not self._pyssa_database_interface.db.isOpen():
            raise pyssa_exception.DatabaseIsClosedError()
        if a_protein_id is None:
            raise pyssa_exception.IllegalArgumentError("a_protein_id", a_protein_id)
        # </editor-fold>

        tmp_sql_query_chain = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                             enums.SQLQueryStatement.GET_CHAINS)
        self._pyssa_database_interface.execute_query(
            tmp_sql_query_chain, params=(a_protein_id,)
        )
        tmp_chain_values = []
        while tmp_sql_query_chain.next():
            tmp_chain_values.append(
                {
                    "id": tmp_sql_query_chain.value(0),
                    "chain_identifier": tmp_sql_query_chain.value(1),
                    "chain_type": tmp_sql_query_chain.value(2),
                    "chain_sequence": tmp_sql_query_chain.value(3)
                }
            )
        return tmp_chain_values

    def _get_all_chains_of_a_protein_as_objects(self,
                                                the_chain_values: list[dict[str, Any]],
                                                a_protein_name: str,
                                                a_protein_id) -> list["chain.Chain"]:
        # <editor-fold desc="Checks">
        if the_chain_values is None:
            raise pyssa_exception.IllegalArgumentError("the_chain_values", the_chain_values)
        if a_protein_name is None:
            raise pyssa_exception.IllegalArgumentError("a_protein_name", a_protein_name)
        if a_protein_id is None:
            raise pyssa_exception.IllegalArgumentError("a_protein_id", a_protein_id)
        # </editor-fold>

        tmp_chains: list = []
        for tmp_chain_key_value_pairs in the_chain_values:
            tmp_seq = sequence.Sequence(a_protein_name, tmp_chain_key_value_pairs["chain_sequence"])
            tmp_chain = chain.Chain(tmp_chain_key_value_pairs["chain_identifier"], tmp_seq, tmp_chain_key_value_pairs["chain_type"])
            tmp_chain.set_id(tmp_chain_key_value_pairs["id"])
            tmp_chain.db_protein_id = a_protein_id
            tmp_chain.pymol_parameters = self._get_data_of_all_pymol_parameters_of_a_chain(tmp_chain.get_id())
            tmp_chains.append(tmp_chain)
        return tmp_chains

    def _get_data_of_all_pymol_parameters_of_a_chain(self, a_chain_id) -> dict[str, Any]:
        """Gets the data of all pymol parameters of a given chain id.

        Notes:
            The return value dict can directly be used for the chain object.
        """
        # <editor-fold desc="Checks">
        if not self._pyssa_database_interface.db.isOpen():
            raise pyssa_exception.DatabaseIsClosedError()
        if a_chain_id is None:
            raise pyssa_exception.IllegalArgumentError("a_chain_id", a_chain_id)
        # </editor-fold>

        tmp_sql_query_pymol_parameters = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                                        enums.SQLQueryStatement.GET_PYMOL_PARAMETERS)
        self._pyssa_database_interface.execute_query(
            tmp_sql_query_pymol_parameters, params=(a_chain_id,)
        )
        tmp_pymol_parameter_values = []
        tmp_sql_query_pymol_parameters.next()
        return {
            enums.PymolParameterEnum.COLOR.value: tmp_sql_query_pymol_parameters.value(0),
            enums.PymolParameterEnum.REPRESENTATION.value: tmp_sql_query_pymol_parameters.value(1),
        }

    def _get_data_of_pymol_selection_of_a_protein(self, a_protein_id) -> dict[str, str]:
        # <editor-fold desc="Checks">
        if not self._pyssa_database_interface.db.isOpen():
            raise pyssa_exception.DatabaseIsClosedError()
        if a_protein_id is None:
            raise pyssa_exception.IllegalArgumentError("a_protein_id", a_protein_id)
        # </editor-fold>

        tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                       enums.SQLQueryStatement.GET_PYMOL_SELECTIONS)
        self._pyssa_database_interface.execute_query(
            tmp_sql_query, params=(a_protein_id,)
        )
        tmp_sql_query.next()
        return {"selection": tmp_sql_query.value(0)}
    # </editor-fold>

    # <editor-fold desc="Protein pairs">
    def _get_protein_pair_data_by_name(self, a_protein_pair_name: str):
        # <editor-fold desc="Checks">
        if not self._pyssa_database_interface.db.isOpen():
            raise pyssa_exception.DatabaseIsClosedError()
        if a_protein_pair_name is None:
            raise pyssa_exception.IllegalArgumentError("a_protein_pair_name", a_protein_pair_name)
        # </editor-fold>

        tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                       enums.SQLQueryStatement.GET_PROTEIN_PAIR_BY_NAME)
        self._pyssa_database_interface.execute_query(
            tmp_sql_query, params=(a_protein_pair_name,)
        )
        tmp_sql_query.next()
        return {
            "id": tmp_sql_query.value(0),
            "protein_1_id": tmp_sql_query.value(1),
            "protein_2_id": tmp_sql_query.value(2),
            "pymol_session": tmp_sql_query.value(3),
            "protein_pair_name": tmp_sql_query.value(4)
        }

    def _get_data_of_all_protein_pairs(self, a_project_id) -> list[dict[str, Any]]:
        # <editor-fold desc="Checks">
        if not self._pyssa_database_interface.db.isOpen():
            raise pyssa_exception.DatabaseIsClosedError()
        if a_project_id is None:
            raise pyssa_exception.IllegalArgumentError("a_project_id", a_project_id)
        # </editor-fold>

        tmp_sql_query_protein_pair = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                                    enums.SQLQueryStatement.GET_PROTEIN_PAIRS)
        self._pyssa_database_interface.execute_query(
            tmp_sql_query_protein_pair, params=(a_project_id,)
        )
        tmp_protein_pair_values = []
        while tmp_sql_query_protein_pair.next():
            tmp_protein_pair_values.append(
                {
                    "id": tmp_sql_query_protein_pair.value(0),
                    "protein_1_id": tmp_sql_query_protein_pair.value(1),
                    "protein_2_id": tmp_sql_query_protein_pair.value(2),
                    "pymol_session": tmp_sql_query_protein_pair.value(3),
                    "protein_pair_name": tmp_sql_query_protein_pair.value(4)
                }
            )
        return tmp_protein_pair_values

    def _get_all_protein_pairs_as_objects(self,
                                          the_protein_pair_values: list[dict[str, Any]],
                                          tmp_project: "project.Project",
                                          the_app_settings: "settings.Settings") -> list["protein_pair.ProteinPair"]:
        # <editor-fold desc="Checks">
        if the_protein_pair_values is None:
            raise pyssa_exception.IllegalArgumentError("the_protein_pair_values", the_protein_pair_values)
        if tmp_project is None:
            raise pyssa_exception.IllegalArgumentError("tmp_project", tmp_project)
        if the_app_settings is None:
            raise pyssa_exception.IllegalArgumentError("the_app_settings", the_app_settings)
        # </editor-fold>

        tmp_protein_pairs: list = []
        for tmp_protein_pair_key_value_pairs in the_protein_pair_values:
            tmp_protein_pair = protein_pair.ProteinPair(
                tmp_project.search_protein(self._get_protein_name_by_id(tmp_protein_pair_key_value_pairs["protein_1_id"])),
                tmp_project.search_protein(self._get_protein_name_by_id(tmp_protein_pair_key_value_pairs["protein_2_id"]))
            )
            tmp_protein_pair.set_id(tmp_protein_pair_key_value_pairs["id"])
            tmp_protein_pair.db_project_id = tmp_project.get_id()
            tmp_protein_pair.name = tmp_protein_pair_key_value_pairs["protein_pair_name"]
            tmp_protein_pair.pymol_session = tmp_protein_pair_key_value_pairs["pymol_session"]
            tmp_protein_pair.distance_analysis = self._get_distance_analysis_as_object(
                tmp_protein_pair_key_value_pairs["id"], the_app_settings
            )
            tmp_protein_pairs.append(tmp_protein_pair)
        return tmp_protein_pairs

    def _get_distance_analysis_as_object(self,
                                         a_protein_pair_id,
                                         the_app_settings) -> "structure_analysis.DistanceAnalysis":
        # <editor-fold desc="Checks">
        if not self._pyssa_database_interface.db.isOpen():
            raise pyssa_exception.DatabaseIsClosedError()
        if a_protein_pair_id is None:
            raise pyssa_exception.IllegalArgumentError("a_protein_pair_id", a_protein_pair_id)
        # </editor-fold>

        tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                       enums.SQLQueryStatement.GET_DISTANCE_ANALYSIS)
        self._pyssa_database_interface.execute_query(
            tmp_sql_query, params=(a_protein_pair_id,)
        )
        tmp_sql_query.next()
        # create distance analysis object
        tmp_distance_analysis_id = tmp_sql_query.value(0)
        tmp_distance_analysis = structure_analysis.DistanceAnalysis(the_app_settings)
        tmp_distance_analysis.name = tmp_sql_query.value(1)
        tmp_distance_analysis.cutoff = tmp_sql_query.value(2)
        tmp_distance_analysis.cycles = tmp_sql_query.value(3)
        tmp_distance_analysis.figure_size = (tmp_sql_query.value(4), tmp_sql_query.value(5))
        tmp_distance_analysis.analysis_results = self._get_distance_analysis_results_as_object(tmp_distance_analysis_id)
        return tmp_distance_analysis

    def _get_distance_analysis_results_as_object(self, a_distance_analysis_id) -> "results.DistanceAnalysisResults":
        # <editor-fold desc="Checks">
        if not self._pyssa_database_interface.db.isOpen():
            raise pyssa_exception.DatabaseIsClosedError()
        if a_distance_analysis_id is None:
            raise pyssa_exception.IllegalArgumentError("a_distance_analysis_id", a_distance_analysis_id)
        # </editor-fold>

        tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                       enums.SQLQueryStatement.GET_DISTANCE_ANALYSIS_RESULTS)
        self._pyssa_database_interface.execute_query(
            tmp_sql_query, params=(a_distance_analysis_id,)
        )
        tmp_sql_query.next()
        tmp_dist_analysis_results_id = tmp_sql_query.value(0)
        tmp_pymol_session = tmp_sql_query.value(1)
        tmp_rmsd = tmp_sql_query.value(2)
        tmp_aligned_aa = tmp_sql_query.value(3)

        tmp_distance_data_records = self._get_all_distance_records_of_a_result(tmp_dist_analysis_results_id)
        return results.DistanceAnalysisResults(tmp_distance_data_records,
                                               tmp_pymol_session,
                                               tmp_rmsd,
                                               tmp_aligned_aa)

    def _get_all_distance_records_of_a_result(self, a_distance_analysis_result_id) -> dict[str, np.ndarray]:
        """Gets all distance records of a specific distance analysis result.

        Notes:
            The return value dict can be directly used.
        """

        # <editor-fold desc="Checks">
        if not self._pyssa_database_interface.db.isOpen():
            raise pyssa_exception.DatabaseIsClosedError()
        if a_distance_analysis_result_id is None:
            raise pyssa_exception.IllegalArgumentError("a_distance_analysis_result_id", a_distance_analysis_result_id)
        # </editor-fold>

        tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                       enums.SQLQueryStatement.GET_DISTANCE_ANALYSIS_RESULT_DATA)
        self._pyssa_database_interface.execute_query(
            tmp_sql_query, params=(a_distance_analysis_result_id,)
        )

        index = []
        prot_1_chains = []
        prot_1_position = []
        prot_1_residue = []
        prot_2_chains = []
        prot_2_position = []
        prot_2_residue = []
        distances = []
        while tmp_sql_query.next():
            index.append(tmp_sql_query.value(0))
            prot_1_chains.append(tmp_sql_query.value(1))
            prot_1_position.append(tmp_sql_query.value(2))
            prot_1_residue.append(tmp_sql_query.value(3))
            prot_2_chains.append(tmp_sql_query.value(4))
            prot_2_position.append(tmp_sql_query.value(5))
            prot_2_residue.append(tmp_sql_query.value(6))
            distances.append(tmp_sql_query.value(7))

        return {
            pyssa_keys.ARRAY_DISTANCE_INDEX: np.array(index),
            pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN: np.array(prot_1_chains),
            pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION: np.array(prot_1_position),
            pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI: np.array(prot_1_residue),
            pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN: np.array(prot_2_chains),
            pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION: np.array(prot_2_position),
            pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI: np.array(prot_2_residue),
            pyssa_keys.ARRAY_DISTANCE_DISTANCES: np.array(distances),
        }

    # </editor-fold>
    # </editor-fold>

    # <editor-fold desc="Public methods">
    # <editor-fold desc="Util ">
    def set_database_filepath(self, a_new_database_filepath: str) -> None:
        """Sets a new database filepath."""
        self._database_filepath = a_new_database_filepath

    def set_application_settings(self, the_app_settings) -> None:
        """Sets the app settings."""
        self._application_settings = the_app_settings

    def build_new_database(self) -> None:
        """Builds a new database for a new project."""
        if self._pyssa_database_interface.connect():
            sql = """-- Project definition
                CREATE TABLE Project (
                    id INTEGER NOT NULL,
                    name TEXT NOT NULL,
                    os TEXT,
                    CONSTRAINT Project_PK PRIMARY KEY (id)
                );
            """
            tmp_sql_query = PyssaSqlQuery.create_sql_query_from_raw_statement_string(
                self._pyssa_database_interface.db, sql
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query,
                params=()
            )
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
            tmp_sql_query = PyssaSqlQuery.create_sql_query_from_raw_statement_string(
                self._pyssa_database_interface.db, sql
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query,
                params=()
            )
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
            tmp_sql_query = PyssaSqlQuery.create_sql_query_from_raw_statement_string(
                self._pyssa_database_interface.db, sql
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query,
                params=()
            )
            sql = """-- "Chain" definition
                CREATE TABLE "Chain" (
                id INTEGER NOT NULL,
                protein_id INTEGER,
                chain_identifier TEXT,
                chain_type TEXT, chain_sequence TEXT,
                CONSTRAINT Chain_PK PRIMARY KEY (id),
                CONSTRAINT Chain_Protein_FK FOREIGN KEY (protein_id) REFERENCES Protein(id));
            """
            tmp_sql_query = PyssaSqlQuery.create_sql_query_from_raw_statement_string(
                self._pyssa_database_interface.db, sql
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query,
                params=()
            )
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
            tmp_sql_query = PyssaSqlQuery.create_sql_query_from_raw_statement_string(
                self._pyssa_database_interface.db, sql
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query,
                params=()
            )
            sql = """-- PyMOLSelection definition
                CREATE TABLE PyMOLSelection (
                    id INTEGER NOT NULL,
                    selection_string TEXT,
                    protein_id INTEGER,
                    CONSTRAINT PyMOLSelection_PK PRIMARY KEY (id),
                    CONSTRAINT PyMOLSelection_Protein_FK FOREIGN KEY (protein_id) REFERENCES Protein(id)
                );
            """
            tmp_sql_query = PyssaSqlQuery.create_sql_query_from_raw_statement_string(
                self._pyssa_database_interface.db, sql
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query,
                params=()
            )
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
            tmp_sql_query = PyssaSqlQuery.create_sql_query_from_raw_statement_string(
                self._pyssa_database_interface.db, sql
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query,
                params=()
            )
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
            tmp_sql_query = PyssaSqlQuery.create_sql_query_from_raw_statement_string(
                self._pyssa_database_interface.db, sql
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query,
                params=()
            )
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
            tmp_sql_query = PyssaSqlQuery.create_sql_query_from_raw_statement_string(
                self._pyssa_database_interface.db, sql
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query,
                params=()
            )
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
            tmp_sql_query = PyssaSqlQuery.create_sql_query_from_raw_statement_string(
                self._pyssa_database_interface.db, sql
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query,
                params=()
            )
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
            tmp_sql_query = PyssaSqlQuery.create_sql_query_from_raw_statement_string(
                self._pyssa_database_interface.db, sql
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query,
                params=()
            )
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
            tmp_sql_query = PyssaSqlQuery.create_sql_query_from_raw_statement_string(
                self._pyssa_database_interface.db, sql
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query,
                params=()
            )
            self._pyssa_database_interface.disconnect()

    def get_last_id(self):
        if not self._pyssa_database_interface.db.isOpen():
            raise pyssa_exception.DatabaseIsClosedError()
        tmp_sql_query = PyssaSqlQuery.create_sql_query(
            self._pyssa_database_interface.db, enums.SQLQueryStatement.GET_LAST_INSERT_ROW_ID
        )
        self._pyssa_database_interface.execute_query(
            tmp_sql_query,
            params=()
        )
        tmp_sql_query.next()
        return tmp_sql_query.value(0)

    def get_next_id_of_protein_table(self):
        # Assuming your_table_name has an auto-incrementing primary key column named 'id'
        if self._pyssa_database_interface.connect():
            tmp_sql_query = PyssaSqlQuery.create_sql_query(
                self._pyssa_database_interface.db, enums.SQLQueryStatement.GET_LATEST_ID_OF_PROTEIN_TABLE
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query,
                params=()
            )
            tmp_sql_query.next()
            latest_id_before_insert = tmp_sql_query.value(0)
            self._pyssa_database_interface.disconnect()
            if latest_id_before_insert is None:
                raise pyssa_exception.IllegalReturnValueError()
            if latest_id_before_insert == "":
                return 0 + 1
            return latest_id_before_insert + 1

    def get_next_id_of_chain_table(self):
        # Assuming your_table_name has an auto-incrementing primary key column named 'id'
        if self._pyssa_database_interface.connect():
            tmp_sql_query = PyssaSqlQuery.create_sql_query(
                self._pyssa_database_interface.db, enums.SQLQueryStatement.GET_LATEST_ID_OF_CHAIN_TABLE
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query,
                params=()
            )
            tmp_sql_query.next()
            latest_id_before_insert = tmp_sql_query.value(0)
            self._pyssa_database_interface.disconnect()
            if latest_id_before_insert is None:
                raise pyssa_exception.IllegalReturnValueError()
            if latest_id_before_insert == "":
                return 0 + 1
            return latest_id_before_insert + 1

    # </editor-fold>

    # <editor-fold desc="Project">
    def get_project_as_object(
            self,
            a_project_name: str,
            a_workspace_path: pathlib.Path,
            the_app_settings,
            the_progress_signal: "custom_signals.ProgressSignal" = custom_signals.ProgressSignal()
    ) -> "project.Project":
        """Creates a project object based on the data from the project database."""
        tmp_project = project.Project(a_project_name, a_workspace_path)
        if self._pyssa_database_interface.connect():
            tmp_sql_query_project_id = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                                      enums.SQLQueryStatement.GET_PROJECT_ID)
            self._pyssa_database_interface.execute_query(
                tmp_sql_query_project_id, params=(a_project_name,)
            )
            tmp_sql_query_project_id.next()
            tmp_project_id = tmp_sql_query_project_id.value(0)
            tmp_project.set_id(tmp_project_id)

            tmp_progress = 10
            the_progress_signal.emit_signal(f"Loading sequence ...", tmp_progress)
            tmp_project.sequences = self._get_all_sequences_as_objects(
                self._get_data_of_all_sequences(tmp_project_id)
            )

            tmp_progress = 25
            the_progress_signal.emit_signal(f"Loading protein ...", tmp_progress)
            tmp_project.proteins = self._get_all_proteins_as_objects(
                self._get_data_of_all_proteins(tmp_project_id)
            )

            tmp_progress = 50
            the_progress_signal.emit_signal(f"Loading protein pair ...", tmp_progress)
            tmp_project.protein_pairs = self._get_all_protein_pairs_as_objects(
                self._get_data_of_all_protein_pairs(tmp_project_id),
                tmp_project,
                the_app_settings
            )

            tmp_progress = 60
            the_progress_signal.emit_signal(f"Loading project data finished.", tmp_progress)
            self._pyssa_database_interface.disconnect()
        return tmp_project

    def insert_new_project(self, a_project_name, an_os) -> int:
        """Writes a new empty project to the database."""
        if self._pyssa_database_interface.connect():
            tmp_sql_query_project_id = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                                      enums.SQLQueryStatement.INSERT_NEW_PROJECT)
            self._pyssa_database_interface.execute_query(
                tmp_sql_query_project_id, params=(a_project_name, an_os)
            )
            tmp_last_id = self.get_last_id()
            self._pyssa_database_interface.disconnect()
            return tmp_last_id

    def update_project_name(self, the_new_project_name: str):
        if self._pyssa_database_interface.connect():
            tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                           enums.SQLQueryStatement.UPDATE_PROJECT_NAME)
            self._pyssa_database_interface.execute_query(
                tmp_sql_query, params=(the_new_project_name, 1)
            )
            self._pyssa_database_interface.disconnect()
    # </editor-fold>

    # <editor-fold desc="Sequences">
    def insert_new_sequence(self, a_seq_record):
        # <editor-fold desc="Checks">
        if a_seq_record is None:
            raise pyssa_exception.IllegalArgumentError("a_seq_record", a_seq_record)
        # </editor-fold>
        if self._pyssa_database_interface.connect():
            tmp_sql_query_project_id = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                                      enums.SQLQueryStatement.INSERT_SEQUENCE)
            tmp_project_id = 1  # fixme: this could lead to problems in the future
            self._pyssa_database_interface.execute_query(
                tmp_sql_query_project_id, params=(a_seq_record.id, a_seq_record.seq, a_seq_record.name, tmp_project_id)
            )
            self._pyssa_database_interface.disconnect()

    def delete_existing_sequence(self, a_seq_record_name):
        if self._pyssa_database_interface.connect():
            tmp_sql_query_project_id = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                                      enums.SQLQueryStatement.DELETE_SEQUENCE)
            self._pyssa_database_interface.execute_query(
                tmp_sql_query_project_id, params=(a_seq_record_name,)
            )
            self._pyssa_database_interface.disconnect()

    def update_sequence_name(self, the_new_seq_name: str, the_old_seq_name: str, the_sequence: str):
        if self._pyssa_database_interface.connect():
            tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                           enums.SQLQueryStatement.UPDATE_SEQUENCE_NAME)
            self._pyssa_database_interface.execute_query(
                tmp_sql_query, params=(the_new_seq_name, the_old_seq_name, the_sequence)
            )
            self._pyssa_database_interface.disconnect()
    # </editor-fold>

    # <editor-fold desc="Proteins">
    def insert_new_protein(self, a_protein: "protein.Protein") -> int:
        """Writes a new protein to the database."""
        tmp_dict = a_protein.get_object_as_dict_for_database()
        tmp_pymol_molecule_object = tmp_dict[enums.DatabaseEnum.PROTEIN_NAME.value]
        tmp_pymol_session = tmp_dict[enums.DatabaseEnum.PROTEIN_PYMOL_SESSION.value]
        tmp_project_id = tmp_dict["project_id"]

        if self._pyssa_database_interface.connect():
            tmp_protein_id = self._get_next_id_of_protein_table()
            tmp_sql_query_protein = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                                      enums.SQLQueryStatement.INSERT_PROTEIN)
            self._pyssa_database_interface.execute_query(
                tmp_sql_query_protein, params=(tmp_pymol_molecule_object, tmp_pymol_session, tmp_project_id)
            )
            if len(a_protein.chains) == 0:
                logger.warning("There are no chains to be inserted into the database.")
            else:
                # chains will be inserted
                for tmp_chain in a_protein.chains:
                    tmp_sql_query_chain = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                                         enums.SQLQueryStatement.INSERT_CHAIN)
                    self._pyssa_database_interface.execute_query(
                        tmp_sql_query_chain,
                        params=(
                            tmp_protein_id,
                            tmp_chain.chain_letter,
                            tmp_chain.chain_type,
                            tmp_chain.chain_sequence.sequence
                        )
                    )
                    tmp_sql_query_pymol_parameter = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                                                   enums.SQLQueryStatement.INSERT_PYMOL_PARAMETER)
                    tmp_chain.set_id(self.get_last_id())
                    self._pyssa_database_interface.execute_query(
                        tmp_sql_query_pymol_parameter,
                        params=(
                            tmp_chain.pymol_parameters[enums.PymolParameterEnum.COLOR.value],
                            tmp_chain.pymol_parameters[enums.PymolParameterEnum.REPRESENTATION.value],
                            tmp_chain.get_id()
                        )
                    )
                tmp_sql_query_pymol_selection = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                                               enums.SQLQueryStatement.INSERT_PYMOL_SELECTION)
                self._pyssa_database_interface.execute_query(
                    tmp_sql_query_pymol_selection,
                    params=(
                        a_protein.pymol_selection.selection_string,
                        tmp_protein_id
                    )
                )
            for tmp_pdb_atom_dict in a_protein.get_pdb_data():
                tmp_sql_query_pdb_atom = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                                        enums.SQLQueryStatement.INSERT_PDB_ATOM)
                self._pyssa_database_interface.execute_query(
                    tmp_sql_query_pdb_atom,
                    params=(
                        tmp_pdb_atom_dict["record_type"],
                        tmp_pdb_atom_dict["atom_number"],
                        tmp_pdb_atom_dict["atom_name"],
                        tmp_pdb_atom_dict["alternate_location_indicator"],
                        tmp_pdb_atom_dict["residue_name"],
                        tmp_pdb_atom_dict["chain_identifier"],
                        tmp_pdb_atom_dict["residue_sequence_number"],
                        tmp_pdb_atom_dict["code_for_insertions_of_residues"],
                        tmp_pdb_atom_dict["x_coord"],
                        tmp_pdb_atom_dict["y_coord"],
                        tmp_pdb_atom_dict["z_coord"],
                        tmp_pdb_atom_dict["occupancy"],
                        tmp_pdb_atom_dict["temperature_factor"],
                        tmp_pdb_atom_dict["segment_identifier"],
                        tmp_pdb_atom_dict["element_symbol"],
                        tmp_pdb_atom_dict["charge"],
                        tmp_protein_id
                    )
                )
            self._pyssa_database_interface.disconnect()
            return tmp_protein_id

    def delete_existing_protein(self, a_protein_id: int):
        if self._pyssa_database_interface.connect():
            tmp_sql_query_protein = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                                   enums.SQLQueryStatement.DELETE_PROTEIN)
            tmp_sql_query_chain = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                                 enums.SQLQueryStatement.DELETE_CHAINS)
            tmp_sql_query_pymol_parameter = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                                           enums.SQLQueryStatement.DELETE_PYMOL_PARAMETER)
            tmp_sql_query_pymol_selection = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                                           enums.SQLQueryStatement.DELETE_PYMOL_SELECTION)
            tmp_sql_query_pdb_atom = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                                    enums.SQLQueryStatement.DELETE_PDB_ATOM)
            self._pyssa_database_interface.execute_query(
                tmp_sql_query_pdb_atom, params=(a_protein_id,)
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query_pymol_selection, params=(a_protein_id,)
            )
            tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                           enums.SQLQueryStatement.GET_CHAINS)
            self._pyssa_database_interface.execute_query(
                tmp_sql_query, params=(a_protein_id,)
            )
            while tmp_sql_query.next():
                tmp_chain_id = tmp_sql_query.value(0)
                self._pyssa_database_interface.execute_query(
                    tmp_sql_query_pymol_parameter, params=(tmp_chain_id,)
                )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query_chain, params=(a_protein_id,)
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query_protein, params=(a_protein_id,)
            )
            self._pyssa_database_interface.disconnect()

    def get_protein_as_object(self, the_pymol_molecule_object: str) -> "protein.Protein":
        # create protein objects
        if self._pyssa_database_interface.connect():
            tmp_protein = self._get_all_proteins_as_objects(
                [self._get_protein_data_by_name(the_pymol_molecule_object)]
            )[0]
            #
            #
            # tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
            #                                                enums.SQLQueryStatement.GET_PROTEIN_BY_NAME)
            # self._pyssa_database_interface.execute_query(
            #     tmp_sql_query, params=(the_pymol_molecule_object,)
            # )
            # # for proteins
            # tmp_sql_query.next()
            # tmp_protein_id = tmp_sql_query.value(0)
            # tmp_protein_name = tmp_sql_query.value(1)
            # tmp_pymol_session = tmp_sql_query.value(2)
            #
            # tmp_protein = protein.Protein(tmp_protein_name)
            # tmp_protein.set_id(tmp_protein_id)
            # tmp_protein.pymol_session = tmp_pymol_session
            #
            # tmp_sql_query_protein = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
            #                                                enums.SQLQueryStatement.GET_CHAINS)
            # self._pyssa_database_interface.execute_query(
            #     tmp_sql_query_protein, params=(tmp_protein_id,)
            # )
            # # for chains
            # while tmp_sql_query_protein.next():
            #     tmp_chain_id = tmp_sql_query_protein.value(0)
            #     tmp_chain_identifier = tmp_sql_query_protein.value(1)
            #     tmp_chain_type = tmp_sql_query_protein.value(2)
            #     tmp_chain_sequence = tmp_sql_query_protein.value(3)
            #
            #     tmp_seq = sequence.Sequence(tmp_protein_name, tmp_chain_sequence)
            #     tmp_chain = chain.Chain(tmp_chain_identifier, tmp_seq, tmp_chain_type)
            #     tmp_chain.set_id(tmp_chain_id)
            #     tmp_chain.db_protein_id = tmp_protein_id
            #
            #     tmp_sql_query_chain = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
            #                                                    enums.SQLQueryStatement.GET_PYMOL_PARAMETERS)
            #     self._pyssa_database_interface.execute_query(
            #         tmp_sql_query_chain, params=(tmp_chain_id,)
            #     )
            #     # for pymol parameters
            #     while tmp_sql_query_chain.next():
            #         tmp_color = tmp_sql_query_chain.value(0)
            #         tmp_representation = tmp_sql_query_chain.value(0)
            #         tmp_chain.pymol_parameters = {
            #             enums.PymolParameterEnum.COLOR.value: tmp_color,
            #             enums.PymolParameterEnum.REPRESENTATION.value: tmp_representation,
            #         }
            #     tmp_protein.chains.append(tmp_chain)
            # tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
            #                                                enums.SQLQueryStatement.GET_PYMOL_SELECTIONS)
            # self._pyssa_database_interface.execute_query(
            #     tmp_sql_query, params=(tmp_protein_id,)
            # )
            # tmp_sql_query.next()
            # tmp_protein.pymol_selection.selection_string = tmp_sql_query.value(0)
            self._pyssa_database_interface.disconnect()
            return tmp_protein
        else:
            raise pyssa_exception.UnableToConnectToDatabaseError()

    def get_pdb_atoms_of_protein(self, the_protein_id: int) -> list[tuple]:
        if self._pyssa_database_interface.connect():
            tmp_sql_query: QtSql.QSqlQuery = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                           enums.SQLQueryStatement.GET_PDB_ATOMS_OF_PROTEIN)
            self._pyssa_database_interface.execute_query(
                tmp_sql_query, params=(the_protein_id,)
            )
            # for proteins
            tmp_structure_infos = []
            while tmp_sql_query.next():

                tmp_structure_infos.append(
                    (
                        tmp_sql_query.value(0),
                        tmp_sql_query.value(1),
                        tmp_sql_query.value(2),
                        tmp_sql_query.value(3),
                        tmp_sql_query.value(4),
                        tmp_sql_query.value(5),
                        tmp_sql_query.value(6),
                        tmp_sql_query.value(7),
                        tmp_sql_query.value(8),
                        tmp_sql_query.value(9),
                        tmp_sql_query.value(10),
                        tmp_sql_query.value(11),
                        tmp_sql_query.value(12),
                        tmp_sql_query.value(13),
                        tmp_sql_query.value(14),
                        tmp_sql_query.value(15),
                    )
                )
            self._pyssa_database_interface.disconnect()
            return tmp_structure_infos

    def update_protein_chain_color(self, a_chain_id: int, a_color: str):
        if self._pyssa_database_interface.connect():
            tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                           enums.SQLQueryStatement.UPDATE_PROTEIN_CHAIN_COLOR)
            self._pyssa_database_interface.execute_query(
                tmp_sql_query, params=(a_color, a_chain_id)
            )
            self._pyssa_database_interface.disconnect()

    def update_protein_name(self, the_new_protein_name: str, the_old_protein_name: str, the_protein_id: int):
        if self._pyssa_database_interface.connect():
            tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                           enums.SQLQueryStatement.UPDATE_PROTEIN_NAME)
            self._pyssa_database_interface.execute_query(
                tmp_sql_query, params=(the_new_protein_name, the_old_protein_name, the_protein_id)
            )
            self._pyssa_database_interface.disconnect()

    def update_protein_pdb_atom_data(self, the_protein_id: int, a_pdb_atom_dict_list: list[dict]):
        if self._pyssa_database_interface.connect():
            tmp_sql_query_delete = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                                   enums.SQLQueryStatement.DELETE_PDB_ATOM)
            tmp_sql_query_insert = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                                   enums.SQLQueryStatement.INSERT_PDB_ATOM)
            self._pyssa_database_interface.execute_query(
                tmp_sql_query_delete, params=(the_protein_id,)
            )
            for tmp_pdb_atom_dict in a_pdb_atom_dict_list:
                self._pyssa_database_interface.execute_query(
                    tmp_sql_query_insert,
                    params=(
                        tmp_pdb_atom_dict["record_type"],
                        tmp_pdb_atom_dict["atom_number"],
                        tmp_pdb_atom_dict["atom_name"],
                        tmp_pdb_atom_dict["alternate_location_indicator"],
                        tmp_pdb_atom_dict["residue_name"],
                        tmp_pdb_atom_dict["chain_identifier"],
                        tmp_pdb_atom_dict["residue_sequence_number"],
                        tmp_pdb_atom_dict["code_for_insertions_of_residues"],
                        tmp_pdb_atom_dict["x_coord"],
                        tmp_pdb_atom_dict["y_coord"],
                        tmp_pdb_atom_dict["z_coord"],
                        tmp_pdb_atom_dict["occupancy"],
                        tmp_pdb_atom_dict["temperature_factor"],
                        tmp_pdb_atom_dict["segment_identifier"],
                        tmp_pdb_atom_dict["element_symbol"],
                        tmp_pdb_atom_dict["charge"],
                        the_protein_id
                    )
                )

    def update_pymol_session_of_protein(self, the_protein_id: int, the_new_pymol_session: str):
        if self._pyssa_database_interface.connect():
            tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                           enums.SQLQueryStatement.UPDATE_PYMOL_SESSION_OF_PROTEIN)
            self._pyssa_database_interface.execute_query(
                tmp_sql_query, params=(str(the_new_pymol_session), int(the_protein_id))
            )
            self._pyssa_database_interface.disconnect()
    # </editor-fold>

    # <editor-fold desc="Protein pairs">
    def insert_new_protein_pair(self, a_protein_pair: "protein_pair.ProteinPair") -> int:
        """Inserts a new protein pair in the project database.

        Note:
            Run this method after distance analysis finished!!!
        """
        if self._pyssa_database_interface.connect():
            tmp_sql_query_protein_pair = PyssaSqlQuery.create_sql_query(
                self._pyssa_database_interface.db, enums.SQLQueryStatement.INSERT_PROTEIN_PAIR
            )
            tmp_sql_query_pymol_parameter = PyssaSqlQuery.create_sql_query(
                self._pyssa_database_interface.db, enums.SQLQueryStatement.INSERT_PYMOL_PARAMETER_PROTEIN_PAIR
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query_protein_pair,
                params=(
                    a_protein_pair.protein_1.get_id(),
                    a_protein_pair.protein_2.get_id(),
                    a_protein_pair.pymol_session,
                    a_protein_pair.db_project_id,
                    a_protein_pair.name
                )
            )
            a_protein_pair.set_id(self.get_last_id())

            for tmp_chain in a_protein_pair.protein_1.chains:
                self._pyssa_database_interface.execute_query(
                    tmp_sql_query_pymol_parameter,
                    params=(
                        a_protein_pair.protein_1.get_id(),
                        tmp_chain.chain_letter,
                        enums.PymolParameterEnum.COLOR.value,
                        tmp_chain.pymol_parameters[enums.PymolParameterEnum.COLOR.value],
                        a_protein_pair.get_id(),
                    )
                )
                self._pyssa_database_interface.execute_query(
                    tmp_sql_query_pymol_parameter,
                    params=(
                        a_protein_pair.protein_1.get_id(),
                        tmp_chain.chain_letter,
                        enums.PymolParameterEnum.REPRESENTATION.value,
                        tmp_chain.pymol_parameters[enums.PymolParameterEnum.REPRESENTATION.value],
                        a_protein_pair.get_id(),
                    )
                )

            for tmp_chain in a_protein_pair.protein_2.chains:
                self._pyssa_database_interface.execute_query(
                    tmp_sql_query_pymol_parameter,
                    params=(
                        a_protein_pair.protein_2.get_id(),
                        tmp_chain.chain_letter,
                        enums.PymolParameterEnum.COLOR.value,
                        tmp_chain.pymol_parameters[enums.PymolParameterEnum.COLOR.value],
                        a_protein_pair.get_id(),
                    )
                )
                self._pyssa_database_interface.execute_query(
                    tmp_sql_query_pymol_parameter,
                    params=(
                        a_protein_pair.protein_2.get_id(),
                        tmp_chain.chain_letter,
                        enums.PymolParameterEnum.REPRESENTATION.value,
                        tmp_chain.pymol_parameters[enums.PymolParameterEnum.REPRESENTATION.value],
                        a_protein_pair.get_id(),
                    )
                )

            tmp_sql_query_distance_analysis = PyssaSqlQuery.create_sql_query(
                self._pyssa_database_interface.db, enums.SQLQueryStatement.INSERT_DISTANCE_ANALYSIS
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query_distance_analysis,
                params=(
                    a_protein_pair.distance_analysis.name,
                    a_protein_pair.distance_analysis.cutoff,
                    a_protein_pair.distance_analysis.cycles,
                    a_protein_pair.get_id(),
                    a_protein_pair.distance_analysis.figure_size[0],
                    a_protein_pair.distance_analysis.figure_size[1],
                )
            )
            tmp_distance_analysis_id = self.get_last_id()
            tmp_sql_query_distance_analysis_results = PyssaSqlQuery.create_sql_query(
                self._pyssa_database_interface.db, enums.SQLQueryStatement.INSERT_DISTANCE_ANALYSIS_RESULTS
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query_distance_analysis_results,
                params=(
                    a_protein_pair.distance_analysis.analysis_results.pymol_session,
                    a_protein_pair.distance_analysis.analysis_results.rmsd,
                    a_protein_pair.distance_analysis.analysis_results.aligned_aa,
                    tmp_distance_analysis_id
                )
            )
            tmp_distance_analysis_results_id = self.get_last_id()
            tmp_distance_data = a_protein_pair.distance_analysis.analysis_results.distance_data
            index = list(tmp_distance_data[pyssa_keys.ARRAY_DISTANCE_INDEX])
            prot_1_chains = list(tmp_distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN])
            prot_1_position = list(tmp_distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION])
            prot_1_residue = list(tmp_distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI])
            prot_2_chains = list(tmp_distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN])
            prot_2_position = list(tmp_distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION])
            prot_2_residue = list(tmp_distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI])
            distances = list(tmp_distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES])
            for i in range(len(index)):
                tmp_params = (
                    int(index[i]),
                    str(prot_1_chains[i]),
                    int(prot_1_position[i]),
                    str(prot_1_residue[i]),
                    str(prot_2_chains[i]),
                    int(prot_2_position[i]),
                    str(prot_2_residue[i]),
                    float(distances[i]),
                    tmp_distance_analysis_results_id,
                )
                tmp_sql_query_insert_distance_data_records = PyssaSqlQuery.create_sql_query(
                    self._pyssa_database_interface.db, enums.SQLQueryStatement.INSERT_DISTANCE_DATA_RECORDS
                )
                self._pyssa_database_interface.execute_query(
                    tmp_sql_query_insert_distance_data_records,
                    params=tmp_params
                )
                i += 1
            self._pyssa_database_interface.disconnect()
            return a_protein_pair.get_id()



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

    def delete_existing_protein_pair(self, a_protein_pair_id: int) -> None:
        if self._pyssa_database_interface.connect():
            tmp_sql_query_protein_pair = PyssaSqlQuery.create_sql_query(
                self._pyssa_database_interface.db, enums.SQLQueryStatement.DELETE_PROTEIN_PAIR
            )
            tmp_sql_query_pymol_parameter = PyssaSqlQuery.create_sql_query(
                self._pyssa_database_interface.db, enums.SQLQueryStatement.DELETE_PYMOL_PARAMETERS_PROTEIN_PAIR
            )
            tmp_sql_query_distance_analysis = PyssaSqlQuery.create_sql_query(
                self._pyssa_database_interface.db, enums.SQLQueryStatement.DELETE_DISTANCE_ANALYSIS
            )
            tmp_sql_query_distance_analysis_results = PyssaSqlQuery.create_sql_query(
                self._pyssa_database_interface.db, enums.SQLQueryStatement.DELETE_DISTANCE_ANALYSIS_RESULTS
            )
            tmp_sql_query_insert_distance_data_records = PyssaSqlQuery.create_sql_query(
                self._pyssa_database_interface.db, enums.SQLQueryStatement.DELETE_DISTANCE_ANALYSIS_RESULT_DATA
            )

            # Get ids of distance analysis related objects
            tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                           enums.SQLQueryStatement.GET_DISTANCE_ANALYSIS)
            self._pyssa_database_interface.execute_query(
                tmp_sql_query, params=(a_protein_pair_id,)
            )
            tmp_sql_query.next()
            tmp_distance_analysis_id = tmp_sql_query.value(0)
            tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                           enums.SQLQueryStatement.GET_DISTANCE_ANALYSIS_RESULTS)
            self._pyssa_database_interface.execute_query(
                tmp_sql_query, params=(tmp_distance_analysis_id,)
            )
            tmp_sql_query.next()
            tmp_dist_analysis_results_id = tmp_sql_query.value(0)

            # delete statements
            self._pyssa_database_interface.execute_query(
                tmp_sql_query_insert_distance_data_records, params=(tmp_dist_analysis_results_id,)
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query_distance_analysis_results, params=(tmp_distance_analysis_id,)
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query_distance_analysis, params=(a_protein_pair_id,)
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query_pymol_parameter, params=(a_protein_pair_id,)
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query_protein_pair, params=(a_protein_pair_id,)
            )

    def get_protein_pair_as_object(self, a_name: str, a_project: "project.Project", the_app_settings):
        # create protein objects
        if self._pyssa_database_interface.connect():
            tmp_protein_pair = self._get_all_protein_pairs_as_objects(
                [self._get_protein_pair_data_by_name(a_name)],
                a_project,
                the_app_settings
            )[0]
            self._pyssa_database_interface.disconnect()
            return tmp_protein_pair
        else:
            raise pyssa_exception.UnableToConnectToDatabaseError()

        # create protein pair objects
        # if self._pyssa_database_interface.connect():
        #     tmp_sql_query_protein_pair = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
        #                                                    enums.SQLQueryStatement.GET_PROTEIN_PAIR_BY_NAME)
        #     self._pyssa_database_interface.execute_query(
        #         tmp_sql_query_protein_pair, params=(a_name,)
        #     )
        #     # for protein pairs
        #     while tmp_sql_query_protein_pair.next():
        #         tmp_protein_pair_id = tmp_sql_query_protein_pair.value(0)
        #         tmp_protein_1_id = tmp_sql_query_protein_pair.value(1)
        #         tmp_protein_2_id = tmp_sql_query_protein_pair.value(2)
        #         tmp_pymol_session = tmp_sql_query_protein_pair.value(3)
        #         tmp_pp_name = tmp_sql_query_protein_pair.value(4)
        #
        #         tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
        #                                                        enums.SQLQueryStatement.GET_PROTEIN_NAME_BY_ID)
        #         self._pyssa_database_interface.execute_query(
        #             tmp_sql_query, params=(tmp_protein_1_id,)
        #         )
        #         tmp_sql_query.next()
        #         tmp_protein_1_name = tmp_sql_query.value(0)
        #         self._pyssa_database_interface.execute_query(
        #             tmp_sql_query, params=(tmp_protein_2_id,)
        #         )
        #         tmp_sql_query.next()
        #         tmp_protein_2_name = tmp_sql_query.value(0)
        #
        #         tmp_protein_pair = protein_pair.ProteinPair(a_project.search_protein(tmp_protein_1_name),
        #                                                     a_project.search_protein(tmp_protein_2_name))
        #         tmp_protein_pair.set_id(tmp_protein_pair_id)
        #         tmp_protein_pair.db_project_id = a_project.get_id()
        #         tmp_protein_pair.name = tmp_pp_name
        #         tmp_protein_pair.pymol_session = tmp_pymol_session
        #
        #         tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
        #                                                        enums.SQLQueryStatement.GET_DISTANCE_ANALYSIS)
        #         self._pyssa_database_interface.execute_query(
        #             tmp_sql_query, params=(tmp_protein_pair_id,)
        #         )
        #         tmp_sql_query.next()
        #         tmp_distance_analysis_id = tmp_sql_query.value(0)
        #         tmp_name = tmp_sql_query.value(1)
        #         tmp_cutoff = tmp_sql_query.value(2)
        #         tmp_cycles = tmp_sql_query.value(3)
        #         tmp_figure_size_x = tmp_sql_query.value(4)
        #         tmp_figure_size_y = tmp_sql_query.value(5)
        #
        #         # create distance analysis object
        #         tmp_distance_analysis = structure_analysis.DistanceAnalysis(the_app_settings)
        #         tmp_distance_analysis.name = tmp_name
        #         tmp_distance_analysis.cutoff = tmp_cutoff
        #         tmp_distance_analysis.cycles = tmp_cycles
        #         tmp_distance_analysis.figure_size = (tmp_figure_size_x, tmp_figure_size_y)
        #
        #         tmp_protein_pair.distance_analysis = tmp_distance_analysis
        #
        #         # create distance analysis results object
        #         tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
        #                                                        enums.SQLQueryStatement.GET_DISTANCE_ANALYSIS_RESULTS)
        #         self._pyssa_database_interface.execute_query(
        #             tmp_sql_query, params=(tmp_distance_analysis_id,)
        #         )
        #         tmp_sql_query.next()
        #         tmp_dist_analysis_results_id = tmp_sql_query.value(0)
        #         tmp_pymol_session = tmp_sql_query.value(1)
        #         tmp_rmsd = tmp_sql_query.value(2)
        #         tmp_aligned_aa = tmp_sql_query.value(3)
        #
        #         index = []
        #         prot_1_chains = []
        #         prot_1_position = []
        #         prot_1_residue = []
        #         prot_2_chains = []
        #         prot_2_position = []
        #         prot_2_residue = []
        #         distances = []
        #
        #         tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
        #                                                        enums.SQLQueryStatement.GET_DISTANCE_ANALYSIS_RESULT_DATA)
        #         self._pyssa_database_interface.execute_query(
        #             tmp_sql_query, params=(tmp_dist_analysis_results_id,)
        #         )
        #     while tmp_sql_query.next():
        #         index.append(tmp_sql_query.value(0))
        #         prot_1_chains.append(tmp_sql_query.value(1))
        #         prot_1_position.append(tmp_sql_query.value(2))
        #         prot_1_residue.append(tmp_sql_query.value(3))
        #         prot_2_chains.append(tmp_sql_query.value(4))
        #         prot_2_position.append(tmp_sql_query.value(5))
        #         prot_2_residue.append(tmp_sql_query.value(6))
        #         distances.append(tmp_sql_query.value(7))
        #
        #     tmp_distance_data_records = {
        #         pyssa_keys.ARRAY_DISTANCE_INDEX: np.array(index),
        #         pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN: np.array(prot_1_chains),
        #         pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION: np.array(prot_1_position),
        #         pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI: np.array(prot_1_residue),
        #         pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN: np.array(prot_2_chains),
        #         pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION: np.array(prot_2_position),
        #         pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI: np.array(prot_2_residue),
        #         pyssa_keys.ARRAY_DISTANCE_DISTANCES: np.array(distances),
        #     }
        #     tmp_dist_analysis_results = results.DistanceAnalysisResults(tmp_distance_data_records,
        #                                                                 tmp_pymol_session,
        #                                                                 tmp_rmsd,
        #                                                                 tmp_aligned_aa)
        #     tmp_protein_pair.distance_analysis.analysis_results = tmp_dist_analysis_results
        #     self._pyssa_database_interface.disconnect()
        #     return tmp_protein_pair

    def get_pymol_parameter_for_certain_protein_chain_in_protein_pair(
            self,
            a_protein_pair_id: int,
            a_protein_id: int,
            a_chain_letter: str,
            a_parameter_name: str) -> tuple:
        if self._pyssa_database_interface.connect():
            tmp_sql_query = PyssaSqlQuery.create_sql_query(
                self._pyssa_database_interface.db,
                enums.SQLQueryStatement.GET_PYMOL_PARAMETER_FOR_PROTEIN_CHAIN_IN_PROTEIN_PAIR
            )
            self._pyssa_database_interface.execute_query(
                tmp_sql_query, params=(
                    a_protein_id, a_chain_letter, a_protein_pair_id, a_parameter_name
                )
            )
            # for protein pairs
            tmp_sql_query.next()
            tmp_query_result = tmp_sql_query.value(0)
            self._pyssa_database_interface.disconnect()
            return tmp_query_result

    def update_pymol_session_of_protein_pair(self, the_protein_pair_id: int, the_new_pymol_session: str):
        if self._pyssa_database_interface.connect():
            tmp_sql_query = PyssaSqlQuery.create_sql_query(self._pyssa_database_interface.db,
                                                           enums.SQLQueryStatement.UPDATE_PYMOL_SESSION_OF_PROTEIN_PAIR)
            self._pyssa_database_interface.execute_query(
                tmp_sql_query, params=(str(the_new_pymol_session), int(the_protein_pair_id))
            )
            self._pyssa_database_interface.disconnect()


    # </editor-fold>
    # </editor-fold>




# class DatabaseManager:
#     """Manages the project database."""
#     _database_filepath: str
#     _connection: sqlite3.Connection
#     _cursor: sqlite3.Cursor
#     _application_settings: "settings.Settings"
#
#     def __init__(self, the_database_filepath: str) -> None:
#         self._database_filepath = the_database_filepath
#         self._connection = None
#
#     # <editor-fold desc="Util methods">
#     def get_last_id(self):
#         sql = """SELECT last_insert_rowid();"""
#         self._cursor.execute(sql)
#         tmp_id, = self._cursor.fetchall()[0]
#         return tmp_id
#
#     def get_database_filepath(self):
#         return self._database_filepath
#
#     def set_application_settings(self, the_app_settings):
#         self._application_settings = the_app_settings
#
#     def get_latest_id_of_protein_table(self):
#         self.open_project_database()
#         # Assuming your_table_name has an auto-incrementing primary key column named 'id'
#         self._cursor.execute("SELECT MAX(id) FROM Protein")
#         latest_id_before_insert = self._cursor.fetchone()[0]
#         self.close_project_database()
#         return latest_id_before_insert if latest_id_before_insert is not None else 0
#
#     def get_latest_id_of_a_specific_table(self, a_table_name):
#         self.open_project_database()
#         # Assuming your_table_name has an auto-incrementing primary key column named 'id'
#         self._cursor.execute(f"SELECT MAX(id) FROM {a_table_name}")
#         latest_id_before_insert = self._cursor.fetchone()[0]
#         return latest_id_before_insert if latest_id_before_insert is not None else 0
#
#     # </editor-fold>
#
#     # <editor-fold desc="Base methods">
#     def build_new_database(self) -> None:
#         """Builds a new database for a new project."""
#         tmp_connection = sqlite3.connect(self._database_filepath)
#         tmp_cursor = tmp_connection.cursor()
#         sql = """-- Project definition
#             CREATE TABLE Project (
#                 id INTEGER NOT NULL,
#                 name TEXT NOT NULL,
#                 os TEXT,
#                 CONSTRAINT Project_PK PRIMARY KEY (id)
#             );
#         """
#         tmp_cursor.execute(sql)
#         sql = """-- SeqRecord definition
#             CREATE TABLE SeqRecord (
#                 id INTEGER NOT NULL,
#                 seq_id TEXT,
#                 seq TEXT,
#                 name TEXT,
#                 project_id INTEGER,
#                 CONSTRAINT SeqRecord_PK PRIMARY KEY (id),
#                 CONSTRAINT SeqRecord_Project_FK FOREIGN KEY (project_id) REFERENCES Project(id)
#             );
#         """
#         tmp_cursor.execute(sql)
#         sql = """-- Protein definition
#             CREATE TABLE Protein (
#                 id INTEGER NOT NULL,
#                 pymol_molecule_object TEXT,
#                 pymol_session TEXT,
#                 project_id INTEGER,
#                 pdb_id INTEGER,
#                 CONSTRAINT Protein_PK PRIMARY KEY (id),
#                 CONSTRAINT Protein_Project_FK FOREIGN KEY (project_id) REFERENCES Project(id)
#             );
#         """
#         tmp_cursor.execute(sql)
#         sql = """-- "Chain" definition
#             CREATE TABLE "Chain" (
#             id INTEGER NOT NULL,
#             protein_id INTEGER,
#             chain_identifier TEXT,
#             chain_type TEXT, chain_sequence TEXT,
#             CONSTRAINT Chain_PK PRIMARY KEY (id),
#             CONSTRAINT Chain_Protein_FK FOREIGN KEY (protein_id) REFERENCES Protein(id));
#         """
#         tmp_cursor.execute(sql)
#         sql = """-- PyMOLParameter definition
#             CREATE TABLE PyMOLParameter (
#                 id INTEGER NOT NULL,
#                 color TEXT,
#                 representation TEXT,
#                 chain_id INTEGER,
#                 CONSTRAINT PyMOLParameter_PK PRIMARY KEY (id),
#                 CONSTRAINT PyMOLParameter_Chain_FK FOREIGN KEY (chain_id) REFERENCES "Chain"(id)
#             );
#         """
#         tmp_cursor.execute(sql)
#         sql = """-- PyMOLSelection definition
#             CREATE TABLE PyMOLSelection (
#                 id INTEGER NOT NULL,
#                 selection_string TEXT,
#                 protein_id INTEGER,
#                 CONSTRAINT PyMOLSelection_PK PRIMARY KEY (id),
#                 CONSTRAINT PyMOLSelection_Protein_FK FOREIGN KEY (protein_id) REFERENCES Protein(id)
#             );
#         """
#         tmp_cursor.execute(sql)
#         sql = """-- PdbAtom definition
#             CREATE TABLE PdbAtom (
#                 id INTEGER NOT NULL,
#                 record_type TEXT(6),
#                 atom_number INTEGER,
#                 atom_name TEXT(4),
#                 alternate_location_indicator TEXT(1),
#                 residue_name TEXT(3),
#                 chain_identifier TEXT(1),
#                 residue_sequence_number INTEGER,
#                 code_for_insertions_of_residues TEXT(1),
#                 x_coord REAL,
#                 y_coord REAL,
#                 z_coord REAL,
#                 occupancy REAL,
#                 temperature_factor REAL,
#                 segment_identifier TEXT(4),
#                 element_symbol TEXT(2),
#                 charge TEXT(2),
#                 protein_id INTEGER,
#                 CONSTRAINT PdbAtom_PK PRIMARY KEY (id),
#                 CONSTRAINT PdbAtom_Protein_FK FOREIGN KEY (protein_id) REFERENCES Protein(id)
#             );
#         """
#         tmp_cursor.execute(sql)
#         sql = """-- ProteinPair definition
#             CREATE TABLE ProteinPair (
#                 id INTEGER NOT NULL,
#                 protein_1_id INTEGER,
#                 protein_2_id INTEGER,
#                 pymol_session_filepath TEXT,
#                 pymol_session TEXT,
#                 project_id INTEGER, name TEXT,
#                 CONSTRAINT ProteinPair_PK PRIMARY KEY (id),
#                 CONSTRAINT ProteinPair_Protein_FK FOREIGN KEY (protein_1_id) REFERENCES Protein(id),
#                 CONSTRAINT ProteinPair_Protein_FK_1 FOREIGN KEY (protein_2_id) REFERENCES Protein(id),
#                 CONSTRAINT ProteinPair_Project_FK FOREIGN KEY (project_id) REFERENCES Project(id)
#             );
#         """
#         tmp_cursor.execute(sql)
#         sql = """-- PyMOLParameterProteinPair definition
#             CREATE TABLE PyMOLParameterProteinPair (
#                 id INTEGER NOT NULL,
#                 protein_id INTEGER,
#                 chain_letter TEXT,
#                 parameter_name TEXT,
#                 parameter_value TEXT,
#                 protein_pair_id INTEGER,
#                 CONSTRAINT PyMOLParameterProteinPair_PK PRIMARY KEY (id),
#                 CONSTRAINT PyMOLParameterProteinPair_ProteinPair_FK FOREIGN KEY (protein_pair_id) REFERENCES ProteinPair(id)
#             );
#         """
#         tmp_cursor.execute(sql)
#         sql = """-- DistanceAnalysis definition
#             CREATE TABLE DistanceAnalysis (
#                 id INTEGER NOT NULL,
#                 name TEXT,
#                 cutoff REAL,
#                 cycles INTEGER,
#                 protein_pair_id INTEGER, figure_size_x REAL, figure_size_y REAL,
#                 CONSTRAINT DistanceAnalysis_PK PRIMARY KEY (id),
#                 CONSTRAINT DistanceAnalysis_ProteinPair_FK FOREIGN KEY (protein_pair_id) REFERENCES ProteinPair(id)
#             );
#         """
#         tmp_cursor.execute(sql)
#         sql = """-- DistanceAnalysisResults definition
#             CREATE TABLE DistanceAnalysisResults (
#                 id INTEGER NOT NULL,
#                 pymol_session TEXT,
#                 rmsd REAL,
#                 aligned_aa TEXT,
#                 distance_analysis_id INTEGER,
#                 CONSTRAINT DistanceAnalysisResults_PK PRIMARY KEY (id),
#                 CONSTRAINT DistanceAnalysisResults_DistanceAnalysis_FK FOREIGN KEY (distance_analysis_id) REFERENCES DistanceAnalysis(id)
#             );
#         """
#         tmp_cursor.execute(sql)
#         sql = """-- DistanceAnalysisResultData definition
#             CREATE TABLE DistanceAnalysisResultData (
#                             id INTEGER NOT NULL,
#                             "my_index" INTEGER,
#                             protein_1_chain TEXT,
#                             protein_1_position INTEGER,
#                             protein_1_residue TEXT,
#                             protein_2_chain TEXT,
#                             protein_2_position INTEGER,
#                             protein_2_residue TEXT,
#                             distances REAL,
#                             distance_analysis_results_id INTEGER,
#                             CONSTRAINT DistanceAnalysisResultData_PK PRIMARY KEY (id),
#                             CONSTRAINT DistanceAnalysisResultData_DistanceAnalysisResults_FK FOREIGN KEY (distance_analysis_results_id) REFERENCES DistanceAnalysisResults(id)
#                         );
#         """
#         tmp_cursor.execute(sql)
#         tmp_connection.close()
#
#     def set_database_filepath(self, a_new_database_filepath: str) -> None:
#         """Sets a new database filepath."""
#         self._database_filepath = a_new_database_filepath
#
#     def open_project_database(self):
#         """Opens a project database"""
#         logger.debug(f"Filepath of the database_manager: {self._database_filepath}")
#         self._connection = sqlite3.connect(self._database_filepath)
#         self._cursor = self._connection.cursor()
#
#     def close_project_database(self):
#         """Closes a project database."""
#         self._cursor.close()
#         self._connection.close()
#
#     def __enter__(self):
#         self._connection = sqlite3.connect(self._database_filepath)
#         self._cursor = self._connection.cursor()
#         return self
#
#     def __exit__(self, exc_type, exc_value, traceback):
#         if self._connection:
#             self._connection.close()
#
#     # </editor-fold>
#
#     # <editor-fold desc="SeqRecord objects inserts">
#     def insert_new_sequence(self, a_seq_record: SeqRecord.SeqRecord):
#         self._insert_sequence(a_seq_record.id, a_seq_record.seq, a_seq_record.name)
#
#     def _insert_sequence(self, a_seq_id, a_seq, a_name):
#         sql = """   INSERT INTO SeqRecord(seq_id, seq, name, project_id)
#                     VALUES (?, ?, ?, ?)
#         """
#         self._cursor.execute(sql, (str(a_seq_id), str(a_seq), str(a_name), 1))  # fixme: not the best solution with the id=1
#         self._connection.commit()
#         return self.get_last_id()
#     # </editor-fold>
#
#     # <editor-fold desc="Delete statements for SeqRecord object">
#     def delete_existing_sequence(self, a_seq_record_name):
#         sql = """
#             DELETE FROM SeqRecord
#             WHERE name = ?
#         """
#         self._cursor.execute(sql, (a_seq_record_name,))
#         self._connection.commit()
#
#     # </editor-fold>
#
#     # <editor-fold desc="Protein object inserts">
#     def insert_new_protein(self, a_protein: "protein.Protein") -> int:
#         """Writes a new protein to the database."""
#         tmp_protein_id = self._insert_protein(a_protein.get_object_as_dict_for_database())
#         a_protein.set_id(tmp_protein_id)
#         if len(a_protein.chains) == 0:
#             logger.warning("There are no chains to be inserted into the database.")
#         else:
#             # chains will be inserted
#             for tmp_chain in a_protein.chains:
#                 tmp_chain_id = self._insert_chain(tmp_protein_id, tmp_chain)
#                 self._insert_pymol_parameter(tmp_chain_id, tmp_chain.pymol_parameters)
#         self._insert_pymol_selection(tmp_protein_id, a_protein.pymol_selection.selection_string)
#         for tmp_pdb_atom_dict in a_protein.get_pdb_data():
#             self._insert_pdb_atom(tmp_protein_id, tmp_pdb_atom_dict)
#         return tmp_protein_id
#
#     def _insert_protein(self, a_protein_obj_dict: dict) -> int:
#         sql = """   INSERT INTO Protein(pymol_molecule_object, pymol_session, project_id)
#                     VALUES (:pymol_molecule_object, :pymol_session, :project_id)
#         """
#         self._cursor.execute(sql, a_protein_obj_dict)
#         self._connection.commit()
#         return self.get_last_id()
#
#     def _insert_chain(self, the_protein_id: int, a_chain) -> int:
#         sql = """   INSERT INTO Chain(protein_id, chain_identifier, chain_type, chain_sequence)
#                     VALUES (?, ?, ?, ?)
#         """
#         tmp_params = (the_protein_id, a_chain.chain_letter, a_chain.chain_type, a_chain.chain_sequence.sequence)
#         self._cursor.execute(sql, tmp_params)
#         self._connection.commit()
#         return self.get_last_id()
#
#     def _insert_pymol_parameter(self, the_chain_id, a_pymol_parameter_dict: dict):
#         sql = """   INSERT INTO PyMOLParameter(color, representation, chain_id)
#                     VALUES (?, ?, ?)
#         """
#         tmp_params = (
#             a_pymol_parameter_dict[enums.PymolParameterEnum.COLOR.value],
#             a_pymol_parameter_dict[enums.PymolParameterEnum.REPRESENTATION.value],
#             the_chain_id
#         )
#         self._cursor.execute(sql, tmp_params)
#         self._connection.commit()
#
#     def _insert_pymol_selection(self, the_protein_id: int, a_selection_string: str):
#         sql = """   INSERT INTO PyMOLSelection(selection_string, protein_id)
#                     VALUES (?, ?)
#         """
#         tmp_params = (
#             a_selection_string,
#             the_protein_id,
#         )
#         self._cursor.execute(sql, tmp_params)
#         self._connection.commit()
#
#     def _insert_pdb_atom(self, the_protein_id: int, a_pdb_atom_dict: dict):
#         sql =
#         tmp_params = (
#             a_pdb_atom_dict["record_type"],"""   INSERT INTO PdbAtom(record_type, atom_number, atom_name, alternate_location_indicator, residue_name,
# #             chain_identifier, residue_sequence_number, code_for_insertions_of_residues,
# #             x_coord, y_coord, z_coord, occupancy, temperature_factor, segment_identifier, element_symbol,
# #             charge, protein_id)
# #                     VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
# #         """
#             a_pdb_atom_dict["atom_number"],
#             a_pdb_atom_dict["atom_name"],
#             a_pdb_atom_dict["alternate_location_indicator"],
#             a_pdb_atom_dict["residue_name"],
#             a_pdb_atom_dict["chain_identifier"],
#             a_pdb_atom_dict["residue_sequence_number"],
#             a_pdb_atom_dict["code_for_insertions_of_residues"],
#             a_pdb_atom_dict["x_coord"],
#             a_pdb_atom_dict["y_coord"],
#             a_pdb_atom_dict["z_coord"],
#             a_pdb_atom_dict["occupancy"],
#             a_pdb_atom_dict["temperature_factor"],
#             a_pdb_atom_dict["segment_identifier"],
#             a_pdb_atom_dict["element_symbol"],
#             a_pdb_atom_dict["charge"],
#             the_protein_id
#         )
#         self._cursor.execute(sql, tmp_params)
#         self._connection.commit()
#
#     # </editor-fold>
#
#     # <editor-fold desc="Delete statements for protein object">
#     def delete_existing_protein(self, a_protein_id: int):
#         self._delete_pdb_atom(a_protein_id)
#         self._delete_pymol_selection(a_protein_id)
#         for tmp_chain_info in self._get_chain(a_protein_id):
#             self._delete_pymol_parameter(tmp_chain_info[0])
#         self._delete_chains(a_protein_id)
#         self._delete_protein(a_protein_id)
#
#     def _delete_protein(self, the_protein_id: int):
#         sql = """
#             DELETE FROM Protein
#             WHERE id = ?
#         """
#         self._cursor.execute(sql, (the_protein_id,))
#         self._connection.commit()
#
#     def _delete_chains(self, the_protein_id: int):
#         sql = """
#             DELETE FROM Chain
#             WHERE protein_id = ?
#         """
#         self._cursor.execute(sql, (the_protein_id,))
#         self._connection.commit()
#
#     def _delete_pymol_parameter(self, the_chain_id):
#         sql = """
#             DELETE FROM PyMOLParameter
#             WHERE chain_id = ?
#         """
#         self._cursor.execute(sql, (the_chain_id,))
#         self._connection.commit()
#
#     def _delete_pymol_selection(self, the_protein_id: int):
#         sql = """
#             DELETE FROM PyMOLSelection
#             WHERE protein_id = ?
#         """
#         self._cursor.execute(sql, (the_protein_id,))
#         self._connection.commit()
#
#     def _delete_pdb_atom(self, the_protein_id: int):
#         sql = """   DELETE FROM PdbAtom
#                     WHERE protein_id = ?
#         """
#         self._cursor.execute(sql, (the_protein_id,))
#         self._connection.commit()
#
#     # </editor-fold>
#
#     # <editor-fold desc="Protein pair object inserts">
#     def insert_new_protein_pair(self, a_protein_pair: "protein_pair.ProteinPair") -> int:
#         """Inserts a new protein pair in the project database.
#
#         Note:
#             Run this method after distance analysis finished!!!
#         """
#         tmp_protein_pair_id = self._insert_protein_pair(a_protein_pair)
#         a_protein_pair.set_id(tmp_protein_pair_id)
#         for tmp_chain in a_protein_pair.protein_1.chains:
#             self._insert_pymol_parameter_protein_pair(
#                 a_protein_pair.protein_1.get_id(),
#                 tmp_chain.chain_letter,
#                 enums.PymolParameterEnum.COLOR.value,
#                 tmp_chain.pymol_parameters[enums.PymolParameterEnum.COLOR.value],
#                 a_protein_pair.get_id(),
#             )
#             self._insert_pymol_parameter_protein_pair(
#                 a_protein_pair.protein_1.get_id(),
#                 tmp_chain.chain_letter,
#                 enums.PymolParameterEnum.REPRESENTATION.value,
#                 tmp_chain.pymol_parameters[enums.PymolParameterEnum.REPRESENTATION.value],
#                 a_protein_pair.get_id(),
#             )
#         for tmp_chain in a_protein_pair.protein_2.chains:
#             self._insert_pymol_parameter_protein_pair(
#                 a_protein_pair.protein_2.get_id(),
#                 tmp_chain.chain_letter,
#                 enums.PymolParameterEnum.COLOR.value,
#                 tmp_chain.pymol_parameters[enums.PymolParameterEnum.COLOR.value],
#                 a_protein_pair.get_id(),
#             )
#             self._insert_pymol_parameter_protein_pair(
#                 a_protein_pair.protein_2.get_id(),
#                 tmp_chain.chain_letter,
#                 enums.PymolParameterEnum.REPRESENTATION.value,
#                 tmp_chain.pymol_parameters[enums.PymolParameterEnum.REPRESENTATION.value],
#                 a_protein_pair.get_id(),
#             )
#         tmp_distance_analysis_id = self._insert_distance_analysis(a_protein_pair.distance_analysis, tmp_protein_pair_id)
#         tmp_distance_analysis_results_id = self._insert_distance_analysis_results(
#             a_protein_pair.distance_analysis.analysis_results, tmp_distance_analysis_id
#         )
#         self._insert_distance_data_records(tmp_distance_analysis_results_id,
#                                            a_protein_pair.distance_analysis.analysis_results.distance_data)
#         return tmp_protein_pair_id
#
#     def _insert_protein_pair(self, a_protein_pair: "protein_pair.ProteinPair") -> int:
#         sql = """   INSERT INTO ProteinPair(protein_1_id, protein_2_id, pymol_session, project_id, name)
#                     VALUES (?, ?, ?, ?, ?)
#         """
#         tmp_params = (a_protein_pair.protein_1.get_id(), a_protein_pair.protein_2.get_id(),
#                       a_protein_pair.pymol_session, a_protein_pair.db_project_id, a_protein_pair.name)
#         self._cursor.execute(sql, tmp_params)
#         self._connection.commit()
#         return self.get_last_id()
#
#     def _insert_pymol_parameter_protein_pair(self, a_protein_id,
#                                              a_chain_letter,
#                                              a_parameter_name,
#                                              a_parameter_value,
#                                              the_protein_pair_id):
#         sql = """   INSERT INTO PyMOLParameterProteinPair(protein_id, chain_letter, parameter_name, parameter_value, protein_pair_id)
#                     VALUES (?, ?, ?, ?, ?)
#         """
#         tmp_params = (a_protein_id, a_chain_letter, a_parameter_name, a_parameter_value, the_protein_pair_id)
#         self._cursor.execute(sql, tmp_params)
#         self._connection.commit()
#
#     def _insert_distance_analysis(self,
#                                   a_distance_analysis: structure_analysis.DistanceAnalysis,
#                                   the_protein_pair_id) -> int:
#         sql = """   INSERT INTO DistanceAnalysis(name, cutoff, cycles, protein_pair_id, figure_size_x, figure_size_y)
#                     VALUES (?, ?, ?, ?, ?, ?)
#         """
#         tmp_params = (
#             a_distance_analysis.name,
#             a_distance_analysis.cutoff,
#             a_distance_analysis.cycles,
#             the_protein_pair_id,
#             a_distance_analysis.figure_size[0],
#             a_distance_analysis.figure_size[1]
#         )
#         self._cursor.execute(sql, tmp_params)
#         self._connection.commit()
#         return self.get_last_id()
#
#     def _insert_distance_analysis_results(self, a_distance_analysis_result: results.DistanceAnalysisResults,
#                                           the_distance_analysis_id: int) -> int:
#         sql = """   INSERT INTO DistanceAnalysisResults(pymol_session, rmsd, aligned_aa, distance_analysis_id)
#                     VALUES (?, ?, ?, ?)
#         """
#         tmp_params = (
#             a_distance_analysis_result.pymol_session,
#             a_distance_analysis_result.rmsd,
#             a_distance_analysis_result.aligned_aa,
#             the_distance_analysis_id
#         )
#         self._cursor.execute(sql, tmp_params)
#         self._connection.commit()
#         return self.get_last_id()
#
#     def _insert_distance_data_records(self, the_distance_analysis_results_id: int, distance_data: dict):
#         index = list(distance_data[pyssa_keys.ARRAY_DISTANCE_INDEX])
#         prot_1_chains = list(distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN])
#         prot_1_position = list(distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION])
#         prot_1_residue = list(distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI])
#         prot_2_chains = list(distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN])
#         prot_2_position = list(distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION])
#         prot_2_residue = list(distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI])
#         distances = list(distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES])
#
#         sql = """   INSERT INTO DistanceAnalysisResultData(my_index, protein_1_chain, protein_1_position, protein_1_residue,
#                                                             protein_2_chain, protein_2_position, protein_2_residue,
#                                                             distances, distance_analysis_results_id)
#                     VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
#         """
#         for i in range(len(index)):
#             tmp_params = (
#                 int(index[i]),
#                 prot_1_chains[i],
#                 int(prot_1_position[i]),
#                 prot_1_residue[i],
#                 prot_2_chains[i],
#                 int(prot_2_position[i]),
#                 prot_2_residue[i],
#                 distances[i],
#                 the_distance_analysis_results_id,
#             )
#             logger.info(tmp_params)
#             self._cursor.execute(sql, tmp_params)
#             i += 1
#         self._connection.commit()
#
#     # </editor-fold>
#
#     # <editor-fold desc="Protein pair delete statements">
#     def delete_existing_protein_pair(self, a_protein_pair_id: int) -> None:
#         # Get ids of distance analysis related objects
#         tmp_distance_analysis_info = self._get_distance_analysis(a_protein_pair_id)
#         tmp_distance_analysis_id, _, _, _, _, _ = tmp_distance_analysis_info
#         tmp_distance_analysis_results_info = self._get_distance_analysis_results(tmp_distance_analysis_id)
#         tmp_dist_analysis_results_id, _, _, _ = tmp_distance_analysis_results_info
#
#         # delete statements
#         self._delete_distance_analysis_result_data(tmp_dist_analysis_results_id)
#         self._delete_distance_analysis_results(tmp_distance_analysis_id)
#         self._delete_distance_analysis(a_protein_pair_id)
#         self._delete_pymol_parameters_protein_pair(a_protein_pair_id)
#         self._delete_protein_pair(a_protein_pair_id)
#
#     def _delete_protein_pair(self, a_protein_pair_id: int) -> None:
#         sql = """
#             DELETE FROM ProteinPair
#             WHERE id = ?
#         """
#         self._cursor.execute(sql, (a_protein_pair_id,))
#         self._connection.commit()
#
#     def _delete_pymol_parameters_protein_pair(self, a_protein_pair_id: int) -> None:
#         sql = """
#             DELETE FROM PyMOLParameterProteinPair
#             WHERE protein_pair_id = ?
#         """
#         self._cursor.execute(sql, (a_protein_pair_id,))
#         self._connection.commit()
#
#     def _delete_distance_analysis(self, a_protein_pair_id: int) -> None:
#         sql = """
#             DELETE FROM DistanceAnalysis
#             WHERE protein_pair_id = ?
#         """
#         self._cursor.execute(sql, (a_protein_pair_id,))
#         self._connection.commit()
#
#     def _delete_distance_analysis_results(self, a_distance_analysis_id: int) -> None:
#         sql = """
#             DELETE FROM DistanceAnalysisResults
#             WHERE distance_analysis_id = ?
#         """
#         self._cursor.execute(sql, (a_distance_analysis_id,))
#         self._connection.commit()
#
#     def _delete_distance_analysis_result_data(self, a_distance_analysis_results_id: int) -> None:
#         sql = """
#             DELETE FROM DistanceAnalysisResultData
#             WHERE distance_analysis_results_id = ?
#         """
#         self._cursor.execute(sql, (a_distance_analysis_results_id,))
#         self._connection.commit()
#     # </editor-fold>
#
#     # <editor-fold desc="Get all data of a specific table">
#     def get_all_protein_table_data(self):
#         """Gets all records from the Protein table.
#
#         Returns:
#             a list of protein records (id, pymol_molecule_object, pdb_filepath, fasta_filepath,
#             export_dirname, pymol_session_filepath, pymol_session)
#         """
#         sql = """
#             SELECT id, pymol_molecule_object, pdb_filepath, fasta_filepath,
#             export_dirname, pymol_session_filepath, pymol_session FROM Protein"""
#         self._cursor.execute(sql)
#         return self._cursor.fetchall()
#
#     def get_all_chain_table_data(self):
#         """Gets all records from the Chain table.
#
#         Returns:
#             a list of chain records (id, protein_id, chain_identifier, chain_type, chain_sequence)
#         """
#         sql = """
#             SELECT id, protein_id, chain_identifier, chain_type, chain_sequence
#             FROM Chain"""
#         self._cursor.execute(sql)
#         return self._cursor.fetchall()
#
#     def get_all_pymol_parameter_table_data(self):
#         """Gets all records from the PyMOLParameter table.
#
#         Returns:
#             a list of pymol parameter records (id, color, representation, chain_id)
#         """
#         sql = """
#             SELECT id, color, representation, chain_id
#             FROM PyMOLParameter"""
#         self._cursor.execute(sql)
#         return self._cursor.fetchall()
#
#     def get_all_pymol_selection_table_data(self):
#         """Gets all records from the PyMOLSelection table.
#
#         Returns:
#             a list of pymol selection records (id, selection_string, protein_id)
#         """
#         sql = """
#             SELECT id, selection_string, protein_id
#             FROM PyMOLSelection"""
#         self._cursor.execute(sql)
#         return self._cursor.fetchall()
#
#     def get_all_pdb_atom_table_data(self):
#         """Gets all records from the PdbAtom table.
#
#         Returns:
#             a list of pdb atom records (id, record_type, atom_number, atom_name, alternate_location_indicator, residue_name,
#             chain_identifier, residue_sequence_number, code_for_insertions_of_residues,
#             x_coord, y_coord, z_coord, occupancy, temperature_factor, segment_identifier, element_symbol,
#             charge, protein_id)
#         """
#         sql = """
#             SELECT id, record_type, atom_number, atom_name, alternate_location_indicator, residue_name,
#             chain_identifier, residue_sequence_number, code_for_insertions_of_residues,
#             x_coord, y_coord, z_coord, occupancy, temperature_factor, segment_identifier, element_symbol,
#             charge, protein_id
#             FROM PdbAtom"""
#         self._cursor.execute(sql)
#         return self._cursor.fetchall()
#
#     # </editor-fold>
#
#     # <editor-fold desc="Project related statements">
#     def insert_new_project(self, a_project_name, an_os) -> int:
#         """Writes a new empty project to the database."""
#         sql = """   INSERT INTO Project(name, os)
#                     VALUES(:name, :os)
#         """
#         self._cursor.execute(sql, {"name": a_project_name, "os": an_os})
#         self._connection.commit()
#         return self.get_last_id()
#
#     def get_project_as_object(self,
#                               a_project_name: str,
#                               a_workspace_path: pathlib.Path,
#                               the_app_settings,
#                               the_progress_signal: "custom_signals.ProgressSignal" = custom_signals.ProgressSignal()) -> "project.Project":
#         """Creates a project object based on the data from the project database."""
#         tmp_project = project.Project(a_project_name, a_workspace_path)
#         tmp_project_id, = self._get_project(a_project_name)
#         tmp_project.set_id(tmp_project_id)
#
#         tmp_progress = 0
#         i = 1
#         for tmp_seq_info in self._get_sequence():
#             tmp_progress = min(tmp_progress + i, 40)
#             the_progress_signal.emit_signal(f"Loading sequence ({i}/{len(self._get_sequence())}) ...", tmp_progress)
#             tmp_seq_id, tmp_seq, tmp_seq_name = tmp_seq_info
#             tmp_sequence = SeqRecord.SeqRecord(tmp_seq, id=tmp_seq_id, name=tmp_seq_name)
#             tmp_project.sequences.append(tmp_sequence)
#             i += 1
#
#         i = 1
#         # create protein objects
#         for tmp_protein_info in self._get_protein(tmp_project_id):
#             tmp_progress = 20
#             tmp_message = f"Loading protein ({i}/{len(self._get_protein(tmp_project_id))}) ..."
#             the_progress_signal.emit_signal(tmp_message, tmp_progress)
#
#             tmp_protein_id, tmp_protein_name, tmp_pymol_session = tmp_protein_info
#             tmp_protein = protein.Protein(tmp_protein_name)
#             tmp_protein.set_id(tmp_protein_id)
#             tmp_protein.pymol_session = tmp_pymol_session
#
#             # Chains
#             j = 5
#             for tmp_chain_info in self._get_chain(tmp_protein_id):
#                 tmp_progress += j
#                 the_progress_signal.emit_signal(tmp_message, min(tmp_progress, 50))
#
#                 tmp_chain_id, tmp_chain_identifier, tmp_chain_type, tmp_chain_sequence = tmp_chain_info
#                 tmp_seq = sequence.Sequence(tmp_protein_name, tmp_chain_sequence)
#                 tmp_chain = chain.Chain(tmp_chain_identifier, tmp_seq, tmp_chain_type)
#                 tmp_chain.set_id(tmp_chain_id)
#                 tmp_chain.db_protein_id = tmp_protein_id
#                 tmp_color, tmp_representation = self._get_pymol_parameter(tmp_chain_id)
#                 tmp_chain.pymol_parameters = {
#                     enums.PymolParameterEnum.COLOR.value: tmp_color,
#                     enums.PymolParameterEnum.REPRESENTATION.value: tmp_representation,
#                 }
#                 tmp_protein.chains.append(tmp_chain)
#                 j += 3
#             # pymol selection
#             tmp_selection, = self._get_pymol_selection(tmp_protein_id)
#             tmp_protein.pymol_selection.selection_string = tmp_selection
#             # add protein to project
#             tmp_project.proteins.append(tmp_protein)
#             i += 1
#
#         i = 1
#         # create protein pair objects
#         for tmp_protein_pair_info in self._get_protein_pair(tmp_project_id):
#             tmp_progress = 20
#             tmp_message = f"Loading protein pair ({i}/{len(self._get_protein_pair(tmp_project_id))}) ..."
#             the_progress_signal.emit_signal(tmp_message, tmp_progress)
#             # create protein pair
#             tmp_protein_pair_id, tmp_protein_1_id, tmp_protein_2_id, tmp_pymol_session, tmp_pp_name = tmp_protein_pair_info
#             tmp_protein_1_name, = self._get_protein_name_by_id(tmp_protein_1_id)
#             tmp_protein_2_name, = self._get_protein_name_by_id(tmp_protein_2_id)
#             tmp_protein_pair = protein_pair.ProteinPair(tmp_project.search_protein(tmp_protein_1_name),
#                                                         tmp_project.search_protein(tmp_protein_2_name))
#             tmp_protein_pair.set_id(tmp_protein_pair_id)
#             tmp_protein_pair.db_project_id = tmp_project_id
#             tmp_protein_pair.name = tmp_pp_name
#             tmp_protein_pair.pymol_session = tmp_pymol_session
#
#             tmp_progress += 20
#             the_progress_signal.emit_signal(tmp_message, tmp_progress)
#             # create distance analysis object
#             tmp_distance_analysis_info = self._get_distance_analysis(tmp_protein_pair_id)
#             tmp_distance_analysis_id, tmp_name, tmp_cutoff, tmp_cycles, tmp_figure_size_x, tmp_figure_size_y = tmp_distance_analysis_info
#             tmp_distance_analysis = structure_analysis.DistanceAnalysis(the_app_settings)
#             tmp_distance_analysis.name = tmp_name
#             tmp_distance_analysis.cutoff = tmp_cutoff
#             tmp_distance_analysis.cycles = tmp_cycles
#             tmp_distance_analysis.figure_size = (tmp_figure_size_x, tmp_figure_size_y)
#
#             tmp_protein_pair.distance_analysis = tmp_distance_analysis
#
#             tmp_progress += 10
#             the_progress_signal.emit_signal(tmp_message, tmp_progress)
#             # create distance analysis results object
#             tmp_distance_analysis_results_info = self._get_distance_analysis_results(tmp_distance_analysis_id)
#             tmp_dist_analysis_results_id, tmp_pymol_session, tmp_rmsd, tmp_aligned_aa = tmp_distance_analysis_results_info
#
#             index = []
#             prot_1_chains = []
#             prot_1_position = []
#             prot_1_residue = []
#             prot_2_chains = []
#             prot_2_position = []
#             prot_2_residue = []
#             distances = []
#             for tmp_distance_data in self._get_distance_analysis_result_data(tmp_dist_analysis_results_id):
#                 tmp_progress += 0.5
#                 print(tmp_distance_data[2])
#                 the_progress_signal.emit_signal(tmp_message, min(round(tmp_progress, 2), 90))
#                 index.append(tmp_distance_data[0])
#                 prot_1_chains.append(tmp_distance_data[1])
#                 prot_1_position.append(tmp_distance_data[2])
#                 prot_1_residue.append(tmp_distance_data[3])
#                 prot_2_chains.append(tmp_distance_data[4])
#                 prot_2_position.append(tmp_distance_data[5])
#                 prot_2_residue.append(tmp_distance_data[6])
#                 distances.append(tmp_distance_data[7])
#
#             tmp_distance_data_records = {
#                 pyssa_keys.ARRAY_DISTANCE_INDEX: np.array(index),
#                 pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN: np.array(prot_1_chains),
#                 pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION: np.array(prot_1_position),
#                 pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI: np.array(prot_1_residue),
#                 pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN: np.array(prot_2_chains),
#                 pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION: np.array(prot_2_position),
#                 pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI: np.array(prot_2_residue),
#                 pyssa_keys.ARRAY_DISTANCE_DISTANCES: np.array(distances),
#             }
#             tmp_dist_analysis_results = results.DistanceAnalysisResults(tmp_distance_data_records,
#                                                                         tmp_pymol_session,
#                                                                         tmp_rmsd,
#                                                                         tmp_aligned_aa)
#             tmp_protein_pair.distance_analysis.analysis_results = tmp_dist_analysis_results
#             tmp_project.protein_pairs.append(tmp_protein_pair)
#             i += 1
#
#         return tmp_project
#
#     def _get_project(self, a_project_name: str) -> tuple:
#         sql = """SELECT id FROM Project WHERE name = ?"""
#         self._cursor.execute(sql, (a_project_name,))
#         return self._cursor.fetchall()[0]
#
#     # </editor-fold>
#
#     # <editor-fold desc="Select statements for protein objects">
#     def get_protein_as_object(self, the_pymol_molecule_object: str) -> "protein.Protein":
#         sql = """SELECT id, pymol_molecule_object, pymol_session FROM Protein WHERE pymol_molecule_object = ?"""
#         self._cursor.execute(sql, (the_pymol_molecule_object,))
#         for tmp_protein_info in self._cursor.fetchall():
#             tmp_protein_id, tmp_protein_name, tmp_pymol_session = tmp_protein_info
#             tmp_protein = protein.Protein(tmp_protein_name)
#             tmp_protein.set_id(tmp_protein_id)
#             tmp_protein.pymol_session = tmp_pymol_session
#             # Chains
#             for tmp_chain_info in self._get_chain(tmp_protein_id):
#                 tmp_chain_id, tmp_chain_identifier, tmp_chain_type, tmp_chain_sequence = tmp_chain_info
#                 tmp_seq = sequence.Sequence(tmp_protein_name, tmp_chain_sequence)
#                 tmp_chain = chain.Chain(tmp_chain_identifier, tmp_seq, tmp_chain_type)
#                 tmp_chain.set_id(tmp_chain_id)
#                 tmp_chain.db_protein_id = tmp_protein_id
#                 tmp_color, tmp_representation = self._get_pymol_parameter(tmp_chain_id)
#                 tmp_chain.pymol_parameters = {
#                     enums.PymolParameterEnum.COLOR.value: tmp_color,
#                     enums.PymolParameterEnum.REPRESENTATION.value: tmp_representation,
#                 }
#                 tmp_protein.chains.append(tmp_chain)
#             # pymol selection
#             tmp_selection, = self._get_pymol_selection(tmp_protein_id)
#             tmp_protein.pymol_selection.selection_string = tmp_selection
#             return tmp_protein
#
#     def _get_protein(self, the_project_id: int) -> list[tuple]:
#         sql = """SELECT id, pymol_molecule_object, pymol_session FROM Protein WHERE project_id = ?"""
#         self._cursor.execute(sql, (the_project_id,))
#         return self._cursor.fetchall()
#
#     def _get_pymol_selection(self, the_protein_id: int) -> tuple:
#         sql = """SELECT selection_string FROM PyMOLSelection WHERE protein_id = ?"""
#         self._cursor.execute(sql, (the_protein_id,))
#         return self._cursor.fetchall()[0]
#
#     def _get_chain(self, the_protein_id: int) -> list[tuple]:
#         sql = """SELECT id, chain_identifier, chain_type, chain_sequence FROM Chain WHERE protein_id = ?"""
#         self._cursor.execute(sql, (the_protein_id,))
#         return self._cursor.fetchall()
#
#     def _get_pymol_parameter(self, the_chain_id: int) -> tuple:
#         sql = """SELECT color, representation FROM PyMOLParameter WHERE chain_id = ?"""
#         self._cursor.execute(sql, (the_chain_id,))
#         return self._cursor.fetchall()[0]
#
#     def _get_protein_name_by_id(self, a_protein_id) -> tuple:
#         sql = """SELECT pymol_molecule_object FROM Protein WHERE id = ?"""
#         self._cursor.execute(sql, (a_protein_id,))
#         return self._cursor.fetchall()[0]
#
#     def get_pdb_atoms_of_protein(self, the_protein_id: int) -> list[tuple]:
#         sql = """
#             SELECT record_type, atom_number, atom_name, alternate_location_indicator, residue_name,
#             chain_identifier, residue_sequence_number, code_for_insertions_of_residues,
#             x_coord, y_coord, z_coord, occupancy, temperature_factor, segment_identifier, element_symbol,
#             charge, protein_id
#             FROM PdbAtom
#             WHERE protein_id = ?"""
#         # tmp_params = (
#         #     a_pdb_atom_dict["record_type"],
#         #     a_pdb_atom_dict["atom_number"],
#         #     a_pdb_atom_dict["atom_name"],
#         #     a_pdb_atom_dict["alternate_location_indicator"],
#         #     a_pdb_atom_dict["residue_name"],
#         #     a_pdb_atom_dict["chain_identifier"],
#         #     a_pdb_atom_dict["residue_sequence_number"],
#         #     a_pdb_atom_dict["code_for_insertions_of_residues"],
#         #     a_pdb_atom_dict["x_coord"],
#         #     a_pdb_atom_dict["y_coord"],
#         #     a_pdb_atom_dict["z_coord"],
#         #     a_pdb_atom_dict["occupancy"],
#         #     a_pdb_atom_dict["temperature_factor"],
#         #     a_pdb_atom_dict["segment_identifier"],
#         #     a_pdb_atom_dict["element_symbol"],
#         #     a_pdb_atom_dict["charge"],
#         #     the_protein_id
#         # )
#         self._cursor.execute(sql, (the_protein_id,))
#         return self._cursor.fetchall()
#
#     # </editor-fold>
#
#     # <editor-fold desc="Select statements for protein pair objects">
#     def get_protein_pair_as_object(self, a_name: str, a_project: "project.Project", the_app_settings) -> "protein_pair.ProteinPair":
#         sql = """SELECT id, protein_1_id, protein_2_id, pymol_session, name FROM ProteinPair WHERE name = ?"""
#         self._cursor.execute(sql, (a_name,))
#         for tmp_protein_pair_info in self._cursor.fetchall():
#             # create protein pair
#             tmp_protein_pair_id, tmp_protein_1_id, tmp_protein_2_id, tmp_pymol_session, tmp_pp_name = tmp_protein_pair_info
#             tmp_protein_1_name, = self._get_protein_name_by_id(tmp_protein_1_id)
#             tmp_protein_2_name, = self._get_protein_name_by_id(tmp_protein_2_id)
#             tmp_protein_pair = protein_pair.ProteinPair(a_project.search_protein(tmp_protein_1_name),
#                                                         a_project.search_protein(tmp_protein_2_name))
#             tmp_protein_pair.set_id(tmp_protein_pair_id)
#             tmp_protein_pair.db_project_id = a_project.get_id()
#             tmp_protein_pair.name = tmp_pp_name
#             tmp_protein_pair.pymol_session = tmp_pymol_session
#
#             # create distance analysis object
#             tmp_distance_analysis_info = self._get_distance_analysis(tmp_protein_pair_id)
#             tmp_distance_analysis_id, tmp_name, tmp_cutoff, tmp_cycles, tmp_figure_size_x, tmp_figure_size_y = tmp_distance_analysis_info
#             tmp_distance_analysis = structure_analysis.DistanceAnalysis(the_app_settings)
#             tmp_distance_analysis.name = tmp_name
#             tmp_distance_analysis.cutoff = tmp_cutoff
#             tmp_distance_analysis.cycles = tmp_cycles
#             tmp_distance_analysis.figure_size = (tmp_figure_size_x, tmp_figure_size_y)
#
#             tmp_protein_pair.distance_analysis = tmp_distance_analysis
#
#             # create distance analysis results object
#             tmp_distance_analysis_results_info = self._get_distance_analysis_results(tmp_distance_analysis_id)
#             tmp_dist_analysis_results_id, tmp_pymol_session, tmp_rmsd, tmp_aligned_aa = tmp_distance_analysis_results_info
#
#             index = []
#             prot_1_chains = []
#             prot_1_position = []
#             prot_1_residue = []
#             prot_2_chains = []
#             prot_2_position = []
#             prot_2_residue = []
#             distances = []
#             for tmp_distance_data in self._get_distance_analysis_result_data(tmp_dist_analysis_results_id):
#                 index.append(int(tmp_distance_data[0]))
#                 prot_1_chains.append(tmp_distance_data[1])
#                 prot_1_position.append(tmp_distance_data[2])
#                 prot_1_residue.append(tmp_distance_data[3])
#                 prot_2_chains.append(tmp_distance_data[4])
#                 prot_2_position.append(tmp_distance_data[5])
#                 prot_2_residue.append(tmp_distance_data[6])
#                 distances.append(tmp_distance_data[7])
#
#             tmp_distance_data_records = {
#                 pyssa_keys.ARRAY_DISTANCE_INDEX: np.array(index),
#                 pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN: np.array(prot_1_chains),
#                 pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION: np.array(prot_1_position),
#                 pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI: np.array(prot_1_residue),
#                 pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN: np.array(prot_2_chains),
#                 pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION: np.array(prot_2_position),
#                 pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI: np.array(prot_2_residue),
#                 pyssa_keys.ARRAY_DISTANCE_DISTANCES: np.array(distances),
#             }
#             tmp_dist_analysis_results = results.DistanceAnalysisResults(tmp_distance_data_records,
#                                                                         tmp_pymol_session,
#                                                                         tmp_rmsd,
#                                                                         tmp_aligned_aa)
#             tmp_protein_pair.distance_analysis.analysis_results = tmp_dist_analysis_results
#             return tmp_protein_pair
#
#     def _get_protein_pair(self, the_project_id) -> list[tuple]:
#         sql = """SELECT id, protein_1_id, protein_2_id, pymol_session, name FROM ProteinPair WHERE project_id = ?"""
#         self._cursor.execute(sql, (the_project_id,))
#         return self._cursor.fetchall()
#
#     def _get_distance_analysis(self, the_protein_pair_id: int) -> tuple:
#         sql = """SELECT id, name, cutoff, cycles, figure_size_x, figure_size_y FROM DistanceAnalysis WHERE protein_pair_id = ?"""
#         self._cursor.execute(sql, (the_protein_pair_id,))
#         return self._cursor.fetchall()[0]
#
#     def _get_distance_analysis_results(self, the_distance_analysis_id: int) -> tuple:
#         sql = """SELECT id, pymol_session, rmsd, aligned_aa FROM DistanceAnalysisResults WHERE distance_analysis_id = ?"""
#         self._cursor.execute(sql, (the_distance_analysis_id,))
#         return self._cursor.fetchall()[0]
#
#     def _get_distance_analysis_result_data(self, the_distance_analysis_results_id) -> list[tuple]:
#         sql = """   SELECT my_index, protein_1_chain, protein_1_position, protein_1_residue,
#                     protein_2_chain, protein_2_position, protein_2_residue, distances
#                     FROM DistanceAnalysisResultData WHERE distance_analysis_results_id = ?"""
#         self._cursor.execute(sql, (the_distance_analysis_results_id,))
#         return self._cursor.fetchall()
#
#     def get_pymol_parameter_for_certain_protein_chain_in_protein_pair(self,
#                                                                       a_protein_pair_id: int,
#                                                                       a_protein_id: int,
#                                                                       a_chain_letter: str,
#                                                                       a_parameter_name: str) -> tuple:
#         sql = """
#             SELECT parameter_value
#             FROM PyMOLParameterProteinPair WHERE protein_id = ? and chain_letter = ? and protein_pair_id = ? and parameter_name = ?"""
#         self._cursor.execute(sql, (a_protein_id, a_chain_letter, a_protein_pair_id, a_parameter_name))
#         return self._cursor.fetchall()[0]
#
#     # </editor-fold>
#
#     # <editor-fold desc="Select statements for sequences">
#     def _get_sequence(self) -> list[tuple]:
#         sql = """SELECT seq_id, seq, name FROM SeqRecord WHERE project_id = ?"""
#         self._cursor.execute(sql, (1,))  # fixme: this is not the best solution
#         return self._cursor.fetchall()
#
#     # </editor-fold>
#
#     # <editor-fold desc="Update statements for protein objects">
#     def update_protein_chain_color(self, a_chain_id: int, a_color: str):
#         sql = """
#             UPDATE PyMOLParameter
#             SET color = ?
#             WHERE chain_id = ?
#         """
#         self._cursor.execute(sql, (a_color, a_chain_id))
#         self._connection.commit()
#
#     def update_protein_chain_representation(self, a_chain_id: int, a_representation: str):
#         sql = """
#             UPDATE PyMOLParameter
#             SET representation = ?
#             WHERE chain_id = ?
#         """
#         self._cursor.execute(sql, (a_representation, a_chain_id))
#         self._connection.commit()
#
#     def update_protein_name(self, the_new_protein_name: str, the_old_protein_name: str, the_protein_id: int):
#         sql = """
#             UPDATE Protein
#             SET pymol_molecule_object = ?
#             WHERE pymol_molecule_object = ? and id = ?
#         """
#         self._cursor.execute(sql, (the_new_protein_name, the_old_protein_name, the_protein_id))
#         self._connection.commit()
#
#     def update_protein_pdb_atom_data(self, the_protein_id: int, a_pdb_atom_dict_list: list[dict]):
#         self._delete_pdb_atom(the_protein_id)
#
#         for tmp_pdb_atom_dict in a_pdb_atom_dict_list:
#             self._insert_pdb_atom(the_protein_id, tmp_pdb_atom_dict)
#
#         # sql = """
#         #     UPDATE PdbAtom
#         #     SET
#         #         record_type = ?,
#         #         atom_number = ?,
#         #         atom_name = ?,
#         #         alternate_location_indicator = ?,
#         #         residue_name = ?,
#         #         chain_identifier = ?,
#         #         residue_sequence_number = ?,
#         #         code_for_insertions_of_residues = ?,
#         #         x_coord = ?,
#         #         y_coord = ?,
#         #         z_coord = ?,
#         #         occupancy = ?,
#         #         temperature_factor = ?,
#         #         segment_identifier = ?,
#         #         element_symbol = ?,
#         #         charge = ?
#         #     WHERE
#         #         protein_id = ?;
#         #     """
#         # tmp_params = (
#         #     a_pdb_atom_dict["record_type"],
#         #     a_pdb_atom_dict["atom_number"],
#         #     a_pdb_atom_dict["atom_name"],
#         #     a_pdb_atom_dict["alternate_location_indicator"],
#         #     a_pdb_atom_dict["residue_name"],
#         #     a_pdb_atom_dict["chain_identifier"],
#         #     a_pdb_atom_dict["residue_sequence_number"],
#         #     a_pdb_atom_dict["code_for_insertions_of_residues"],
#         #     a_pdb_atom_dict["x_coord"],
#         #     a_pdb_atom_dict["y_coord"],
#         #     a_pdb_atom_dict["z_coord"],
#         #     a_pdb_atom_dict["occupancy"],
#         #     a_pdb_atom_dict["temperature_factor"],
#         #     a_pdb_atom_dict["segment_identifier"],
#         #     a_pdb_atom_dict["element_symbol"],
#         #     a_pdb_atom_dict["charge"],
#         #     the_protein_id
#         # )
#         # self._cursor.execute(sql, tmp_params)
#         # self._connection.commit()
#
#     # </editor-fold>
#
#     def update_pymol_parameter_for_certain_protein_chain_in_protein_pair(self,
#                                                                          a_protein_pair_id: int,
#                                                                          a_protein_id: int,
#                                                                          a_chain_letter: str,
#                                                                          the_parameter_name: str,
#                                                                          the_new_parameter_value: str):
#         sql = """
#                     UPDATE PyMOLParameterProteinPair
#                     SET parameter_value = ?
#                     WHERE protein_id = ? and chain_letter = ? and protein_pair_id = ? and parameter_name = ?
#                 """
#         self._cursor.execute(sql, (the_new_parameter_value, a_protein_id, a_chain_letter, a_protein_pair_id, the_parameter_name))
#         self._connection.commit()
#         # sql = """
#         #     UPDATE PyMOLParameterProteinPair
#         #     FROM PyMOLParameterProteinPair WHERE protein_id = ? and chain_letter = ? and protein_pair_id = ? and parameter_name = ?"""
#         # self._cursor.execute(sql, (a_protein_id, a_chain_letter, a_protein_pair_id, a_parameter_name))
#         # return self._cursor.fetchall()[0]
#
#     def update_pymol_session_of_protein_pair(self, the_protein_pair_id: int, the_new_pymol_session: str):
#         sql = """
#                             UPDATE ProteinPair
#                             SET pymol_session = ?
#                             WHERE id = ?
#                         """
#         self._cursor.execute(sql, (str(the_new_pymol_session), int(the_protein_pair_id)))
#         self._connection.commit()
#
#     def update_pymol_session_of_protein(self, the_protein_id: int, the_new_pymol_session: str):
#         sql = """
#                             UPDATE Protein
#                             SET pymol_session = ?
#                             WHERE id = ?
#                         """
#         self._cursor.execute(sql, (str(the_new_pymol_session), int(the_protein_id)))
#         self._connection.commit()
#
#     def update_sequence_name(self, the_new_seq_name: str, the_old_seq_name: str, the_sequence: str):
#         sql = """
#             UPDATE SeqRecord
#             SET name = ?
#             WHERE name = ? and seq = ?
#         """
#         self._cursor.execute(sql, (the_new_seq_name, the_old_seq_name, the_sequence))
#         self._connection.commit()
#
#     def update_project_name(self, the_new_project_name: str):
#         sql = """
#             UPDATE Project
#             SET name = ?
#             WHERE id = ?
#         """
#         self._cursor.execute(sql, (the_new_project_name, 1))
#         self._connection.commit()
