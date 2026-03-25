# project_database.py
"""
ProjectDatabase: the single public API for all persistence on one project.

The write_queue attribute exposes a ProjectWriteQueue for fire-and-forget
writes from any thread.  Direct read/write methods on this class are
synchronous and safe to call from any thread (each gets its own connection
via ProjectConnectionPool).
"""
from __future__ import annotations

import logging
from contextlib import contextmanager
from typing import Generator

from src.pyssa.gui.qt import QtSql

from src.pyssa.io_pyssa.db_pyssa.connection import ProjectConnectionPool
from src.pyssa.io_pyssa.db_pyssa.repositories import (
    ProjectRepository, SequenceRepository, ProteinRepository,
    ChainRepository, PyMOLParameterRepository, PyMOLSelectionRepository,
    PdbAtomRepository, ProteinPairRepository, ProteinPairPyMOLParamRepository,
    DistanceAnalysisRepository, DistanceAnalysisResultsRepository,
    DistanceDataRepository,
)
from src.pyssa.io_pyssa.db_pyssa.service import ProjectService
from src.pyssa.io_pyssa.db_pyssa.write_queue import ProjectWriteQueue

logger = logging.getLogger(__name__)

_SCHEMA_SQL = """
PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS Project (
    id   INTEGER PRIMARY KEY,
    name TEXT NOT NULL,
    os   TEXT
);
CREATE TABLE IF NOT EXISTS SeqRecord (
    id INTEGER PRIMARY KEY, seq_id TEXT, seq TEXT, name TEXT,
    project_id INTEGER REFERENCES Project(id)
);
CREATE TABLE IF NOT EXISTS Protein (
    id INTEGER PRIMARY KEY, pymol_molecule_object TEXT,
    pymol_session TEXT, project_id INTEGER REFERENCES Project(id), pdb_id INTEGER
);
CREATE TABLE IF NOT EXISTS Chain (
    id INTEGER PRIMARY KEY, protein_id INTEGER REFERENCES Protein(id),
    chain_identifier TEXT, chain_type TEXT, chain_sequence TEXT
);
CREATE TABLE IF NOT EXISTS PyMOLParameter (
    id INTEGER PRIMARY KEY, color TEXT, representation TEXT,
    chain_id INTEGER REFERENCES Chain(id)
);
CREATE TABLE IF NOT EXISTS PyMOLSelection (
    id INTEGER PRIMARY KEY, selection_string TEXT,
    protein_id INTEGER REFERENCES Protein(id)
);
CREATE TABLE IF NOT EXISTS PdbAtom (
    id INTEGER PRIMARY KEY,
    record_type TEXT(6), atom_number INTEGER, atom_name TEXT(4),
    alternate_location_indicator TEXT(1), residue_name TEXT(3),
    chain_identifier TEXT(1), residue_sequence_number INTEGER,
    code_for_insertions_of_residues TEXT(1),
    x_coord REAL, y_coord REAL, z_coord REAL,
    occupancy REAL, temperature_factor REAL,
    segment_identifier TEXT(4), element_symbol TEXT(2), charge TEXT(2),
    protein_id INTEGER REFERENCES Protein(id)
);
CREATE TABLE IF NOT EXISTS ProteinPair (
    id INTEGER PRIMARY KEY,
    protein_1_id INTEGER REFERENCES Protein(id),
    protein_2_id INTEGER REFERENCES Protein(id),
    pymol_session TEXT, pymol_session_filepath TEXT,
    project_id INTEGER REFERENCES Project(id), name TEXT
);
CREATE TABLE IF NOT EXISTS PyMOLParameterProteinPair (
    id INTEGER PRIMARY KEY, protein_id INTEGER, chain_letter TEXT,
    parameter_name TEXT, parameter_value TEXT,
    protein_pair_id INTEGER REFERENCES ProteinPair(id)
);
CREATE TABLE IF NOT EXISTS DistanceAnalysis (
    id INTEGER PRIMARY KEY, name TEXT, cutoff REAL, cycles INTEGER,
    protein_pair_id INTEGER REFERENCES ProteinPair(id),
    figure_size_x REAL, figure_size_y REAL
);
CREATE TABLE IF NOT EXISTS DistanceAnalysisResults (
    id INTEGER PRIMARY KEY, pymol_session TEXT, rmsd REAL, aligned_aa TEXT,
    distance_analysis_id INTEGER REFERENCES DistanceAnalysis(id)
);
CREATE TABLE IF NOT EXISTS DistanceAnalysisResultData (
    id INTEGER PRIMARY KEY, my_index INTEGER,
    protein_1_chain TEXT, protein_1_position INTEGER, protein_1_residue TEXT,
    protein_2_chain TEXT, protein_2_position INTEGER, protein_2_residue TEXT,
    distances REAL,
    distance_analysis_results_id INTEGER REFERENCES DistanceAnalysisResults(id)
);
"""


class ProjectDatabase:
    """All persistence for a single PySSA project, safe from any thread.

    Args:
        db_path:    Absolute path to the SQLite file.
        project_id: Short stable string identifying the project; used to
                    namespace Qt connection names so multiple open projects
                    never clash (e.g. the project name or a UUID).

    Attributes:
        write_queue:  Fire-and-forget write API backed by ThreadRuntime.
        service:      High-level operations (full project load, atomic
                      multi-table inserts/deletes).
    """

    def __init__(self, db_path: str, project_id: str) -> None:
        if not db_path:
            raise ValueError("db_path cannot be empty")
        if not project_id:
            raise ValueError("project_id cannot be empty")

        # Validate that db_path is a string path, not necessarily existing yet
        # (it may be created later via initialise_schema)
        import pathlib
        try:
            self._db_path = str(pathlib.Path(db_path).resolve())
        except (TypeError, ValueError) as e:
            raise ValueError(f"Invalid db_path: {e}") from e

        self._project_id = project_id
        self._pool = ProjectConnectionPool(self._db_path, project_id)

        # Stateless repository singletons
        self._repo_project      = ProjectRepository()
        self._repo_sequences    = SequenceRepository()
        self._repo_proteins     = ProteinRepository()
        self._repo_chains       = ChainRepository()
        self._repo_pymol_params = PyMOLParameterRepository()
        self._repo_pymol_sel    = PyMOLSelectionRepository()
        self._repo_pdb_atoms    = PdbAtomRepository()
        self._repo_pairs        = ProteinPairRepository()
        self._repo_pair_params  = ProteinPairPyMOLParamRepository()
        self._repo_dist_anal    = DistanceAnalysisRepository()
        self._repo_dist_res     = DistanceAnalysisResultsRepository()
        self._repo_dist_data    = DistanceDataRepository()

        # Public facades
        self.service     = ProjectService(self)
        self.write_queue = ProjectWriteQueue(self)

        logger.info("ProjectDatabase initialized for '%s' at '%s'", project_id, self._db_path)

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def initialise_schema(self) -> None:
        """Create all tables if they do not exist yet (idempotent)."""
        with self._conn() as db:
            for statement in _SCHEMA_SQL.strip().split(";"):
                stmt = statement.strip()
                if stmt:
                    from .queries import run
                    run(db, stmt)
        logger.info("Schema initialised for project.")

    def close(self, drain_timeout: float = 10.0) -> None:
        """Drain pending writes then release all thread-local connections.

        Always call this before discarding a ProjectDatabase instance.
        """
        self.write_queue.drain(timeout=drain_timeout)
        self._pool.close_all()

    # ------------------------------------------------------------------
    # Transaction helper
    # ------------------------------------------------------------------

    @contextmanager
    def transaction(self) -> Generator[QtSql.QSqlDatabase, None, None]:
        """Wrap the calling thread's connection in BEGIN / COMMIT / ROLLBACK.

        Yields the QSqlDatabase so repositories can be called directly
        without an extra pool round-trip.

        Example::

            with db.transaction() as tx:
                protein_id = db._repo_proteins.insert(tx, ...)
                db._repo_chains.insert(tx, protein_id, ...)
        """
        with self._pool.connection() as conn:
            conn.transaction()
            try:
                yield conn
            except Exception:
                conn.rollback()
                raise
            else:
                conn.commit()

    # ------------------------------------------------------------------
    # Private connection helper
    # ------------------------------------------------------------------

    @contextmanager
    def _conn(self) -> Generator[QtSql.QSqlDatabase, None, None]:
        with self._pool.connection() as db:
            yield db

    # ------------------------------------------------------------------
    # Project
    # ------------------------------------------------------------------

    def get_project_id(self, name: str) -> int | None:
        """Retrieves the unique identifier of a project based on its name.

        Args:
            name: The name of the project whose ID is being retrieved. Must be a non-empty string.

        Returns:
            The unique identifier of the project if it exists, or None if no project with
            the provided name is found.

        Raises:
        ValueError
            If the provided project name is an empty string.
        """
        if not name:
            raise ValueError("Project name cannot be empty")
        with self._conn() as db:
            return self._repo_project.get_id(db, name)

    def insert_project(self, name: str, os: str) -> int:
        if not name:
            raise ValueError("Project name cannot be empty")
        if not os:
            raise ValueError("OS cannot be empty")
        with self._conn() as db:
            project_id = self._repo_project.insert(db, name, os)
            logger.info("Inserted project '%s' with id=%d", name, project_id)
            return project_id

    def update_project_name(self, new_name: str, project_id: int) -> None:
        if not new_name:
            raise ValueError("New project name cannot be empty")
        if project_id is None or project_id < 1:
            raise ValueError("Invalid project_id")
        with self._conn() as db:
            self._repo_project.update_name(db, new_name, project_id)
            logger.info("Updated project id=%d to name='%s'", project_id, new_name)

    # ------------------------------------------------------------------
    # Sequences
    # ------------------------------------------------------------------

    def get_all_sequences(self, project_id: int) -> list[dict]:
        with self._conn() as db:
            return self._repo_sequences.get_all(db, project_id)

    def insert_sequence(self, seq_id: str, seq: str, name: str, project_id: int) -> None:
        with self._conn() as db:
            self._repo_sequences.insert(db, seq_id, seq, name, project_id)

    def delete_sequence(self, name: str) -> None:
        with self._conn() as db:
            self._repo_sequences.delete(db, name)

    def update_sequence_name(self, new_name: str, old_name: str, sequence: str) -> None:
        with self._conn() as db:
            self._repo_sequences.update_name(db, new_name, old_name, sequence)

    # ------------------------------------------------------------------
    # Proteins
    # ------------------------------------------------------------------

    def get_all_proteins(self, project_id: int) -> list[dict]:
        with self._conn() as db:
            return self._repo_proteins.get_all(db, project_id)

    def get_protein_by_name(self, name: str) -> dict | None:
        with self._conn() as db:
            return self._repo_proteins.get_by_name(db, name)

    def get_protein_name_by_id(self, protein_id: int) -> str | None:
        with self._conn() as db:
            return self._repo_proteins.get_name_by_id(db, protein_id)

    def insert_protein_full(self, protein_obj: object) -> int:
        if protein_obj is None:
            raise ValueError("protein_obj cannot be None")
        protein_id = self.service.insert_protein(protein_obj)
        logger.info("Inserted protein with id=%d", protein_id)
        return protein_id

    def delete_protein_full(self, protein_id: int) -> None:
        if protein_id is None or protein_id < 1:
            raise ValueError("Invalid protein_id")
        self.service.delete_protein(protein_id)
        logger.info("Deleted protein with id=%d", protein_id)

    def delete_chain(self, protein_id: int, chain_id: int) -> None:
        with self._conn() as db:
            self._repo_chains.delete_one(db, protein_id, chain_id)

    def get_pdb_atoms(self, protein_id: int) -> list[tuple]:
        with self._conn() as db:
            return self._repo_pdb_atoms.get_all(db, protein_id)

    def replace_pdb_atoms(self, protein_id: int, atom_dicts: list[dict]) -> None:
        with self.transaction() as tx:
            self._repo_pdb_atoms.replace_all(tx, protein_id, atom_dicts)

    def update_protein_name(self, new_name: str, old_name: str, protein_id: int) -> None:
        with self._conn() as db:
            self._repo_proteins.update_name(db, new_name, old_name, protein_id)

    def update_protein_session(self, protein_id: int, pymol_session: str) -> None:
        with self._conn() as db:
            self._repo_proteins.update_session(db, protein_id, pymol_session)

    def update_chain_color(self, chain_id: int, color: str) -> None:
        with self._conn() as db:
            self._repo_pymol_params.update_color(db, chain_id, color)

    def get_chain_color(self, chain_id: int) -> str | None:
        with self._conn() as db:
            return self._repo_pymol_params.get(db, chain_id).get("color")

    # ------------------------------------------------------------------
    # Protein pairs
    # ------------------------------------------------------------------

    def get_all_protein_pairs(self, project_id: int) -> list[dict]:
        with self._conn() as db:
            return self._repo_pairs.get_all(db, project_id)

    def get_protein_pair_by_name(self, name: str) -> dict | None:
        with self._conn() as db:
            return self._repo_pairs.get_by_name(db, name)

    def insert_protein_pair_full(self, pair_obj: object) -> int:
        if pair_obj is None:
            raise ValueError("pair_obj cannot be None")
        pair_id = self.service.insert_protein_pair(pair_obj)
        logger.info("Inserted protein pair with id=%d", pair_id)
        return pair_id

    def delete_protein_pair_full(self, pair_id: int) -> None:
        if pair_id is None or pair_id < 1:
            raise ValueError("Invalid pair_id")
        self.service.delete_protein_pair(pair_id)
        logger.info("Deleted protein pair with id=%d", pair_id)

    def update_protein_pair_session(self, pair_id: int, pymol_session: str) -> None:
        with self._conn() as db:
            self._repo_pairs.update_session(db, pair_id, pymol_session)

    def get_pair_pymol_param(
            self, pair_id: int, protein_id: int, chain_letter: str, param_name: str,
    ) -> object:
        with self._conn() as db:
            return self._repo_pair_params.get(
                db, protein_id, chain_letter, pair_id, param_name
            )

    def update_pair_chain_color(
            self, protein_id: int, chain_letter: str, pair_id: int, color: str,
    ) -> None:
        from src.pyssa.util.enums import PymolParameterEnum
        with self._conn() as db:
            self._repo_pair_params.update_color(
                db, color, protein_id, chain_letter, pair_id,
                PymolParameterEnum.COLOR.value,
            )
