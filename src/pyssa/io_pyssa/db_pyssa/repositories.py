# repositories.py
"""
Repository classes: one class per database entity, all SQL centralised here.

Each repository is *stateless* – it receives a db connection on every call,
which makes the hot/cold project model trivial: just pass the right
connection pool's db handle.

Repositories return plain data (dicts, namedtuples, numpy arrays) so the
callers can build domain objects without the repository knowing about them.
"""
from __future__ import annotations

import logging
from typing import Any

import numpy as np
from src.pyssa.gui.qt import QtSql

from src.pyssa.io_pyssa.db_pyssa.queries import run, scalar, rows, last_insert_id, prepare, execute

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# SQL constants  (all SQL lives here – no scattered string literals elsewhere)
# ---------------------------------------------------------------------------

class _SQL:
    # Project
    GET_PROJECT_ID          = "SELECT id FROM Project WHERE name = ?"
    INSERT_PROJECT          = "INSERT INTO Project (name, os) VALUES (?, ?)"
    UPDATE_PROJECT_NAME     = "UPDATE Project SET name = ? WHERE id = ?"

    # SeqRecord
    GET_SEQUENCES           = "SELECT id, seq, name FROM SeqRecord WHERE project_id = ?"
    INSERT_SEQUENCE         = "INSERT INTO SeqRecord (seq_id, seq, name, project_id) VALUES (?, ?, ?, ?)"
    DELETE_SEQUENCE         = "DELETE FROM SeqRecord WHERE name = ?"
    UPDATE_SEQUENCE_NAME    = "UPDATE SeqRecord SET name = ? WHERE name = ? AND seq = ?"

    # Protein
    GET_PROTEINS            = "SELECT id, pymol_molecule_object, pymol_session FROM Protein WHERE project_id = ?"
    GET_PROTEIN_BY_NAME     = "SELECT id, pymol_molecule_object, pymol_session FROM Protein WHERE pymol_molecule_object = ?"
    GET_PROTEIN_NAME_BY_ID  = "SELECT pymol_molecule_object FROM Protein WHERE id = ?"
    GET_LATEST_PROTEIN_ID   = "SELECT MAX(id) FROM Protein"
    INSERT_PROTEIN          = "INSERT INTO Protein (pymol_molecule_object, pymol_session, project_id) VALUES (?, ?, ?)"
    DELETE_PROTEIN          = "DELETE FROM Protein WHERE id = ?"
    UPDATE_PROTEIN_NAME     = "UPDATE Protein SET pymol_molecule_object = ? WHERE pymol_molecule_object = ? AND id = ?"
    UPDATE_PROTEIN_SESSION  = "UPDATE Protein SET pymol_session = ? WHERE id = ?"

    # Chain
    GET_CHAINS              = "SELECT id, chain_identifier, chain_type, chain_sequence FROM Chain WHERE protein_id = ?"
    GET_LATEST_CHAIN_ID     = "SELECT MAX(id) FROM Chain"
    INSERT_CHAIN            = "INSERT INTO Chain (protein_id, chain_identifier, chain_type, chain_sequence) VALUES (?, ?, ?, ?)"
    DELETE_CHAINS           = "DELETE FROM Chain WHERE protein_id = ?"
    DELETE_SPECIFIC_CHAIN   = "DELETE FROM Chain WHERE protein_id = ? AND id = ?"

    # PyMOLParameter (per chain, for plain proteins)
    GET_PYMOL_PARAMS        = "SELECT color, representation FROM PyMOLParameter WHERE chain_id = ?"
    INSERT_PYMOL_PARAM      = "INSERT INTO PyMOLParameter (color, representation, chain_id) VALUES (?, ?, ?)"
    DELETE_PYMOL_PARAM      = "DELETE FROM PyMOLParameter WHERE chain_id = ?"
    UPDATE_CHAIN_COLOR      = "UPDATE PyMOLParameter SET color = ? WHERE chain_id = ?"

    # PyMOLSelection
    GET_PYMOL_SELECTION     = "SELECT selection_string FROM PyMOLSelection WHERE protein_id = ?"
    INSERT_PYMOL_SELECTION  = "INSERT INTO PyMOLSelection (selection_string, protein_id) VALUES (?, ?)"
    DELETE_PYMOL_SELECTION  = "DELETE FROM PyMOLSelection WHERE protein_id = ?"

    # PdbAtom
    GET_PDB_ATOMS           = """
        SELECT record_type, atom_number, atom_name, alternate_location_indicator,
               residue_name, chain_identifier, residue_sequence_number,
               code_for_insertions_of_residues, x_coord, y_coord, z_coord,
               occupancy, temperature_factor, segment_identifier,
               element_symbol, charge
        FROM PdbAtom WHERE protein_id = ?
    """
    INSERT_PDB_ATOM         = """
        INSERT INTO PdbAtom (
            record_type, atom_number, atom_name, alternate_location_indicator,
            residue_name, chain_identifier, residue_sequence_number,
            code_for_insertions_of_residues, x_coord, y_coord, z_coord,
            occupancy, temperature_factor, segment_identifier,
            element_symbol, charge, protein_id
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """
    DELETE_PDB_ATOMS        = "DELETE FROM PdbAtom WHERE protein_id = ?"

    # ProteinPair
    GET_PROTEIN_PAIRS       = """
        SELECT id, protein_1_id, protein_2_id, pymol_session, name
        FROM ProteinPair WHERE project_id = ?
    """
    GET_PROTEIN_PAIR_BY_NAME = """
        SELECT id, protein_1_id, protein_2_id, pymol_session, name
        FROM ProteinPair WHERE name = ?
    """
    INSERT_PROTEIN_PAIR     = """
        INSERT INTO ProteinPair (protein_1_id, protein_2_id, pymol_session, project_id, name)
        VALUES (?, ?, ?, ?, ?)
    """
    DELETE_PROTEIN_PAIR     = "DELETE FROM ProteinPair WHERE id = ?"
    UPDATE_PAIR_SESSION     = "UPDATE ProteinPair SET pymol_session = ? WHERE id = ?"

    # PyMOLParameterProteinPair
    GET_PAIR_PYMOL_PARAM    = """
        SELECT parameter_value FROM PyMOLParameterProteinPair
        WHERE protein_id = ? AND chain_letter = ? AND protein_pair_id = ? AND parameter_name = ?
    """
    INSERT_PAIR_PYMOL_PARAM = """
        INSERT INTO PyMOLParameterProteinPair
            (protein_id, chain_letter, parameter_name, parameter_value, protein_pair_id)
        VALUES (?, ?, ?, ?, ?)
    """
    DELETE_PAIR_PYMOL_PARAMS = "DELETE FROM PyMOLParameterProteinPair WHERE protein_pair_id = ?"
    UPDATE_PAIR_CHAIN_COLOR  = """
        UPDATE PyMOLParameterProteinPair
        SET parameter_value = ?
        WHERE protein_id = ? AND chain_letter = ? AND protein_pair_id = ? AND parameter_name = ?
    """

    # DistanceAnalysis
    GET_DISTANCE_ANALYSIS   = """
        SELECT id, name, cutoff, cycles, figure_size_x, figure_size_y
        FROM DistanceAnalysis WHERE protein_pair_id = ?
    """
    INSERT_DISTANCE_ANALYSIS = """
        INSERT INTO DistanceAnalysis (name, cutoff, cycles, protein_pair_id, figure_size_x, figure_size_y)
        VALUES (?, ?, ?, ?, ?, ?)
    """
    DELETE_DISTANCE_ANALYSIS = "DELETE FROM DistanceAnalysis WHERE protein_pair_id = ?"

    # DistanceAnalysisResults
    GET_DISTANCE_RESULTS    = """
        SELECT id, pymol_session, rmsd, aligned_aa
        FROM DistanceAnalysisResults WHERE distance_analysis_id = ?
    """
    INSERT_DISTANCE_RESULTS = """
        INSERT INTO DistanceAnalysisResults (pymol_session, rmsd, aligned_aa, distance_analysis_id)
        VALUES (?, ?, ?, ?)
    """
    DELETE_DISTANCE_RESULTS = "DELETE FROM DistanceAnalysisResults WHERE distance_analysis_id = ?"

    # DistanceAnalysisResultData
    GET_DISTANCE_DATA       = """
        SELECT my_index, protein_1_chain, protein_1_position, protein_1_residue,
               protein_2_chain, protein_2_position, protein_2_residue, distances
        FROM DistanceAnalysisResultData WHERE distance_analysis_results_id = ?
    """
    INSERT_DISTANCE_DATA    = """
        INSERT INTO DistanceAnalysisResultData (
            my_index, protein_1_chain, protein_1_position, protein_1_residue,
            protein_2_chain, protein_2_position, protein_2_residue,
            distances, distance_analysis_results_id
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
    """
    DELETE_DISTANCE_DATA    = "DELETE FROM DistanceAnalysisResultData WHERE distance_analysis_results_id = ?"


# ---------------------------------------------------------------------------
# Repository classes
# ---------------------------------------------------------------------------

class ProjectRepository:
    def get_id(self, db: QtSql.QSqlDatabase, name: str) -> int | None:
        logger.debug("Getting project ID for name='%s'", name)
        result = scalar(db, _SQL.GET_PROJECT_ID, name)
        logger.debug("Project ID for '%s': %s", name, result)
        return result

    def insert(self, db: QtSql.QSqlDatabase, name: str, os: str) -> int:
        logger.debug("Inserting project: name='%s', os='%s'", name, os)
        run(db, _SQL.INSERT_PROJECT, name, os)
        project_id = last_insert_id(db)
        logger.info("Inserted project '%s' with id=%d", name, project_id)
        return project_id

    def update_name(self, db: QtSql.QSqlDatabase, new_name: str, project_id: int) -> None:
        logger.debug("Updating project id=%d to name='%s'", project_id, new_name)
        run(db, _SQL.UPDATE_PROJECT_NAME, new_name, project_id)
        logger.info("Updated project id=%d name", project_id)


class SequenceRepository:
    def get_all(self, db: QtSql.QSqlDatabase, project_id: int) -> list[dict]:
        return [
            {"id": r[0], "seq": r[1], "name": r[2]}
            for r in rows(db, _SQL.GET_SEQUENCES, project_id)
        ]

    def insert(self, db: QtSql.QSqlDatabase, seq_id: str, seq: str, name: str, project_id: int) -> None:
        run(db, _SQL.INSERT_SEQUENCE, seq_id, seq, name, project_id)

    def delete(self, db: QtSql.QSqlDatabase, name: str) -> None:
        run(db, _SQL.DELETE_SEQUENCE, name)

    def update_name(self, db: QtSql.QSqlDatabase, new_name: str, old_name: str, sequence: str) -> None:
        run(db, _SQL.UPDATE_SEQUENCE_NAME, new_name, old_name, sequence)


class ProteinRepository:
    def get_all(self, db: QtSql.QSqlDatabase, project_id: int) -> list[dict]:
        logger.debug("Getting all proteins for project_id=%d", project_id)
        proteins = [
            {"id": r[0], "name": r[1], "pymol_session": r[2]}
            for r in rows(db, _SQL.GET_PROTEINS, project_id)
        ]
        logger.debug("Found %d proteins for project_id=%d", len(proteins), project_id)
        return proteins

    def get_by_name(self, db: QtSql.QSqlDatabase, name: str) -> dict | None:
        logger.debug("Getting protein by name='%s'", name)
        result = rows(db, _SQL.GET_PROTEIN_BY_NAME, name)
        if not result:
            logger.debug("Protein '%s' not found", name)
            return None
        r = result[0]
        return {"id": r[0], "name": r[1], "pymol_session": r[2]}

    def get_name_by_id(self, db: QtSql.QSqlDatabase, protein_id: int) -> str | None:
        logger.debug("Getting protein name for id=%d", protein_id)
        return scalar(db, _SQL.GET_PROTEIN_NAME_BY_ID, protein_id)

    def next_id(self, db: QtSql.QSqlDatabase) -> int:
        latest = scalar(db, _SQL.GET_LATEST_PROTEIN_ID)
        return (latest or 0) + 1

    def insert(self, db: QtSql.QSqlDatabase, name: str, pymol_session: str, project_id: int) -> int:
        logger.debug("Inserting protein '%s' for project_id=%d", name, project_id)
        run(db, _SQL.INSERT_PROTEIN, name, pymol_session, project_id)
        protein_id = last_insert_id(db)
        logger.info("Inserted protein '%s' with id=%d", name, protein_id)
        return protein_id

    def delete(self, db: QtSql.QSqlDatabase, protein_id: int) -> None:
        logger.debug("Deleting protein id=%d", protein_id)
        run(db, _SQL.DELETE_PROTEIN, protein_id)
        logger.info("Deleted protein id=%d", protein_id)

    def update_name(self, db: QtSql.QSqlDatabase, new_name: str, old_name: str, protein_id: int) -> None:
        logger.debug("Updating protein id=%d from '%s' to '%s'", protein_id, old_name, new_name)
        run(db, _SQL.UPDATE_PROTEIN_NAME, new_name, old_name, protein_id)
        logger.info("Updated protein id=%d name", protein_id)

    def update_session(self, db: QtSql.QSqlDatabase, protein_id: int, pymol_session: str) -> None:
        session_size = len(pymol_session) if pymol_session else 0
        logger.debug("Updating PyMOL session for protein id=%d (size=%d bytes)", protein_id, session_size)
        run(db, _SQL.UPDATE_PROTEIN_SESSION, pymol_session, protein_id)
        logger.info("Updated PyMOL session for protein id=%d", protein_id)


class ChainRepository:
    def get_all(self, db: QtSql.QSqlDatabase, protein_id: int) -> list[dict]:
        return [
            {
                "id": r[0],
                "chain_identifier": r[1],
                "chain_type": r[2],
                "chain_sequence": r[3],
            }
            for r in rows(db, _SQL.GET_CHAINS, protein_id)
        ]

    def insert(
        self, db: QtSql.QSqlDatabase,
        protein_id: int, identifier: str, chain_type: str, sequence: str,
    ) -> int:
        run(db, _SQL.INSERT_CHAIN, protein_id, identifier, chain_type, sequence)
        return last_insert_id(db)

    def delete_all(self, db: QtSql.QSqlDatabase, protein_id: int) -> None:
        run(db, _SQL.DELETE_CHAINS, protein_id)

    def delete_one(self, db: QtSql.QSqlDatabase, protein_id: int, chain_id: int) -> None:
        run(db, _SQL.DELETE_SPECIFIC_CHAIN, protein_id, chain_id)


class PyMOLParameterRepository:
    """PyMOL parameters attached to individual protein chains."""

    def get(self, db: QtSql.QSqlDatabase, chain_id: int) -> dict:
        result = rows(db, _SQL.GET_PYMOL_PARAMS, chain_id)
        if not result:
            return {}
        return {"color": result[0][0], "representation": result[0][1]}

    def insert(self, db: QtSql.QSqlDatabase, color: str, representation: str, chain_id: int) -> None:
        run(db, _SQL.INSERT_PYMOL_PARAM, color, representation, chain_id)

    def delete(self, db: QtSql.QSqlDatabase, chain_id: int) -> None:
        run(db, _SQL.DELETE_PYMOL_PARAM, chain_id)

    def update_color(self, db: QtSql.QSqlDatabase, chain_id: int, color: str) -> None:
        run(db, _SQL.UPDATE_CHAIN_COLOR, color, chain_id)


class PyMOLSelectionRepository:
    def get(self, db: QtSql.QSqlDatabase, protein_id: int) -> str | None:
        return scalar(db, _SQL.GET_PYMOL_SELECTION, protein_id)

    def insert(self, db: QtSql.QSqlDatabase, selection_string: str, protein_id: int) -> None:
        run(db, _SQL.INSERT_PYMOL_SELECTION, selection_string, protein_id)

    def delete(self, db: QtSql.QSqlDatabase, protein_id: int) -> None:
        run(db, _SQL.DELETE_PYMOL_SELECTION, protein_id)


class PdbAtomRepository:
    def get_all(self, db: QtSql.QSqlDatabase, protein_id: int) -> list[tuple]:
        return rows(db, _SQL.GET_PDB_ATOMS, protein_id)

    def insert_many(self, db: QtSql.QSqlDatabase, protein_id: int, atom_dicts: list[dict]) -> None:
        q = prepare(db, _SQL.INSERT_PDB_ATOM)
        for a in atom_dicts:
            execute(
                q,
                a["record_type"], a["atom_number"], a["atom_name"],
                a["alternate_location_indicator"], a["residue_name"],
                a["chain_identifier"], a["residue_sequence_number"],
                a["code_for_insertions_of_residues"],
                a["x_coord"], a["y_coord"], a["z_coord"],
                a["occupancy"], a["temperature_factor"],
                a["segment_identifier"], a["element_symbol"],
                a["charge"], protein_id,
            )

    def delete_all(self, db: QtSql.QSqlDatabase, protein_id: int) -> None:
        run(db, _SQL.DELETE_PDB_ATOMS, protein_id)

    def replace_all(self, db: QtSql.QSqlDatabase, protein_id: int, atom_dicts: list[dict]) -> None:
        self.delete_all(db, protein_id)
        self.insert_many(db, protein_id, atom_dicts)


class ProteinPairRepository:
    def get_all(self, db: QtSql.QSqlDatabase, project_id: int) -> list[dict]:
        logger.debug("Getting all protein pairs for project_id=%d", project_id)
        pairs = [
            {
                "id": r[0], "protein_1_id": r[1], "protein_2_id": r[2],
                "pymol_session": r[3], "name": r[4],
            }
            for r in rows(db, _SQL.GET_PROTEIN_PAIRS, project_id)
        ]
        logger.debug("Found %d protein pairs for project_id=%d", len(pairs), project_id)
        return pairs

    def get_by_name(self, db: QtSql.QSqlDatabase, name: str) -> dict | None:
        logger.debug("Getting protein pair by name='%s'", name)
        result = rows(db, _SQL.GET_PROTEIN_PAIR_BY_NAME, name)
        if not result:
            logger.debug("Protein pair '%s' not found", name)
            return None
        r = result[0]
        return {
            "id": r[0], "protein_1_id": r[1], "protein_2_id": r[2],
            "pymol_session": r[3], "name": r[4],
        }

    def insert(
        self, db: QtSql.QSqlDatabase,
        protein_1_id: int, protein_2_id: int,
        pymol_session: str, project_id: int, name: str,
    ) -> int:
        logger.debug(
            "Inserting protein pair '%s' (protein_1=%d, protein_2=%d) for project_id=%d",
            name, protein_1_id, protein_2_id, project_id
        )
        run(db, _SQL.INSERT_PROTEIN_PAIR, protein_1_id, protein_2_id, pymol_session, project_id, name)
        pair_id = last_insert_id(db)
        logger.info("Inserted protein pair '%s' with id=%d", name, pair_id)
        return pair_id

    def delete(self, db: QtSql.QSqlDatabase, pair_id: int) -> None:
        logger.debug("Deleting protein pair id=%d", pair_id)
        run(db, _SQL.DELETE_PROTEIN_PAIR, pair_id)
        logger.info("Deleted protein pair id=%d", pair_id)

    def update_session(self, db: QtSql.QSqlDatabase, pair_id: int, pymol_session: str) -> None:
        session_size = len(pymol_session) if pymol_session else 0
        logger.debug("Updating PyMOL session for pair id=%d (size=%d bytes)", pair_id, session_size)
        run(db, _SQL.UPDATE_PAIR_SESSION, pymol_session, pair_id)
        logger.info("Updated PyMOL session for pair id=%d", pair_id)


class ProteinPairPyMOLParamRepository:
    def get(
        self, db: QtSql.QSqlDatabase,
        protein_id: int, chain_letter: str, pair_id: int, param_name: str,
    ) -> Any:
        return scalar(db, _SQL.GET_PAIR_PYMOL_PARAM, protein_id, chain_letter, pair_id, param_name)

    def insert(
        self, db: QtSql.QSqlDatabase,
        protein_id: int, chain_letter: str, param_name: str, param_value: str, pair_id: int,
    ) -> None:
        run(db, _SQL.INSERT_PAIR_PYMOL_PARAM, protein_id, chain_letter, param_name, param_value, pair_id)

    def delete_all_for_pair(self, db: QtSql.QSqlDatabase, pair_id: int) -> None:
        run(db, _SQL.DELETE_PAIR_PYMOL_PARAMS, pair_id)

    def update_color(
        self, db: QtSql.QSqlDatabase,
        color: str, protein_id: int, chain_letter: str, pair_id: int, param_name: str,
    ) -> None:
        run(db, _SQL.UPDATE_PAIR_CHAIN_COLOR, color, protein_id, chain_letter, pair_id, param_name)


class DistanceAnalysisRepository:
    def get(self, db: QtSql.QSqlDatabase, pair_id: int) -> dict | None:
        result = rows(db, _SQL.GET_DISTANCE_ANALYSIS, pair_id)
        if not result:
            return None
        r = result[0]
        return {
            "id": r[0], "name": r[1], "cutoff": r[2],
            "cycles": r[3], "figure_size_x": r[4], "figure_size_y": r[5],
        }

    def insert(
        self, db: QtSql.QSqlDatabase,
        name: str, cutoff: float, cycles: int, pair_id: int,
        figure_size_x: float, figure_size_y: float,
    ) -> int:
        run(db, _SQL.INSERT_DISTANCE_ANALYSIS, name, cutoff, cycles, pair_id, figure_size_x, figure_size_y)
        return last_insert_id(db)

    def delete(self, db: QtSql.QSqlDatabase, pair_id: int) -> None:
        run(db, _SQL.DELETE_DISTANCE_ANALYSIS, pair_id)


class DistanceAnalysisResultsRepository:
    def get(self, db: QtSql.QSqlDatabase, analysis_id: int) -> dict | None:
        result = rows(db, _SQL.GET_DISTANCE_RESULTS, analysis_id)
        if not result:
            return None
        r = result[0]
        return {"id": r[0], "pymol_session": r[1], "rmsd": r[2], "aligned_aa": r[3]}

    def insert(
        self, db: QtSql.QSqlDatabase,
        pymol_session: str, rmsd: float, aligned_aa: str, analysis_id: int,
    ) -> int:
        run(db, _SQL.INSERT_DISTANCE_RESULTS, pymol_session, rmsd, aligned_aa, analysis_id)
        return last_insert_id(db)

    def delete(self, db: QtSql.QSqlDatabase, analysis_id: int) -> None:
        run(db, _SQL.DELETE_DISTANCE_RESULTS, analysis_id)


class DistanceDataRepository:
    def get_as_arrays(self, db: QtSql.QSqlDatabase, results_id: int) -> dict[str, np.ndarray]:
        data = rows(db, _SQL.GET_DISTANCE_DATA, results_id)
        if not data:
            empty = np.array([])
            return {k: empty for k in (
                "index", "prot_1_chain", "prot_1_position", "prot_1_residue",
                "prot_2_chain", "prot_2_position", "prot_2_residue", "distances",
            )}
        transposed = list(zip(*data))
        keys = (
            "index", "prot_1_chain", "prot_1_position", "prot_1_residue",
            "prot_2_chain", "prot_2_position", "prot_2_residue", "distances",
        )
        return {k: np.array(v) for k, v in zip(keys, transposed)}

    def insert_many(self, db: QtSql.QSqlDatabase, results_id: int, data: dict[str, np.ndarray]) -> None:
        q = prepare(db, _SQL.INSERT_DISTANCE_DATA)
        indices       = data["index"].tolist()
        prot1_chains  = data["prot_1_chain"].tolist()
        prot1_pos     = data["prot_1_position"].tolist()
        prot1_res     = data["prot_1_residue"].tolist()
        prot2_chains  = data["prot_2_chain"].tolist()
        prot2_pos     = data["prot_2_position"].tolist()
        prot2_res     = data["prot_2_residue"].tolist()
        distances     = data["distances"].tolist()
        for row in zip(indices, prot1_chains, prot1_pos, prot1_res,
                       prot2_chains, prot2_pos, prot2_res, distances):
            execute(q, *row, results_id)

    def delete(self, db: QtSql.QSqlDatabase, results_id: int) -> None:
        run(db, _SQL.DELETE_DISTANCE_DATA, results_id)
