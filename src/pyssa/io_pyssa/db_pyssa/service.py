# service.py
"""
ProjectService: orchestrates multi-repository operations that must execute
in a single transaction, and handles the mapping between domain objects and
raw repository data.

This layer knows about domain objects (Protein, ProteinPair, etc.) so that
ProjectDatabase doesn't have to import them.
"""
from __future__ import annotations

import logging
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .project_database import ProjectDatabase
    from src.pyssa.internal.data_structures import (
        protein, protein_pair, structure_analysis, settings, project
    )
    from src.pyssa.internal.thread.async_pyssa import custom_signals

from src.pyssa.util import enums, pyssa_keys
from src.pyssa.internal.data_structures import (
    chain as chain_module,
    sequence as sequence_module,
    results as results_module,
    structure_analysis as analysis_module,
)
from Bio import SeqRecord
from Bio.Seq import Seq

logger = logging.getLogger(__name__)


class ProjectService:
    """Coordinates multi-repository transactions and object assembly.

    Receives a ProjectDatabase (not the pool directly) so it can use the
    existing transaction() helper and call repo methods through db.
    """

    def __init__(self, db: "ProjectDatabase") -> None:
        self._db = db

    # ------------------------------------------------------------------
    # Full project load
    # ------------------------------------------------------------------

    def load_project(
        self,
        project_name: str,
        workspace_path: object,
        app_settings: "settings.Settings",
        progress_signal: "custom_signals.ProgressSignal | None" = None,
    ) -> "project.Project":
        """Reconstruct a full Project domain object from the database."""
        from src.pyssa.internal.data_structures import project as project_module

        def _emit(msg: str, pct: int) -> None:
            if progress_signal is not None:
                progress_signal.emit_signal(msg, pct)

        tmp_project = project_module.Project(project_name, workspace_path)

        project_id = self._db.get_project_id(project_name)
        if project_id is None:
            raise ValueError(f"Project '{project_name}' not found in database.")
        tmp_project.set_id(project_id)

        _emit("Loading sequences …", 10)
        tmp_project.sequences = self._load_sequences(project_id)

        _emit("Loading proteins …", 25)
        tmp_project.proteins = self._load_proteins(project_id)

        _emit("Loading protein pairs …", 50)
        tmp_project.protein_pairs = self._load_protein_pairs(
            project_id, tmp_project, app_settings
        )

        _emit("Loading project data finished.", 60)
        return tmp_project

    # ------------------------------------------------------------------
    # Sequence assembly
    # ------------------------------------------------------------------

    def _load_sequences(self, project_id: int) -> list[SeqRecord.SeqRecord]:
        """Load all sequences for a project from the database.

        Args:
            project_id: The database ID of the project.

        Returns:
            A list of fully constructed ``SeqRecord`` objects, with the
            sequence wrapped in a ``Bio.Seq.Seq`` instance and the
            sequence name used as both ``id`` and ``name``.
        """
        raw = self._db.get_all_sequences(project_id)
        return [
            SeqRecord.SeqRecord(Seq(r["seq"]), id=r["name"], name=r["name"])
            for r in raw
        ]

    # ------------------------------------------------------------------
    # Protein assembly
    # ------------------------------------------------------------------

    def _load_proteins(self, project_id: int) -> list["protein.Protein"]:
        from src.pyssa.internal.data_structures import protein as protein_module
        proteins = []
        for raw in self._db.get_all_proteins(project_id):
            proteins.append(self._build_protein(raw))
        return proteins

    def _build_protein(self, raw: dict) -> "protein.Protein":
        """Reconstruct a single Protein domain object from raw DB data.

        Args:
            raw: A dict with keys ``id``, ``name``, and ``pymol_session``
                 as returned by :class:`ProteinRepository`.

        Returns:
            A fully populated :class:`Protein` instance, including chains,
            PyMOL selection string, and PDB atom data loaded from the DB.
        """
        from src.pyssa.internal.data_structures import protein as protein_module
        p = protein_module.Protein(raw["name"])
        p.set_id(raw["id"])
        p.pymol_session = raw["pymol_session"]
        p.chains = self._load_chains(raw["id"], raw["name"])
        p.pymol_selection.selection_string = self._load_selection(raw["id"])
        p.set_pdb_data(self._load_pdb_atoms(raw["id"]))
        return p

    def _load_chains(self, protein_id: int, protein_name: str) -> list:
        with self._db._conn() as db:
            raw_chains = self._db._repo_chains.get_all(db, protein_id)
            chains = []
            for rc in raw_chains:
                seq = sequence_module.Sequence(protein_name, rc["chain_sequence"])
                c = chain_module.Chain(rc["chain_identifier"], seq, rc["chain_type"])
                c.set_id(rc["id"])
                c.db_protein_id = protein_id
                c.pymol_parameters = self._db._repo_pymol_params.get(db, rc["id"])
                chains.append(c)
        return chains

    def _load_pdb_atoms(self, protein_id: int) -> list[dict]:
        """Load and convert PDB atom rows for one protein from the database.

        The repository returns raw tuples whose column order matches the
        ``SELECT`` column list in :data:`_SQL.GET_PDB_ATOMS`:
        ``record_type, atom_number, atom_name, alternate_location_indicator,
        residue_name, chain_identifier, residue_sequence_number,
        code_for_insertions_of_residues, x_coord, y_coord, z_coord,
        occupancy, temperature_factor, segment_identifier,
        element_symbol, charge``.

        Args:
            protein_id: The database ID of the protein.

        Returns:
            A list of atom dicts ready for use with :meth:`Protein.set_pdb_data`.
        """
        _FIELDS = (
            "record_type", "atom_number", "atom_name",
            "alternate_location_indicator", "residue_name",
            "chain_identifier", "residue_sequence_number",
            "code_for_insertions_of_residues",
            "x_coord", "y_coord", "z_coord",
            "occupancy", "temperature_factor",
            "segment_identifier", "element_symbol", "charge",
        )
        with self._db._conn() as db:
            raw_rows = self._db._repo_pdb_atoms.get_all(db, protein_id)
        return [dict(zip(_FIELDS, row)) for row in raw_rows]

    def _load_selection(self, protein_id: int) -> str:
        with self._db._conn() as db:
            sel_str = self._db._repo_pymol_sel.get(db, protein_id)
        return sel_str if sel_str is not None else ""

    # ------------------------------------------------------------------
    # Protein pair assembly
    # ------------------------------------------------------------------

    def _load_protein_pairs(
        self,
        project_id: int,
        project_obj: "project.Project",
        app_settings: "settings.Settings",
    ) -> list["protein_pair.ProteinPair"]:
        from src.pyssa.internal.data_structures import protein_pair as pp_module
        pairs = []
        for raw in self._db.get_all_protein_pairs(project_id):
            prot1_name = self._db.get_protein_name_by_id(raw["protein_1_id"])
            prot2_name = self._db.get_protein_name_by_id(raw["protein_2_id"])
            pair = pp_module.ProteinPair(
                project_obj.search_protein(prot1_name),
                project_obj.search_protein(prot2_name),
            )
            pair.set_id(raw["id"])
            pair.db_project_id = project_obj.get_id()
            pair.name = raw["name"]
            pair.pymol_session = raw["pymol_session"]
            pair.distance_analysis = self._load_distance_analysis(raw["id"], app_settings)
            pairs.append(pair)
        return pairs

    def _load_distance_analysis(
        self, pair_id: int, app_settings: "settings.Settings"
    ) -> "structure_analysis.DistanceAnalysis":
        with self._db._conn() as db:
            raw = self._db._repo_dist_anal.get(db, pair_id)
        if raw is None:
            return analysis_module.DistanceAnalysis(app_settings)

        da = analysis_module.DistanceAnalysis(app_settings)
        da.name = raw["name"]
        da.cutoff = raw["cutoff"]
        da.cycles = raw["cycles"]
        da.figure_size = (raw["figure_size_x"], raw["figure_size_y"])
        da.analysis_results = self._load_distance_results(raw["id"])
        return da

    def _load_distance_results(self, analysis_id: int) -> "results_module.DistanceAnalysisResults":
        with self._db._conn() as db:
            raw = self._db._repo_dist_res.get(db, analysis_id)
            if raw is None:
                return results_module.DistanceAnalysisResults(None, None, None, None)

            distance_data = self._db._repo_dist_data.get_as_arrays(db, raw["id"])

        return results_module.DistanceAnalysisResults(
            distance_data, raw["pymol_session"], raw["rmsd"], raw["aligned_aa"]
        )

    # ------------------------------------------------------------------
    # Full protein insert (transactional)
    # ------------------------------------------------------------------

    def insert_protein(self, protein_obj: "protein.Protein") -> int:
        d = protein_obj.get_object_as_dict_for_database()
        with self._db.transaction() as tx:
            protein_id = self._db._repo_proteins.insert(
                tx,
                d[enums.DatabaseEnum.PROTEIN_NAME.value],
                d[enums.DatabaseEnum.PROTEIN_PYMOL_SESSION.value],
                d["project_id"],
            )
            for c in protein_obj.chains:
                chain_id = self._db._repo_chains.insert(
                    tx,
                    protein_id,
                    c.chain_letter,
                    c.chain_type,
                    c.chain_sequence.sequence,
                )
                c.set_id(chain_id)
                self._db._repo_pymol_params.insert(
                    tx,
                    c.pymol_parameters[enums.PymolParameterEnum.COLOR.value],
                    c.pymol_parameters[enums.PymolParameterEnum.REPRESENTATION.value],
                    chain_id,
                )
            if protein_obj.chains:
                self._db._repo_pymol_sel.insert(
                    tx,
                    protein_obj.pymol_selection.selection_string,
                    protein_id,
                )
            self._db._repo_pdb_atoms.insert_many(tx, protein_id, protein_obj.get_pdb_data())
        return protein_id

    # ------------------------------------------------------------------
    # Full protein delete (transactional)
    # ------------------------------------------------------------------

    def delete_protein(self, protein_id: int) -> None:
        with self._db.transaction() as tx:
            # Delete PyMOL parameters for each chain before deleting chains
            for rc in self._db._repo_chains.get_all(tx, protein_id):
                self._db._repo_pymol_params.delete(tx, rc["id"])
            self._db._repo_pdb_atoms.delete_all(tx, protein_id)
            self._db._repo_pymol_sel.delete(tx, protein_id)
            self._db._repo_chains.delete_all(tx, protein_id)
            self._db._repo_proteins.delete(tx, protein_id)

    # ------------------------------------------------------------------
    # Full protein pair insert (transactional)
    # ------------------------------------------------------------------

    def insert_protein_pair(self, pair_obj: "protein_pair.ProteinPair") -> int:
        with self._db.transaction() as tx:
            pair_id = self._db._repo_pairs.insert(
                tx,
                pair_obj.protein_1.get_id(),
                pair_obj.protein_2.get_id(),
                pair_obj.pymol_session,
                pair_obj.db_project_id,
                pair_obj.name,
            )
            pair_obj.set_id(pair_id)

            for protein, chains in (
                (pair_obj.protein_1, pair_obj.protein_1.chains),
                (pair_obj.protein_2, pair_obj.protein_2.chains),
            ):
                for c in chains:
                    for param in (
                        enums.PymolParameterEnum.COLOR,
                        enums.PymolParameterEnum.REPRESENTATION,
                    ):
                        self._db._repo_pair_params.insert(
                            tx,
                            protein.get_id(),
                            c.chain_letter,
                            param.value,
                            # We pass an empty string as the value because these
                            # are not being used anymore but still exist in the
                            # DB for backward compatibility
                            "",
                            pair_id,
                        )

            da = pair_obj.distance_analysis
            analysis_id = self._db._repo_dist_anal.insert(
                tx,
                da.name, da.cutoff, da.cycles, pair_id,
                da.figure_size[0], da.figure_size[1],
            )
            results_id = self._db._repo_dist_res.insert(
                tx,
                da.analysis_results.pymol_session,
                da.analysis_results.rmsd,
                da.analysis_results.aligned_aa,
                analysis_id,
            )
            self._db._repo_dist_data.insert_many(
                tx, results_id, da.analysis_results.distance_data
            )
        return pair_id

    # ------------------------------------------------------------------
    # Full protein pair delete (transactional)
    # ------------------------------------------------------------------

    def delete_protein_pair(self, pair_id: int) -> None:
        with self._db.transaction() as tx:
            raw_da = self._db._repo_dist_anal.get(tx, pair_id)
            if raw_da:
                raw_res = self._db._repo_dist_res.get(tx, raw_da["id"])
                if raw_res:
                    self._db._repo_dist_data.delete(tx, raw_res["id"])
                self._db._repo_dist_res.delete(tx, raw_da["id"])
            self._db._repo_dist_anal.delete(tx, pair_id)
            self._db._repo_pair_params.delete_all_for_pair(tx, pair_id)
            self._db._repo_pairs.delete(tx, pair_id)

    # ------------------------------------------------------------------
    # Update operations
    # ------------------------------------------------------------------

    def update_protein_name(
        self, protein_id: int, new_name: str, old_name: str
    ) -> None:
        """Update the name of a protein.

        Args:
            protein_id: The database ID of the protein.
            new_name: The new name for the protein.
            old_name: The old name (used for validation).
        """
        with self._db._conn() as db:
            self._db._repo_proteins.update_name(db, new_name, old_name, protein_id)

    def update_protein_chain_color(self, chain_id: int, color: str) -> None:
        """Update the color of a specific chain.

        Args:
            chain_id: The database ID of the chain.
            color: The new color value.
        """
        with self._db._conn() as db:
            self._db._repo_pymol_params.update_color(db, chain_id, color)

    def update_protein_session(self, protein_id: int, pymol_session: str) -> None:
        """Update the PyMOL session for a protein.

        Args:
            protein_id: The database ID of the protein.
            pymol_session: The PyMOL session data.
        """
        with self._db._conn() as db:
            self._db._repo_proteins.update_session(db, protein_id, pymol_session)

    def update_protein_pair_session(
        self, pair_id: int, pymol_session: str
    ) -> None:
        """Update the PyMOL session for a protein pair.

        Args:
            pair_id: The database ID of the protein pair.
            pymol_session: The PyMOL session data.
        """
        with self._db._conn() as db:
            self._db._repo_pairs.update_session(db, pair_id, pymol_session)

    def update_sequence_name(
        self, new_name: str, old_name: str, sequence: str
    ) -> None:
        """Update the name of a sequence.

        Args:
            new_name: The new name for the sequence.
            old_name: The old name (used to find the sequence).
            sequence: The sequence string (additional validation).
        """
        with self._db._conn() as db:
            self._db._repo_sequences.update_name(db, new_name, old_name, sequence)

    def delete_specific_chain(self, protein_id: int, chain_id: int) -> None:
        """Delete a specific chain from a protein.

        Args:
            protein_id: The database ID of the protein.
            chain_id: The database ID of the chain to delete.
        """
        with self._db.transaction() as tx:
            # Delete PyMOL parameters first (foreign key dependency)
            self._db._repo_pymol_params.delete(tx, chain_id)
            # Then delete the chain
            self._db._repo_chains.delete_one(tx, protein_id, chain_id)

    def update_protein_pdb_atoms(
        self, protein_id: int, atom_dicts: list[dict]
    ) -> None:
        """Replace all PDB atom data for a protein.

        Args:
            protein_id: The database ID of the protein.
            atom_dicts: List of dictionaries containing atom data.
        """
        with self._db.transaction() as tx:
            self._db._repo_pdb_atoms.replace_all(tx, protein_id, atom_dicts)
