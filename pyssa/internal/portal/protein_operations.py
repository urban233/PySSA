#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2022
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/urban233/PySSA>
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
"""Module for protein operations in pymol."""
import logging
import pathlib
import pymol
from typing import TYPE_CHECKING
from pymol import cmd

from pyssa.internal.portal import pymol_safeguard
from pyssa.internal.portal import pymol_io
from pyssa.internal.data_structures import chain
from pyssa.io_pyssa import safeguard
from pyssa.util import constants
from pyssa.util import protein_util
from pyssa.logging_pyssa import log_handlers
from pyssa.internal.data_structures import sequence

if TYPE_CHECKING:
    from pyssa.internal.data_structures import protein

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def remove_solvent_molecules_in_protein() -> None:
    """Removes solvent molecules in a protein."""
    if not pymol_safeguard.PymolSafeguard.check_if_protein_in_session():
        raise pymol.CmdException("No protein is in pymol session.")
    try:
        cmd.remove("solvent")
    except pymol.CmdException:
        print("No solvent molecules needs to be removed.")


def remove_organic_molecules_in_protein() -> None:
    """Removes organic molecules in a protein."""
    if not pymol_safeguard.PymolSafeguard.check_if_protein_in_session():
        raise pymol.CmdException("No protein is in pymol session.")
    try:
        cmd.remove("organic")
    except pymol.CmdException:
        print("No organic molecules needs to be removed.")


def get_protein_chains(molecule_object: str, dirname: pathlib.Path, basename: str) -> list[chain.Chain]:
    """This function divides the chains from a protein, into protein and non-protein chains.

    Args:
        molecule_object:
            the name of the protein which is also used within pymol
        dirname:
             the filepath where the pdb file is stored
        basename:
             the name of the file with extension

    Returns:
        a list of all chains from the protein, divided into protein and non-protein chains
    """
    # <editor-fold desc="Checks">
    if not safeguard.Safeguard.check_if_value_is_not_none(molecule_object) or molecule_object == "":
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")
    if not safeguard.Safeguard.check_if_value_is_not_none(dirname):
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")
    if not safeguard.Safeguard.check_filepath(dirname):
        logger.error("The directory does not exist.")
        raise NotADirectoryError("The directory does not exist.")
    if not safeguard.Safeguard.check_if_value_is_not_none(basename) or basename == "":
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")

    # </editor-fold>

    cmd.reinitialize()
    pymol_io.load_protein(dirname, basename, molecule_object)
    tmp_chains: list[str] = cmd.get_chains()
    i = 0
    chains_of_protein: list[chain.Chain] = []
    for tmp_chain in tmp_chains:
        sequence_of_chain = cmd.get_model(f"chain {tmp_chain}")
        if sequence_of_chain.atom[0].resn in constants.AMINO_ACID_CODE:
            # TODO: can this be better?
            fasta_sequence_of_chain = cmd.get_fastastr(f"chain {tmp_chain}")
            fasta_sequence_of_chain_without_header = fasta_sequence_of_chain[fasta_sequence_of_chain.find("\n") :]
            complete_sequence_of_chain = sequence.Sequence(
                molecule_object,
                fasta_sequence_of_chain_without_header.replace("\n", ""),
            )
            chains_of_protein.append(chain.Chain(tmp_chain, complete_sequence_of_chain, constants.CHAIN_TYPE_PROTEIN))
        else:
            # TODO: this should produce a sequence with the non-protein atoms
            fasta_sequence_of_chain = cmd.get_fastastr(f"chain {tmp_chain}")
            fasta_sequence_of_chain_without_header = fasta_sequence_of_chain[fasta_sequence_of_chain.find("\n") :]
            complete_sequence_of_chain = sequence.Sequence(
                molecule_object,
                fasta_sequence_of_chain_without_header.replace("\n", ""),
            )
            chains_of_protein.append(
                chain.Chain(tmp_chain, complete_sequence_of_chain, constants.CHAIN_TYPE_NON_PROTEIN),
            )
        i += 1
    return chains_of_protein


def get_protein_sequences_from_protein(molecule_object: str, chains: list[chain.Chain]) -> list["sequence.Sequence"]:
    """This function gets all sequences from protein chains only.

    Args:
        molecule_object:
            the name of the protein which is also used within pymol
        chains:
            a list of chains which occur in the protein

    Returns:
        a protein sequence object with amino acid sequences
    """
    # <editor-fold desc="Checks">
    if not safeguard.Safeguard.check_if_value_is_not_none(molecule_object) or molecule_object == "":
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")
    if not safeguard.Safeguard.check_if_value_is_not_none(chains) or not safeguard.Safeguard.check_if_list_is_empty(
        chains,
    ):
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")

    # </editor-fold>
    protein_sequences: list[sequence.Sequence] = []
    for tmp_chain in protein_util.filter_chains_for_protein_chains(chains):
        protein_sequences.append(tmp_chain.chain_sequence)
    return protein_sequences


def get_protein_sequence_length_from_protein(a_protein: "protein.Protein") -> int:
    """Gets the length of the protein sequence.

    Args:
        a_protein: a protein object of which the sequence length is calculated.
    """
    fasta_prot_1 = cmd.get_fastastr(a_protein.pymol_selection.selection_string)
    return len(fasta_prot_1[fasta_prot_1.find("\n") :])
