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
"""Module for protein operations in pymol"""
import logging
import pathlib

import pymol
from pymol import cmd
from pyssa.internal.portal import pymol_safeguard
from pyssa.internal.portal import pymol_io
from pyssa.internal.data_structures import chain
from pyssa.internal.data_structures import sequence
from pyssa.io_pyssa import safeguard
from pyssa.util import constants
from pyssa.util import types
from pyssa.util import protein_util
from pyssa.logging_pyssa import log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def remove_solvent_molecules_in_protein():
    """This function removes solvent molecules in a protein.

    """
    if not pymol_safeguard.PymolSafeguard.check_if_protein_in_session():
        raise pymol.CmdException("No protein is in pymol session.")
    try:
        cmd.remove("solvent")
    except pymol.CmdException:
        print("No solvent molecules needs to be removed.")


def remove_organic_molecules_in_protein():
    """This function removes organic molecules in a protein.

    """
    if not pymol_safeguard.PymolSafeguard.check_if_protein_in_session():
        raise pymol.CmdException("No protein is in pymol session.")
    try:
        cmd.remove("organic")
    except pymol.CmdException:
        print("No organic molecules needs to be removed.")


def get_protein_chains(molecule_object: str, filepath: pathlib.Path, filename: str) -> list[types.CHAIN]:
    """This function divides the chains from a protein, into protein and non-protein chains.

    Args:
        molecule_object:
            the name of the protein which is also used within pymol
        filepath:
             the filepath where the pdb file is stored
        filename:
             the name of the file with extension

    Returns:
        a list of all chains from the protein, divided into protein and non-protein chains
    """
    # <editor-fold desc="Checks">
    if not safeguard.Safeguard.check_if_value_is_not_none(molecule_object) or molecule_object == "":
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")
    if not safeguard.Safeguard.check_if_value_is_not_none(filepath):
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")
    if not safeguard.Safeguard.check_filepath(filepath):
        logger.error("The directory does not exist.")
        raise NotADirectoryError("The directory does not exist.")
    if not safeguard.Safeguard.check_if_value_is_not_none(filename) or filename == "":
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")

    # </editor-fold>

    pymol_io.load_protein(filepath, filename, molecule_object)
    tmp_chains: list[str] = cmd.get_chains()
    i = 0
    chains_of_protein: list[types.CHAIN] = []
    for tmp_chain in tmp_chains:
        molecules_of_chain = cmd.get_model(f"chain {tmp_chain}")
        if molecules_of_chain.atom[0].resn in constants.AMINO_ACID_CODE:
            chains_of_protein.append(chain.Chain(tmp_chain, molecules_of_chain, constants.CHAIN_TYPE_PROTEIN))
        else:
            chains_of_protein.append(chain.Chain(tmp_chain, molecules_of_chain, constants.CHAIN_TYPE_NON_PROTEIN))
        i += 1
    return chains_of_protein


def get_protein_sequences_from_protein(molecule_object, chains: list[types.CHAIN]) -> types.PROTEIN_SEQUENCE:
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
    if not safeguard.Safeguard.check_if_value_is_not_none(chains) or not safeguard.Safeguard.check_if_list_is_empty(chains):
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")

    # </editor-fold>
    protein_sequences = []
    for tmp_chain in protein_util.filter_chains_for_protein_chains(chains):
        protein_sequences.append(tmp_chain.sequence)
    return sequence.ProteinSequence(molecule_object, protein_sequences)
