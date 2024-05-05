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
"""Module contains helper functions for the protein class."""
import logging
from pyssa.io_pyssa import safeguard
from pyssa.internal.data_structures import chain
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pyssa.internal.data_structures import sequence
    from pyssa.internal.data_structures import protein


logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def check_if_protein_is_from_file_or_id(molecule_object: str) -> tuple:
    """This function checks if a protein is from a .pdb file or from a PDB ID.

    Args:
        molecule_object:
            the name of the protein which is also used within pymol
    Returns:
        a tuple with the molecule object and the basename
    """
    # <editor-fold desc="Checks">
    safeguard.Safeguard.check_if_value_is_not_none(molecule_object, logger)
    if molecule_object == "":
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")

    # </editor-fold>

    if molecule_object.find(".pdb") != -1:
        # pdb file is given
        tmp_molecule_object = molecule_object.replace(".pdb", "")
        tmp_basename = molecule_object
    else:
        # PDB ID is given
        tmp_basename = f"{molecule_object}.pdb"
        tmp_molecule_object = molecule_object

    return tmp_molecule_object, tmp_basename


def filter_chains_for_protein_chains(chains: list["chain.Chain"]) -> list["chain.Chain"]:
    """Filters the chains for protein chains only.

    Args:
        chains: a list of chains which occur in the protein
    Returns: a list of protein chains only
    """
    # <editor-fold desc="Checks">
    safeguard.Safeguard.check_if_value_is_not_none(chains, logger)
    if not safeguard.Safeguard.check_if_list_is_empty(
        chains,
    ):
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")
    # </editor-fold>

    protein_chains = []
    for tmp_chain in chains:
        if tmp_chain.chain_type == constants.CHAIN_TYPE_PROTEIN:
            protein_chains.append(tmp_chain)
    return protein_chains


def get_chains_as_list_of_tuples(chains: list["chain.Chain"]) -> list[tuple[str, "sequence.Sequence", str]]:
    """Gets the chains as a list of tuples containing the chain letter, sequence and type."""
    chains_information = []
    for tmp_chain in chains:
        chains_information.append((tmp_chain.chain_letter, tmp_chain.chain_sequence, tmp_chain.chain_type))
    return chains_information


def create_chains_from_list_of_tuples(
        chains_as_list_of_tuples: list[tuple[str, "sequence.Sequence", str]]
) -> list["chain.Chain"]:
    """Creates chain objects from a list of tuples containing the chain letter, sequence and type."""
    chains: list["chain.Chain"] = []
    for tmp_chain_information in chains_as_list_of_tuples:
        chains.append(chain.Chain(tmp_chain_information[0], tmp_chain_information[1], tmp_chain_information[2]))
    return chains


def get_chains_from_list_of_chain_names(a_protein: "protein.Protein", chain_names: list) -> list["chain.Chain"]:
    """Gets the chains from a list of chain names.

    Args:
        a_protein: the protein from which the chains are used.
        chain_names: a list of chain names.
    """
    chains: list["chain.Chain"] = []
    for tmp_chain in a_protein.chains:
        for tmp_chain_name in chain_names:
            if tmp_chain.chain_letter == tmp_chain_name:
                chains.append(tmp_chain)
                break
    return chains
