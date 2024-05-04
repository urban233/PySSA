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
from typing import TYPE_CHECKING

from auxiliary_pymol import auxiliary_pymol_client
from pyssa.internal.data_structures import chain, job
from pyssa.io_pyssa import safeguard
from pyssa.util import constants, enums
from pyssa.util import protein_util
from pyssa.logging_pyssa import log_handlers
from pyssa.internal.data_structures import sequence

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def get_protein_chains(a_pdb_filepath, the_main_socket, a_socket) -> list[chain.Chain]:
    """Divides the chains from a protein, into protein and non-protein chains.

    Returns:
        a list of all chains from the protein, divided into protein and non-protein chains
    """
    tmp_job_description = job.GeneralPurposeJobDescription(enums.JobShortDescription.GET_ALL_CHAINS_OF_GIVEN_PROTEIN)
    tmp_job_description.setup_dict({enums.JobDescriptionKeys.PDB_FILEPATH.value: str(a_pdb_filepath)})
    tmp_reply = auxiliary_pymol_client.send_request_to_auxiliary_pymol(
        the_main_socket, a_socket, tmp_job_description
    )
    tmp_data: list[tuple] = tmp_reply["data"]
    tmp_chains_of_protein: list[chain.Chain] = []
    for tmp_chain_object_values in tmp_data:
        tmp_chain_letter, tmp_sequence_object_values, tmp_chain_type = tmp_chain_object_values
        tmp_sequence = sequence.Sequence(tmp_sequence_object_values[0], tmp_sequence_object_values[1])
        tmp_chains_of_protein.append(chain.Chain(tmp_chain_letter, tmp_sequence, tmp_chain_type))
    # cmd.reinitialize()
    # pymol_io.load_protein(dirname, basename, molecule_object)
    # tmp_chains: list[str] = cmd.get_chains()
    # i = 0
    # chains_of_protein: list[chain.Chain] = []
    # for tmp_chain in tmp_chains:
    #     sequence_of_chain = cmd.get_model(f"chain {tmp_chain}")
    #     if sequence_of_chain.atom[0].resn in constants.AMINO_ACID_CODE:
    #         # TODO: can this be better?
    #         fasta_sequence_of_chain = cmd.get_fastastr(f"chain {tmp_chain}")
    #         fasta_sequence_of_chain_without_header = fasta_sequence_of_chain[fasta_sequence_of_chain.find("\n") :]
    #         complete_sequence_of_chain = sequence.Sequence(
    #             molecule_object,
    #             fasta_sequence_of_chain_without_header.replace("\n", ""),
    #         )
    #         chains_of_protein.append(chain.Chain(tmp_chain, complete_sequence_of_chain, constants.CHAIN_TYPE_PROTEIN))
    #     else:
    #         # TODO: this should produce a sequence with the non-protein atoms
    #         fasta_sequence_of_chain = cmd.get_fastastr(f"chain {tmp_chain}")
    #         fasta_sequence_of_chain_without_header = fasta_sequence_of_chain[fasta_sequence_of_chain.find("\n") :]
    #         complete_sequence_of_chain = sequence.Sequence(
    #             molecule_object,
    #             fasta_sequence_of_chain_without_header.replace("\n", ""),
    #         )
    #         chains_of_protein.append(
    #             chain.Chain(tmp_chain, complete_sequence_of_chain, constants.CHAIN_TYPE_NON_PROTEIN),
    #         )
    #     i += 1
    return tmp_chains_of_protein


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
    safeguard.Safeguard.check_if_value_is_not_none(molecule_object, logger)
    if molecule_object == "":
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")
    safeguard.Safeguard.check_if_value_is_not_none(chains, logger)
    if not safeguard.Safeguard.check_if_list_is_empty(
        chains,
    ):
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")

    # </editor-fold>
    protein_sequences: list[sequence.Sequence] = []
    for tmp_chain in protein_util.filter_chains_for_protein_chains(chains):
        protein_sequences.append(tmp_chain.chain_sequence)
    return protein_sequences


def get_chain_letter_of_first_protein_sequence(chains: list[chain.Chain]):
    tmp_protein_chains = protein_util.filter_chains_for_protein_chains(chains)
    return tmp_protein_chains[0].chain_letter
