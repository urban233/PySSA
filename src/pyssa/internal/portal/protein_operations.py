#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
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

import zmq

from src.auxiliary_pymol import auxiliary_pymol_client
from src.pyssa.internal.data_structures import chain, job
from src.pyssa.util import enums, exception
from src.pyssa.util import protein_util
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.internal.data_structures import sequence

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


def get_protein_chains(
    a_pdb_filepath: str, the_main_socket: zmq.Socket, a_socket: zmq.Socket
) -> list["chain.Chain"]:
  """Divides the chains from a protein, into protein and non-protein chains.

  Args:
      a_pdb_filepath: A string representing the path to the PDB file.
      the_main_socket: A ZMQ socket object used for communication with the main server.
      a_socket: A ZMQ socket object used for communication with the auxiliary PyMOL server.

  Returns:
      A list of Chain objects representing the protein chains obtained from the PDB file.

  Raises:
      exception.IllegalArgumentError: If any of the args are None or `a_pdb_filepath` is an empty string.
  """
  # <editor-fold desc="Checks">
  if a_pdb_filepath is None or a_pdb_filepath == "":
    logger.error("a_pdb_filepath is either None or an empty string.")
    raise exception.IllegalArgumentError(
        "a_pdb_filepath is either None or an empty string."
    )
  if the_main_socket is None:
    logger.error("the_main_socket is None.")
    raise exception.IllegalArgumentError("the_main_socket is None.")
  if a_socket is None:
    logger.error("a_socket is None.")
    raise exception.IllegalArgumentError("a_socket is None.")

  # </editor-fold>

  tmp_job_description = job.GeneralPurposeJobDescription(
      enums.JobShortDescription.GET_ALL_CHAINS_OF_GIVEN_PROTEIN
  )
  tmp_job_description.setup_dict(
      {enums.JobDescriptionKeys.PDB_FILEPATH.value: str(a_pdb_filepath)}
  )
  tmp_reply = auxiliary_pymol_client.send_request_to_auxiliary_pymol(
      the_main_socket,
      a_socket,
      tmp_job_description,
  )
  tmp_data: list[tuple] = tmp_reply["data"]
  tmp_chains_of_protein: list[chain.Chain] = []
  for tmp_chain_object_values in tmp_data:
    tmp_chain_letter, tmp_sequence_object_values, tmp_chain_type = (
        tmp_chain_object_values
    )
    tmp_sequence = sequence.Sequence(
        tmp_sequence_object_values[0], tmp_sequence_object_values[1]
    )
    tmp_chains_of_protein.append(
        chain.Chain(tmp_chain_letter, tmp_sequence, tmp_chain_type)
    )
  return tmp_chains_of_protein


def get_protein_sequences_from_protein(
    molecule_object: str, chains: list["chain.Chain"]
) -> list["sequence.Sequence"]:
  """Gets all sequences from protein chains only.

  Args:
      molecule_object (str): The name of the protein which is also used within PyMOL.
      chains (list[chain.Chain]): A list of chains which occur in the protein.

  Returns:
      A protein sequence object with amino acid sequences.

  Raises:
      exception.IllegalArgumentError: If `molecule_object` is either None or an empty string.
      exception.IllegalArgumentError: If `chains` is either None or an empty list.
  """
  # <editor-fold desc="Checks">
  if molecule_object is None or molecule_object == "":
    logger.error("molecule_object is either None or an empty string.")
    raise exception.IllegalArgumentError(
        "molecule_object is either None or an empty string."
    )
  if chains is None or len(chains) == 0:
    logger.error("chains is either None or an empty list.")
    raise exception.IllegalArgumentError(
        "chains is either None or an empty list."
    )

  # </editor-fold>

  protein_sequences: list[sequence.Sequence] = []
  for tmp_chain in protein_util.filter_chains_for_protein_chains(chains):
    protein_sequences.append(tmp_chain.chain_sequence)
  return protein_sequences


def get_chain_letter_of_first_protein_sequence(
    chains: list["chain.Chain"],
) -> str:
  """Gets the chain letter of the first protein sequence.

  Args:
      chains (list[chain.Chain]): A list of objects of type "chain.Chain" representing protein chains.

  Returns:
      The chain letter of the first protein sequence found in the given list of chains.

  Raises:
      exception.IllegalArgumentError: If `chains` is either None or an empty list.
  """
  # <editor-fold desc="Checks">
  if chains is None or len(chains) == 0:
    logger.error("chains is either None or an empty list.")
    raise exception.IllegalArgumentError(
        "chains is either None or an empty list."
    )

  # </editor-fold>

  tmp_protein_chains = protein_util.filter_chains_for_protein_chains(chains)
  return tmp_protein_chains[0].chain_letter
