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
from pyssa.util import constants, exception
from typing import TYPE_CHECKING

if TYPE_CHECKING:
  from pyssa.internal.data_structures import sequence
  from pyssa.internal.data_structures import protein


logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


def check_if_protein_is_from_file_or_id(molecule_object: str) -> tuple:
  """Checks if a protein is from a .pdb file or from a PDB ID.

  Args:
      molecule_object (str): The name of the protein which is also used within pymol.

  Returns:
      A tuple with the molecule object and the basename.

  Raises:
      exception.IllegalArgumentError: If molecule_object is either None or an empty string.
  """
  # <editor-fold desc="Checks">
  if molecule_object is None or molecule_object == "":
    logger.error("molecule_object is either None or an empty string.")
    raise exception.IllegalArgumentError(
        "molecule_object is either None or an empty string."
    )

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


def filter_chains_for_protein_chains(
    chains: list["chain.Chain"],
) -> list["chain.Chain"]:
  """Filters the chains for protein chains only.

  Args:
      chains (list[chain.Chain]): A list of chains which occur in the protein.

  Returns:
      A list of protein chains only.

  Raises:
      exception.IllegalArgumentError: If chains is either None or an empty list.
  """
  # <editor-fold desc="Checks">
  if chains is None or len(chains) == 0:
    logger.error("chains is either None or an empty list.")
    raise exception.IllegalArgumentError(
        "chains is either None or an empty list."
    )

  # </editor-fold>

  protein_chains = []
  for tmp_chain in chains:
    if tmp_chain.chain_type == constants.CHAIN_TYPE_PROTEIN:
      protein_chains.append(tmp_chain)
  return protein_chains


def get_chains_as_list_of_tuples(
    chains: list["chain.Chain"],
) -> list[tuple[str, "sequence.Sequence", str]]:
  """Gets the chains as a list of tuples containing the chain letter, sequence and type.

  Args:
      chains (list[chain.Chain]): A list of chains which occur in the protein.

  Returns:
      A list of tuples containing the chain letter, sequence and type.

  Raises:
      exception.IllegalArgumentError: If chains is either None or an empty list.
  """
  # <editor-fold desc="Checks">
  if chains is None or len(chains) == 0:
    logger.error("chains is either None or an empty list.")
    raise exception.IllegalArgumentError(
        "chains is either None or an empty list."
    )

  # </editor-fold>

  chains_information = []
  for tmp_chain in chains:
    chains_information.append(
        (tmp_chain.chain_letter, tmp_chain.chain_sequence, tmp_chain.chain_type)
    )
  return chains_information


def create_chains_from_list_of_tuples(
    chains_as_list_of_tuples: list[tuple[str, "sequence.Sequence", str]],
) -> list["chain.Chain"]:
  """Creates chain objects from a list of tuples containing the chain letter, sequence and type.

  Args:
      chains_as_list_of_tuples (list[tuple[str, "sequence.Sequence", str]]): A list of tuples containing the chain letter, sequence and type.

  Returns:
      A list of chain objects.

  Raises:
      exception.IllegalArgumentError: If chains_as_list_of_tuples is None.
  """
  # <editor-fold desc="Checks">
  if chains_as_list_of_tuples is None:
    logger.error("chains_as_list_of_tuples is None.")
    raise exception.IllegalArgumentError("chains_as_list_of_tuples is None.")

  # </editor-fold>

  chains: list["chain.Chain"] = []
  for tmp_chain_information in chains_as_list_of_tuples:
    chains.append(
        chain.Chain(
            tmp_chain_information[0],
            tmp_chain_information[1],
            tmp_chain_information[2],
        )
    )
  return chains


def get_chains_from_list_of_chain_names(
    a_protein: "protein.Protein", chain_names: list
) -> list["chain.Chain"]:
  """Gets the chains from a list of chain names.

  Args:
      a_protein (protein.Protein): the protein from which the chains are used.
      chain_names (list): a list of chain names.

  Returns:
      A list of chain objects.

  Raises:
      exception.IllegalArgumentError: If a_protein is None.
      exception.IllegalArgumentError: If chain_names is either None or an empty list.
  """
  # <editor-fold desc="Checks">
  if a_protein is None:
    logger.error("a_protein is None.")
    raise exception.IllegalArgumentError("a_protein is None.")
  if chain_names is None or len(chain_names) == 0:
    logger.error("chain_names is either None or an empty list.")
    raise exception.IllegalArgumentError(
        "chain_names is either None or an empty list."
    )

  # </editor-fold>

  chains: list["chain.Chain"] = []
  for tmp_chain in a_protein.chains:
    for tmp_chain_name in chain_names:
      if tmp_chain.chain_letter == tmp_chain_name:
        chains.append(tmp_chain)
        break
  return chains
