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
"""Module contains the selection class."""
import logging
from typing import TYPE_CHECKING
from src.pyssa.io_pyssa import safeguard
from src.pyssa.logging_pyssa import log_handlers

from src.pyssa.util import exception

if TYPE_CHECKING:
  from src.pyssa.internal.data_structures import chain

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class Selection:
  """Contains information about a selection."""

  # <editor-fold desc="Class attribute">
  selection_string: str
  """A pymol conform selection string."""

  molecule_object: str
  """The name of the protein which is also used within pymol."""

  selection_chain_letters: list[str]
  """A list of all current chains from the selection."""

  # </editor-fold>

  def __init__(self, molecule_object: str) -> None:
    """Constructor.

    Args:
        molecule_object: the name of the protein which is also used within pymol

    Raises:
        exception.IllegalArgumentError: If `molecule_object` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if molecule_object is None or molecule_object == "":
      logger.error("molecule_object is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "molecule_object is either None or an empty string."
      )

    # </editor-fold>

    self.molecule_object = molecule_object

  def set_selections_from_chains_ca(self, chains: list["chain.Chain"]) -> None:
    """Sets a selection based on the chains of the protein. The selection selects only the alpha-C's.

    Args:
        chains: A list of chains.

    Raises:
        exception.IllegalArgumentError: If `chains` is None.
    """
    # <editor-fold desc="Checks">
    if chains is None:
      logger.error("chains is None.")
      raise exception.IllegalArgumentError("chains is None.")

    # </editor-fold>

    seperator = ", "
    tmp_list = []
    self.selection_chain_letters = []
    for tmp_chain in chains:
      tmp_selection = f"/{self.molecule_object}//{tmp_chain.chain_letter}//CA"
      self.selection_chain_letters.append(tmp_chain.chain_letter)
      tmp_list.append(tmp_selection)
      self.selection_string = seperator.join(tmp_list)

  def set_selections_without_chains_ca(self) -> None:
    """Sets a selection without any chains of the protein.

    Notes:
        The selection selects only the alpha-C's.
    """
    self.selection_string = f"/{self.molecule_object}////CA"

  def set_single_selection(
      self, segi: str, chain: str, resi: str, atom_name: str
  ) -> None:
    """Creates a single pymol selection with only one chain and one resi.

    Args:
        segi: A segment identifier.
        chain: A chain identifier.
        resi: A residue name or position.
        atom_name: A type of atom.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if segi is None:
      logger.error("segi is None.")
      raise exception.IllegalArgumentError("segi is None.")
    if chain is None:
      logger.error("chain is None.")
      raise exception.IllegalArgumentError("chain is None.")
    if resi is None:
      logger.error("resi is None.")
      raise exception.IllegalArgumentError("resi is None.")
    if atom_name is None:
      logger.error("atom_name is None.")
      raise exception.IllegalArgumentError("atom_name is None.")

    # </editor-fold>

    self.selection_string = (
        f"/{self.molecule_object}/{segi}/{chain}/{resi}/{atom_name}"
    )

  def set_selection_for_a_single_chain(self, a_chain_letter: str) -> None:
    """Sets a selection for a single chain.

    Args:
        a_chain_letter: A letter of a chain.

    Raises:
        exception.IllegalArgumentError: If `a_chain_letter` is None.
    """
    # <editor-fold desc="Checks">
    if a_chain_letter is None:
      logger.error("a_chain_letter is None.")
      raise exception.IllegalArgumentError("a_chain_letter is None.")

    # </editor-fold>

    self.selection_string = f"/{self.molecule_object}//{a_chain_letter}"

  def set_selection_for_first_ca_atom_in_a_given_chain(
      self, a_chain_letter: str
  ) -> None:
    """Sets a selection for alpha C atoms of a single chain.

    Args:
        a_chain_letter (str): The letter identifier for the chain.

    Raises:
        exception.IllegalArgumentError: If `a_chain_letter` is None.
    """
    # <editor-fold desc="Checks">
    if a_chain_letter is None:
      logger.error("a_chain_letter is None.")
      raise exception.IllegalArgumentError("a_chain_letter is None.")

    # </editor-fold>

    self.selection_string = f"first chain {a_chain_letter} and name CA"

  def set_custom_selection(self, a_sele_string: str) -> None:
    """Sets a custom selection as selection string.

    Args:
        a_sele_string: a custom pymol selection string.

    Raises:
        exception.IllegalArgumentError: If `a_sele_string` is None.
    """
    # <editor-fold desc="Checks">
    if a_sele_string is None:
      logger.error("a_sele_string is None.")
      raise exception.IllegalArgumentError("a_sele_string is None.")

    # </editor-fold>

    self.selection_string = a_sele_string
