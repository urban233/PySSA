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
"""Module contains the chain class."""
import logging

from src.pyssa.controller import pymol_session_manager
from src.pyssa.internal.portal import graphic_operations
from src.pyssa.util import enums
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.internal.data_structures import sequence
from src.pyssa.util import exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class Chain:
  """This class contains information about a single chain."""

  # <editor-fold desc="Class attributes">
  _id: int
  """The id of the chain."""

  chain_letter: str
  """A letter of the chain."""

  chain_sequence: sequence.Sequence
  """A sequence of the chain."""

  chain_type: str
  """A type of the chain, whether it is a protein, or nuclein acid chain or something different."""

  pymol_parameters: dict
  """A dict of parameters that can be changed in pymol."""

  db_protein_id: int
  """The protein id from the project database."""

  # </editor-fold>

  def __init__(
      self,
      chain_letter: str,
      chain_sequence: "sequence.Sequence",
      chain_type: str,
  ) -> None:
    """Constructor.

    Args:
        chain_letter (str): The letter representing the chain.
        chain_sequence (sequence.Sequence): The sequence of the chain.
        chain_type (str): The type of the chain.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if chain_letter is None:
      logger.error("chain_letter is None.")
      raise exception.IllegalArgumentError("chain_letter is None.")
    if chain_sequence is None:
      logger.error("chain_sequence is None.")
      raise exception.IllegalArgumentError("chain_sequence is None.")
    if chain_type is None:
      logger.error("chain_type is None.")
      raise exception.IllegalArgumentError("chain_type is None.")

    # </editor-fold>

    self.chain_letter = chain_letter
    self.chain_sequence = chain_sequence
    self.chain_type = chain_type
    self.pymol_parameters = {
        enums.PymolParameterEnum.COLOR.value: "green",
        enums.PymolParameterEnum.REPRESENTATION.value: "cartoon",
    }

  def get_id(self) -> int:
    """Returns the id of the object.

    Returns:
        The id value.
    """
    return self._id

  def set_id(self, value: int) -> None:
    """Sets the id of the object.

    Args:
        value: The value to be set as the ID.

    Raises:
        exception.IllegalArgumentError: If `value` is either None or a negative integer.
    """
    # <editor-fold desc="Checks">
    if value is None or value < 0:
      logger.error("value is either None or a negative integer.")
      raise exception.IllegalArgumentError(
          "value is either None or a negative integer."
      )

    # </editor-fold>

    self._id = value

  def get_color(
      self,
      a_selection_string: str,
      the_pymol_session_manager: "pymol_session_manager.PymolSessionManager",
  ) -> tuple[str, bool]:
    """Gets the color for a given selection.

    Args:
        a_selection_string (str): The selection string used to identify the objects in Pymol.
        the_pymol_session_manager (pymol_session_manager.PymolSessionManager): An instance of the PymolSessionManager class.

    Returns:
        A tuple containing the color string and a boolean indicating whether the objects are colored by elements.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_selection_string is None:
      logger.error("a_selection_string is None.")
      raise exception.IllegalArgumentError("a_selection_string is None.")
    if the_pymol_session_manager is None:
      logger.error("the_pymol_session_manager is None.")
      raise exception.IllegalArgumentError("the_pymol_session_manager is None.")

    # </editor-fold>

    (
        self.pymol_parameters[enums.PymolParameterEnum.COLOR.value],
        _,
        tmp_is_colored_by_elements,
    ) = the_pymol_session_manager.get_chain_color(
        a_selection_string, self.chain_letter
    )
    return (
        self.pymol_parameters[enums.PymolParameterEnum.COLOR.value],
        tmp_is_colored_by_elements,
    )

  def get_representation_state(self, a_selection_string: str) -> dict:
    """Gets the representation state for a given selection.

    Notes:
        DO NOT USE THIS METHOD ANYMORE!

    Args:
        a_selection_string (str): A string representing selected elements in a graphical representation.

    Returns:
        A dictionary containing the representation state of the selected elements.

    Raises:
        exception.IllegalArgumentError: If `a_selection_string` is None.
    """
    # <editor-fold desc="Checks">
    if a_selection_string is None:
      logger.error("a_selection_string is None.")
      raise exception.IllegalArgumentError("a_selection_string is None.")

    # </editor-fold>

    return graphic_operations.get_chain_repr_state(
        a_selection_string, self.chain_letter
    )
