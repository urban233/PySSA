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
"""Module contains the ProteinSequence as a class."""
import pathlib
import logging
from pyssa.io_pyssa import safeguard
from pyssa.io_pyssa import filesystem_io
from pyssa.logging_pyssa import log_handlers
from pyssa.util import exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class Sequence:
  """Contains information about a protein sequence."""

  # <editor-fold desc="Class attributes">
  name: str
  """The name of the protein."""

  sequence: str
  """The sequence of the protein."""

  # </editor-fold>

  def __init__(self, protein_name: str, single_sequence: str) -> None:
    """Constructor.

    Args:
        protein_name (str): the name of the protein
        single_sequence (str): a sequence as string of the protein

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if protein_name is None:
      logger.error("protein_name is None.")
      raise exception.IllegalArgumentError("protein_name is None.")
    if single_sequence is None:
      logger.error("single_sequence is None.")
      raise exception.IllegalArgumentError("single_sequence is None.")

    # </editor-fold>

    self.name: str = protein_name
    self.sequence: str = single_sequence
