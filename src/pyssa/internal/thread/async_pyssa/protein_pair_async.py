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
"""Module for all asynchronous functions that are related to the protein pair object."""
import logging
from typing import TYPE_CHECKING, Optional
from src.pyssa.controller import pymol_session_manager
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"

if TYPE_CHECKING:
  from src.pyssa.internal.data_structures import protein_pair


def color_protein_pair_by_rmsd_value(
    a_protein_pair: "protein_pair.ProteinPair",
    the_pymol_session_manager: "pymol_session_manager.PymolSessionManager",
) -> tuple[str, Optional["protein_pair.ProteinPair"]]:
  """Colors protein pair by RMSD value.

  Args:
      a_protein_pair (protein_pair.ProteinPair): The protein pair to be colored.
      the_pymol_session_manager (pymol_session_manager.PymolSessionManager): The Pymol session manager object.

  Returns:
      A tuple containing the result and the protein pair.
      If an exception occurs during coloring, an empty string and None will be returned.
      Otherwise, the string "result" and the protein pair will be returned.
  """
  # <editor-fold desc="Checks">
  if a_protein_pair is None:
    logger.error("a_protein_pair is None.")
    return "", None
  if the_pymol_session_manager is None:
    logger.error("the_pymol_session_manager is None.")
    return "", None

  # </editor-fold>

  try:
    the_pymol_session_manager.color_protein_pair_by_rmsd(a_protein_pair)
  except Exception as e:
    logger.error(e)
    return "", None
  else:
    return "result", a_protein_pair
