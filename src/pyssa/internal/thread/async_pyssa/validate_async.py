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
"""Module for validation functions that are used by async functions."""
import logging
import os
import pathlib
from typing import Optional

from src.pyssa.io_pyssa import bio_data
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import constants
from src.pyssa.util import exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


def validate_add_protein_view_input(
    the_entered_text: str, placeholder: int
) -> tuple[int, bool, Optional[str]]:
  """Validates the input for the add_protein_view.

  Args:
      the_entered_text: The text entered by the user.
      placeholder: The placeholder parameter so that the function can be used with LegacyTask.

  Returns:
      A tuple containing an integer, a boolean, and a string.

      The first value in the tuple represents the validation status:
          - 1: Validated as a PDB ID.
          - 2: Validated as a file path.
          - -1: Validation failed due to an exception.

      The second value in the tuple represents the validation result:
          - True: Valid.
          - False: Invalid.

      The third value in the tuple depends on the validation status:
          - If the validation status is 1 or 2, this value is the corresponding ID (PDB ID or file name without extension).
          - If the validation status is -1, this value is None.
  """
  # <editor-fold desc="Checks">
  if the_entered_text is None:
    logger.error("the_entered_text is None.")
    return -1, False, None

  # </editor-fold>

  try:
    # checks if a pdb id was entered
    if len(the_entered_text) == 4:
      # pdb id is used
      pdb_id = the_entered_text.upper()
      tmp_filepath: str = str(
          pathlib.Path(f"{constants.SCRATCH_DIR}/{pdb_id}.pdb")
      )
      bio_data.download_pdb_file(pdb_id, tmp_filepath)
      if os.path.exists(tmp_filepath):
        os.remove(tmp_filepath)
        return 1, True, pdb_id
      else:
        return 1, False, pdb_id
    else:
      if os.path.exists(the_entered_text):
        return 2, True, pathlib.Path(the_entered_text).name.replace(".pdb", "")
      else:
        return 2, False, pathlib.Path(the_entered_text).name.replace(".pdb", "")
  except Exception as e:
    logger.error(e)
    return -1, False, None
