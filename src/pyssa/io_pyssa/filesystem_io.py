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
"""Module for handling filesystem operations."""
import os
import pathlib
import shutil
import logging
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import constants

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class FilesystemCleaner:
  """Class for cleaning up the filesystem."""

  def __init__(self) -> None:
    """Empty constructor."""
    pass

  @staticmethod
  def clean_prediction_scratch_folder() -> None:
    """Deletes the scratch folder for fasta and pdb files and creates new ones."""
    if pathlib.Path(constants.PREDICTION_FASTA_DIR).exists():
      shutil.rmtree(constants.PREDICTION_FASTA_DIR)
    pathlib.Path(constants.PREDICTION_FASTA_DIR).mkdir(parents=True)
    if pathlib.Path(constants.PREDICTION_PDB_DIR).exists():
      shutil.rmtree(constants.PREDICTION_PDB_DIR)
    pathlib.Path(constants.PREDICTION_PDB_DIR).mkdir(parents=True)
