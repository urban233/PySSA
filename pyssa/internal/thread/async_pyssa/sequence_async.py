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
"""Module for all asynchronous functions that are related to the sequence/SeqRecord object."""
import logging
from typing import TYPE_CHECKING
from Bio import SeqIO

from pyssa.logging_pyssa import log_handlers
from pyssa.util import exit_codes, exception

if TYPE_CHECKING:
  from Bio import SeqRecord

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


def save_selected_protein_sequence_as_fasta_file(
    a_seq_record: "SeqRecord.SeqRecord",
    a_filepath: str,
) -> tuple[int, str]:
  """Saves a given protein sequence to a fasta file.

  Args:
      a_seq_record (SeqRecord.SeqRecord): The protein sequence to be saved as a fasta file.
      a_filepath (str): The file path where the fasta file will be saved.

  Returns:
      A tuple containing the exit code and message indicating the success or failure of saving the sequence as a fasta file.
  """
  # <editor-fold desc="Checks">
  if a_seq_record is None:
    logger.error("a_seq_record is None.")
    return (
        exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0],
        exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[1],
    )
  if a_filepath is None or a_filepath == "":
    logger.error("a_filepath is either None or an empty string.")
    return (
        exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0],
        exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[1],
    )

  # </editor-fold>

  try:
    print(a_seq_record.id)
    if a_seq_record.id == "<unknown id>":
      a_seq_record.id = a_seq_record.name
    with open(a_filepath, "w") as file_handler:
      SeqIO.write(a_seq_record, file_handler, "fasta")
  except Exception as e:
    logger.error(f"Saving sequence to fasta file ended with error: {e}")
    return (
        exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0],
        exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[1],
    )
  else:
    return exit_codes.EXIT_CODE_ZERO[0], exit_codes.EXIT_CODE_ZERO[1]
