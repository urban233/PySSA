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

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class Sequence:
    """This class contains information about a protein sequence."""

    # <editor-fold desc="Class attributes">
    """
    name of the protein
    """
    name: str
    """
    sequence of the protein
    """
    sequence: str

    # </editor-fold>

    def __init__(self, protein_name: str, single_sequence: str) -> None:
        """Constructor.

        Args:
            protein_name: the name of the protein
            single_sequence: a sequence as string of the protein
        Raises:
            ValueError: raised if an argument is illegal
        """
        # <editor-fold desc="Checks">
        # if len(protein_sequence[0]) == 1:
        #     if not safeguard.Safeguard.check_if_protein_sequence_is_valid_one_letter(protein_sequence):
        #         logger.error("An argument is illegal.")
        #         raise ValueError("An argument is illegal.")
        # elif len(protein_sequence[0]) == 3:
        #     if not safeguard.Safeguard.check_if_protein_sequence_is_valid_three_letter(protein_sequence):
        #         logger.error("An argument is illegal.")
        #         raise ValueError("An argument is illegal.")
        # else:
        #     logger.error("An argument is illegal.")
        #     raise ValueError("An argument is illegal.")

        # </editor-fold>

        self.name: str = protein_name
        self.sequence: str = single_sequence
