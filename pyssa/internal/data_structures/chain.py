#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2022
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
"""This module contains the chain class."""
import logging
import pathlib
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants
from pyssa.util import types

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class Chain:
    """This class contains information about a single chain."""

    # <editor-fold desc="Class attributes">
    """
    letter of the chain
    """
    chain: str = ""
    """
    sequence of the chain
    """
    sequence: str = ""
    """
    type of the chain, whether it is a protein, or nuclein acid chain or something different
    """
    chain_type: str = ""

    # </editor-fold>

    def __init__(self, chain: str, sequence: str, chain_type: str) -> None:
        """Constructor

        Args:
            chain:
                letter of the chain
            sequence:
                sequence of the chain
            chain_type:
                type of the chain, whether it is a protein, or nuclein acid chain or something different
        """
        self.chain = chain
        self.sequence = sequence
        self.chain_type = chain_type
