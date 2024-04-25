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
"""Module for protein_pair operations in pymol."""
import pymol
import logging
from pymol import cmd
from pyssa.io_pyssa import safeguard
from pyssa.io_pyssa import filesystem_io
from pyssa.util import constants
from pyssa.logging_pyssa import log_handlers
from pyssa.internal.portal import pymol_io

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def save_session_of_protein_pair(name_of_protein_pair: str) -> str:
    """Saves the pymol session of the protein pair.

    Args:
        name_of_protein_pair (str): Name of the protein pair.
    """
    return pymol_io.convert_pymol_session_to_base64_string(name_of_protein_pair)
