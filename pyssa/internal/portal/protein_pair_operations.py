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
"""Module for protein_pair operations in pymol"""
import pymol
import logging
from pymol import cmd
from pyssa.io_pyssa import safeguard
from pyssa.io_pyssa import filesystem_io
from pyssa.util import constants
from pyssa.logging_pyssa import log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def align_protein_pair(target_selection, mobile_selection, alignment_filename) -> tuple:
    """This function aligns two pymol selections and creates an alignment file.

    Args:
        target_selection:
            first protein selection
        mobile_selection:
            second protein selection
        alignment_filename:
            name of the alignment file

    Returns:
        a tuple with rmsd and aligned amino acids
    """
    # <editor-fold desc="Checks">
    if not safeguard.Safeguard.check_if_value_is_not_none(target_selection):
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")
    if not safeguard.Safeguard.check_if_value_is_not_none(mobile_selection):
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")
    if not safeguard.Safeguard.check_if_value_is_not_none(alignment_filename):
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")

    # </editor-fold>

    tmp_settings = filesystem_io.ObjectDeserializer(constants.SETTINGS_DIR, constants.SETTINGS_FILENAME).deserialize_settings()
    results = cmd.align(target=target_selection, mobile=mobile_selection,
                        cutoff=tmp_settings.cutoff, cycles=tmp_settings.cycles,
                        object=alignment_filename, quiet=0)
    return results


