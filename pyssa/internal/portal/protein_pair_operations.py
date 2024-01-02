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
from pyssa.util import globals
from pyssa.logging_pyssa import log_handlers
from pyssa.internal.portal import pymol_io

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def align_protein_pair(target_selection: str, mobile_selection: str, alignment_filename: str) -> tuple:
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

    tmp_settings = filesystem_io.ObjectDeserializer(
        constants.SETTINGS_DIR,
        constants.SETTINGS_FILE_NAME,
    ).deserialize_settings()
    results = cmd.align(
        target=target_selection,
        mobile=mobile_selection,
        cutoff=tmp_settings.cutoff,
        cycles=tmp_settings.cycles,
        object=alignment_filename,
        quiet=0,
    )
    return results  # noqa: RET504


def color_protein_pair(
    pymol_molecule_object_ref: str,
    pymol_molecule_object_model: str,
    color_ref: str = "green",
    color_model: str = "blue",
) -> None:
    """This function colors both the reference and the model Protein.

    Note:
        Only the official colors from PyMOL are supported. These can
        be looked up under the `color values`_ page.

    Args:
        pymol_molecule_object_ref: the molecule object string of the reference protein
        pymol_molecule_object_model: the molecule object string of the model protein
        color_ref (str, optional):
            defines color for the reference Protein
        color_model (str, optional):
            defines color for the model Protein

    Raises:
        pymol.CmdException:
            Exception is raised if one or both proteins
            does not exist as pymol objects.

    .. _color values:
        https://pymolwiki.org/index.php/Color_Values
    """
    # argument test
    # checks if either the reference or the model is an actual object in the memory
    if not safeguard.Safeguard.check_if_protein_is_in_pymol(pymol_molecule_object_ref):
        raise pymol.CmdException("The reference is not in the pymol session as an object.")
    if not safeguard.Safeguard.check_if_protein_is_in_pymol(pymol_molecule_object_model):
        raise pymol.CmdException("The model is not in the pymol session as an object.")
    # actual color cmd command
    if globals.g_settings.color_vision_mode == constants.CVM_NORMAL:
        color_ref = constants.CVM_NORMAL_PROT_1_COLOR
        color_model = constants.CVM_NORMAL_PROT_2_COLOR
    elif globals.g_settings.color_vision_mode == constants.CVM_DEUTERANOPIA:
        color_ref = constants.CVM_DEUTERANOPIA_PROT_1_COLOR
        color_model = constants.CVM_DEUTERANOPIA_PROT_2_COLOR
    elif globals.g_settings.color_vision_mode == constants.CVM_PROTANOPIA:
        color_ref = constants.CVM_PROTANOPIA_PROT_1_COLOR
        color_model = constants.CVM_PROTANOPIA_PROT_2_COLOR
    elif globals.g_settings.color_vision_mode == constants.CVM_TRITANOPIA:
        color_ref = constants.CVM_TRITANOPIA_PROT_1_COLOR
        color_model = constants.CVM_TRITANOPIA_PROT_2_COLOR
    cmd.color(color_ref, pymol_molecule_object_ref)
    cmd.color(color_model, pymol_molecule_object_model)


def save_session_of_protein_pair(name_of_protein_pair: str) -> str:
    """Saves the pymol session of the protein pair.

    Args:
        name_of_protein_pair (str): Name of the protein pair.
    """
    return pymol_io.convert_pymol_session_to_base64_string(name_of_protein_pair)
