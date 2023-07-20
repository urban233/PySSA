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
"""Module for in- and output processes in pymol"""
import pathlib
import pymol
import os
from pymol import cmd

from pyssa.internal.portal import graphic_operations
from pyssa.io_pyssa import safeguard
from pyssa.io_pyssa import binary_data
from pyssa.io_pyssa import path_util
from pyssa.util import constants
from pyssa.util import globals


def load_protein(filepath: pathlib.Path, basename: str, molecule_object: str) -> None:
    """This function loads a protein in pymol through a protein object.

    Args:
        protein_obj:
            object which is in instance of the protein class

    """
    if not safeguard.Safeguard.check_filepath(f"{filepath}/{basename}"):
        raise FileNotFoundError
    cmd.load(f"{filepath}/{basename}", object=molecule_object)
    if globals.g_settings.color_vision_mode == constants.CVM_NORMAL:
        color_prot_1 = constants.CVM_NORMAL_PROT_1_COLOR
    elif globals.g_settings.color_vision_mode == constants.CVM_DEUTERANOPIA:
        color_prot_1 = constants.CVM_DEUTERANOPIA_PROT_1_COLOR
    elif globals.g_settings.color_vision_mode == constants.CVM_PROTANOPIA:
        color_prot_1 = constants.CVM_PROTANOPIA_PROT_1_COLOR
    elif globals.g_settings.color_vision_mode == constants.CVM_TRITANOPIA:
        color_prot_1 = constants.CVM_TRITANOPIA_PROT_1_COLOR
    else:
        color_prot_1 = "green"
    cmd.color(color_prot_1, molecule_object)


def fetch_protein_from_pdb(filepath: pathlib.Path, filename: str, molecule_object: str) -> None:
    """This function fetches a protein in pymol from the PDB.

    Args:
        protein_obj:
            object which is in instance of the protein class
    Raises:
        FileNotFoundError: If file not found.
        ValueError: If PDB ID couldn't be found in PDB.

    """
    if not safeguard.Safeguard.check_filepath(f"{filepath}/{filename}"):
        raise FileNotFoundError
    try:
        cmd.fetch(code=molecule_object, type="pdb", path=filepath, file=filename)
    except pymol.CmdException:
        raise ValueError("PDB ID is invalid.")


def save_protein_to_pdb_file(export_filepath: pathlib.Path, molecule_object: str) -> None:
    """This function saves a protein from the current pymol session as a .pdb file.

    Args:
        protein_obj:
            object which is in instance of the protein class
    Raises:
        NotADirectoryError: If directory is not found.

    """
    if not safeguard.Safeguard.check_filepath(f"{export_filepath}"):
        raise NotADirectoryError(f"The filepath {export_filepath} does not exists.")
    # save the pdb file under the path (export_data_dir)
    cmd.save(f"{export_filepath}/{molecule_object}.pdb")


def load_pymol_session(pymol_session_file):
    """This function loads a pymol session file into the current pymol session.

    Args:
        pymol_session_file:
            filepath of the session file
    Raises:
        FileNotFoundError: If file not found.

    """
    if not safeguard.Safeguard.check_filepath(f"{pymol_session_file}"):
        raise FileNotFoundError
    cmd.load(pymol_session_file)
    graphic_operations.setup_default_session_graphic_settings()


def convert_pymol_session_to_base64_string(pymol_molecule_object) -> str:
    session_filepath = pathlib.Path(f"{constants.SCRATCH_DIR}/{pymol_molecule_object}_session.pse")
    cmd.save(session_filepath)
    base64_string = binary_data.create_base64_string_from_file(path_util.FilePath(session_filepath))
    os.remove(session_filepath)
    return base64_string
