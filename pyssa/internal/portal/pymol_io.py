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
"""Module for in- and output processes in pymol."""
import pathlib
import pymol
import os
from pymol import cmd
from pyssa.util import exception
from pyssa.internal.data_structures import protein
from pyssa.internal.portal import graphic_operations
from pyssa.io_pyssa import safeguard
from pyssa.io_pyssa import binary_data
from pyssa.io_pyssa import path_util
from pyssa.util import constants, tools


def load_protein(filepath: pathlib.Path, basename: str, molecule_object: str) -> None:
    """Loads a protein in pymol through a protein object.

    Args:
        filepath: the filepath to the protein
        basename: the name of the protein file
        molecule_object: the name of the protein in PyMOL
    """
    if not safeguard.Safeguard.check_filepath(pathlib.Path(f"{filepath}/{basename}")):
        raise FileNotFoundError
    cmd.load(f"{filepath}/{basename}", object=molecule_object)
    # fixme: this is code for the color vision feature
    # if globals.g_settings.color_vision_mode == constants.CVM_NORMAL:
    #     color_prot_1 = constants.CVM_NORMAL_PROT_1_COLOR
    # elif globals.g_settings.color_vision_mode == constants.CVM_DEUTERANOPIA:
    #     color_prot_1 = constants.CVM_DEUTERANOPIA_PROT_1_COLOR
    # elif globals.g_settings.color_vision_mode == constants.CVM_PROTANOPIA:
    #     color_prot_1 = constants.CVM_PROTANOPIA_PROT_1_COLOR
    # elif globals.g_settings.color_vision_mode == constants.CVM_TRITANOPIA:
    #     color_prot_1 = constants.CVM_TRITANOPIA_PROT_1_COLOR
    # else:
    #     color_prot_1 = "green"
    color_prot_1 = "green"
    graphic_operations.setup_default_session_graphic_settings()
    cmd.color(color_prot_1, molecule_object)


def fetch_protein_from_pdb(filepath: pathlib.Path, filename: str, molecule_object: str) -> None:
    """This function fetches a protein in pymol from the PDB.

    Args:
        filepath: the filepath to the protein
        filename: the name of the protein file with extension
        molecule_object: the name of the protein in PyMOL
    Raises:
        FileNotFoundError: If file not found.
        ValueError: If PDB ID couldn't be found in PDB.

    """
    if not safeguard.Safeguard.check_filepath(pathlib.Path(f"{filepath}/{filename}")):
        raise FileNotFoundError
    try:
        cmd.fetch(code=molecule_object, type="pdb", path=filepath, file=filename)
    except pymol.CmdException:
        raise ValueError("PDB ID is invalid.")


def get_protein_from_pdb(pdb_id: str) -> "protein.Protein":
    """Fetches a protein from the PDB and creates a protein object."""
    if len(pdb_id) != 4:
        raise exception.IllegalArgumentError("PDB ID is invalid!")
    try:
        # PDB ID as input: the pdb file gets saved in a scratch directory where it gets deleted immediately
        cmd.fetch(pdb_id, type="pdb", path=constants.SCRATCH_DIR)
        graphic_operations.setup_default_session_graphic_settings()
        tmp_protein = protein.Protein(
            molecule_object=pdb_id,
            pdb_filepath=path_util.FilePath(pathlib.Path(f"{constants.SCRATCH_DIR}/{pdb_id}.pdb")),
        )
        return tmp_protein  # noqa: RET504
    except pymol.CmdException:
        tools.clean_scratch_folder()
        # TODO: add message that fetching the reference failed
    except FileNotFoundError:
        print("File could not be found.")


def save_protein_to_pdb_file(export_filepath: pathlib.Path, molecule_object: str) -> tuple[bool, pathlib.Path]:
    """Saves a protein from the current pymol session as a .pdb file.

    Args:
        export_filepath: the filepath to save the protein as pdb file.
        molecule_object: the name of the protein in PyMOL.

    Raises:
        NotADirectoryError: If directory is not found.
    """
    if not safeguard.Safeguard.check_filepath(export_filepath.parent):
        raise NotADirectoryError(f"The filepath {export_filepath.parent} does not exists.")
    # save the pdb file under the path (export_data_dir)
    tmp_filepath: pathlib.Path = pathlib.Path(f"{export_filepath.parent}/{molecule_object}.pdb")
    cmd.save(str(tmp_filepath))
    if os.path.exists(tmp_filepath):
        return True, tmp_filepath
    raise FileNotFoundError(f"The filepath {tmp_filepath} was not found!")


def load_pymol_session(pymol_session_filepath: pathlib.Path) -> None:
    """This function loads a pymol session file into the current pymol session.

    Args:
        pymol_session_filepath: the filepath of the session file
    Raises:
        FileNotFoundError: If file not found.

    """
    if not safeguard.Safeguard.check_filepath(pymol_session_filepath):
        raise FileNotFoundError
    cmd.load(str(pymol_session_filepath))
    graphic_operations.setup_default_session_graphic_settings()


def convert_pymol_session_to_base64_string(pymol_molecule_object: str) -> str:
    """This function converts a pymol session file into a base64 string.

    Args:
        pymol_molecule_object (str): PyMOL molecule object to be converted.
    """
    session_filepath = pathlib.Path(f"{constants.SCRATCH_DIR}/{pymol_molecule_object}_session.pse")
    cmd.save(session_filepath)
    base64_string = binary_data.create_base64_string_from_file(path_util.FilePath(session_filepath))
    os.remove(session_filepath)
    return base64_string
