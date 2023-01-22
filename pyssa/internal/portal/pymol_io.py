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
import pymol
from pymol import cmd
from pyssa.io_pyssa import safeguard
from pyssa.internal.data_structures import protein


def load_protein(protein_obj: protein.Protein) -> None:
    """This function loads a protein in pymol through a protein object.

    Args:
        protein_obj:
            object which is in instance of the protein class

    """
    if not safeguard.Safeguard.check_filepath(f"{protein_obj.filepath}/{protein_obj.filename}"):
        raise FileNotFoundError
    cmd.load(f"{protein_obj.filepath}/{protein_obj.filename}", object=protein_obj.molecule_object)


def fetch_protein_from_pdb(protein_obj: protein.Protein) -> None:
    """This function fetches a protein in pymol from the PDB.

    Args:
        protein_obj:
            object which is in instance of the protein class
    Raises:
        FileNotFoundError: If file not found.
        ValueError: If PDB ID couldn't be found in PDB.

    """
    if not safeguard.Safeguard.check_filepath(f"{protein_obj.filepath}/{protein_obj.filename}"):
        raise FileNotFoundError
    try:
        cmd.fetch(code=protein_obj.molecule_object, type="pdb", path=protein_obj.filepath, file=protein_obj.filename)
    except pymol.CmdException:
        raise ValueError("PDB ID is invalid.")


def save_protein_to_pdb_file(protein_obj: protein.Protein) -> None:
    """This function saves a protein from the current pymol session as a .pdb file.

    Args:
        protein_obj:
            object which is in instance of the protein class
    Raises:
        NotADirectoryError: If directory is not found.

    """
    if not safeguard.Safeguard.check_filepath(f"{protein_obj.export_filepath}"):
        raise NotADirectoryError(f"The filepath {protein_obj.export_filepath} does not exists.")
    # save the pdb file under the path (export_data_dir)
    cmd.save(f"{protein_obj.export_filepath}/{protein_obj.molecule_object}.pdb")
