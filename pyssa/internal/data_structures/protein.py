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
"""Module for the protein class"""
import json
import os
import pathlib
import logging
import pymol
from pymol import cmd
from pyssa.io_pyssa import safeguard
from pyssa.internal.portal import pymol_io
from pyssa.internal.portal import protein_operations
from pyssa.internal.portal import graphic_operations
from pyssa.internal.data_structures import selection
from pyssa.util import protein_util
from pyssa.util import types
from pyssa.logging_pyssa import log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class Protein:
    """This class stores one protein in a PyMOL compatible form"""

    # <editor-fold desc="Class attributes">
    """
    the name of the protein which is also used within pymol
    """
    molecule_object: str
    """
    a pymol conform selection 
    """
    selection: types.SELECTION
    """
    the primary sequence of the protein
    """
    sequence: types.PROTEIN_SEQUENCE
    """
    a list of chains which occur in the protein
    """
    chains: list[types.CHAIN] = []
    """
    the filepath where the pdb file is stored
    """
    filepath: pathlib.Path
    """
    the name of the file with extension
    """
    filename: str
    """
    a directory where all results related to the protein will be stored
    """
    export_filepath: pathlib.Path

    # </editor-fold>

    def __init__(self, molecule_object: str, filepath: pathlib.Path = "", export_filepath: pathlib.Path = "") -> None:
        """Constructor.

        Args:
            molecule_object (str):
                the name of the protein which is also used within pymol
            filepath (Path, optional):
                directory where the pdb files of both model and
                reference are stored
            export_filepath (Path, optional):
                directory where all results related to the Protein
                will be stored.
                All subdirectories like ``images``, ``alignment_files``
                and ``distances`` will be created automatically.
                The export_data_dir will then function as parent directory.

        Raises:
            NotADirectoryError: If directory not found.
            FileNotFoundError: If file not found.
        """
        # <editor-fold desc="Checks">
        if not safeguard.Safeguard.check_if_value_is_not_none(molecule_object) or molecule_object == "":
            logger.error("An argument is illegal.")
            raise ValueError("An argument is illegal.")
        if filepath != "":
            if not safeguard.Safeguard.check_if_value_is_not_none(filepath):
                logger.error("An argument is illegal.")
                raise ValueError("An argument is illegal.")
            if not safeguard.Safeguard.check_filepath(filepath):
                logger.error("The directory does not exist.")
                raise NotADirectoryError("The directory does not exist.")
        if export_filepath != "":
            if not safeguard.Safeguard.check_if_value_is_not_none(export_filepath):
                logger.error("An argument is illegal.")
                raise ValueError("An argument is illegal.")
            if not safeguard.Safeguard.check_filepath(export_filepath):
                logger.error("The directory does not exist.")
                raise NotADirectoryError("The directory does not exist.")

        # </editor-fold>

        self.molecule_object, self.filename = protein_util.check_if_protein_is_from_file_or_id(molecule_object)
        # check if pdb file exists
        if not safeguard.Safeguard.check_filepath(pathlib.Path(f"{filepath}/{self.filename}")):
            logger.error("PDB file was not found.")
            raise FileNotFoundError("PDB file was not found.")
        self.filepath = filepath
        self.export_filepath = export_filepath
        self.chains = protein_operations.get_protein_chains(self.molecule_object, self.filepath, self.filename)
        self.sequence = protein_operations.get_protein_sequences_from_protein(self.molecule_object, self.chains)
        self.selection = selection.Selection(self.molecule_object)

    def set_selections_from_chains_ca(self):
        """This function sets a selection based on the chains of the protein. The selection selects only the alpha-C's.

        """
        self.selection.set_selections_from_chains_ca(protein_util.filter_chains_for_protein_chains(self.chains))

    def set_selection_without_chains_ca(self):
        """This function sets a selection without any chains of the protein. The selection selects only the alpha-C's.

        """
        self.selection.set_selections_without_chains_ca()

    def clean_pdb_file(self) -> None:
        """This function cleans a pdb file from the PDB

        Raises:
            AttributeError: If no export directory is given
        """
        # argument test
        if self.export_filepath == "":
            raise AttributeError("A export directory must be defined!")

        pymol_io.fetch_protein_from_pdb(self)
        protein_operations.remove_solvent_molecules_in_protein()
        protein_operations.remove_organic_molecules_in_protein()
        # check if path exists where the data will be exported,
        # if not the directory will be created
        if not os.path.exists(f"{self.export_filepath}"):
            os.mkdir(f"{self.export_filepath}")
        # save the pdb file under the path (export_data_dir)
        pymol_io.save_protein_to_pdb_file(self)

    def clean_protein(self, new_protein=False):
        cmd.reinitialize()
        self.load_protein()
        if new_protein is False:
            try:
                protein_operations.remove_solvent_molecules_in_protein()
                protein_operations.remove_organic_molecules_in_protein()
            except pymol.CmdException:
                return
            os.remove(f"{self.filepath}/{self.filename}")
            pymol_io.save_protein_to_pdb_file(self)
        else:
            clean_prot = self.duplicate_protein()
            cmd.reinitialize()
            clean_prot.load_protein()
            protein_operations.remove_solvent_molecules_in_protein()
            protein_operations.remove_organic_molecules_in_protein()
            clean_prot.molecule_object = f"{clean_prot.molecule_object}_cleaned"
            clean_prot.filename = f"{clean_prot.molecule_object}.pdb"
            pymol_io.save_protein_to_pdb_file(self)
            return clean_prot

    def load_protein(self) -> None:
        pymol_io.load_protein(self.filepath, self.filename, self.molecule_object)

    def show_resi_as_balls_and_sticks(self) -> None:
        graphic_operations.show_protein_selection_as_balls_and_sticks(self.selection)

    def hide_resi_as_balls_and_sticks(self) -> None:
        graphic_operations.hide_protein_selection_as_balls_and_sticks(self.selection)

    def zoom_resi_protein_position(self) -> None:
        graphic_operations.zoom_to_residue_in_protein_position(self.selection)

    def serialize_protein(self, filepath, filename) -> None:
        """This function serialize the protein object

        """
        if not os.path.exists(filepath):
            print(f"The filepath: {filepath} does not exists!")
            return
        if filename.find(".pdb"):
            filename = filename.replace(".pdb", "")
        protein_dict = self.__dict__
        update = {
            "filepath": str(self.filepath),
            "export_data_dir": str(self.export_filepath),
        }
        protein_dict.update(update)
        protein_file = open(pathlib.Path(f"{filepath}/{filename}.json"), "w", encoding="utf-8")
        json.dump(protein_dict, protein_file, indent=4)

    @staticmethod
    def deserialize_protein(protein_obj_json_file):
        """This function constructs the protein object from
        the json file

        Returns:
            a complete protein object deserialized from a json file
        """
        try:
            protein_obj_file = open(protein_obj_json_file, "r", encoding="utf-8")
        except FileNotFoundError:
            print(f"There is no valid protein json file under: {protein_obj_json_file}")
            return
        protein_dict = json.load(protein_obj_file)
        if protein_dict.get("export_data_dir") == "None":
            update = {"export_data_dir": None}
            protein_dict.update(update)
        tmp_protein = Protein(protein_dict.get("molecule_object"),
                              protein_dict.get("filepath"),
                              protein_dict.get("export_data_dir"))
        tmp_protein.filename = protein_dict.get("filename")
        tmp_protein.set_sequence(protein_dict.get("sequence"))
        tmp_protein.set_selection(protein_dict.get("selection"))
        tmp_protein.set_chains(protein_dict.get("chains"))
        return tmp_protein

    def duplicate_protein(self):
        tmp_protein = Protein(self.molecule_object, self.filepath, self.export_filepath)
        tmp_protein.set_sequence(self.sequence)
        tmp_protein.set_selection(self.selection)
        tmp_protein.set_chains(self.chains)
        return tmp_protein
