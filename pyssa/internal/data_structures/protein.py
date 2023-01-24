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

import pymol
from pymol import cmd
from pyssa.io_pyssa import safeguard
from pyssa.internal.portal import pymol_io
from pyssa.internal.portal import protein_operations
import pyssa.internal.portal as graphic_operations
# from pyssa.internal.portal import graphic_operations


class Protein:
    """This class stores one protein in a PyMOL compatible form

    Attributes:
        molecule_object:
            The name of the protein which is also used within pymol.
        filepath:
            The filepath where the pdb file is stored.
        export_filepath:
            The filepath where results are saved. Default is None.
        selection:
            A pymol selection string which needs to conform with the selection algebra
            from pymol.
        sequence:
            The primary sequence of the protein.
        chains:
            A list of chains which occur in the protein.
        filename:
            The name of the file with extension
    """

    def __init__(self, molecule_object: str, filepath: pathlib.Path = None,
                 export_filepath: pathlib.Path = None):
        """Constructor.

        Args:
            molecule_object (str):
                name of the reference Protein in the pymol session
            filepath (str, optional):
                directory where the pdb files of both model and
                reference are stored
            export_filepath (str, optional):
                directory where all results related to the Protein
                will be stored.
                All subdirectories like ``images``, ``alignment_files``
                and ``distances`` will be created automatically.
                The export_data_dir will then function as parent directory.

        Raises:
            NotADirectoryError: If directory not found.
            FileNotFoundError: If file not found.
        """
        self.molecule_object = molecule_object
        self.filepath: pathlib.Path = filepath
        self.export_filepath: pathlib.Path = export_filepath
        self.filename = ""
        self.selection: str = ""
        self.sequence: str = ""
        self.chains: list[str] = []

        # argument test
        if molecule_object.find(".pdb") != -1:
            # pdb file is given
            self.molecule_object = molecule_object.replace(".pdb", "")
            self.filename = molecule_object
        else:
            # PDB ID is given
            self.filename = f"{self.molecule_object}.pdb"

        if filepath == "None":
            filepath = None
        if filepath is not None:
            if not os.path.exists(f"{self.filepath}"):
                raise NotADirectoryError(f"The path {filepath} was not "
                                         f"found.")
        if export_filepath == "None":
            export_filepath = None
        if export_filepath is not None:
            export_filepath = pathlib.Path(export_filepath)
            if not os.path.exists(f"{self.export_filepath}"):
                raise NotADirectoryError(f"The path {export_filepath} was not "
                                         f"found.")
        if not os.path.exists(pathlib.Path(f"{filepath}/{self.filename}")):
            path = pathlib.Path(f"{filepath}/{self.filename}")
            raise FileNotFoundError(f"Path {path} does not exists")

    def set_selection(self, selection: str) -> None:
        """This function sets a selection for the Protein object.

        Args:
            selection (str):
                A pymol conform selection as a string.

        Example:
            This is a pymol conform selection::

                "pymol_6omn////CA"
        """
        self.selection = selection

    def set_sequence(self, sequence=None):
        """This function sets the sequence for the protein, either automatically with the help
        of the get_fastastr method from PyMOL or directly through a function argument

        Args:
            sequence:
                string of the sequence which gets set for the protein
        """
        if sequence is None:
            cmd.reinitialize()
            if self.filepath is None:
                # this means a PDB ID was given
                pymol_io.fetch_protein_from_pdb(self)
                self.sequence = cmd.get_fastastr('all')
            else:
                # a .pdb file was given
                pymol_io.load_protein(self)
                self.sequence = cmd.get_fastastr('all')
        else:
            self.sequence = sequence

    def set_chains(self, chains=None):
        """This function sets the chains for the protein, either automatically with the help
        of the get_chains method from PyMOL or directly through a function argument

        """
        if chains is None:
            cmd.reinitialize()
            if self.filepath is None:
                # this means a PDB ID was given
                pymol_io.fetch_protein_from_pdb(self)
                self.chains = cmd.get_chains(self.molecule_object)
            else:
                # a .pdb file was given
                pymol_io.load_protein(self)
                self.chains = cmd.get_chains(self.molecule_object)
        else:
            self.chains = chains

    def create_selection_from_chains(self) -> str:
        """This function creates a selection with given chains

        Returns:
            selection (str):
                the selection string for pymol

        Raises:
            ValueError: If there is an empty chains list.
            ValueError: If the selection remains empty.
        """
        if len(self.chains) == 0:
            raise ValueError("No chains were added to the protein!")
        selection = ""
        seperator = ", "
        tmp_list = []
        for chain in self.chains:
            tmp_selection = f"/{self.molecule_object}//{chain}//CA"
            tmp_list.append(tmp_selection)
            selection = seperator.join(tmp_list)
        if selection == "":
            raise ValueError("The selection is still empty!")
        return selection

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
        pymol_io.load_protein(self)
        # if not safeguard.Safeguard.check_filepath(f"{self.filepath}/{self.filename}"):
        #     raise FileNotFoundError
        # cmd.load(f"{self.filepath}/{self.filename}", object=self.molecule_object)

    def show_resi_as_balls_and_sticks(self) -> None:
        graphic_operations.show_protein_selection_as_balls_and_sticks(self)

    def hide_resi_as_balls_and_sticks(self) -> None:
        graphic_operations.hide_protein_selection_as_balls_and_sticks(self)

    def zoom_resi_protein_position(self) -> None:
        graphic_operations.zoom_to_residue_in_protein_position(self)

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
