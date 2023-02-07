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
import shutil

import pymol
from pymol import cmd
from pyssa.io_pyssa import safeguard
from pyssa.internal.portal import pymol_io
from pyssa.internal.portal import protein_operations
from pyssa.internal.portal import graphic_operations
from pyssa.internal.data_structures import selection
from pyssa.util import protein_util
from pyssa.util import pyssa_keys
from pyssa.io_pyssa import filesystem_io
from pyssa.logging_pyssa import log_handlers
from pyssa.io_pyssa import path_util
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pyssa.internal.data_structures import sequence
    from pyssa.internal.data_structures import chain


logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class Protein:
    """This class stores one protein in a PyMOL compatible form"""

    # <editor-fold desc="Class attributes">
    """
    the name of the protein which is also used within pymol
    """
    _pymol_molecule_object: str
    """
    a pymol conform selection 
    """
    pymol_selection: selection.Selection
    """
    a list of chains which occur in the protein
    """
    chains: list['chain.Chain']
    """
    the full filepath where the pdb file is stored
    """
    pdb_filepath: path_util.FilePath
    """
    a directory where all results related to the protein will be stored
    """
    export_dirname: pathlib.Path
    """
    the full filepath where the session file is stored
    """
    pymol_session_filepath: path_util.FilePath

    # </editor-fold>

    def __init__(self, molecule_object: str, proteins_dirname: pathlib.Path, pdb_filepath: path_util.FilePath) -> None:
        """Constructor.

        Args:
            molecule_object (str):
                the name of the protein which is also used within pymol
            proteins_dirname (Path):
                the directory where the protein subdirs will be created
            pdb_filepath (FilePath, optional):
                the filepath of the pdb file

        Raises:
            NotADirectoryError: If directory not found.
            FileNotFoundError: If file not found.
        """
        # <editor-fold desc="Checks">
        if not safeguard.Safeguard.check_if_value_is_not_none(molecule_object) or molecule_object == "":
            logger.error("An argument is illegal.")
            raise ValueError("An argument is illegal.")
        if not safeguard.Safeguard.check_if_value_is_not_none(proteins_dirname):
            logger.error("An argument is illegal.")
            raise ValueError("An argument is illegal.")
        # if not safeguard.Safeguard.check_filepath(proteins_dirname):
        #     logger.error("The directory does not exist.")
        #     raise NotADirectoryError("The directory does not exist.")
        if not safeguard.Safeguard.check_if_value_is_not_none(pdb_filepath.get_filepath()):
            logger.error("An argument is illegal.")
            raise ValueError("An argument is illegal.")
        if not safeguard.Safeguard.check_filepath(pdb_filepath.get_filepath()):
            logger.error("PDB file was not found.")
            raise NotADirectoryError("PDB file was not found.")

        # </editor-fold>

        self._pymol_molecule_object = pdb_filepath.get_filename()
        protein_dirname = pathlib.Path(f"{proteins_dirname}/{self._pymol_molecule_object}")
        self.protein_subdirs = {
            pyssa_keys.PROTEINS_SUBDIR: pathlib.Path(proteins_dirname),
            pyssa_keys.PROTEIN_SUBDIR: pathlib.Path(f"{protein_dirname}"),
            pyssa_keys.PROTEIN_SEQUENCE_SUBDIR: pathlib.Path(f"{protein_dirname}/sequence"),
            pyssa_keys.PROTEIN_PDB_SUBDIR: pathlib.Path(f"{protein_dirname}/pdb"),
            pyssa_keys.PROTEIN_SESSION_SUBDIR: pathlib.Path(f"{protein_dirname}/session"),
            pyssa_keys.PROTEIN_RESULTS_SUBDIR: pathlib.Path(f"{protein_dirname}/results"),
            pyssa_keys.PROTEIN_OBJECTS_SUBDIR: pathlib.Path(f"{protein_dirname}/.objects"),
        }
        for key in self.protein_subdirs:
            if not os.path.exists(self.protein_subdirs.get(key)):
                os.mkdir(self.protein_subdirs.get(key))
        self.export_dirname = self.protein_subdirs.get(pyssa_keys.PROTEIN_RESULTS_SUBDIR)
        try:
            shutil.move(src=pdb_filepath.get_filepath(), dst=self.protein_subdirs.get(pyssa_keys.PROTEIN_PDB_SUBDIR))
        except shutil.Error as e:
            print(f"An error occurred while moving the file: {e}")

        self.pdb_filepath = path_util.FilePath(pathlib.Path(f"{self.protein_subdirs.get(pyssa_keys.PROTEIN_PDB_SUBDIR)}/{pdb_filepath.get_basename()}"))
        self.chains = protein_operations.get_protein_chains(self._pymol_molecule_object,
                                                            self.pdb_filepath.get_dirname(),
                                                            self.pdb_filepath.get_basename())
        self.pymol_selection = selection.Selection(self._pymol_molecule_object)
        self.pymol_selection.selection_string = ""
        self.load_protein_in_pymol()
        cmd.save(pathlib.Path(f"{self.protein_subdirs.get(pyssa_keys.PROTEIN_SESSION_SUBDIR)}/{self._pymol_molecule_object}_session.pse"))
        self.pymol_session_filepath = path_util.FilePath(pathlib.Path(f"{self.protein_subdirs.get(pyssa_keys.PROTEIN_SESSION_SUBDIR)}/{self._pymol_molecule_object}_session.pse"))

    def get_molecule_object(self) -> str:
        """This function gets the molecule object.

        Returns:
            the pymol_molecule_object
        """
        return self._pymol_molecule_object

    def set_molecule_object(self, value) -> None:
        """This function sets the molecule object.

        Args:
            value:
                a new molecule object
        """
        self._pymol_molecule_object = value
        for tmp_chain in self.chains:
            tmp_chain.chain_sequence.name = value
        self.pymol_selection.molecule_object = value

    def get_all_sequences(self) -> list['sequence.ProteinSequence']:
        tmp_sequences = []
        for tmp_chain in self.chains:
            tmp_sequences.append(tmp_chain.chain_sequence)
        return tmp_sequences

    def get_protein_sequences(self) -> list['sequence.ProteinSequence']:
        return protein_operations.get_protein_sequences_from_protein(self._pymol_molecule_object, self.chains)

    def load_protein_in_pymol(self) -> None:
        src = self.pdb_filepath.get_dirname()
        src2 = self.pdb_filepath.get_filename()
        pymol_io.load_protein(self.pdb_filepath.get_dirname(), self.pdb_filepath.get_basename(), self._pymol_molecule_object)

    def color_protein_in_pymol(self) -> None:
        # TODO: needs to be implemented
        pass

    def set_selections_from_chains_ca(self) -> None:
        """This function sets a selection based on the chains of the protein. The selection selects only the alpha-C's.

        """
        self.pymol_selection.set_selections_from_chains_ca(protein_util.filter_chains_for_protein_chains(self.chains))

    def set_selection_without_chains_ca(self) -> None:
        """This function sets a selection without any chains of the protein. The selection selects only the alpha-C's.

        """
        self.pymol_selection.set_selections_without_chains_ca()

    def clean_pdb_file(self) -> None:
        """This function cleans a pdb file from the PDB

        Raises:
            AttributeError: If no export directory is given
        """
        # argument test
        if self.export_dirname == "":
            raise AttributeError("A export directory must be defined!")

        pymol_io.fetch_protein_from_pdb(self.pdb_filepath.get_dirname(), self.pdb_filepath.get_filename(), self._pymol_molecule_object)
        protein_operations.remove_solvent_molecules_in_protein()
        protein_operations.remove_organic_molecules_in_protein()
        # check if path exists where the data will be exported,
        # if not the directory will be created
        if not os.path.exists(f"{self.export_dirname}"):
            os.mkdir(f"{self.export_dirname}")
        # save the pdb file under the path (export_data_dir)
        pymol_io.save_protein_to_pdb_file(self.export_dirname, self._pymol_molecule_object)

    def clean_protein(self, new_protein=False):
        cmd.reinitialize()
        self.load_protein_in_pymol()
        if new_protein is False:
            try:
                protein_operations.remove_solvent_molecules_in_protein()
                protein_operations.remove_organic_molecules_in_protein()
            except pymol.CmdException:
                return
            os.remove(self.pdb_filepath.get_filepath())
            pymol_io.save_protein_to_pdb_file(self.export_dirname, self._pymol_molecule_object)
        else:
            clean_prot = self.duplicate_protein()
            cmd.reinitialize()
            clean_prot.load_protein_in_pymol()
            protein_operations.remove_solvent_molecules_in_protein()
            protein_operations.remove_organic_molecules_in_protein()
            clean_prot.molecule_object = f"{clean_prot.molecule_object}_cleaned"
            clean_prot.filename = f"{clean_prot.molecule_object}.pdb"
            pymol_io.save_protein_to_pdb_file(clean_prot.export_dirname, clean_prot._pymol_molecule_object)
            return clean_prot

    def show_resi_as_balls_and_sticks(self) -> None:
        graphic_operations.show_protein_selection_as_balls_and_sticks(self.pymol_selection.selection_string)

    def hide_resi_as_balls_and_sticks(self) -> None:
        graphic_operations.hide_protein_selection_as_balls_and_sticks(self.pymol_selection.selection_string)

    def zoom_resi_protein_position(self) -> None:
        graphic_operations.zoom_to_residue_in_protein_position(self.pymol_selection.selection_string)

    def serialize_protein(self) -> None:
        """This function serialize the protein object

        """
        protein_serializer = filesystem_io.ObjectSerializer(self, self.protein_subdirs.get(pyssa_keys.PROTEIN_OBJECTS_SUBDIR), self._pymol_molecule_object)
        protein_dict = {
            "molecule_object": self.get_molecule_object(),
            "proteins_dirname": str(self.protein_subdirs.get(pyssa_keys.PROTEINS_SUBDIR)),
            "import_data_dir": str(self.pdb_filepath.get_filepath()),
            "export_data_dir": str(self.export_dirname),
            "filename": str(self.pdb_filepath.get_filename()),
            "selection": self.pymol_selection.selection_string,
            "chains": protein_util.get_chains_as_list_of_tuples(self.chains),
            "pymol_session_file": str(self.pymol_session_filepath.get_filepath())
        }
        protein_serializer.set_custom_object_dict(protein_dict)
        protein_serializer.serialize_object()

    @staticmethod
    def deserialize_protein(protein_obj_json_file: path_util.FilePath):
        """This function constructs the protein object from
        the json file

        Returns:
            a complete protein object deserialized from a json file
        """
        return filesystem_io.ObjectDeserializer(protein_obj_json_file.get_dirname(),
                                                protein_obj_json_file.get_filename()).deserialize_protein()

    def duplicate_protein(self):
        tmp_protein = Protein(molecule_object=self._pymol_molecule_object, proteins_dirname=self.protein_subdirs.get(pyssa_keys.PROTEIN_SUBDIR), pdb_filepath=self.pdb_filepath.get_filepath())
        # tmp_protein.pymol_selection.selection_string = self.pymol_selection.selection_string
        # tmp_protein.chains = self.chains
        # # TODO: create new session file for duplicate
        return tmp_protein
