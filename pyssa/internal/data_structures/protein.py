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
"""Module for the protein class."""
import os
import pathlib
import logging
import uuid
import pymol
from Bio import PDB
from PyQt5 import QtCore

from pyssa.controller import database_manager
from pyssa.io_pyssa import safeguard
from pyssa.internal.portal import pymol_io
from pyssa.internal.portal import protein_operations
from pyssa.internal.portal import graphic_operations
from pyssa.internal.data_structures import selection
from pyssa.util import protein_util, enums
from pyssa.util import exception
from pyssa.io_pyssa.xml_pyssa import element_names
from pyssa.io_pyssa.xml_pyssa import attribute_names
from pyssa.io_pyssa import binary_data
from pyssa.io_pyssa import bio_data
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants
from pyssa.io_pyssa import path_util
from typing import TYPE_CHECKING, TextIO, Union, Any
from xml.etree import ElementTree
from pyssa.internal.data_structures import chain


if TYPE_CHECKING:
    from pyssa.internal.data_structures import sequence


logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class Protein:
    """Stores one protein in a PyMOL compatible form."""

    # <editor-fold desc="Class attributes">
    """
    the unique identifier of the protein
    """
    _id: int
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
    chains: list["chain.Chain"] = []
    """
    a base64 string of the pymol session
    """
    pymol_session: str = ""
    """
    a list of pdb information
    """
    _pdb_data: list[dict]
    """
    the project id from the database
    """
    db_project_id: int

    # </editor-fold>

    def __init__(
        self,
        molecule_object: str
    ) -> None:
        """Constructor.

        Args:
            molecule_object (str): the name of the protein which is also used within pymol
            the_project_id (int): the id of the project in the database

        Raises:
            NotADirectoryError: If directory not found.
            FileNotFoundError: If file not found.
        """
        # <editor-fold desc="Checks">
        safeguard.Safeguard.check_if_value_is_not_none(molecule_object, logger)
        if molecule_object == "":
            logger.error("An argument is illegal.")
            raise ValueError("An argument is illegal.")
        if molecule_object.find(" "):
            molecule_object = molecule_object.replace(" ", "_")

        # </editor-fold>

        self._pymol_molecule_object = molecule_object
        self.pymol_selection = selection.Selection(self._pymol_molecule_object)
        self.pymol_selection.selection_string = ""
        self.chains: list["chain.Chain"] = []
        self._pdb_data = []

        #self._id = uuid.uuid4()
        #self.pdb_cache_path = pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{self._id}.pdb")
        #if pdb_data is not None:
        #    self._pdb_data = pdb_data
        # if pdb_filepath == "" and pdb_xml_string != "":
        #     self._pymol_molecule_object = molecule_object
        #     if not os.path.exists(constants.CACHE_PROTEIN_DIR):
        #         os.mkdir(constants.CACHE_PROTEIN_DIR)
        #     bio_data.convert_xml_string_to_pdb_file(pdb_xml_string, self.pdb_cache_path)
        #     self.chains = protein_operations.get_protein_chains(
        #         molecule_object,
        #         constants.CACHE_PROTEIN_DIR,
        #         f"{self._id}.pdb",
        #     )
        #     self._pdb_data = bio_data.convert_pdb_xml_string_to_list(pdb_xml_string)
        #     os.remove(f"{constants.CACHE_PROTEIN_DIR}/{self._id}.pdb")
        #     self.pymol_selection = selection.Selection(self._pymol_molecule_object)
        #     self.pymol_selection.selection_string = ""
        #     self.load_protein_in_pymol()
        #     if protein_operations.count_states_of_molecule_object(self._pymol_molecule_object) > 1:
        #         logger.info(f"Protein {molecule_object} has more than one state. Trying to consolidate to one state.")
        #         try:
        #             protein_operations.consolidate_molecule_object_to_first_state(self._pymol_molecule_object)
        #             tmp_pdb_cache_filepath = pathlib.Path(
        #                 f"{constants.CACHE_PROTEIN_DIR}/{self._pymol_molecule_object}.pdb",
        #             )
        #             pymol_io.save_protein_to_pdb_file(
        #                 tmp_pdb_cache_filepath,
        #                 self._pymol_molecule_object,
        #             )
        #             self._pdb_data = bio_data.convert_pdb_xml_string_to_list(
        #                 bio_data.convert_pdb_file_into_xml_element(path_util.FilePath(tmp_pdb_cache_filepath)),
        #             )
        #             self.load_protein_in_pymol()
        #         except Exception as e:
        #             logger.error(f"Protein states could not be consolidated! Ran into error: {e}")
        #         else:
        #             logger.info(f"Consolidation of protein {molecule_object} finished without any errors.")
        #     # saves pymol session into a base64 string
        #     self.pymol_session = pymol_io.convert_pymol_session_to_base64_string(self._pymol_molecule_object)
        # elif pdb_filepath != "" and pdb_xml_string == "":
        #     self._pymol_molecule_object = pdb_filepath.get_filename().replace(" ", "_")
        #     self.chains = protein_operations.get_protein_chains(
        #         molecule_object,
        #         pdb_filepath.get_dirname(),
        #         pdb_filepath.get_basename(),
        #     )
        #     self._pdb_data = bio_data.convert_pdb_xml_string_to_list(
        #         bio_data.convert_pdb_file_into_xml_element(pdb_filepath),
        #     )
        #     self.pymol_selection = selection.Selection(self._pymol_molecule_object)
        #     self.pymol_selection.selection_string = ""
        #     self.load_protein_in_pymol()
        #     if protein_operations.count_states_of_molecule_object(self._pymol_molecule_object) > 1:
        #         logger.info(f"Protein {molecule_object} has more than one state. Trying to consolidate to one state.")
        #         try:
        #             protein_operations.consolidate_molecule_object_to_first_state(self._pymol_molecule_object)
        #             tmp_pdb_cache_filepath = pathlib.Path(
        #                 f"{constants.CACHE_PROTEIN_DIR}/{self._pymol_molecule_object}.pdb",
        #             )
        #             pymol_io.save_protein_to_pdb_file(
        #                 tmp_pdb_cache_filepath,
        #                 self._pymol_molecule_object,
        #             )
        #             self._pdb_data = bio_data.convert_pdb_xml_string_to_list(
        #                 bio_data.convert_pdb_file_into_xml_element(path_util.FilePath(tmp_pdb_cache_filepath)),
        #             )
        #         except Exception as e:
        #             logger.error(f"Protein states could not be consolidated! Ran into error: {e}")
        #         else:
        #             logger.info(f"Consolidation of protein {molecule_object} finished without any errors.")
        #
        #     # saves pymol session into a base64 string
        #     self.pymol_session = pymol_io.convert_pymol_session_to_base64_string(self._pymol_molecule_object)
        # elif pdb_filepath == "" and pdb_xml_string == "":
        #     self._pymol_molecule_object = molecule_object
        #     self.pymol_selection = selection.Selection(self._pymol_molecule_object)
        #     self.pymol_selection.selection_string = ""
        # else:
        #     raise ValueError("Function has too many arguments.")

    def add_protein_structure_data_from_pdb_db(self, a_pdb_id) -> None:
        """Adds protein structure data based on a protein from the PDB database."""
        tmp_pdb_filepath = pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{a_pdb_id}.pdb")
        bio_data.download_pdb_file(a_pdb_id, tmp_pdb_filepath)
        self.chains = protein_operations.get_protein_chains(
            self._pymol_molecule_object,
            tmp_pdb_filepath.parent,
            tmp_pdb_filepath.name,
        )
        self._pdb_data = bio_data.parse_pdb_file(tmp_pdb_filepath)
        self.check_states_and_reduce_to_one_state_if_necessary()
        try:
            os.remove(str(pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{a_pdb_id}.pdb")))
            os.remove(str(pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{self._pymol_molecule_object}.pdb")))
        except Exception as e:
            logger.error(f"Could not delete pdb file! Ran into error: {e}")

    def add_protein_structure_data_from_local_pdb_file(self, a_filepath: pathlib.Path) -> None:
        """Adds protein structure data based on a protein from the local filesystem."""
        self.chains = protein_operations.get_protein_chains(
            self._pymol_molecule_object,
            a_filepath.parent,
            a_filepath.name,
        )
        self._pdb_data = bio_data.parse_pdb_file(a_filepath)
        self.check_states_and_reduce_to_one_state_if_necessary()
        try:
            os.remove(f"{constants.CACHE_PROTEIN_DIR}/{self._pymol_molecule_object}.pdb")
        except Exception as e:
            logger.error(f"Could not delete pdb file! Ran into error: {e}")

    def add_id_to_all_chains(self, the_last_chain_id: int) -> None:
        i = the_last_chain_id + 1
        for tmp_chain in self.chains:
            tmp_chain.set_id(i)
            i += 1

    def create_new_pymol_session(self) -> None:
        """Creates a new pymol session by loading the protein into pymol."""
        self.load_protein_in_pymol()

    def save_pymol_session_as_base64_string(self) -> None:
        """Sets the active pymol session as base64 string into the protein object."""
        self.pymol_session = pymol_io.convert_pymol_session_to_base64_string(self._pymol_molecule_object)

    def check_states_and_reduce_to_one_state_if_necessary(self) -> None:
        """Checks the state of the protein in pymol and deletes all states which are not the first if necessary."""
        if protein_operations.count_states_of_molecule_object(self._pymol_molecule_object) > 1:
            logger.info(f"Protein {self._pymol_molecule_object} has more than one state. Trying to consolidate to one state.")
            try:
                protein_operations.consolidate_molecule_object_to_first_state(self._pymol_molecule_object)
                tmp_pdb_cache_filepath = pathlib.Path(
                    f"{constants.CACHE_PROTEIN_DIR}/{self._pymol_molecule_object}.pdb",
                )
                pymol_io.save_protein_to_pdb_file(
                    tmp_pdb_cache_filepath,
                    self._pymol_molecule_object,
                )
                self._pdb_data = bio_data.parse_pdb_file(tmp_pdb_cache_filepath)
            except Exception as e:
                logger.error(f"Protein states could not be consolidated! Ran into error: {e}")
            else:
                logger.info(f"Consolidation of protein {self._pymol_molecule_object} finished without any errors.")

    def set_id(self, value):
        self._id = value

    def get_molecule_object(self) -> str:
        """This function gets the molecule object.

        Returns:
            the pymol_molecule_object
        """
        return self._pymol_molecule_object

    def set_molecule_object(self, value: str) -> None:
        """This function sets the molecule object.

        Args:
            value:
                a new molecule object
        """
        self._pymol_molecule_object = value
        for tmp_chain in self.chains:
            if tmp_chain.chain_type != "non_protein_chain":
                tmp_chain.chain_sequence.name = value
        self.pymol_selection.molecule_object = value

    def get_all_sequences(self) -> list["sequence.Sequence"]:
        """Gets all sequences of the protein as list of sequences."""
        tmp_sequences = []
        for tmp_chain in self.chains:
            tmp_sequences.append(tmp_chain.chain_sequence)
        return tmp_sequences

    def get_protein_sequences(self) -> list["sequence.Sequence"]:
        """Gets all protein sequences of the protein as list of sequences."""
        return protein_operations.get_protein_sequences_from_protein(self._pymol_molecule_object, self.chains)

    def get_sequences_based_on_chain_letter(self, chain_letter: str) -> "sequence.Sequence":
        """Gets the sequence based on the provided chain letter.

        Args:
            chain_letter: the letter of which the sequence is to be needed.
        """
        for tmp_chain in self.chains:
            if chain_letter == tmp_chain.chain_letter:
                return tmp_chain.chain_sequence
        return None

    def get_pdb_data(self) -> list:
        """Gets the pdb information of the protein."""
        return self._pdb_data

    def set_pdb_data(self, pdb_data: list) -> None:
        self._pdb_data = pdb_data

    def get_id(self) -> int:
        """Gets the id of the protein."""
        return self._id

    def write_fasta_file(self, a_filepath: pathlib.Path) -> None:
        """Writes a fasta file the specified filepath.

        Raises:
            IllegalArgumentError: If the file is None or the file could not be found.
            DirectoryDoesNotExistError: If the parent directory of the filepath could not be found.
        """
        # <editor-fold desc="Checks">
        if a_filepath is None or not os.path.exists(a_filepath):
            exception.IllegalArgumentError("")
        if not os.path.exists(a_filepath.parent):
            raise exception.DirectoryDoesNotExistError("")

        # </editor-fold>

        try:
            fasta_file: TextIO = open(f"{a_filepath}/{self._pymol_molecule_object}.fasta", "w")
        except OSError:
            raise exception.FastaFilesNotCreatedError("")

        # Writes fasta file header
        fasta_file.write(f">{self._pymol_molecule_object}\n")

        # Writes fasta file sequence content
        i: int = 0
        seq_objs: list[sequence.Sequence] = self.get_protein_sequences()
        for tmp_sequence in seq_objs:
            if i == len(self.get_protein_sequences()) - 1:
                # should be the last entry
                fasta_file.write(tmp_sequence.sequence)
            else:
                fasta_file.write(f"{tmp_sequence.sequence}:")
            i += 1

        logger.info(f"Fasta file for sequence {self._pymol_molecule_object} written.")
        fasta_file.close()

    def append_chain(self, chain_name: str, chain_sequence: "sequence.Sequence", chain_type: str) -> None:
        """Appends a chain to the protein.

        Args:
            chain_name: name of the chain to append to the protein.
            chain_sequence: sequence of the chain to append to the protein.
            chain_type: type of chain to append to the protein.
        """
        # TODO: add check if chain or any type of chain information already exists
        self.chains.append(chain.Chain(chain_letter=chain_name, chain_sequence=chain_sequence, chain_type=chain_type))

    def add_chain_names_to_chains(self) -> None:
        """Adds the chain names to the individual chains."""
        i = 0
        for tmp_chain in self.chains:
            if tmp_chain.chain_letter == "":
                tmp_chain.chain_letter = constants.chain_dict[i]
                i += 1
            else:
                raise ValueError("Chain name exists.")

    def get_chain_by_letter(self, a_letter: str) -> "chain.Chain":
        """Returns the chain of a certain letter"""
        for tmp_chain in self.chains:
            if tmp_chain.chain_letter == a_letter:
                return tmp_chain
        return chain.Chain("", "", "")

    def load_protein_in_pymol(self) -> None:
        """Load a protein in PyMOL.

        Raises:
            UnableToCreatePdbFileError:If the pdb file cannot be created.
            UnableToOpenPdbFileError:If the pdb file cannot be opened.
            UnableToLoadPdbFileError:If the pdb file cannot be loaded.
        """
        # <editor-fold desc="Checks">
        pdb_filepath = pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{self._pymol_molecule_object}.pdb")
        if not os.path.exists(constants.CACHE_PROTEIN_DIR):
            os.mkdir(constants.CACHE_PROTEIN_DIR)

        # </editor-fold>

        if len(self._pdb_data) == 0:
            raise ValueError("No pdb data in current object!")

        try:
            bio_data.build_pdb_file(self._pdb_data, pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{self._pymol_molecule_object}.pdb"))
        except exception.IllegalArgumentError:
            logger.error(f"The argument pdb data is not usable: {self._pdb_data}.")
            raise exception.UnableToCreatePdbFileError("")
        except exception.DirectoryNotFoundError:
            logger.error(f"The argument pdb_filepath is illegal: {pdb_filepath}!")
            raise exception.UnableToCreatePdbFileError("")
        except PermissionError:
            logger.error(f"The argument pdb_filepath is illegal: {pdb_filepath}!")
            raise exception.UnableToCreatePdbFileError("")
        except exception.UnableToOpenFileError:
            logger.error("pdb file could not be opened for writing.")
            raise exception.UnableToOpenFileError("")

        try:
            pymol_io.load_protein(constants.CACHE_PROTEIN_DIR, f"{self._pymol_molecule_object}.pdb", self._pymol_molecule_object)
        except exception.UnableToLoadProteinError:
            logger.error("Protein can not be loaded in PyMOL!")
            raise exception.UnableToLoadProteinError

    def load_protein_pymol_session(self) -> None:
        """Loads the protein in the pymol session based on the base64 data."""
        tmp_session_path = pathlib.Path(
            f"{constants.CACHE_PYMOL_SESSION_DIR}/{self._pymol_molecule_object}_session.pse",
        )
        binary_data.write_binary_file_from_base64_string(tmp_session_path, self.pymol_session)
        pymol_io.load_pymol_session(tmp_session_path)

    def color_protein_in_pymol(self, color: str, a_selection: str) -> None:
        """This function colors the protein, in the given color.

        Args:
            color:
                a PyMOL color string
            a_selection:
                a PyMOL conform selection string
        Notes:
            The correct selection string needs to be in a_protein before running this function!!!

        """
        self.pymol_selection.set_custom_selection(a_selection)
        graphic_operations.color_protein(color, self.pymol_selection.selection_string)

    def set_selections_from_chains_ca(self) -> None:
        """Sets a selection based on the chains of the protein.

        Notes:
            The selection selects only the alpha-C's.
        """
        self.pymol_selection.set_selections_from_chains_ca(protein_util.filter_chains_for_protein_chains(self.chains))

    def set_selection_without_chains_ca(self) -> None:
        """Sets a selection without any chains of the protein.

        Notes:
            The selection selects only the alpha-C's.
        """
        self.pymol_selection.set_selections_without_chains_ca()

    def clean_pdb_file(self) -> None:
        """This function cleans a pdb file from the PDB.

        Raises:
            AttributeError: If no export directory is given
        """
        # argument test
        if self.export_dirname == "":
            raise AttributeError("A export directory must be defined!")

        pymol_io.fetch_protein_from_pdb(
            self.pdb_filepath.get_dirname(),
            self.pdb_filepath.get_filename(),
            self._pymol_molecule_object,
        )
        protein_operations.remove_solvent_molecules_in_protein()
        protein_operations.remove_organic_molecules_in_protein()
        # check if path exists where the data will be exported,
        # if not the directory will be created
        if not os.path.exists(f"{self.export_dirname}"):
            os.mkdir(f"{self.export_dirname}")
        # save the pdb file under the path (export_data_dir)
        pymol_io.save_protein_to_pdb_file(self.export_dirname, self._pymol_molecule_object)

    def clean_protein(self, new_protein: bool = False):  # noqa: ANN201 #TODO: needs to be redone
        """Cleans the protein from all sugar and solvent molecules.

        Args:
            new_protein: a flag to determine whether to create a new protein object.
        """
        if new_protein is False:
            try:
                self.load_protein_pymol_session()
                protein_operations.remove_solvent_molecules_in_protein()
                protein_operations.remove_organic_molecules_in_protein()
            except pymol.CmdException:
                return  # noqa: RET502 #TODO: needs more thoughts
            # tmp_full_pdb_path = pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{self._id}.pdb")
            tmp_was_successful, tmp_pdb_filepath = pymol_io.save_protein_to_pdb_file(constants.CACHE_PROTEIN_DIR, str(self._id))
            if tmp_was_successful:
                self._pdb_data = bio_data.parse_pdb_file(tmp_pdb_filepath)
                logger.debug(self._pdb_data)
            else:
                logger.error("The protein could not be cleaned, because the new pdb file could not be found!")
                raise RuntimeError("The protein could not be cleaned, because the new pdb file could not be found!")
        else:
            clean_prot = self.duplicate_protein()
            # clean_prot.load_protein_in_pymol()
            clean_prot.load_protein_pymol_session()
            protein_operations.remove_solvent_molecules_in_protein()
            protein_operations.remove_organic_molecules_in_protein()
            clean_prot_molecule_object = f"{clean_prot.get_molecule_object()}_cleaned"
            clean_prot.set_molecule_object(clean_prot_molecule_object)
            tmp_full_pdb_path = pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{clean_prot.get_id()}.pdb")
            pymol_io.save_protein_to_pdb_file(constants.CACHE_PROTEIN_DIR, str(clean_prot.get_id()))
            cleaned_prot = Protein(
                molecule_object=clean_prot_molecule_object,
            )
            return cleaned_prot  # noqa: RET504 #TODO: needs more thoughts

    def show_resi_as_balls_and_sticks(self) -> None:
        """Shows the residues of the selection as sticks representation."""
        graphic_operations.show_protein_selection_as_balls_and_sticks(self.pymol_selection.selection_string)

    def hide_resi_as_balls_and_sticks(self) -> None:
        """Hides the residues of the selection as sticks representation."""
        graphic_operations.hide_protein_selection_as_balls_and_sticks(self.pymol_selection.selection_string)

    def zoom_resi_protein_position(self) -> None:
        """Zoomes to the residues of the selection."""
        graphic_operations.zoom_to_residue_in_protein_position(self.pymol_selection.selection_string)

    def serialize_protein(self, xml_proteins_element: ElementTree.Element) -> None:
        """Serializes the protein object."""
        tmp_protein = ElementTree.SubElement(xml_proteins_element, element_names.PROTEIN)
        tmp_protein.set(attribute_names.ID, str(self._id))
        tmp_protein.set(attribute_names.PROTEIN_MOLECULE_OBJECT, self._pymol_molecule_object)
        tmp_protein.set(attribute_names.PROTEIN_SELECTION, self.pymol_selection.selection_string)
        tmp_protein.append(bio_data.convert_pdb_data_list_to_xml_string(self._pdb_data))
        tmp_session_data = ElementTree.SubElement(tmp_protein, element_names.PROTEIN_SESSION)
        tmp_session_data.set(attribute_names.PROTEIN_SESSION, self.pymol_session)

    def get_object_as_dict_for_database(self) -> dict:
        """Serializes the protein object to a dict which can be inserted into a database."""
        return {
            enums.DatabaseEnum.PROTEIN_NAME.value: self._pymol_molecule_object,
            enums.DatabaseEnum.PROTEIN_PYMOL_SESSION.value: self.pymol_session,
            "project_id": self.db_project_id,
        }

    def write_protein_to_xml_structure(self, an_xml_writer: QtCore.QXmlStreamWriter):
        an_xml_writer.writeStartElement('protein')
        an_xml_writer.writeAttribute('id', str(self._id))
        an_xml_writer.writeAttribute('pymol_molecule_object', str(self._pymol_molecule_object))
        an_xml_writer.writeAttribute('pymol_selection', str(self.pymol_selection.selection_string))
        # Chains
        an_xml_writer.writeStartElement('chains')
        for tmp_chain in self.chains:
            tmp_chain.serialize(an_xml_writer)
        an_xml_writer.writeEndElement()
        # PDB data
        an_xml_writer.writeStartElement('pdb_data')
        for tmp_atom in self._pdb_data:
            an_xml_writer.writeStartElement('atom')
            an_xml_writer.writeCharacters(tmp_atom)
            an_xml_writer.writeEndElement()
        an_xml_writer.writeEndElement()
        # Session data
        an_xml_writer.writeStartElement('session_data')
        an_xml_writer.writeAttribute('session', self.pymol_session)
        an_xml_writer.writeEndElement()
        an_xml_writer.writeEndElement()

    def duplicate_protein(self):  # noqa: ANN201
        """Duplicates the protein object."""
        tmp_protein = Protein(
            molecule_object=self._pymol_molecule_object,
        )
        logger.debug(tmp_protein.chains[0])
        # tmp_protein.pymol_selection.selection_string = self.pymol_selection.selection_string
        # tmp_protein.chains = self.chains
        # # TODO: create new session file for duplicate
        return tmp_protein

    def set_all_attributes(self, attrib_dict: dict, pdb_data: list, pymol_session: str) -> None:
        """Sets all attributes of the protein."""
        """The unique identifier of the protein."""
        self._id = attrib_dict[attribute_names.ID]
        self._pymol_molecule_object = attrib_dict[attribute_names.PROTEIN_MOLECULE_OBJECT]
        self.pymol_selection.selection_string = attrib_dict[attribute_names.PROTEIN_SELECTION]
        self._pdb_data = pdb_data
        self.pymol_session = pymol_session

    # def create_plain_text_memory_mirror(self) -> list[Union[tuple[str, str], tuple[str, Any]]]:
    #     """Creates a plain text memory mirror of the protein object."""
    #     mirror = [
    #         ("_id", str(self._id)),
    #         ("_pymol_molecule_object", str(self._pymol_molecule_object)),
    #         ("pymol_selection.molecule_object", str(self.pymol_selection.molecule_object)),
    #         ("pymol_selection.selection_string", str(self.pymol_selection.selection_string)),
    #     ]
    #     for tmp_chain in self.chains:
    #         mirror.append(("tmp_chain.chain", str(tmp_chain.chain_letter)))
    #         mirror.append(("tmp_chain.chain_type", str(tmp_chain.chain_type)))
    #         mirror.append(("tmp_chain.chain_sequence", str(tmp_chain.chain_sequence)))
    #
    #     for tmp_pdb_line in self._pdb_data:
    #         mirror.append(("_pdb_data", tmp_pdb_line))
    #
    #     mirror.append(("pdb_filepath", str(self.pdb_filepath)))
    #     mirror.append(("fasta_filepath", str(self.fasta_filepath)))
    #     mirror.append(("export_dirname", str(self.export_dirname)))
    #     mirror.append(("pymol_session_filepath", str(self.pymol_session_filepath)))
    #     mirror.append(("pymol_session", str(self.pymol_session)))
    #     return mirror
