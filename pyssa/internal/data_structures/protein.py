#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/zielesny/PySSA>
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

import zmq

from auxiliary_pymol import auxiliary_pymol_client
from pyssa.internal.portal import protein_operations
from pyssa.internal.data_structures import selection, job
from pyssa.util import enums
from pyssa.util import exception
from pyssa.io_pyssa import bio_data
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants
from typing import TYPE_CHECKING, TextIO, Optional
from pyssa.internal.data_structures import chain


if TYPE_CHECKING:
    from pyssa.internal.data_structures import sequence


logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class Protein:
    """Stores one protein in a PyMOL compatible form."""

    # <editor-fold desc="Class attributes">
    _id: int
    """The unique identifier of the protein."""
    
    _pymol_molecule_object: str
    """The name of the protein which is also used within pymol."""
    
    pymol_selection: selection.Selection
    """A pymol conform selection."""
    
    chains: list["chain.Chain"] = []
    """A list of chains which occur in the protein."""
    
    pymol_session: str = ""
    """A base64 string of the pymol session."""
    
    _pdb_data: list[dict]
    """A list of pdb information."""
    
    db_project_id: int
    """The project id from the database."""

    # </editor-fold>

    def __init__(
        self,
        molecule_object: str,
    ) -> None:
        """Constructor.

        Args:
            molecule_object (str): The name of the protein which is also used within pymol.

        Raises:
            exception.IllegalArgumentError: If `molecule_object` is either None or an empty string.
        """
        # <editor-fold desc="Checks">
        if molecule_object is None or molecule_object == "":
            logger.error("molecule_object is either None or an empty string.")
            raise exception.IllegalArgumentError("molecule_object is either None or an empty string.")
        if molecule_object.find(" "):
            molecule_object = molecule_object.replace(" ", "_")

        # </editor-fold>

        self._pymol_molecule_object = molecule_object
        self.pymol_selection = selection.Selection(self._pymol_molecule_object)
        self.pymol_selection.selection_string = ""
        self.chains: list["chain.Chain"] = []
        self._pdb_data = []

    def add_protein_structure_data_from_pdb_db(self, a_pdb_id: str, the_main_socket: zmq.Socket, a_socket: zmq.Socket) -> None:
        """Adds protein structure data based on a protein from the PDB database.

        Args:
            a_pdb_id (str): The ID of the protein structure in the PDB database.
            the_main_socket (zmq.Socket): The main ZeroMQ socket for communication.
            a_socket (zmq.Socket): An auxiliary ZeroMQ socket for communication.
        
        Raises:
            exception.IllegalArgumentError: If any of the arguments are None.
        """
        # <editor-fold desc="Checks">
        if a_pdb_id is None:
            logger.error("a_pdb_id is None.")
            raise exception.IllegalArgumentError("a_pdb_id is None.")
        if the_main_socket is None:
            logger.error("the_main_socket is None.")
            raise exception.IllegalArgumentError("the_main_socket is None.")
        if a_socket is None:
            logger.error("a_socket is None.")
            raise exception.IllegalArgumentError("a_socket is None.")
        
        # </editor-fold>
        
        tmp_pdb_filepath = pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{a_pdb_id}.pdb")
        bio_data.download_pdb_file(a_pdb_id, tmp_pdb_filepath)
        self.chains = protein_operations.get_protein_chains(
            tmp_pdb_filepath, the_main_socket, a_socket,
        )

        tmp_job_description = job.GeneralPurposeJobDescription(
            enums.JobShortDescription.CONSOLIDATE_MOLECULE_OBJECT_TO_FIRST_STATE,
        )
        tmp_job_description.setup_dict({enums.JobDescriptionKeys.PDB_FILEPATH.value: str(tmp_pdb_filepath)})
        print(tmp_job_description.job_information)
        tmp_reply = auxiliary_pymol_client.send_request_to_auxiliary_pymol(
            the_main_socket, a_socket, tmp_job_description,
        )
        print(tmp_reply)
        self._pdb_data = bio_data.parse_pdb_file(tmp_reply["data"])
        #self.check_states_and_reduce_to_one_state_if_necessary()
        try:
            os.remove(str(pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{a_pdb_id}.pdb")))
        except Exception as e:
            logger.error(f"Could not delete pdb file! Ran into error: {e}")

    def add_protein_structure_data_from_local_pdb_file(self, a_filepath: pathlib.Path, the_main_socket: zmq.Socket, a_socket: zmq.Socket) -> None:
        """Adds protein structure data based on a protein from the local filesystem.

        Args:
            a_filepath (pathlib.Path): The filepath of the PDB file to be processed.
            the_main_socket (zmq.Socket): The main socket for communication.
            a_socket (zmq.Socket): A secondary socket for communication.

        Raises:
            exception.IllegalArgumentError: If any of the arguments are None.
            FileNotFoundError: If the PDB file specified by `a_filepath` cannot be found.
        """
        # <editor-fold desc="Checks">
        if a_filepath is None:
            logger.error("a_filepath is None.")
            raise exception.IllegalArgumentError("a_filepath is None.")
        if the_main_socket is None:
            logger.error("the_main_socket is None.")
            raise exception.IllegalArgumentError("the_main_socket is None.")
        if a_socket is None:
            logger.error("a_socket is None.")
            raise exception.IllegalArgumentError("a_socket is None.")
        
        # </editor-fold>
        
        self.chains = protein_operations.get_protein_chains(
            a_filepath, the_main_socket, a_socket,
        )
        tmp_job_description = job.GeneralPurposeJobDescription(
            enums.JobShortDescription.CONSOLIDATE_MOLECULE_OBJECT_TO_FIRST_STATE,
        )
        tmp_job_description.setup_dict({enums.JobDescriptionKeys.PDB_FILEPATH.value: str(a_filepath)})
        tmp_reply = auxiliary_pymol_client.send_request_to_auxiliary_pymol(
            the_main_socket, a_socket, tmp_job_description,
        )
        if tmp_reply["data"] == "":
            logger.error(f"Could not find pdb file! Ran into error: {tmp_reply}")
            raise FileNotFoundError()
        self._pdb_data = bio_data.parse_pdb_file(tmp_reply["data"])
        try:
            os.remove(pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{self._pymol_molecule_object}.pdb"))
        except Exception as e:
            logger.error(f"Could not delete pdb file! Ran into error: {e}")

    def add_id_to_all_chains(self, the_last_chain_id: int) -> None:
        """Adds an ID to each chain in the list of chains.

        Args:
            the_last_chain_id: int. The ID to start assigning to each chain.
        
        Raises:
            exception.IllegalArgumentError: If the_last_chain_id is None.
        """
        # <editor-fold desc="Checks">
        if the_last_chain_id is None:
            logger.error("the_last_chain_id is None.")
            raise exception.IllegalArgumentError("the_last_chain_id is None.")
        
        # </editor-fold>
        
        i = the_last_chain_id
        for tmp_chain in self.chains:
            tmp_chain.set_id(i)
            i += 1

    def create_new_pymol_session(self, the_main_socket: zmq.Socket, the_general_purpose_socket: zmq.Socket) -> None:
        """Creates a new pymol session by loading the protein into PyMOL.

        Args:
            the_main_socket (zmq.Socket): The main ZMQ socket used for communication with the auxiliary PyMOL process.
            the_general_purpose_socket (zmq.Socket): The general purpose ZMQ socket used for communication with the auxiliary PyMOL process.

        Raises:
            exception.IllegalArgumentError: If the_main_socket or the_general_purpose_socket is None.
            ValueError: If there is no pdb data in the current object.
            UnableToCreatePdbFileError: If an error occurs while building the pdb file.
            UnableToOpenFileError: If the pdb file could not be opened for writing.

        """
        # <editor-fold desc="Checks">
        if the_main_socket is None:
            logger.error("the_main_socket is None.")
            raise exception.IllegalArgumentError("the_main_socket is None.")
        if the_general_purpose_socket is None:
            logger.error("the_general_purpose_socket is None.")
            raise exception.IllegalArgumentError("the_general_purpose_socket is None.")
        pdb_filepath = pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{self._pymol_molecule_object}.pdb")
        if not os.path.exists(constants.CACHE_PROTEIN_DIR):
            os.mkdir(constants.CACHE_PROTEIN_DIR)
        if len(self._pdb_data) == 0:
            raise ValueError("No pdb data in current object!")
        
        # </editor-fold>

        try:
            bio_data.build_pdb_file(self._pdb_data,
                                    pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{self._pymol_molecule_object}.pdb"))
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

        tmp_job_description = job.GeneralPurposeJobDescription(
            enums.JobShortDescription.CREATE_NEW_PROTEIN_PYMOL_SESSION)
        tmp_job_description.setup_dict({enums.JobDescriptionKeys.PDB_FILEPATH.value: str(
            pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{self._pymol_molecule_object}.pdb"))})
        tmp_reply = auxiliary_pymol_client.send_request_to_auxiliary_pymol(
            the_main_socket, the_general_purpose_socket, tmp_job_description,
        )
        self.pymol_session = tmp_reply["data"][0]

    def set_id(self, a_value: int) -> None:
        """Sets the id of the protein object.

        Args:
            a_value (int): The new value to set for the ID.
            
        Raises:
            exception.IllegalArgumentError: If a_value is either None or has a value less than 0.
        """
        # <editor-fold desc="Checks">
        if a_value is None or a_value < 0:
            logger.error("a_value is either None or has a value less than 0.")
            raise exception.IllegalArgumentError("a_value is either None or has a value less than 0.")
        
        # </editor-fold>
        
        self._id = a_value

    def get_molecule_object(self) -> str:
        """Gets the molecule object.

        Returns:
            the pymol_molecule_object
        """
        return self._pymol_molecule_object

    def set_molecule_object(self, a_value: str) -> None:
        """Sets the molecule object.

        Args:
            a_value (str): A new molecule object.
        
        Raises:
            exception.IllegalArgumentError: If a_value is either None or has a value less than 0.
        """
        # <editor-fold desc="Checks">
        if a_value is None or a_value == "":
            logger.error("a_value_to_check is either None or an empty string.")
            raise exception.IllegalArgumentError("a_value_to_check is either None or an empty string.")
        
        # </editor-fold>
        
        self._pymol_molecule_object = a_value
        for tmp_chain in self.chains:
            if tmp_chain.chain_type != "non_protein_chain":
                tmp_chain.chain_sequence.name = a_value
        self.pymol_selection.molecule_object = a_value

    def get_all_sequences(self) -> list["sequence.Sequence"]:
        """Gets all sequences of the protein as list of sequences.
        
        Returns:
            A list of sequence objects.
        """
        tmp_sequences: list["sequence.Sequence"] = []
        for tmp_chain in self.chains:
            tmp_sequences.append(tmp_chain.chain_sequence)
        return tmp_sequences

    def get_protein_sequences(self) -> list["sequence.Sequence"]:
        """Gets all protein sequences of the protein as list of sequences.
        
        Returns:
            A list of sequence objects.
        """
        return protein_operations.get_protein_sequences_from_protein(self._pymol_molecule_object, self.chains)

    def get_sequences_based_on_chain_letter(self, chain_letter: str) -> Optional["sequence.Sequence"]:
        """Gets the sequence based on the provided chain letter.

        Args:
            chain_letter: the letter of which the sequence is to be needed.
        
        Raises:
            exception.IllegalArgumentError: If chain_letter is either None or an empty string.
        """
        # <editor-fold desc="Checks">
        if chain_letter is None or chain_letter == "":
            logger.error("chain_letter is either None or an empty string.")
            raise exception.IllegalArgumentError("chain_letter is either None or an empty string.")
        
        # </editor-fold>
        
        for tmp_chain in self.chains:
            if chain_letter == tmp_chain.chain_letter:
                return tmp_chain.chain_sequence
        return None

    def get_pdb_data(self) -> list:
        """Gets the pdb information of the protein.
        
        Returns:
            A list with the pdb data.
        """
        return self._pdb_data

    def set_pdb_data(self, pdb_data: list) -> None:
        """Sets the pdb information of the protein.
        
        Args:
            pdb_data: A list containing the pdb data.
        
        Raises:
            exception.IllegalArgumentError: If pdb_data is None.
        """
        # <editor-fold desc="Checks">
        if pdb_data is None:
            logger.error("pdb_data is None.")
            raise exception.IllegalArgumentError("pdb_data is None.")
        
        # </editor-fold>
        
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
            chain_name (str): name of the chain to append to the protein.
            chain_sequence (sequence.Sequence): sequence of the chain to append to the protein.
            chain_type (str): type of chain to append to the protein.
        
        Raises:
            exception.IllegalArgumentError: If chain_name is either None or an empty string.
            exception.IllegalArgumentError: If chain_sequence is None.
            exception.IllegalArgumentError: If chain_type is either None or an empty string.
        """
        # <editor-fold desc="Checks">
        if chain_name is None or chain_name == "":
            logger.error("chain_name is either None or an empty string.")
            raise exception.IllegalArgumentError("chain_name is either None or an empty string.")
        if chain_sequence is None:
            logger.error("chain_sequence is None.")
            raise exception.IllegalArgumentError("chain_sequence is None.")
        if chain_type is None or chain_type == "":
            logger.error("chain_type is either None or an empty string.")
            raise exception.IllegalArgumentError("chain_type is either None or an empty string.")
        # TODO: add check if chain or any type of chain information already exists
        
        # </editor-fold>
        
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

    def get_chain_by_letter(self, a_letter: str) -> Optional["chain.Chain"]:
        """Returns the chain of a certain letter.
        
        Args:
            a_letter (str): The letter of the chain to be returned.
        
        Returns:
            The chain of a certain letter or None if the letter could not be found.
        """
        # <editor-fold desc="Checks">
        if a_letter is None or a_letter == "":
            logger.error("a_letter is either None or an empty string.")
            raise exception.IllegalArgumentError("a_letter is either None or an empty string.")
        
        # </editor-fold>
        
        for tmp_chain in self.chains:
            if tmp_chain.chain_letter == a_letter:
                return tmp_chain
        return None

    # def set_selections_from_chains_ca(self) -> None:
    #     """Sets a selection based on the chains of the protein.
    # 
    #     Notes:
    #         The selection selects only the alpha-C's.
    #     """
    #     self.pymol_selection.set_selections_from_chains_ca(protein_util.filter_chains_for_protein_chains(self.chains))

    # def set_selection_without_chains_ca(self) -> None:
    #     """Sets a selection without any chains of the protein.
    # 
    #     Notes:
    #         The selection selects only the alpha-C's.
    #     """
    #     self.pymol_selection.set_selections_without_chains_ca()

    def clean_protein(self, the_main_socket: zmq.Socket, the_general_purpose_socket: zmq.Socket) -> None:
        """Cleans the protein objects and sets the new pdb data into the instance.

        Args:
            the_main_socket (zmq.Socket): The main ZeroMQ socket for communication with the auxiliary pymol client.
            the_general_purpose_socket (zmq.Socket): The general purpose ZeroMQ socket for communication with the auxiliary pymol client.
        
        Raises:
            exception.IllegalArgumentError: If any of the arguments are None.
        """
        # <editor-fold desc="Checks">
        if the_main_socket is None:
            logger.error("the_main_socket is None.")
            raise exception.IllegalArgumentError("the_main_socket is None.")
        if the_general_purpose_socket is None:
            logger.error("the_general_purpose_socket is None.")
            raise exception.IllegalArgumentError("the_general_purpose_socket is None.")
        
        # </editor-fold>
        
        tmp_job_description = job.GeneralPurposeJobDescription(
            enums.JobShortDescription.CLEAN_PROTEIN_UPDATE_STRUCTURE)
        tmp_job_description.setup_dict(
            {
                enums.JobDescriptionKeys.PYMOL_SESSION.value: str(self.pymol_session),
                enums.JobDescriptionKeys.PROTEIN_NAME.value: str(self._pymol_molecule_object),
            },
        )
        tmp_reply = auxiliary_pymol_client.send_request_to_auxiliary_pymol(
            the_main_socket, the_general_purpose_socket, tmp_job_description,
        )
        if tmp_reply["data"][0].find("ERROR") != -1:
            raise ValueError(f"{tmp_reply['data'][0]}; {tmp_reply['data'][1]}")

        self.pymol_session, tmp_pdb_filepath = tmp_reply["data"]
        tmp_pdb_data = bio_data.parse_pdb_file(tmp_pdb_filepath)
        self._pdb_data = tmp_pdb_data

    def get_object_as_dict_for_database(self) -> dict:
        """Serializes the protein object to a dict which can be inserted into a database.
        
        Returns:
            A dict with the molecule_object, the pymol_session and the project id.
        """
        return {
            enums.DatabaseEnum.PROTEIN_NAME.value: self._pymol_molecule_object,
            enums.DatabaseEnum.PROTEIN_PYMOL_SESSION.value: self.pymol_session,
            "project_id": self.db_project_id,
        }

    # def duplicate_protein(self):  # noqa: ANN201
    #     """Duplicates the protein object."""
    #     tmp_protein = Protein(
    #         molecule_object=self._pymol_molecule_object,
    #     )
    #     logger.debug(tmp_protein.chains[0])
    #     # tmp_protein.pymol_selection.selection_string = self.pymol_selection.selection_string
    #     # tmp_protein.chains = self.chains
    #     # # TODO: create new session file for duplicate
    #     return tmp_protein

    # def set_all_attributes(self, attrib_dict: dict, pdb_data: list, pymol_session: str) -> None:
    #     """Sets all attributes of the protein."""
    #     """The unique identifier of the protein."""
    #     self._id = attrib_dict[attribute_names.ID]
    #     self._pymol_molecule_object = attrib_dict[attribute_names.PROTEIN_MOLECULE_OBJECT]
    #     self.pymol_selection.selection_string = attrib_dict[attribute_names.PROTEIN_SELECTION]
    #     self._pdb_data = pdb_data
    #     self.pymol_session = pymol_session
