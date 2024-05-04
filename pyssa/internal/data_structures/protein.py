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

from auxiliary_pymol import auxiliary_pymol_client
from pyssa.controller import database_manager
from pyssa.io_pyssa import safeguard
from pyssa.internal.portal import protein_operations
from pyssa.internal.data_structures import selection, job
from pyssa.util import protein_util, enums
from pyssa.util import exception
from pyssa.io_pyssa.xml_pyssa import attribute_names
from pyssa.io_pyssa import bio_data
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants
from typing import TYPE_CHECKING, TextIO, Union, Any
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

    def add_protein_structure_data_from_pdb_db(self, a_pdb_id, the_main_socket, a_socket) -> None:
        """Adds protein structure data based on a protein from the PDB database."""
        tmp_pdb_filepath = pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{a_pdb_id}.pdb")
        bio_data.download_pdb_file(a_pdb_id, tmp_pdb_filepath)
        self.chains = protein_operations.get_protein_chains(
            tmp_pdb_filepath, the_main_socket, a_socket
        )

        tmp_job_description = job.GeneralPurposeJobDescription(
            enums.JobShortDescription.CONSOLIDATE_MOLECULE_OBJECT_TO_FIRST_STATE
        )
        tmp_job_description.setup_dict({enums.JobDescriptionKeys.PDB_FILEPATH.value: str(tmp_pdb_filepath)})
        print(tmp_job_description.job_information)
        tmp_reply = auxiliary_pymol_client.send_request_to_auxiliary_pymol(
            the_main_socket, a_socket, tmp_job_description
        )
        print(tmp_reply)
        self._pdb_data = bio_data.parse_pdb_file(tmp_reply["data"])
        #self.check_states_and_reduce_to_one_state_if_necessary()
        try:
            os.remove(str(pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{a_pdb_id}.pdb")))
        except Exception as e:
            logger.error(f"Could not delete pdb file! Ran into error: {e}")

    def add_protein_structure_data_from_local_pdb_file(self, a_filepath: pathlib.Path, the_main_socket, a_socket) -> None:
        """Adds protein structure data based on a protein from the local filesystem."""
        self.chains = protein_operations.get_protein_chains(
            a_filepath, the_main_socket, a_socket
        )
        tmp_job_description = job.GeneralPurposeJobDescription(
            enums.JobShortDescription.CONSOLIDATE_MOLECULE_OBJECT_TO_FIRST_STATE
        )
        tmp_job_description.setup_dict({enums.JobDescriptionKeys.PDB_FILEPATH.value: str(a_filepath)})
        tmp_reply = auxiliary_pymol_client.send_request_to_auxiliary_pymol(
            the_main_socket, a_socket, tmp_job_description
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
        i = the_last_chain_id
        for tmp_chain in self.chains:
            tmp_chain.set_id(i)
            i += 1

    def create_new_pymol_session(self, the_main_socket, the_general_purpose_socket) -> None:
        """Creates a new pymol session by loading the protein into PyMOL.

        Raises:
            UnableToCreatePdbFileError:If the pdb file cannot be created.
            UnableToOpenPdbFileError:If the pdb file cannot be opened.
            UnableToLoadPdbFileError:If the pdb file cannot be loaded.
        """
        # <editor-fold desc="Checks">
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
            the_main_socket, the_general_purpose_socket, tmp_job_description
        )
        self.pymol_session = tmp_reply["data"][0]

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

    def clean_protein(self, the_main_socket, the_general_purpose_socket):
        tmp_job_description = job.GeneralPurposeJobDescription(
            enums.JobShortDescription.CLEAN_PROTEIN_UPDATE_STRUCTURE)
        tmp_job_description.setup_dict(
            {
                enums.JobDescriptionKeys.PYMOL_SESSION.value: str(self.pymol_session),
                enums.JobDescriptionKeys.PROTEIN_NAME.value: str(self._pymol_molecule_object)
            }
        )
        tmp_reply = auxiliary_pymol_client.send_request_to_auxiliary_pymol(
            the_main_socket, the_general_purpose_socket, tmp_job_description
        )
        if tmp_reply["data"][0].find("ERROR") != -1:
            raise ValueError(f"{tmp_reply['data'][0]}; {tmp_reply['data'][1]}")

        self.pymol_session, tmp_pdb_filepath = tmp_reply["data"]
        tmp_pdb_data = bio_data.parse_pdb_file(tmp_pdb_filepath)
        self._pdb_data = tmp_pdb_data

    def get_object_as_dict_for_database(self) -> dict:
        """Serializes the protein object to a dict which can be inserted into a database."""
        return {
            enums.DatabaseEnum.PROTEIN_NAME.value: self._pymol_molecule_object,
            enums.DatabaseEnum.PROTEIN_PYMOL_SESSION.value: self.pymol_session,
            "project_id": self.db_project_id,
        }

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
