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
"""Module for the project class."""
import collections
import logging
import pathlib
import platform
from datetime import datetime
from xml.etree import ElementTree
from xml.dom import minidom
from typing import TYPE_CHECKING
from PyQt5 import QtCore
import numpy as np
from Bio import SeqRecord

from pyssa.logging_pyssa import log_handlers
from pyssa.io_pyssa import filesystem_io
from pyssa.io_pyssa import safeguard
from pyssa.io_pyssa.xml_pyssa import element_names
from pyssa.io_pyssa.xml_pyssa import attribute_names
from pyssa.internal.data_structures.data_classes import basic_protein_info

if TYPE_CHECKING:
    from pyssa.internal.data_structures.data_classes import current_session
    from pyssa.internal.data_structures import protein
    from pyssa.internal.data_structures import protein_pair
    from pyssa.internal.data_structures import settings

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class Project:
    """Class for the projects used in the plugin."""

    # <editor-fold desc="Class attributes">
    """
    the id of the project in the database
    """
    _id: int
    """
    the name of the project
    """
    _project_name: str
    """
    the absolute path of the current workspace
    """
    _workspace: pathlib.Path
    """
    the used OS
    """
    _operating_system = platform.system()
    # """
    # all top layer folder paths of the project
    # """
    # folder_paths: dict[str, pathlib.Path]
    # _session_file_name: str = "session_file_model_s.pse"
    """
    a list of all protein objects of the project
    """
    proteins: list["protein.Protein"]
    """
    a list of all protein_pair objects of the project
    """
    protein_pairs: list["protein_pair.ProteinPair"]
    """
    a list of all sequence objects of the project
    """
    sequences: list[SeqRecord.SeqRecord]
    """
    the xml writer to create the project.xml file
    """
    xml_writer: QtCore.QXmlStreamWriter
    # </editor-fold>

    def __init__(self, a_project_name: str = "", a_workspace_path: pathlib.Path = pathlib.Path("")) -> None:
        """Constructor.

        Args:
            a_project_name (str): the name of the project
            a_workspace_path (pathlib.Path): the path of the workspace

        Raises:
            exception.IllegalArgumentError: if one of the arguments are None
            exception.DirectoryDoesNotExistError: if the workspace folder does not exist
        """
        # <editor-fold desc="Checks">
        safeguard.Safeguard.check_if_value_is_not_none(a_project_name, logger)
        safeguard.Safeguard.check_if_value_is_not_none(a_workspace_path, logger)

        # </editor-fold>

        self._project_name: str = a_project_name
        self._workspace: pathlib.Path = a_workspace_path
        self.proteins: list["protein.Protein"] = []
        self.protein_pairs: list["protein_pair.ProteinPair"] = []
        self.sequences: list[SeqRecord.SeqRecord] = []
        self.xml_writer: QtCore.QXmlStreamWriter = QtCore.QXmlStreamWriter()
        self.xml_writer.setAutoFormatting(True)

    def set_id(self, an_id: int) -> None:
        self._id = an_id

    def get_id(self):
        return self._id

    def set_workspace_path(self, a_workspace_path: pathlib.Path) -> None:
        """Setter for the workspace path.

        Args:
            a_workspace_path (pathlib.Path): the new workspace path

        Raises:
            ValueError: if one of the arguments are None or path does not exist
        """
        # <editor-fold desc="Checks">
        safeguard.Safeguard.check_if_value_is_not_none(a_workspace_path, logger)

        if not safeguard.Safeguard.check_filepath(a_workspace_path):
            msg = "The given workspace path does not exists!"
            logger.error(msg)
            raise ValueError(msg)

        # </editor-fold>

        self._workspace = a_workspace_path

    def add_existing_protein(self, a_protein: "protein.Protein") -> None:
        """Adds an existing protein object to the project.

        Args:
            a_protein (protein.Protein): an existing protein object

        Raises:
            ValueError: if one of the arguments are None
        """
        # <editor-fold desc="Checks">
        safeguard.Safeguard.check_if_value_is_not_none(a_protein, logger)

        # </editor-fold>

        self.proteins.append(a_protein)

    def add_protein_pair(self, a_protein_pair: "protein_pair.ProteinPair") -> None:
        """Adds an existing protein_pair object to the project.

        Args:
            a_protein_pair (protein_pair.ProteinPair): An existing protein pair object.

        Raises:
            ValueError: If one of the arguments are None
        """
        # <editor-fold desc="Checks">
        safeguard.Safeguard.check_if_value_is_not_none(a_protein_pair, logger)

        # </editor-fold>

        self.protein_pairs.append(a_protein_pair)

    def get_number_of_proteins(self) -> int:
        """Returns the number of proteins in the project.

        Returns:
            the number of proteins in the project.
        """
        return len(self.proteins)

    def get_project_name(self) -> str:
        """Getter for the project name.

        Returns:
            the project name
        """
        return self._project_name

    def get_project_xml_path(self) -> pathlib.Path:
        """Returns the path of the projects xml file.

        Returns:
            the path of the projects xml file
        """
        return pathlib.Path(f"{self._workspace}/{self.get_project_name()}.xml")

    def get_database_filepath(self) -> pathlib.Path:
        return pathlib.Path(f"{self._workspace}/{self.get_project_name()}.db")

    def serialize_project(self, a_filepath: pathlib.Path) -> None:
        """Serializes the project to a xml file.

        Args:
            a_filepath (pathlib.Path): the filepath of the project xml

        Raises:
            ValueError: if one of the arguments are None or path not exist
        """
        # <editor-fold desc="Checks">
        safeguard.Safeguard.check_if_value_is_not_none(a_filepath, logger)

        # </editor-fold>

        file = QtCore.QFile(str(a_filepath))
        if not file.open(QtCore.QFile.WriteOnly | QtCore.QFile.Text):
            print("Error: Cannot open file for writing")
            return
        self.xml_writer.setDevice(file)
        self.xml_writer.writeStartDocument()
        # Project itself
        self.xml_writer.writeStartElement("project")
        self.xml_writer.writeStartElement("project_info")
        self.xml_writer.writeAttribute("name", self._project_name)
        self.xml_writer.writeAttribute("workspace_path", str(self._workspace))
        self.xml_writer.writeAttribute("os", self._operating_system)
        self.xml_writer.writeEndElement()
        # Sequences
        if len(self.sequences) > 0:
            self.xml_writer.writeStartElement("sequences")
            for tmp_seq in self.sequences:
                self.xml_writer.writeStartElement("sequence")
                self.xml_writer.writeAttribute("name", tmp_seq.name)
                self.xml_writer.writeAttribute("seq", tmp_seq.seq)
                self.xml_writer.writeEndElement()
            self.xml_writer.writeEndElement()
        # Proteins
        if len(self.proteins) > 0:
            self.xml_writer.writeStartElement("proteins")
            for tmp_protein in self.proteins:
                tmp_protein.write_protein_to_xml_structure(self.xml_writer)
            self.xml_writer.writeEndElement()
        # Protein pairs
        if len(self.protein_pairs) > 0:
            self.xml_writer.writeStartElement("protein_pairs")
            for tmp_protein_pair in self.protein_pairs:
                tmp_protein_pair.serialize(self.xml_writer)
            self.xml_writer.writeEndElement()
        self.xml_writer.writeEndElement()
        self.xml_writer.writeEndDocument()
        file.close()

        # project_root = ElementTree.Element(element_names.PROJECT)
        # # setup project information tree
        # project_info = ElementTree.SubElement(project_root, element_names.PROJECT_INFO)
        # project_info.set(attribute_names.PROJECT_NAME, self._project_name)
        # project_info.set(attribute_names.PROJECT_WORKSPACE_PATH, str(self._workspace))
        # project_info.set(attribute_names.PROJECT_CREATION_DATE, self._date)
        # project_info.set(attribute_names.PROJECT_OS, self._operating_system)
        #
        # xml_proteins_element = ElementTree.SubElement(project_root, element_names.PROTEINS)
        # for tmp_protein in self.proteins:
        #     tmp_protein.serialize_protein(xml_proteins_element)
        # xml_protein_pairs_element = ElementTree.SubElement(project_root, element_names.PROTEIN_PAIRS)
        # for tmp_protein_pair in self.protein_pairs:
        #     tmp_protein_pair.serialize_protein_pair(xml_protein_pairs_element)
        # xml_sequences_element = ElementTree.SubElement(project_root, element_names.SEQUENCES)
        # for tmp_sequence in self.sequences:
        #     tmp_sequence_xml = ElementTree.SubElement(xml_sequences_element, element_names.SEQUENCE)
        #     tmp_sequence_xml.set(attribute_names.SEQUENCE_NAME, str(tmp_sequence.name))
        #     tmp_sequence_xml.set(
        #         attribute_names.SEQUENCE_SEQ,
        #         str(tmp_sequence.seq),
        #     )
        #
        # xml_string = minidom.parseString(ElementTree.tostring(project_root)).toprettyxml(indent="   ")
        # with open(a_filepath, "w", encoding="utf-8") as tmp_file:
        #     tmp_file.write(xml_string)

    @staticmethod
    def deserialize_project(a_filepath: pathlib.Path, app_settings: "settings.Settings") -> "Project":
        """Constructs the protein object from the xml file.

        Args:
            a_filepath (pathlib.Path): the filepath of the project xml file
            app_settings (settings.Settings): the settings object of the main window

        Raises:
            ValueError: if one of the arguments are None or path not exist

        Returns:
            a complete project object deserialized from a xml file
        """
        # <editor-fold desc="Checks">
        safeguard.Safeguard.check_if_value_is_not_none(a_filepath, logger)
        safeguard.Safeguard.check_filepath(a_filepath)
        safeguard.Safeguard.check_if_value_is_not_none(app_settings, logger)

        # </editor-fold>

        return filesystem_io.XmlDeserializer.deserialize_project(a_filepath, app_settings)

    def search_protein(self, a_protein_name: str) -> "protein.Protein":
        """Searches the project for a specific protein name.

        Args:
            a_protein_name (str): the name of the protein to search

        Raises:
            IllegalArgumentError: If the argument is None.
        """
        # <editor-fold desc="Checks">
        safeguard.Safeguard.check_if_value_is_not_none(a_protein_name, logger)

        # </editor-fold>

        for tmp_protein in self.proteins:
            if tmp_protein.get_molecule_object() == a_protein_name:
                return tmp_protein
        print(f"No matching protein with the name {a_protein_name} found.")  # noqa: RET503

    def search_protein_pair(self, a_protein_pair_name: str) -> "protein_pair.ProteinPair":
        """Searches all protein_pairs within the project and returns true if the project contains the pair.

        Args:
            a_protein_pair_name: the name of the protein pair to search

        Raises:
            ValueError: if one of the arguments are None

        Returns:
            if the protein pair is found, the protein pair
        """
        # <editor-fold desc="Checks">
        safeguard.Safeguard.check_if_value_is_not_none(a_protein_pair_name, logger)

        # </editor-fold>

        for tmp_protein_pair in self.protein_pairs:
            if tmp_protein_pair.name == a_protein_pair_name:
                return tmp_protein_pair
        print(f"No matching protein with the name {a_protein_pair_name} found.")  # noqa: RET503

    def save_pymol_session(self, a_current_session: "current_session.CurrentSession") -> None:
        """Saves the pymol session.

        Args:
            a_current_session (current_session.CurrentSession): the CurrentSession object which should be saved

        Raises:
            ValueError: if one of the arguments are None
        """
        # <editor-fold desc="Checks">
        safeguard.Safeguard.check_if_value_is_not_none(a_current_session, logger)

        # </editor-fold>

        if a_current_session.type == "protein":
            tmp_protein = self.search_protein(a_current_session.name)
            tmp_protein.pymol_session = a_current_session.session
        elif a_current_session.type == "protein_pair":
            tmp_protein_pair = self.search_protein_pair(a_current_session.name)
            tmp_protein_pair.pymol_session = a_current_session.session

    def check_if_protein_is_in_any_protein_pair(self, a_protein_name: str) -> bool:
        """Checks if a certain protein is part of an existing protein pair.

        Args:
            a_protein_name (str): the name of the protein to check

        Raises:
            ValueError: if one of the arguments are None
        """
        # <editor-fold desc="Checks">
        safeguard.Safeguard.check_if_value_is_not_none(a_protein_name, logger)

        # </editor-fold>

        protein_obj = self.search_protein(a_protein_name)
        for tmp_protein_pair in self.protein_pairs:
            if protein_obj == tmp_protein_pair.protein_1 or tmp_protein_pair.protein_2:
                return True
            if protein_obj.get_molecule_object() == tmp_protein_pair.protein_1.get_molecule_object().replace("_1", ""):
                return True
        return False

    def delete_specific_protein(self, a_protein_name: str) -> None:
        """Deletes a certain protein from the project based on the protein name.

        Args:
            a_protein_name (str): the name of the protein to delete

        Raises:
            ValueError: if one of the arguments are None
        """
        # <editor-fold desc="Checks">
        safeguard.Safeguard.check_if_value_is_not_none(a_protein_name, logger)

        # </editor-fold>

        protein_obj = self.search_protein(a_protein_name)
        if protein_obj in self.proteins:
            self.proteins.remove(protein_obj)
        else:
            raise ValueError("An argument is not in the list.")

    def delete_specific_protein_pair(self, a_protein_pair_name: str) -> None:
        """Deletes a certain protein from the project based on the protein name.

        Args:
            a_protein_pair_name (str): the name of the protein pair to delete

        Raises:
            ValueError: if one of the arguments are None
        """
        # <editor-fold desc="Checks">
        safeguard.Safeguard.check_if_value_is_not_none(a_protein_pair_name, logger)

        # </editor-fold>

        protein_pair_obj = self.search_protein_pair(a_protein_pair_name)
        if protein_pair_obj in self.protein_pairs:
            self.protein_pairs.remove(protein_pair_obj)
        else:
            raise ValueError("An argument is not in the list.")

    def convert_list_of_proteins_to_list_of_protein_infos(self) -> np.ndarray:
        """Converts the list of proteins into an array of basic_protein_info objects."""
        tmp_protein_infos: collections.deque = collections.deque()
        for tmp_protein in self.proteins:
            tmp_protein_infos.append(
                basic_protein_info.BasicProteinInfo(
                    tmp_protein.get_molecule_object(),
                    tmp_protein.get_id(),
                    self._project_name,
                ),
            )
        return np.array(list(tmp_protein_infos))
