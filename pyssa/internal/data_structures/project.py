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
"""Module for the project class"""
import json
import logging
import os
import pathlib
import platform
from datetime import datetime
from pathlib import Path
from util import gui_utils
from internal.data_structures import protein, protein_pair
from pyssa.internal.analysis_types import distance_analysis
from pyssa.logging_pyssa import log_handlers
from pyssa.io_pyssa import filesystem_io
from pyssa.io_pyssa.xml_pyssa import element_names
from pyssa.io_pyssa.xml_pyssa import attribute_names
from pyssa.util import constants
from xml.etree import ElementTree
from xml.dom import minidom
from pyssa.io_pyssa import binary_data


logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class Project:
    """This class is for the projects used in the plugin"""

    # <editor-fold desc="Class attributes">
    """
    the name of the project
    """
    _project_name: str
    """
    the absolute path of the current workspace
    """
    _workspace: pathlib.Path
    """
    the current date
    """
    _date = str(datetime.now())
    """
    the used OS
    """
    _operating_system = platform.system()
    """
    a list which contains all ATOM and HETATM lines of the .pdb file
    """
    _pdb_data: list[str]
    """
    all top layer folder paths of the project
    """
    folder_paths: dict[str, pathlib.Path]
    # _session_file_name: str = "session_file_model_s.pse"
    """
    a list of all protein objects of the project
    """
    proteins: list[protein.Protein]
    """
    a list of all protein_pair objects of the project
    """
    protein_pairs: list[protein_pair.ProteinPair]

    # </editor-fold>

    def __init__(self, project_name: str, workspace_path: pathlib.Path) -> None:
        """Constructor

        Args:
            project_name:
                name of the project
            workspace_path:
                path of the workspace
        """
        self._project_name: str = project_name
        self._workspace: pathlib.Path = workspace_path
        self.proteins: list[protein.Protein] = []
        self.protein_pairs: list[protein_pair.ProteinPair] = []
        self.create_folder_paths()

    def add_existing_protein(self, value_protein: protein.Protein) -> None:
        """This function adds an existing protein object to the project.

        Args:
            value_protein:
                name of the existing protein object
        """
        self.proteins.append(value_protein)

    def add_protein_pair(self, value_protein_pair: protein_pair.ProteinPair) -> None:
        """This function adds an existing protein_pair object to the project.

        Args:
            value_protein_pair:
                name of the existing protein_pairs object
        """
        self.protein_pairs.append(value_protein_pair)

    def create_folder_paths(self):
        self.folder_paths = {
            "project": Path(f"{self._workspace}/{self._project_name}"),
            "proteins": Path(f"{self._workspace}/{self._project_name}/proteins"),
            "protein_pairs": Path(f"{self._workspace}/{self._project_name}/protein_pairs"),
            "analysis": Path(f"{self._workspace}/{self._project_name}/protein_pairs/analysis"),
            "distance_analysis": Path(f"{self._workspace}/{self._project_name}/protein_pairs/analysis/distance_analysis"),
        }
        # Path(f"{self._workspace}/{self._project_name}/pdb"),
        # Path(f"{self._workspace}/{self._project_name}/results"),
        # Path(f"{self._workspace}/{self._project_name}/.objects/"),
        # Path(f"{self._workspace}/{self._project_name}/.objects/proteins"),
        # Path(f"{self._workspace}/{self._project_name}/.objects/protein_pairs"),
        # Path(f"{workspace_path}/{self._project_name}/results/alignment_files"),
        # Path(f"{workspace_path}/{self._project_name}/results/distance_csv"),
        # Path(f"{workspace_path}/{self._project_name}/results/images"),
        # Path(f"{workspace_path}/{self._project_name}/results/images/interesting_regions"),
        # Path(f"{workspace_path}/{self._project_name}/results/plots"),
        # Path(f"{workspace_path}/{self._project_name}/results/plots/distance_histogram"),
        # Path(f"{workspace_path}/{self._project_name}/results/plots/distance_plot"),
        # Path(f"{workspace_path}/{self._project_name}/results/sessions/"),

    def create_project_tree(self):
        for key in self.folder_paths:
            if not os.path.exists(self.folder_paths[key]):
                os.mkdir(self.folder_paths[key])
            else:
                raise IsADirectoryError

    def get_number_of_proteins(self):
        return len(self.proteins)

    def get_project_name(self):
        return self._project_name

    def get_project_xml_path(self):
        return pathlib.Path(f"{self._workspace}/{self.get_project_name()}.xml")

    # def serialize_project(self, filepath, filename) -> None:
    #     """This function serialize the protein object
    #
    #     """
    #     project_serializer = filesystem_io.ObjectSerializer(self, filepath, filename)
    #     project_dict = {
    #         '_project_name': str(self._project_name),
    #         '_workspace': str(self._workspace),
    #         '_date': str(self._date),
    #         '_operating_system': str(self._operating_system),
    #     }
    #     project_serializer.set_custom_object_dict(project_dict)
    #     project_serializer.serialize_object()
    #
    #     for tmp_protein in self.proteins:
    #         tmp_protein.serialize_protein()
    #     for tmp_protein_pair in self.protein_pairs:
    #         tmp_protein_pair.serialize_protein_pair()
    #     for tmp_distance_analysis in self.distance_analysis:
    #         tmp_distance_analysis.serialize_distance_analysis()

    def serialize_project(self, filepath) -> None:
        project_root = ElementTree.Element(element_names.PROJECT)
        # setup project information tree
        project_info = ElementTree.SubElement(project_root, element_names.PROJECT_INFO)
        project_info.set(attribute_names.PROJECT_NAME, self._project_name)
        project_info.set(attribute_names.PROJECT_WORKSPACE_PATH, str(self._workspace))
        project_info.set(attribute_names.PROJECT_CREATION_DATE, self._date)
        project_info.set(attribute_names.PROJECT_OS, self._operating_system)

        xml_proteins_element = ElementTree.SubElement(project_root, element_names.PROTEINS)
        for tmp_protein in self.proteins:
            tmp_protein.serialize_protein(xml_proteins_element)
        xml_protein_pairs_element = ElementTree.SubElement(project_root, element_names.PROTEIN_PAIRS)
        for tmp_protein_pair in self.protein_pairs:
            tmp_protein_pair.serialize_protein_pair(xml_protein_pairs_element)
        # for i in range(len(self.protein_pairs)):
        #     logger.debug(i)
        #     logger.debug(self.distance_analysis)
        #     tmp_protein_pair = self.protein_pairs[i]
        #     xml_distance_analysis_element = tmp_protein_pair.serialize_protein_pair(xml_protein_pairs_element)
        #     self.distance_analysis[i].serialize_distance_analysis(xml_distance_analysis_element)

        # for tmp_protein_pair in self.protein_pairs:
        #     xml_distance_analysis_element = tmp_protein_pair.serialize_protein_pair(xml_protein_pairs_element)

        # xml_distance_analysis_element = ElementTree.SubElement()
        # for tmp_distance_analysis in self.distance_analysis:
        #     tmp_distance_analysis.serialize_distance_analysis()

        xmlstr = minidom.parseString(ElementTree.tostring(project_root)).toprettyxml(indent="   ")
        with open(filepath, "w") as f:
            f.write(xmlstr)

    @staticmethod
    def deserialize_project(filepath, app_settings):
        """This function constructs the protein object from
        the json file

        Returns:
            a complete project object deserialized from a json file
        """
        return filesystem_io.XmlDeserializer(filepath).deserialize_project(app_settings)

    def search_protein(self, protein_name):
        for tmp_protein in self.proteins:
            if tmp_protein.get_molecule_object() == protein_name:
                return tmp_protein
        print(f"No matching protein with the name {protein_name} found.")

    def search_protein_pair(self, protein_pair_name) -> 'protein_pair.ProteinPair':
        """This function searches all protein_pairs within the project and returns true if the project contains the pair

        Args:
            protein_pair_name:
                name of the pair to search
        Returns:
            True: is in the project
            False: is NOT in the project
        """
        for tmp_protein_pair in self.protein_pairs:
            if tmp_protein_pair.name == protein_pair_name:
                return tmp_protein_pair
        print(f"No matching protein with the name {protein_pair_name} found.")

    def get_specific_protein_pair(self, protein_pair_name):
        """This function gets a specific protein_pair by name from the project

        Args:
            protein_pair_name:
                name of the pair to get
        Returns:
            the actual protein pair
        """
        for tmp_protein_pair in self.protein_pairs:
            if tmp_protein_pair.name == protein_pair_name:
                return tmp_protein_pair
        print(f"No matching protein with the name {protein_pair_name} found.")

    def get_specific_protein_pair_tuple(self, protein_pair_name):
        """This function gets a specific protein_pair by name from the project

        Args:
            protein_pair_name:
                name of the pair to get
        Returns:
            the actual protein pair
        """
        for tmp_protein_pair in self.protein_pairs:
            if tmp_protein_pair[2].name == protein_pair_name:
                return tmp_protein_pair
        print(f"No matching protein with the name {protein_pair_name} found.")

    def delete_specific_protein(self, protein_name):
        protein_obj = self.search_protein(protein_name)
        if protein_obj in self.proteins:
            self.proteins.remove(protein_obj)
            self.serialize_project(self.get_project_xml_path())
        else:
            raise ValueError("An argument is not in the list.")

    def dump_project_to_file(self):
        current_time = datetime.now()
        filename = f"project_dump-{self._project_name}-{current_time.year}-{current_time.month:02d}-{current_time.day:02d}_{current_time.hour:02d}-{current_time.minute:02d}.ascii"
        memory_dump_file = open(f"{constants.SCRATCH_DIR}/{filename}", "w")
        memory_dump_file.write("----- Project information \n")
        memory_dump_file.write(f"_project_name: {self._project_name} \n")
        memory_dump_file.write(f"_workspace: {self._workspace} \n")
        memory_dump_file.write(f"_date: {self._date} \n")
        memory_dump_file.write(f"_operating_system: {self._operating_system} \n")
        for tmp_path in self.folder_paths:
            memory_dump_file.write(f"folder_paths: {tmp_path} \n")

        memory_dump_file.write("----- ----- ----- ----- \n")
        memory_dump_file.write("----- Protein information \n")
        for tmp_protein in self.proteins:
            for tmp_info in tmp_protein.create_plain_text_memory_mirror():
                memory_dump_file.write(f"{tmp_info} \n")

        memory_dump_file.write("----- ----- ----- ----- \n")
        memory_dump_file.write("----- Protein pair information \n")
        for tmp_protein_pair in self.protein_pairs:
            for tmp_info in tmp_protein_pair.create_plain_text_memory_mirror():
                memory_dump_file.write(f"{tmp_info} \n")
        memory_dump_file.close()
