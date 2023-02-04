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
import fnmatch
import os
import json
import pathlib
import shutil
import logging
from PyQt5 import QtWidgets
from pyssa.logging_pyssa import log_handlers
from pyssa.internal.data_structures.data_classes import protein_info
from pyssa.internal.data_structures import protein
from pyssa.internal.data_structures import protein_pair
from pyssa.internal.data_structures import project
from pyssa.internal.data_structures import settings
from pyssa.internal.data_structures import sequence
from pyssa.io_pyssa import safeguard
from pyssa.util import constants
from pyssa.util import tools

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class ObjectSerializer:
    object_dict: dict

    def __init__(self, object_to_serialize, filepath, filename):
        """Constructor

        Args:
            object_to_serialize:
                objects which gets serialized
            filepath:
                filepath of the json file, without the file name
            filename:
                name of the json file WITHOUT .json extension!

        """
        self.object = object_to_serialize
        self.filepath = filepath
        self.filename = filename

    def create_standard_object_dict(self):
        self.object_dict = self.object.__dict__
    
    def set_custom_object_dict(self, custom_dict: dict):
        self.object_dict = custom_dict
    
    def update_object_dict(self, dict_with_updates: dict):
        self.object_dict.update(dict_with_updates)
    
    def serialize_object(self):
        tmp_object_file = open(pathlib.Path(f"{self.filepath}/{self.filename}.json"), "w", encoding="utf-8")
        json.dump(self.object_dict, tmp_object_file, indent=4)
        tmp_object_file.close()


class ObjectDeserializer:
    
    def __init__(self, filepath, filename):
        tmp_object_file = open(f"{filepath}/{filename}.json", "r", encoding="utf-8")
        self.object_dict = json.load(tmp_object_file)
    
    def deserialize_protein(self) -> protein.Protein:
        if self.object_dict.get("export_data_dir") == "None":
            update = {"export_data_dir": None}
            self.object_dict.update(update)
        tmp_protein = protein.Protein(self.object_dict.get("molecule_object"),
                                      self.object_dict.get("filepath"),
                                      self.object_dict.get("export_data_dir")
                                      )
        tmp_protein.filename = self.object_dict.get("filename")
        tmp_protein.set_sequence(self.object_dict.get("sequence"))
        tmp_protein.set_selection(self.object_dict.get("selection"))
        tmp_protein.set_chains(self.object_dict.get("chains"))
        return tmp_protein
    
    def deserialize_protein_pair(self):
        tmp_protein_pair = protein_pair.ProteinPair(self.object_dict.get("prot_1_molecule_object"),
                                                    self.object_dict.get("prot_2_molecule_object"),
                                                    pathlib.Path(self.object_dict.get("results_dir")),
                                                    )
        tmp_protein_pair.cutoff = self.object_dict.get("cutoff")
        tmp_protein_pair.name = self.object_dict.get("name")

        if self.object_dict.get("prot_1_export_data_dir") == "None":
            export_data_dir = None
        else:
            export_data_dir = self.object_dict.get("prot_1_export_data_dir")
        tmp_protein_1 = protein.Protein(self.object_dict.get("prot_1_filename"),
                                        self.object_dict.get("prot_1_import_data_dir"),
                                        export_data_dir=export_data_dir,
                                        )
        tmp_protein_1.molecule_object = self.object_dict.get("prot_1_molecule_object")
        tmp_protein_1.set_sequence(self.object_dict.get("prot_1_sequence"))
        tmp_protein_1.set_selection(self.object_dict.get("prot_1_selection"))
        tmp_protein_1.set_chains(self.object_dict.get("prot_1_chains"))
        tmp_protein_pair.ref_obj = tmp_protein_1

        if self.object_dict.get("prot_2_export_data_dir") == "None":
            export_data_dir = None
        else:
            export_data_dir = self.object_dict.get("prot_2_export_data_dir")
        tmp_protein_2 = protein.Protein(self.object_dict.get("prot_2_filename"),
                                        self.object_dict.get("prot_2_import_data_dir"),
                                        export_data_dir=export_data_dir,
                                        )
        tmp_protein_2.molecule_object = self.object_dict.get("prot_2_molecule_object")
        tmp_protein_2.set_sequence(self.object_dict.get("prot_2_sequence"))
        tmp_protein_2.set_selection(self.object_dict.get("prot_2_selection"))
        tmp_protein_2.set_chains(self.object_dict.get("prot_2_chains"))
        tmp_protein_pair.model_obj = tmp_protein_2

        return tmp_protein_pair
    
    def deserialize_project(self):
        tmp_project: project.Project = project.Project(self.object_dict.get("_project_name"), self.object_dict.get("_workspace"))
        tmp_project.set_folder_paths(self.object_dict.get("_folder_paths"))
        tmp_project.project_path = self.object_dict.get("project_path")
        # adding proteins to projects object
        if len(os.listdir(tmp_project.get_objects_proteins_path())) > 0:
            for tmp_protein_file in os.listdir(tmp_project.get_objects_proteins_path()):
                tmp_protein = protein.Protein.deserialize_protein(
                    pathlib.Path(f"{tmp_project.get_objects_proteins_path()}/{tmp_protein_file}"))
                tmp_project.add_existing_protein(tmp_protein)
        # adding protein pairs to projects object
        if len(os.listdir(tmp_project.get_objects_protein_pairs_path())) > 0:
            for tmp_protein_pair_file in os.listdir(tmp_project.get_objects_protein_pairs_path()):
                tmp_protein_pair = protein_pair.ProteinPair.deserialize_protein_pair(
                    pathlib.Path(f"{tmp_project.get_objects_protein_pairs_path()}/{tmp_protein_pair_file}"))
                tmp_project.add_protein_pair(tmp_protein_pair)
        return tmp_project
    
    def deserialize_settings(self):
        tmp_settings: settings.Settings = settings.Settings(self.object_dict.get("dir_settings"), self.object_dict.get("filename"))
        if safeguard.Safeguard.check_filepath(self.object_dict.get("workspace_path")):
            tmp_settings.workspace_path = self.object_dict.get("workspace_path")
        else:
            raise ValueError
        if safeguard.Safeguard.check_filepath(self.object_dict.get("prediction_path")):
            tmp_settings.prediction_path = self.object_dict.get("prediction_path")
        else:
            raise ValueError
        if safeguard.Safeguard.check_if_number_is_positive(int(self.object_dict.get("cycles"))):
            tmp_settings.cycles = self.object_dict.get("cycles")
        else:
            raise ValueError
        if safeguard.Safeguard.check_if_number_is_positive(float(self.object_dict.get("cutoff"))):
            tmp_settings.cutoff = self.object_dict.get("cutoff")
        else:
            raise ValueError
        if int(self.object_dict.get("app_launch")) == 0 or int(self.object_dict.get("app_launch")) == 1:
            tmp_settings.app_launch = self.object_dict.get("app_launch")
        else:
            raise ValueError
        if not safeguard.Safeguard.check_filepath(self.object_dict.get("dir_settings")):
            raise ValueError
        if not self.object_dict.get("filename") == "settings.json":
            raise ValueError
        if int(self.object_dict.get("wsl_install")) == 0 or int(self.object_dict.get("wsl_install")) == 1:
            tmp_settings.wsl_install = int(self.object_dict.get("wsl_install"))
        else:
            raise ValueError
        if int(self.object_dict.get("local_colabfold")) == 0 or int(self.object_dict.get("local_colabfold")) == 1:
            tmp_settings.local_colabfold = int(self.object_dict.get("local_colabfold"))
        else:
            raise ValueError
        tmp_settings.wsl_username = self.object_dict.get("wsl_username")
        return tmp_settings

    def deserialize_sequence(self):
        """This function deserializes the sequence object from a .json file.

        Returns:
            a ProteinSequence object
        """
        return sequence.ProteinSequence(self.object_dict.get("name"), self.object_dict.get("sequence"))


class ProjectScanner:

    def __init__(self, project_obj: project.Project):
        self.project = project_obj

    def scan_project_for_valid_proteins(self, list_view_project_proteins=None) -> list[str]:
        """This function scans the project pdb path and optionally fills a list view

        Args:
            list_view_project_proteins:

        Returns:

        """
        project_proteins: list[str] = os.listdir(self.project.get_pdb_path())
        pattern = "*.pdb"
        # iterates over possible project directories
        if list_view_project_proteins is not None:
            for tmp_protein in project_proteins:
                if fnmatch.fnmatch(tmp_protein, pattern):
                    list_view_project_proteins.addItem(tmp_protein)
        return project_proteins


class WorkspaceScanner:

    def __init__(self, workspace_path):
        self.workspace_path = workspace_path

    def scan_workspace_for_valid_projects(self, list_new_projects):
        workspace_projects: list[str] = os.listdir(self.workspace_path)
        valid_directories = []
        # iterates over possible project directories
        for directory in workspace_projects:
            try:
                directory_content = os.listdir(f"{self.workspace_path}/{directory}")
                # iterates over the content in a single project directory
                for content in directory_content:
                    if content == "project.json":
                        valid_directories.append(directory)
            except NotADirectoryError:
                print(f"This: {directory} is not a directory.")

        valid_directories.sort()
        for project in valid_directories:
            list_new_projects.addItem(project)
        return valid_directories

    def scan_workspace_for_non_duplicate_proteins(self, valid_projects: list, current_project_name: str,
                                                  list_widget: QtWidgets.QListWidget) -> tuple[dict, list]:
        """This function scans the workspace directory for protein structures and eliminates all duplicates

        Args:
            valid_projects (list):
                a list of all projects within the workspace
            current_project_name (str):
                name of the currently loaded project
            list_widget (Qt.QtWidgets.QListWidget)
                list widget which is needed to temporarily store the results from the function "scan_project_for_valid_proteins"

        Returns:
            dict which contains all proteins without duplicates
        """
        """Var: workspace_proteins is a list which contains all proteins from all projects in the workspace"""
        workspace_proteins = []
        protein_names = []
        protein_tuples_notation = []
        for valid_project in valid_projects:
            # if valid_project != current_project_name: # I don't know why this if-statement should be important
            """Var: project_proteins is a list which contains all proteins from a single project"""
            project_proteins = tools.scan_project_for_valid_proteins(pathlib.Path(f"{self.workspace_path}/{valid_project}"),
                                                                     list_widget)
            list_widget.clear()
            for protein in project_proteins:
                tmp_protein = protein_info.ProteinInfo(protein,
                                                       pathlib.Path(f"{self.workspace_path}/{valid_project}/pdb/{protein}"))
                workspace_proteins.append(tmp_protein)
                if tmp_protein.name not in protein_names:
                    protein_names.append(tmp_protein.name)
        # this for-loop is necessary for the creation of the protein dictionary
        for protein in workspace_proteins:
            protein_tuples_notation.append(protein.get_tuple_notation())
        protein_dict = dict(protein_tuples_notation)
        return protein_dict, protein_names


class FilesystemCleaner:

    def __int__(self):
        pass

    @staticmethod
    def clean_prediction_scratch_folder():
        shutil.rmtree(constants.PREDICTION_FASTA_DIR)
        os.mkdir(constants.PREDICTION_FASTA_DIR)
        shutil.rmtree(constants.PREDICTION_PDB_DIR)
        os.mkdir(constants.PREDICTION_PDB_DIR)
