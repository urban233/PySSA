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

import os
from pathlib import Path
from xml.dom import minidom
from pymol import Qt
from utils import gui_utils
from utils import tools


class Project:
    """This class is for the projects used in the plugin
    Var:
        job_name:
            name of the job (project-only as default)
        project_name:
            name of the project
        pdb_file:
            path of the loaded pdb file
        pdb_id:
            name of the PDB ID
        pdb_models:
            path of the models which are used for the analysis
        ref_chains:
            chains which are used for the reference protein
        model_chains:
            chains which are used for the model protein
        results_path:
            path where all results are stored
    """
    _job_name: str = "project_only"
    _project_name: str = ""
    _pdb_file: str = ""
    _pdb_id: str = ""
    _pdb_model: str = ""
    _ref_chains: str = ""
    _model_chains: str = ""
    _folder_paths: list[Path] = []
    _session_file_name: str = "session_file_model_s.pse"

    def __init__(self, project_name, workspace_path) -> None:
        """Constructor

        Args:
            project_name:
                name of the project
            workspace_path:
                path of the workspace
        """
        project_name_with_underscores = project_name.replace(" ", "_")
        self._project_name = project_name_with_underscores
        # project_folder_path = Path(f"{workspace_path}/{project_name_with_underscores}")

        self._folder_paths = [
            Path(f"{workspace_path}/{project_name_with_underscores}"),
            Path(f"{workspace_path}/{project_name_with_underscores}/pdb"),
            Path(f"{workspace_path}/{project_name_with_underscores}/results"),
            Path(f"{workspace_path}/{project_name_with_underscores}/results/alignment_files"),
            Path(f"{workspace_path}/{project_name_with_underscores}/results/distance_csv"),
            Path(f"{workspace_path}/{project_name_with_underscores}/results/images"),
            Path(f"{workspace_path}/{project_name_with_underscores}/results/images/interesting_regions"),
            Path(f"{workspace_path}/{project_name_with_underscores}/results/plots"),
            Path(f"{workspace_path}/{project_name_with_underscores}/results/plots/distance_histogram"),
            Path(f"{workspace_path}/{project_name_with_underscores}/results/plots/distance_plot"),
            Path(f"{workspace_path}/{project_name_with_underscores}/results/sessions/"),
        ]

    def create_project_tree(self) -> bool:
        """This function creates the directory structure for a new project.

        """
        # check if the project folder already exists
        if os.path.exists(self.get_project_path()):
            detailed_message = f"The project is located under: {self.get_project_path()}."
            flag = gui_utils.warning_message_project_exists(self._project_name, detailed_message,
                                                            self.get_project_path())
            if flag is False:
                return False

        for folder in self._folder_paths:
            os.mkdir(folder)

    def get_session_file(self) -> str:
        """This function gets the value of the session_file variable

        Returns (str):
            session_file as complete path with file name
        """
        return f"{self._folder_paths[10]}/{self._session_file_name}"

    def get_job_name(self) -> str:
        """This function gets the value of the job_name variable

        Returns (str):
            job_name
        """
        return self._job_name

    def set_job_name(self, value: str) -> None:
        """This function gets the value of the job_name variable

        """
        self._job_name = value

    def get_project_name(self) -> str:
        """This function gets the value of the project_name variable

        Returns (str):
            project_name
        """
        return self._project_name

    def set_project_name(self, value: str) -> None:
        """This function gets the value of the project_name variable

        """
        self._project_name = value

    def get_pdb_file(self) -> str:
        """This function gets the value of the pdb_file variable

        Returns (str):
            project_name
        """
        return self._pdb_file

    def set_pdb_file(self, value: str) -> None:
        """This function gets the value of the pdb_file variable

        """
        self._pdb_file = value

    def get_pdb_id(self) -> str:
        """This function gets the value of the pdb_id variable

        Returns (str):
            project_name
        """
        return self._pdb_id

    def set_pdb_id(self, value: str) -> None:
        """This function gets the value of the pdb_id variable

        """
        self._pdb_id = value

    def get_pdb_models(self) -> str:
        """This function gets the value of the pdb_models variable

        Returns (str):
            project_name
        """
        return self._pdb_model

    def set_pdb_model(self, value: str) -> None:
        """This function gets the value of the pdb_models variable

        """
        self._pdb_model = value

    def get_ref_chains(self) -> str:
        """This function gets the value of the ref_chains variable

        Returns (str):
            project_name
        """
        return self._ref_chains

    def set_ref_chains(self, value: str) -> None:
        """This function gets the value of the ref_chains variable

        """
        self._ref_chains = value

    def get_model_chains(self) -> str:
        """This function gets the value of the model_chains variable

        Returns (str):
            project_name
        """
        return self._model_chains

    def set_model_chains(self, value: str) -> None:
        """This function gets the value of the model_chains variable

        """
        self._model_chains = value

    def get_results_path(self) -> str:
        """This function gets the value of the results_path variable

        Returns (str):
            project_name
        """
        return str(self._folder_paths[2])

    def set_results_path(self, value: str) -> None:
        """This function gets the value of the results_path variable

        """
        self._results_path = value

    def get_all_properties(self) -> list:
        """This functions gets the value of all properties in this class

        Returns:
            properties_list which contains the values of all properties
        """
        properties_list = [
            self._job_name,
            self._project_name,
            self._pdb_file,
            self._pdb_id,
            self._pdb_model,
            self._ref_chains,
            self._model_chains,
            self._results_path,
        ]
        return properties_list

    def set_all_properties(self, value_list: list):
        """This function sets all properties in this class

        Args:
            value_list (list):
                list which contains the values of all properties

        Notes:
            if a property should be empty use "" as value

        """
        if len(value_list) != 8:
            raise ValueError
        self._job_name = value_list[0]
        self._project_name = value_list[1]
        self._pdb_file = value_list[2]
        self._pdb_id = value_list[3]
        self._pdb_model = value_list[4]
        self._ref_chains = value_list[5]
        self._model_chains = value_list[6]
        self._results_path = value_list[7]

    def get_model_filename(self) -> str:
        """This function returns the model filename

        Returns:
            model filename
        """
        return Qt.QtCore.QFileInfo(self._pdb_model).baseName()

    def get_project_path(self) -> str:
        """This function returns the project path of the project

        Returns:
            the project path
        """
        return str(self._folder_paths[0])

    def get_pdb_path(self) -> str:
        """This function returns the pdb path of the project

        Returns:
            the pdb path
        """
        return str(self._folder_paths[1])

    def create_xml_file(self, attribute="value") -> None:
        """This function create the settings xml with the format:

        Args:
            attribute (optional str):
                defines the attribute to access the elements

        <?xml version="1.0" ?>
        <root>
            <project_name DEFAULT_ATTRIBUTE=project_name_value/>
            <pdb_file DEFAULT_ATTRIBUTE=pdb_file_path>
            <pdb_id DEFAULT_ATTRIBUTE=pdb_id_value/>
            <pdb_models DEFAULT_ATTRIBUTE=pdb_models_value/>
            <ref_chains DEFAULT_ATTRIBUTE=ref_chains_value/>
            <model_chains DEFAULT_ATTRIBUTE=model_chains_value/>
            <results_path DEFAULT_ATTRIBUTE=results_path_value/>
        </root>
        """
        root = minidom.Document()
        root_node = root.createElement("root")
        root.appendChild(root_node)

        job_name_path_node = root.createElement("job_name")
        job_name_path_node.setAttribute(attribute, self._job_name)
        root_node.appendChild(job_name_path_node)

        project_name_path_node = root.createElement("project_name")
        project_name_path_node.setAttribute(attribute, self._project_name)
        root_node.appendChild(project_name_path_node)

        pdb_file_path_node = root.createElement("pdb_file")
        pdb_file_path_node.setAttribute(attribute, self._pdb_file)
        root_node.appendChild(pdb_file_path_node)

        pdb_id_path_node = root.createElement("pdb_id")
        pdb_id_path_node.setAttribute(attribute, self._pdb_id)
        root_node.appendChild(pdb_id_path_node)

        pdb_models_path_node = root.createElement("pdb_model")
        pdb_models_path_node.setAttribute(attribute, self.get_pdb_models())
        root_node.appendChild(pdb_models_path_node)

        ref_chains_path_node = root.createElement("ref_chains")
        ref_chains_path_node.setAttribute(attribute, self._ref_chains)
        root_node.appendChild(ref_chains_path_node)

        model_chains_path_node = root.createElement("model_chains")
        model_chains_path_node.setAttribute(attribute, self._model_chains)
        root_node.appendChild(model_chains_path_node)

        results_path_path_node = root.createElement("results_path")
        results_path_path_node.setAttribute(attribute, self.get_results_path())
        root_node.appendChild(results_path_path_node)

        with open(f"{self.get_project_path()}/{self.get_project_name()}.xml", "w", encoding="utf-8") as file:
            file.write(root.toprettyxml())

    def create_project(self, workspace_path, model, MODEL_OBJ_NAME, txt_box_project_name,
                       txt_box_load_reference, txt_box_chain_ref, txt_box_chain_model, job="") -> tuple[str, str]:
        """This function creates a project folder and project xml file

        Args:
            workspace_path:
                path of the current workspace
            model:
                full filepath of the model
            MODEL_OBJ_NAME:
                only the name of the model
            txt_box_project_name:
                text box of the project name
            txt_box_load_reference:
                text box where the path of the reference is stored
            txt_box_chain_ref:
                text box of the chains which should be used for the reference
            txt_box_chain_model:
                text box of the chains which should be used for the model
            job (optional):
                if the project is part of a job request, set the job name

        Returns:
            project_pdb_path:
                which is the path of the pdb directory within the project folder
            project_results_path:
                which is the path of the results directory within the project
                folder
        """
        if txt_box_project_name != "":
            project_name = txt_box_project_name.text()
        else:
            project_name = f"project_of_model_{MODEL_OBJ_NAME}"
        project_name_with_underscores = project_name.replace(" ", "_")

        if job == "":
            project_folder_path = Path(f"{workspace_path}/{project_name_with_underscores}")
            # check if the project folder already exists
            if os.path.exists(project_folder_path):
                detailed_message = f"The project is located under: {project_folder_path}."
                gui_utils.warning_message_project_exists(project_name, detailed_message,
                                                         project_folder_path)
                return "", ""
            os.mkdir(project_folder_path)
        else:
            # check if the project folder already exists
            project_folder_path = Path(f"{workspace_path}/{job}/{project_name_with_underscores}")
            if os.path.exists(project_folder_path):
                detailed_message = f"The project is located under: {project_folder_path}."
                gui_utils.warning_message_project_exists(project_name, detailed_message,
                                                         project_folder_path)
                raise IsADirectoryError
            os.mkdir(project_folder_path)

        # setting values for the project xml file
        if job != "":
            self._job_name = job
        self._project_name = project_name_with_underscores
        if len(txt_box_load_reference.text()) == 4:
            self._pdb_id = txt_box_load_reference.text()
        self._pdb_file = txt_box_load_reference.text()
        if txt_box_chain_ref.text() != "":
            # TODO: * implement input check with QValidator and Regex
            self._ref_chains = txt_box_chain_ref.text()
        if txt_box_chain_model.text() != "":
            self._model_chains = txt_box_chain_model.text()
        self._results_path = f"{project_folder_path}/results"
        # sets the filepath of the model in the project xml file
        self._pdb_model = model
        self.create_xml_file(f"{str(project_folder_path)}/{project_name_with_underscores}.xml")

        # creates a pdb folder
        tools.create_directory(Path(f"{project_folder_path}"), "pdb")
        project_pdb_path = f"{project_folder_path}/pdb"
        tools.create_directory(Path(f"{project_folder_path}"), "results")
        project_results_path = f"{project_folder_path}/results"
        return project_pdb_path, project_results_path
