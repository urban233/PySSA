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
"""Module for the job class"""

import os
from pathlib import Path
from xml.dom import minidom
from pyssa.gui.data_structures import project


class Job:
    """This class is for the jobs used in the plugin

    Var:
        job_name:
            name of the job
        job_path:
            path of the job
        project_list:
            list of Project objects
    """
    _job_name: str = ""
    _job_path: Path = ""
    _project_list: [project_utils.Project] = []

    def __init__(self, job_name, workspace_path) -> None:
        """Constructor

        Args:
            job_name:
                name of the job
            workspace_path:
                path of the workspace
        """
        job_name_with_underscores = job_name.replace(" ", "_")
        self._job_name = job_name_with_underscores
        self._job_path = Path(f"{workspace_path}/{self._job_name}")
        os.mkdir(self._job_path)

    def get_job_path(self) -> Path:
        """This function gets the value of the job_path variable

        Returns (pathlib.Path):
            job_path
        """
        return self._job_path

    def get_job_name(self) -> str:
        """This function gets the value of the job_name variable

        Returns (str):
            job_name
        """
        return self._job_name

    def add_project_to_job(self, project: project_utils.Project) -> None:
        """This function appends a Project object to the project_list

        Args:
            project:
                single project as Project object
        """
        self._project_list.append(project)

    def create_xml_file(self, attribute="value") -> None:
        """This function creates the job xml file

        Args:
            attribute (optional):
                attribute in the xml file
        """
        root = minidom.Document()
        root_node = root.createElement("root")
        root.appendChild(root_node)

        job_name_path_node = root.createElement("job_name")
        job_name_path_node.setAttribute(attribute, self._job_name)
        root_node.appendChild(job_name_path_node)

        for index in range(len(self._project_list)):
            xml_tag = root.createElement(f"project_{index}")
            xml_tag.setAttribute(attribute, self._project_list[index].get_project_path())
            root_node.appendChild(xml_tag)
            index += 1

        # save xml file to filesystem
        with open(f"{self.get_job_path()}/{self.get_job_name()}.xml", "w", encoding="utf-8") as file:
            file.write(root.toprettyxml())
