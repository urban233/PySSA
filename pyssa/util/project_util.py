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
"""Module contains helper function for the project class."""
import logging
import os
import pathlib
from pyssa.io_pyssa import path_util
from pyssa.logging_pyssa import log_handlers
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pyssa.internal.data_structures import project


logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def get_all_filepaths_from_project(
    app_project: "project.Project",
    subfolder: str,
    extension: str,
) -> list["path_util.FilePath"]:
    """Gets all filepaths from a given project.

    Args:
        app_project: the app project to get the filepaths from.
        subfolder: the subfolder to search in.
        extension: the extension to search for.
    """
    tmp_dirnames = []
    folder = os.path.join(app_project.folder_paths["project"], subfolder)
    if len(os.listdir(folder)) > 0:
        for tmp_folder in os.listdir(folder):
            tmp_dirnames.append(pathlib.Path(f"{folder}/{tmp_folder}"))
    else:
        logger.warning(f"The subfolder {subfolder} is empty. An empty list will be returned.")
        return []
    filepaths = []
    for tmp_dirname in tmp_dirnames:
        for basename in os.listdir(f"{tmp_dirname}/{extension}"):
            filepaths.append(path_util.FilePath(f"{tmp_dirname}/{extension}/{basename}"))
    return filepaths


def get_all_protein_json_filepaths_from_project(app_project: "project.Project") -> list["path_util.FilePath"]:
    """Gets all protein json filepaths from a given project."""
    return get_all_filepaths_from_project(app_project, "proteins", ".objects")


def get_all_pdb_filepaths_from_project(app_project: "project.Project") -> list["path_util.FilePath"]:
    """Gets all protein pdb filepaths from a given project."""
    return get_all_filepaths_from_project(app_project, "proteins", "pdb")


def get_all_protein_pair_json_filepaths_from_project(app_project: "project.Project") -> list["path_util.FilePath"]:
    """Gets all protein pair json filepaths from a given project."""
    return get_all_filepaths_from_project(app_project, "protein_pairs", ".objects")


def get_all_distance_analysis_json_filepaths_from_project(app_project: "project.Project") -> list["path_util.FilePath"]:
    """Gets all distance analysis json filepaths from a given project."""
    return get_all_filepaths_from_project(app_project, str(pathlib.Path("analysis/distance_analysis")), ".objects")
