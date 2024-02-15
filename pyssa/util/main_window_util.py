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
"""Module for the main window utility functions."""
import sys
import pathlib
import os
import requests
from typing import TYPE_CHECKING

from PyQt5 import QtWidgets
from pyssa.gui.ui.dialogs import dialog_settings_global
from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.internal.portal import pymol_io, graphic_operations
from pyssa.internal.data_structures import protein
from pyssa.io_pyssa import path_util
from pyssa.util import constants, gui_utils
from pyssa.util import globals

if TYPE_CHECKING:
    from pyssa.internal.data_structures import settings
    from pyssa.internal.data_structures import project


def check_operating_system() -> None:
    """Checks operating system and sets plugin path based upon operating system."""
    if sys.platform.startswith("win32"):
        # Windows path
        tmp_pymol_root_path: str = (
            "C:\\ProgramData\\pyssa\\mambaforge_pyssa\\pyssa-mamba-env\\Lib\\site-packages\\pymol"  # noqa: E501
        )
        globals.g_os = "win32"
        globals.g_plugin_path = pathlib.Path(
            f"{tmp_pymol_root_path}\\pymol_path\\data\\startup\\{constants.PLUGIN_NAME}",
        )
    elif sys.platform.startswith("linux"):
        # Linux path
        tmp_linux_pymol_plugin_root_path: str = f"/home/{os.getlogin()}/.local/pyssa/pyssa-mamba-env/lib/python3.10/site-packages/pmg_tk/startup"  # noqa: E501
        globals.g_os = "linux"
        globals.g_plugin_path = f"{tmp_linux_pymol_plugin_root_path}/{constants.PLUGIN_NAME}"

    constants.PYSSA_LOGGER.info(f"Started on platform: {globals.g_os}")


def check_version_number() -> None:
    """Check version of pyssa against the remote."""
    constants.PYSSA_LOGGER.info("Checking version of localhost with remote ...")
    url = "https://w-hs.sciebo.de/s/0kR11n8PkLDo1gB/download"
    try:
        response = requests.get(url)

        if response.status_code == 503:
            constants.PYSSA_LOGGER.warning(
                "The connection to the update servers failed. If could not determine if a new version is avaliable.",
            )
            basic_boxes.ok(
                "No connection",
                "The connection to the update servers failed. "
                "If could not determine if a new version is avaliable. You can start the PySSA now.",
                QtWidgets.QMessageBox.Information,
            )
        else:
            print(f"Latest version: {response.text}")
            if response.text != constants.VERSION_NUMBER[1:]:
                constants.PYSSA_LOGGER.info("here is a new version of PySSA available.")
                basic_boxes.ok(
                    "New version",
                    "There is a new version for your PySSA!\n"
                    f"To install the latest version {response.text}, open the PySSA Installer and click on update.",
                    QtWidgets.QMessageBox.Information,
                )
    except requests.exceptions.RequestException as e:
        constants.PYSSA_LOGGER.error(f"Downloading the version file of the remote failed! Error: {e}")
    constants.PYSSA_LOGGER.info("Checking version of localhost with remote finished.")


def setup_app_settings(the_app_settings: "settings.Settings") -> "settings.Settings":
    """Sets up application settings."""
    try:
        tmp_settings: "settings.Settings" = the_app_settings.deserialize_settings()
    except Exception as e:
        constants.PYSSA_LOGGER.warning("The settings file is damaged or outdated.")
        gui_utils.error_dialog_settings(
            "The settings file is damaged or outdated. You have to restore the settings to use PySSA!",
            "",
            the_app_settings,
        )
        tmp_settings: "settings.Settings" = the_app_settings
    return tmp_settings


def add_protein_to_project(the_protein_information: tuple, a_project: "project.Project") -> "project.Project":
    """Adds a protein based on the filepath/ PDB id to a project.

    Args:
        the_protein_information (tuple[str, int]): a pair of either filepath or id and the length of the first one.
        a_project: a project where the protein should be added.

    Note:
        If a PDB id is used than the first element of the tuple contains the id and the second the length = 4.
        If a filepath is used than the first element of the tuple contains the id and the second the length of
        the entire filepath.
    """
    # TODO: checks are needed
    if the_protein_information[1] == 4:
        # PDB ID is used
        tmp_protein = pymol_io.get_protein_from_pdb(the_protein_information[0])
    else:
        # Filepath is used
        pdb_filepath: "path_util.FilePath" = path_util.FilePath(pathlib.Path(the_protein_information[0]))
        graphic_operations.setup_default_session_graphic_settings()
        tmp_protein_name: str = pdb_filepath.get_filename().replace(" ", "_")
        tmp_protein = protein.Protein(
            molecule_object=tmp_protein_name,
            pdb_filepath=pdb_filepath,
        )
    a_project.add_existing_protein(tmp_protein)
    return a_project
