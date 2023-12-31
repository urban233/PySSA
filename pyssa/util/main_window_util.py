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
from pyssa.util import constants, gui_utils
from pyssa.util import globals

if TYPE_CHECKING:
    from pyssa.internal.data_structures import settings


def check_operating_system() -> None:
    """Checks operating system and sets plugin path based upon operating system."""
    if sys.platform.startswith("win32"):
        # Windows path
        globals.g_os = "win32"
        globals.g_plugin_path = pathlib.Path(
            f"C:\\ProgramData\\pyssa\\mambaforge_pyssa\\pyssa-mamba-env\\Lib\\site-packages\\pymol\\pymol_path\\data\\startup\\{constants.PLUGIN_NAME}",
        )
    elif sys.platform.startswith("linux"):
        # Linux path
        globals.g_os = "linux"
        globals.g_plugin_path = f"/home/{os.getlogin()}/.local/pyssa/pyssa-mamba-env/lib/python3.10/site-packages/pmg_tk/startup/{constants.PLUGIN_NAME}"

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
    except ValueError:
        constants.PYSSA_LOGGER.warning("The settings file is damaged or outdated.")
        gui_utils.error_dialog_settings(
            "The settings file is damaged or outdated. You have to restore the settings to use PySSA!",
            "",
            the_app_settings,
        )
        tmp_settings: "settings.Settings" = the_app_settings

    if globals.g_os == "win32":
        constants.PYSSA_LOGGER.info("Checking if WSL2 is installed ...")
        if dialog_settings_global.is_wsl2_installed():
            tmp_settings.wsl_install = 1
            constants.PYSSA_LOGGER.info("WSL2 is installed.")
        else:
            tmp_settings.wsl_install = 0
            constants.PYSSA_LOGGER.warning("WSL2 is NOT installed.")
    else:
        tmp_settings.wsl_install = 1

    constants.PYSSA_LOGGER.info("Checking if Local Colabfold is installed ...")
    if dialog_settings_global.is_local_colabfold_installed():
        tmp_settings.local_colabfold = 1
        constants.PYSSA_LOGGER.info("Local Colabfold is installed.")
    else:
        tmp_settings.local_colabfold = 0
        constants.PYSSA_LOGGER.warning("Local Colabfold is NOT installed.")
    return tmp_settings
