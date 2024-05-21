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
"""Module for the main window utility functions."""
import logging
import sys
import pathlib
import os
import requests
from typing import TYPE_CHECKING

from pyssa.gui.ui.custom_dialogs import custom_message_box
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants, tools, exception
from pyssa.util import globals

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"

if TYPE_CHECKING:
  from pyssa.internal.data_structures import settings


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
    tmp_linux_pymol_plugin_root_path: str = (
        f"/home/{os.getlogin()}/.local/pyssa/pyssa-mamba-env/lib/python3.10/site-packages/pmg_tk/startup"  # noqa: E501
    )
    globals.g_os = "linux"
    globals.g_plugin_path = (
        f"{tmp_linux_pymol_plugin_root_path}/{constants.PLUGIN_NAME}"
    )

  logger.info(f"Started on platform: {globals.g_os}")


def check_version_number() -> None:
  """Check version of pyssa against the remote."""
  logger.info("Checking version of localhost with remote ...")
  url = "https://w-hs.sciebo.de/s/0kR11n8PkLDo1gB/download"
  try:
    response = requests.get(url)

    if response.status_code == 503:
      logger.warning(
          "The connection to the refresh servers failed. If could not determine if a new version is available.",
      )
      tmp_dialog = custom_message_box.CustomMessageBoxOk(
          "The connection to the refresh servers failed. "
          "It could not determine if a new version is available. You can start the PySSA now.",
          "No Internet Connection",
          custom_message_box.CustomMessageBoxIcons.INFORMATION.value,
      )
      tmp_dialog.exec_()
    else:
      print(f"Latest version: {response.text}")
      if response.text != constants.VERSION_NUMBER[1:]:
        logger.info("here is a new version of PySSA available.")
        tmp_dialog = custom_message_box.CustomMessageBoxOk(
            "There is a new version for your PySSA!\n"
            f"To install the latest version {response.text}, open the PySSA Installer and click on Update.",
            "New Version!",
            custom_message_box.CustomMessageBoxIcons.INFORMATION.value,
        )
        tmp_dialog.exec_()
  except requests.exceptions.RequestException as e:
    logger.error(
        f"Downloading the version file of the remote failed! Error: {e}"
    )
  logger.info("Checking version of localhost with remote finished.")


def setup_app_settings(
    the_app_settings: "settings.Settings",
) -> "settings.Settings":
  """Sets up application settings.

  Args:
      the_app_settings: The settings object for the application.

  Returns:
      The deserialized settings object.

  Raises:
      exception.IllegalArgumentError: If the_app_settings is None.
  """
  # <editor-fold desc="Checks">
  if the_app_settings is None:
    logger.error("the_app_settings is None.")
    raise exception.IllegalArgumentError("the_app_settings is None.")

  # </editor-fold>

  try:
    tmp_settings: "settings.Settings" = the_app_settings.deserialize_settings()
  except Exception as e:
    logger.warning("The settings file is damaged or outdated.")
    tmp_dialog = custom_message_box.CustomMessageBoxOk(
        "The settings file is damaged or outdated. You have to restore the settings to use PySSA!",
        "Restore Settings",
        custom_message_box.CustomMessageBoxIcons.WARNING.value,
    )
    tmp_dialog.exec_()
    tools.restore_default_settings(the_app_settings)
    tmp_settings: "settings.Settings" = the_app_settings.deserialize_settings()
  return tmp_settings
