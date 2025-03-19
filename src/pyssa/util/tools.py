#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
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
"""Module for functions which can be used across the entire project."""
import pathlib
import json
from typing import Union
from typing import Optional

import requests

from src.pyssa.internal.data_structures import settings
from src.pyssa.util import constants, exception


def check_internet_connectivity() -> bool:
  """Checks the connection to the internet.

  Returns:
      A boolean indicating whether the internet is available or not.
  """
  timeout: float = 6
  try:
    response = requests.get("https://www.google.com", timeout=timeout)
    return response.status_code == 200
  except requests.RequestException as e:
    print(f"Could not connect to internet: {e}")
    return False


def download_file(url: str, filepath: str) -> None:
  """Downloads a file from the given URL and saves it to the specified filepath."""
  try:
    response = requests.get(url, stream=True)
    response.raise_for_status()
    with open(filepath, 'wb') as file:
      for chunk in response.iter_content(chunk_size=8192):
        file.write(chunk)
  except requests.RequestException as e:
    print(f"Failed to download file: {e}")


def restore_default_settings(settings_obj: "settings.Settings") -> None:
  """Creates a settings.json file which is filled with the pre-defined values.

  Args:
      settings_obj: The settings object where the default values should be stored.
  """
  if settings_obj is None:
    raise exception.IllegalArgumentError("settings_obj is None.")
  settings_obj.restore_settings(
      constants.SETTINGS_DIR, constants.SETTINGS_FILENAME
  )


def get_latest_version(json_file_path: Union[str, pathlib.Path]) -> Optional[str]:
  """Gets the latest version of the remote version history.

  Args:
    json_file_path: Filepath to the version history JSON file.
  """
  with open(json_file_path, 'r') as file:
    data = json.load(file)

  version_history = data['versionHistory']

  if not version_history:
    return None  # Return None if there are no versions

  # Convert version strings to tuples of integers for proper comparison
  latest_entry = max(version_history,
                     key=lambda x: tuple(map(int, x['version'].split('.'))))

  return latest_entry['version']


def get_latest_release(json_file_path: Union[str, pathlib.Path]) -> Optional[dict]:
  """Gets the latest version of the remote version history.

  Args:
    json_file_path: Filepath to the version history JSON file.
  """
  with open(json_file_path, 'r') as file:
    data = json.load(file)

  version_history = data['versionHistory']

  if not version_history:
    return None

  # Find the entry with the highest semantic version
  latest_entry = max(version_history,
                     key=lambda x: tuple(map(int, x['version'].split('.'))))

  return {
    'version': latest_entry['version'],
    'releaseDate': latest_entry['releaseDate'],
    'releaseUrl': latest_entry['releaseUrl']
  }
