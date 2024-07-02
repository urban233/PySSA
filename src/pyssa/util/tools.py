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
import socket
from urllib.error import URLError
from urllib.request import urlopen
from src.pyssa.internal.data_structures import settings
from src.pyssa.util import constants, exception


def check_internet_connectivity() -> bool:
  """Checks the connection to the internet.

  Returns:
      A boolean indicating whether the internet is available or not.
  """
  timeout: float = 6
  try:
    urlopen("https://www.google.com", timeout=timeout)
  except URLError:
    return False
  except socket.timeout:
    return False
  except Exception as e:
    constants.PYSSA_LOGGER.error(
        f"Could not connect to internet, due to an unknown error: {e}"
    )
    return False
  else:
    return True


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
