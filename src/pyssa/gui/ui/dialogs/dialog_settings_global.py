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
"""Module for the global settings dialog."""
import logging
import subprocess

# setup logger
logging.basicConfig(level=logging.DEBUG)


def is_wsl2_installed() -> bool:
  """Checks if the WSL2 is installed.

  Returns:
      True if the WSL2 is installed, False otherwise.
  """
  output = subprocess.run("wsl --list --verbose")
  if output.returncode == 0:
    return True
  return False


def is_local_colabfold_installed() -> bool:
  """Checks if the local colabfold is installed.

  Returns:
      True if the local colabfold is installed, False otherwise.
  """
  powershell_results = subprocess.run(["wsl", "-d", "almaColabfold9", "ls"])
  if powershell_results.returncode == 0:
    subprocess.run(["wsl", "--shutdown"])
    return True
  subprocess.run(["wsl", "--shutdown"])
  return False
