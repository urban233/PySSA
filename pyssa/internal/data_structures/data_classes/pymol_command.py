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
"""Module contains the CurrentSession dataclass."""
from dataclasses import dataclass

from pyssa_pymol import pymol_enums

__docformat__ = "google"


@dataclass
class PyMOLCommand:
  """Holds information about a pymol command."""

  # <editor-fold desc="Class attributes">
  command: "pymol_enums.CommandEnum"
  """The pymol command."""

  arguments: tuple
  """The command arguments."""

  # </editor-fold>

  def get_command(self) -> dict:
    """Returns a dictionary with the command and its arguments.

    Returns:
        dict: A dictionary containing the command and its arguments.
    """
    return {
        "command": self.command.value,
        "arguments": self.arguments,
    }
