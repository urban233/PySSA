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
"""Module contains ResidueColor dataclass."""
from dataclasses import dataclass

__docformat__ = "google"


@dataclass
class ResidueColorConfig:
  """Holds information about a residue's color configuration."""

  # <editor-fold desc="Class attributes">
  carbon_color: str
  """The color of the carbon atoms."""

  nitrogen_color: str
  """The color of the nitrogen atoms."""

  oxygen_color: str
  """The color of the oxygen atoms."""

  # </editor-fold>

  def atoms_are_colored_by_elements(self) -> bool:
    """Checks if atoms are colored correctly based on their elements.

    Returns:
        True if atoms are colored by their elements, False otherwise.
    """
    if (
        self.carbon_color == "grey70"
        and self.nitrogen_color == "N-blue"
        and self.oxygen_color == "O-red"
    ):
      return True
    return False
