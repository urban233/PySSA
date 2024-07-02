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
"""Module contains the AnalysisRunInfo dataclass."""
from dataclasses import dataclass


@dataclass
class AnalysisRunInfo:
  """Contains all information about an analysis run."""

  # <editor-fold desc="Class attributes">
  _protein_name_1: str
  """The first protein name."""

  protein_chains_1: list
  """The chains of the first protein."""

  _protein_name_2: str
  """The second protein name."""

  protein_chains_2: list
  """The chains of the second protein."""

  analysis_name: str
  """The name of the analysis run."""

  # </editor-fold>

  def get_protein_name_1(self) -> str:
    """Gets the name of the first protein.

    Returns:
        The name of the first protein.
    """
    return self._protein_name_1.replace(".pdb", "")

  def get_protein_name_2(self) -> str:
    """Gets the name of the second protein.

    Returns:
        The name of the first protein.
    """
    return self._protein_name_2.replace(".pdb", "")

  def are_protein_names_identical(self) -> bool:
    """Checks if the two protein names are identical.

    Returns:
        A boolean indicating if the two protein names are identical.
    """
    if self._protein_name_1 == self._protein_name_2:
      return True
    return False
