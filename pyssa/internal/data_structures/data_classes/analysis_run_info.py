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
"""Module contains the AnalysisRunInfo dataclass."""
from dataclasses import dataclass


@dataclass
class AnalysisRunInfo:
    """Class which contains all information about an analysis run."""

    _protein_name_1: str
    protein_chains_1: list
    _protein_name_2: str
    protein_chains_2: list
    analysis_name: str

    def get_protein_name_1(self) -> str:
        """Gets the name of the first protein."""
        return self._protein_name_1.replace(".pdb", "")

    def get_protein_name_2(self) -> str:
        """Gets the name of the second protein."""
        return self._protein_name_2.replace(".pdb", "")

    def are_protein_names_identical(self) -> bool:
        """Checks if the two protein names are identical."""
        if self._protein_name_1 == self._protein_name_2:
            return True
        return False
