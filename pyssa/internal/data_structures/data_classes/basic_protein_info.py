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
"""Module contains the BasicProteinInfo dataclass."""
from dataclasses import dataclass


@dataclass
class BasicProteinInfo:
    """Holds basic information about a protein."""

    # <editor-fold desc="Class attributes">
    name: str
    """Name of the protein."""
    
    id: str  # noqa: A003
    """Id of the protein."""
    
    project_name: str
    """Name of the project that this protein belongs to."""
    
    # </editor-fold>

    def __eq__(self, other) -> bool:  # noqa: ANN001 #TODO: needs to be checked
        """Checks if two basic protein information objects are equal.
        
        Returns:
            A boolean indicating equality.
        """
        if not isinstance(other, BasicProteinInfo):
            return False
        return (self.name, self.id, self.project_name) == (other.name, other.id, other.project_name)

    def __hash__(self) -> int:
        """Combines the hash values of the attributes to create a unique hash.
        
        Returns:
            An hash value.
        """
        return hash((self.name, self.id, self.project_name))
