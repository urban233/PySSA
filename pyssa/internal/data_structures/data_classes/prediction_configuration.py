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
"""Module contains the PredictionConfiguration dataclass."""
from dataclasses import dataclass

__docformat__ = "google"


@dataclass
class PredictionConfiguration:
    """Holds information about the prediction configuration."""

    # <editor-fold desc="Class attributes">
    amber_force_field: bool
    """A boolean indicating if the amber force field should be used."""
    
    templates: str
    """A string indicating which templates should be used."""
    
    # </editor-fold>

    def get_tuple_notation(self) -> tuple[bool, str]:
        """Retrieves the tuple notation of the current instance.

        Returns:
            A tuple containing a boolean value indicating the status of the amber force field
            and a string representing the templates.
        """
        tuple_notation: tuple[bool, str] = (self.amber_force_field, self.templates)
        return tuple_notation
