    #
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2022
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
"""Module for safeguards in pymol"""
from pymol import cmd


class PymolSafeguard:
    """This class is used to collect all safeguard functions in one place for pymol.

    """
    def __init__(self):
        """Constructor

        """
        pass

    @staticmethod
    def check_if_protein_in_session() -> bool:
        """This function checks if the protein exists in current pymol session or not.

        Returns:
            True: if protein exists
            False: if protein does NOT exist
        """
        if len(cmd.get_chains()) == 0:
            return False
        else:
            return True
