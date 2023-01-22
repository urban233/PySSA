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
from dataclasses import dataclass


@dataclass
class PyMOLSelection:
    obj_name: str
    segi: str
    chain: list[str]
    resi: list[int]
    atom_name: str

    def get_single_selection(self) -> str:
        """This function creates a single pymol selection with only one chain and one resi

        Returns:
            one pymol selection string
        """
        if len(self.chain) > 1 or len(self.resi) > 1:
            raise ValueError
        return f"/{self.obj_name}/{self.segi}/{self.chain[0]}/{self.resi[0]}/{self.atom_name}"

    def get_complex_selection(self) -> str:
        """This function creates one pymol selection string but with multiple chains and resi's

        Returns:
            one pymol selection string
        """
        # TODO: implement function
        pass

    def get_all_single_selections(self) -> list[str]:
        """This function creates multiple pymol selection strings with multiple chains and resi's

        Returns:
            a list of pymol selection strings
        """
        # TODO: implement function
        pass

    def get_custom_selection(self, chain, resi) -> str:
        """This function creates a single pymol selection strings based on the args chain and resi.

        Notes:
            only ONE chain and ONE resi is supported with this function!

        Returns:
            one pymol selection string
        """
        return f"/{self.obj_name}/{self.segi}/{chain}/{resi}/{self.atom_name}"
