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
"""Module that implements a safeguard in form of static methods."""
import os
import pathlib

from pymol import cmd


class Safeguard:
    """This class is used to collect all safeguard functions in one place."""

    def __init__(self) -> None:
        """Constructor."""
        pass

    @staticmethod
    def check_filepath(a_filepath: pathlib.Path) -> bool:
        """This function checks if the filepath exists or not.

        Args:
            a_filepath (pathlib.Path): the filepath which should be checked

        Returns:
            True: if path exists
            False: if path does NOT exist
        """
        if not os.path.exists(a_filepath):
            return False
        return True

    @staticmethod
    def check_if_number_is_positive(value) -> bool:  # noqa: ANN001
        """This function checks if a number is positive.

        Returns:
            True: if number is positive (0 is also positive)
            False: if number is negative
        """
        if value >= 0:
            return True
        return False

    @staticmethod
    def check_if_protein_is_in_pymol(molecule_object: str) -> bool:
        """Checks if a protein is loaded into pymol."""
        if cmd.get_model(molecule_object) is not None:
            return True
        return False

    @staticmethod
    def check_if_list_is_empty(value: list) -> bool:
        """Checks if a list is empty."""
        if len(value) != 0:
            return True
        return False

    @staticmethod
    def check_if_value_is_in_table_v_header(value, table) -> bool:  # noqa: ANN001
        """Checks if a value is in a vertical header of a QTableWidget."""
        for i in range(table.rowCount()):
            header = table.verticalHeaderItem(i).text()
            new_value = value
            if header == new_value:
                return True
        return False

    @staticmethod
    def check_if_value_is_not_none(value) -> bool:  # noqa: ANN001
        """This function checks if a value is None or not.

        Args:
            value:
                any kind of variable
        Returns:
            True: if NOT None
            False: if None
        """
        if value is None:
            return False
        return True
