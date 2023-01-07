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
import logging
import os
from pyssa.gui.data_structures import settings


class Safeguard:
    """This class is used to collect all safeguard functions in one place.

    """
    def __init__(self):
        """Constructor

        """
        pass

    @staticmethod
    def check_filepath(filepath) -> bool:
        """This function checks if the filepath exists or not.

        Returns:
            True: if path exists
            False: if path does NOT exist
        """
        if not os.path.exists(filepath):
            logging.critical(f"The path {filepath} does NOT exist!")
            return False
        elif os.path.exists(filepath):
            # logging.info(f"The path {filepath} does exist!")
            return True

    @staticmethod
    def check_if_number_is_positive(value) -> bool:
        """This function checks if a number is positive

        Returns:
            True: if number is positive (0 is also positive)
            False: if number is negative
        """
        if value >= 0:
            return True
        else:
            return False
