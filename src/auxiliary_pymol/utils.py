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
"""Contains basic utility functions."""
import base64
import os
import pathlib
__docformat__ = "google"


def create_base64_string_from_file(a_filepath: str) -> str:
    """Creates a base64 string from a binary file.

    Args:
        a_filepath (str): a filepath to a binary file.

    Returns:
        A base64 encoded string or an empty string if a_filepath is None or could not be found.
    """
    # <editor-fold desc="Checks">
    if a_filepath is None or not os.path.exists(a_filepath):
        return ""
    
    # </editor-fold>
        
    with open(a_filepath, "rb") as binary_file:
        binary_data = binary_file.read()
        binary_file.close()
    # Encode binary data to base64 and decode to utf-8 string
    return base64.b64encode(binary_data).decode("utf-8")


def write_binary_file_from_base64_string(filepath: pathlib.Path, base64_data: str) -> bool:
    """Writes base64 data to a binary file.

    Args:
        filepath (pathlib.Path): a filepath to a binary file.
        base64_data (str): a base64 string to write.
    
    Returns:
        A boolean value indicating success or failure.
    """
    # <editor-fold desc="Checks">
    if filepath is None:
        return False
    if base64_data is None or base64_data == "":
        return False
    
    # </editor-fold>
    
    try:
        binary_data_export = base64.b64decode(base64_data)
        directory = os.path.dirname(filepath)
        if not os.path.exists(directory):
            os.makedirs(directory)
        with open(filepath, "wb") as f:
            f.write(binary_data_export)
    except Exception:
        return False
    else:
        return True
