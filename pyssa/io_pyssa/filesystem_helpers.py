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
import fnmatch
import os
import pathlib
from PyQt5 import QtCore
from pyssa.io_pyssa import safeguard
from pyssa.logging_pyssa import log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def get_file_path_and_name(full_filepath) -> tuple[pathlib.Path, str]:
    """This function gets the file name and path from a full filepath

    Args:
        full_filepath (str):
            the full filepath with file extension

    Returns:
        path:
            absolute path where the file is stored
        file_name:
            name of the file with file extension
    """
    if not safeguard.Safeguard.check_filepath(full_filepath):
        raise NotADirectoryError
    file_info = QtCore.QFileInfo(str(full_filepath))
    file_name = file_info.baseName()
    path = pathlib.Path(file_info.canonicalPath())
    return path, file_name


def get_filepath_from_full_filepath(full_filepath) -> pathlib.Path:
    """This function gets the file path from a full filepath

    Args:
        full_filepath (str):
            the full filepath with file extension

    Returns:
        path:
            absolute path where the file is stored
    """
    if not safeguard.Safeguard.check_filepath(full_filepath):
        raise NotADirectoryError
    return pathlib.Path(QtCore.QFileInfo(str(full_filepath)).canonicalPath())


def get_filename_from_full_filepath(full_filepath) -> str:
    """This function gets the file name from a full filepath

    Args:
        full_filepath (str):
            the full filepath with file extension

    Returns:
        file_name:
            name of the file with file extension
    """
    if not safeguard.Safeguard.check_filepath(full_filepath):
        raise NotADirectoryError
    return QtCore.QFileInfo(str(full_filepath)).baseName()


def filter_directory_for_filetype(path, filetype: str) -> list[str]:
    """This function filters a directory for a given filetype.

    Args:
        path:
            full path of the directory to filter
        filetype:
            file extension which should get filtered, WITHOUT .

    Returns:
        a list with all files which match the given extension
    """
    if not safeguard.Safeguard.check_filepath(path):
        raise NotADirectoryError
    valid_files: list[str] = []
    pattern = f"*.{filetype}"
    # iterates over possible project directories
    for tmp_file in os.listdir(path):
        if fnmatch.fnmatch(tmp_file, pattern):
            valid_files.append(tmp_file)
    return valid_files


def create_generic_dictionary_from_directory(path) -> dict:
    """This function creates a dictionary which is based on the content of a directory. This should increase
    search speed.

    Args:
        path:
            full path of the directory

    Raises:
        NotADirectoryError: if directory does not exist

    Returns:
        dictionary which consists of a generic key (token_x) and a file or directory as value
    """
    if not safeguard.Safeguard.check_filepath(path):
        logger.error("The directory does not exists!")
        raise NotADirectoryError("The directory does not exist!")
    i = 0
    tmp_generic_dict: dict = {}
    for item in os.listdir(path):
        tmp_generic_dict.update({f"token_{i}": item})
        i += 1
    return tmp_generic_dict


def create_directory(a_directory_path: 'pathlib.Path') -> None:
    """Creates a directory if it doesn't already exist."""
    if not os.path.exists(str(a_directory_path)):
        os.mkdir(str(a_directory_path))
