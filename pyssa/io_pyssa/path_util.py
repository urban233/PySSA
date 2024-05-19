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
"""Module containing a specific FilePath class."""
import logging
import os.path
import pathlib
from PyQt5 import QtCore
from pyssa.logging_pyssa import log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class FilePath:
    """This class is a hybrid class of pathlib.Path and QtCore.QFileInfo."""

    # <editor-fold desc="Class attributes">
    _filepath: pathlib.Path
    """
    A complete absolute path to a specific file
    """
    _dirname: pathlib.Path
    """
    A complete absolute path to a specific directory
    """
    _basename: str
    """
    The name of the file with the file extension
    """
    _filename: str
    """
    The name of the file WITHOUT file extension
    """
    _extension: str
    """
    The file extension of the file
    """

    # </editor-fold>

    def __init__(self, filepath: pathlib.Path) -> None:
        """Constructor.

        Args:
            filepath: An existing filepath.

        Raises:
            FileNotFoundError: If given filepath does not exist.
        """
        if os.path.exists(str(filepath)):
            # complete file path was given
            self._filepath = filepath
            self._dirname, self._basename = self.split_filepath()
            self._filename, self._extension = self.split_file_extension_from_name()
        else:
            raise FileNotFoundError("You need to pass an existing filepath!")

    def split_filepath(self) -> tuple[pathlib.Path, str]:
        """Splits the filepath in path to the parent directory and the actual filepath.
        
        Returns:
            A tuple with the path and the filename.
        """
        string_version = str(self._filepath)
        file_info = QtCore.QFileInfo(string_version)
        return pathlib.Path(file_info.canonicalPath()), f"{file_info.baseName()}.{file_info.suffix()}"

    def split_file_extension_from_name(self) -> tuple[str, str]:
        """Splits the extension of the filepath from its name.
        
        Returns:
            A tuple with the filename and the extension.
        """
        file_info = QtCore.QFileInfo(str(self._filepath))
        return file_info.baseName(), f".{file_info.suffix()}"

    def get_filepath(self) -> pathlib.Path:
        """Gets the filepath.
        
        Returns:
            The `_filepath` of the instance.
        """
        return self._filepath

    def get_dirname(self) -> pathlib.Path:
        """Gets the directory path of the filepath.
        
        Returns:
            The `_dirname` of the instance.
        """
        return self._dirname

    def get_basename(self) -> str:
        """Gets the basename of the filepath.
        
        Returns:
            The `_basename` of the instance.
        """
        return self._basename

    def get_filename(self) -> str:
        """Gets the filename of the filepath.
        
        Returns:
            The `_filename` of the instance.
        """
        return self._filename

    def get_file_extension(self) -> str:
        """Gets the file extension of the filepath as a string.
        
        Returns:
            The `_extension` of the instance.
        """
        return self._extension
