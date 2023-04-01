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
import os.path
import pathlib
from PyQt5 import QtCore
from pyssa.logging_pyssa import log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


# TODO: create docstrings and checks
class FilePath:
    """This class is a hybrid class of pathlib.Path and QtCore.QFileInfo."""

    # <editor-fold desc="Class attributes">
    """
    a complete absolute path to a specific file
    """
    _filepath: pathlib.Path
    """
    a complete absolute path to a specific directory
    """
    _dirname: pathlib.Path
    """
    the name of the file with the file extension
    """
    _basename: str
    """
    the name of the file WITHOUT file extension
    """
    _filename: str
    """
    the file extension of the file
    """
    _extension: str

    # </editor-fold>

    def __init__(self, filepath):
        if os.path.exists(filepath):
            # complete file path was given
            self._filepath = filepath
            self._dirname, self._basename = self.split_filepath()
            self._filename, self._extension = self.split_file_extension_from_name()
        else:
            raise FileNotFoundError("You need to pass an existing filepath!")

    def split_filepath(self) -> tuple[pathlib.Path, str]:
        string_version = str(self._filepath)
        file_info = QtCore.QFileInfo(string_version)
        src = pathlib.Path(file_info.canonicalPath())
        return pathlib.Path(file_info.canonicalPath()), f"{file_info.baseName()}.{file_info.suffix()}"

    def merge_filepath_to_filename(self) -> pathlib.Path:
        return pathlib.Path(str(self._filepath)+str(self._filename))

    def split_file_extension_from_name(self) -> tuple[str, str]:
        file_info = QtCore.QFileInfo(str(self._filepath))
        return file_info.baseName(), f".{file_info.suffix()}"

    def get_filepath(self):
        return self._filepath

    def get_dirname(self):
        return self._dirname

    def get_basename(self):
        return self._basename

    def get_filename(self):
        return self._filename

    def get_file_extension(self):
        return self._extension

    def check_if_path_exists(self) -> bool:
        if os.path.exists(self._filepath):
            return True
        else:
            return False
