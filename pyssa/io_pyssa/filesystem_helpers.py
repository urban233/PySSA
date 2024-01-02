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
"""Module for filesystem helper functions."""
import logging
import os
import pathlib
import shutil

from pyssa.logging_pyssa import log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def create_directory(a_directory_path: "pathlib.Path") -> None:
    """Creates a directory if it doesn't already exist."""
    if not os.path.exists(str(a_directory_path)):
        os.mkdir(str(a_directory_path))


def delete_directory(a_directory_path: pathlib.Path) -> None:
    """Deletes a directory and all its subdirectories."""
    if os.path.exists(a_directory_path):
        shutil.rmtree(a_directory_path)
