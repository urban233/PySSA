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
"""Module for in- and output processes in pymol."""
import logging
import pathlib
import pymol
import os
from pymol import cmd

from pyssa.logging_pyssa import log_handlers
from pyssa.util import exception
from pyssa.internal.data_structures import protein
from pyssa.internal.portal import graphic_operations
from pyssa.io_pyssa import safeguard
from pyssa.io_pyssa import binary_data
from pyssa.io_pyssa import path_util
from pyssa.util import constants, tools

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def load_pymol_session(pymol_session_filepath: pathlib.Path) -> None:
    """This function loads a pymol session file into the current pymol session.

    Args:
        pymol_session_filepath: the filepath of the session file
    Raises:
        FileNotFoundError: If file not found.

    """
    if not safeguard.Safeguard.check_filepath(pymol_session_filepath):
        raise FileNotFoundError
    logger.debug("Starting to load pymol session ...")
    cmd.load(str(pymol_session_filepath))
    logger.debug("Finished loading pymol session.")


def save_pymol_session_as_base64_string(pymol_molecule_object: str) -> pathlib.Path:
    session_filepath = pathlib.Path(f"{constants.SCRATCH_DIR}/{pymol_molecule_object}_session.pse")
    cmd.save(session_filepath)
    return session_filepath


def convert_pymol_session_to_base64_string(pymol_molecule_object: str) -> str:
    """This function converts a pymol session file into a base64 string.

    Args:
        pymol_molecule_object (str): PyMOL molecule object to be converted.
    """
    session_filepath = save_pymol_session_as_base64_string(pymol_molecule_object)
    base64_string = binary_data.create_base64_string_from_file(path_util.FilePath(session_filepath))
    os.remove(session_filepath)
    return base64_string


def get_all_scenes_from_pymol_session() -> list[str]:
    return cmd.get_scene_list()
