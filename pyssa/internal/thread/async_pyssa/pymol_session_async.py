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
"""Module for all asynchronous functions that are related to the pymol session."""
import logging

from pyssa.controller import database_manager
from pyssa.logging_pyssa import log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def load_protein_pymol_session(a_protein, the_pymol_session_manager, needs_to_be_reinitialized_flag: bool = False) -> tuple:
    if needs_to_be_reinitialized_flag:
        logger.info("The current session is not empty. Reinitialize session now ...")
        the_pymol_session_manager.reinitialize_session()
        logger.info("Reinitializing session finished.")
    try:
        logger.info(f"Loading session of {a_protein.get_molecule_object()}")
        the_pymol_session_manager.load_protein_session(a_protein)
        the_pymol_session_manager.current_scene_name = "base"
        the_pymol_session_manager.load_current_scene()
        the_pymol_session_manager.get_all_scenes_in_current_session()
    except RuntimeError:
        logger.error("Loading the session failed due to a RuntimeError!")
        return the_pymol_session_manager, False
    except Exception as e:
        logger.error(f"Loading the session failed because this error was raised: {e}")
        return the_pymol_session_manager, False
    else:
        logger.info("Loading the session finished without errors.")
        return the_pymol_session_manager, True


def load_protein_pair_pymol_session(a_protein_pair, the_pymol_session_manager, needs_to_be_reinitialized_flag: bool = False) -> tuple:
    if needs_to_be_reinitialized_flag:
        logger.info("The current session is not empty. Reinitialize session now ...")
        the_pymol_session_manager.reinitialize_session()
        logger.info("Reinitializing session finished.")
    try:
        logger.info(f"Loading session of {a_protein_pair.name}")
        the_pymol_session_manager.load_protein_pair_session(a_protein_pair)
    except Exception as e:
        logger.error(f"Loading the session failed because this error was raised: {e}")
        return 0, False
    else:
        logger.info("Loading the session finished without errors.")
        return 0, True


def save_protein_pymol_session_to_database(
        the_interface_manager: "interface_manager.InterfaceManager",
        placeholder: int) -> tuple:
    try:
        with database_manager.DatabaseManager(
                str(the_interface_manager.get_current_project().get_database_filepath())) as db_manager:
            tmp_protein = the_interface_manager.get_current_active_protein_object()
            tmp_protein.pymol_session = the_interface_manager.pymol_session_manager.save_current_session_as_base64()
            db_manager.update_pymol_session_of_protein(
                tmp_protein.get_id(),
                tmp_protein.pymol_session
            )
    except Exception as e:
        logger.error(f"Unexpected error occurred. Exception: {e}")
        return 0, False
    else:
        return 0, True


def save_protein_pair_pymol_session_to_database(
        the_interface_manager: "interface_manager.InterfaceManager",
        placeholder: int) -> tuple:
    try:
        with database_manager.DatabaseManager(
                str(the_interface_manager.get_current_project().get_database_filepath())) as db_manager:
            tmp_protein_pair = the_interface_manager.get_current_active_protein_pair_object()
            tmp_protein_pair.pymol_session = the_interface_manager.pymol_session_manager.save_current_session_as_base64()
            db_manager.update_pymol_session_of_protein_pair(
                tmp_protein_pair.get_id(),
                tmp_protein_pair.pymol_session
            )
    except Exception as e:
        logger.error(f"Unexpected error occurred. Exception: {e}")
        return 0, False
    else:
        return 0, True

