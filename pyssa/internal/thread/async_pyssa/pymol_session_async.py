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
"""Module for all asynchronous functions that are related to the pymol session."""
import logging
import time
from typing import Optional

from pyssa.controller import database_manager, pymol_session_manager
from pyssa.internal.data_structures import protein
from pyssa.internal.data_structures.data_classes import residue_color_config
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)

__docformat__ = "google"


def load_protein_pymol_session(
        a_protein: "protein.Protein",
        the_pymol_session_manager: "pymol_session_manager.PymolSessionManager",
        needs_to_be_reinitialized_flag: bool = False
) -> tuple[Optional["pymol_session_manager.PymolSessionManager"], bool]:
    """
    Loads a protein pymol session.

    Args:
        a_protein (protein.Protein): The protein object to load into the PyMOL session.
        the_pymol_session_manager (pymol_session_manager.PymolSessionManager): The PyMOL session manager object.
        needs_to_be_reinitialized_flag (bool): Flag indicating whether the PyMOL session needs to be reinitialized. Default is False.

    Returns:
        A tuple containing the updated PyMOL session manager object and a boolean indicating whether the loading of the session was successful.

    Raises:
        RuntimeError: If loading the session fails with a RuntimeError.
    """
    # <editor-fold desc="Checks">
    if a_protein is None:
        logger.error("a_protein is None.")
        return None, False
    if the_pymol_session_manager is None:
        logger.error("the_pymol_session_manager is None.")
        return None, False
    if needs_to_be_reinitialized_flag is None:
        logger.error("needs_to_be_reinitialized_flag is None.")
        return None, False

    # </editor-fold>

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
    try:
        if needs_to_be_reinitialized_flag:
            logger.info("The current session is not empty. Reinitialize session now ...")
            the_pymol_session_manager.reinitialize_session()
            logger.info("Reinitializing session finished.")
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


def load_scene(the_pymol_session_manager: "pymol_session_manager.PymolSessionManager", a_scene_name: str) -> tuple[bool]:
    """
    Loads a scene in the current PyMOL session.

    Args:
        the_pymol_session_manager (pymol_session_manager.PymolSessionManager): The PyMOL session manager instance.
        a_scene_name (str): The name of the scene to be loaded.

    Returns:
        A tuple containing a single boolean value indicating the success of the operation.
    """
    # <editor-fold desc="Checks">
    if the_pymol_session_manager is None:
        logger.error("the_pymol_session_manager is None.")
        return (False,)
    if a_scene_name is None or a_scene_name == "":
        logger.error("a_scene_name is either None or an empty string.")
        return (False,)

    # </editor-fold>

    time.sleep(0.5)  # This should increase the stability because the scene command has more time in between command calls
    try:
        the_pymol_session_manager.load_scene(a_scene_name)
    except Exception as e:
        logger.error(f"Unexpected error occurred. Exception: {e}")
        return (False,)
    else:
        return (True,)


def color_pymol_selection(
        a_color: str,
        a_pymol_selection: str,
        the_pymol_session_manager: "pymol_session_manager.PymolSessionManager"
) -> tuple[bool, str]:
    """
    Colors a given pymol selection with a given color.

    Args:
        a_color (str): The color to be applied to the protein chain.
        a_pymol_selection (str): The selection of the protein chain(s) to be colored.
        the_pymol_session_manager: The PymolSessionManager object.

    Returns:
        A tuple containing a single boolean value indicating the success of the operation and the chain color.
    """
    # <editor-fold desc="Checks">
    if a_color is None or a_color == "":
        logger.error("a_color is either None or an empty string!")
        return False, ""
    if a_color not in constants.PYMOL_COLORS_WITH_INDICES.values():
        logger.error("a_color is not part of the PYMOL_COLORS_WITH_INDICES dict!")
        return False, ""
    if a_pymol_selection is None or a_pymol_selection == "":
        logger.error("a_pymol_selection is either None or an empty string.")
        return False, ""
    if the_pymol_session_manager is None:
        logger.error("the_pymol_session_manager is None.")
        return False, ""

    # </editor-fold>

    try:
        the_pymol_session_manager.color_protein(
            a_color, a_pymol_selection
        )
    except Exception as e:
        logger.error(f"Unexpected error occurred. Exception: {e}")
        return False, ""
    else:
        return True, a_color


def reset_color_pymol_selection_atoms_by_element(
        a_protein_name: str,
        a_chain_letter: str,
        the_current_active_chain_color_of_protein: str,
        a_pymol_selection: str,
        the_pymol_session_manager: "pymol_session_manager.PymolSessionManager"
) -> tuple:
    """
    Reset the color of protein chain atoms by element for a given protein.

    Args:
        a_protein_name (str): The name of the protein.
        a_chain_letter (str): The chain letter of the protein.
        the_current_active_chain_color_of_protein (str): The current active chain color of the protein.
        a_pymol_selection (str): The selection in PyMOL where the protein chain atoms are located.
        the_pymol_session_manager (pymol_session_manager.PymolSessionManager): The PyMOL session manager.

    Returns:
        tuple: A tuple containing two values: success status (bool) and the updated chain color (str).
    """
    # <editor-fold desc="Checks">
    if a_protein_name is None or a_protein_name == "":
        logger.error("a_protein_name is either None or an empty string.")
        return False, ""
    if a_chain_letter is None or a_chain_letter == "":
        logger.error("a_chain_letter is either None or an empty string.")
        return False, ""
    if the_current_active_chain_color_of_protein is None or the_current_active_chain_color_of_protein == "":
        logger.error("the_current_active_chain_color_of_protein is either None or an empty string.")
        return False, ""
    if a_pymol_selection is None or a_pymol_selection == "":
        logger.error("a_pymol_selection is either None or an empty string.")
        return False, ""
    if the_pymol_session_manager is None:
        logger.error("the_pymol_session_manager is None.")
        return False, ""

    # </editor-fold>

    try:
        tmp_residue_config: "residue_color_config.ResidueColorConfig" = the_pymol_session_manager.get_residue_color_config_of_a_given_selection(
            a_protein_name, a_chain_letter
        )
        if tmp_residue_config.atoms_are_colored_by_elements():
            tmp_chain_color = the_current_active_chain_color_of_protein
        else:
            tmp_chain_color = tmp_residue_config.carbon_color
        # This if-statement is needed for the case that a pymol session is loaded,
        # that contains an atom coloring
        if tmp_chain_color == "By Element":
            tmp_chain_color = "green"

        the_pymol_session_manager.color_protein(tmp_chain_color, f"{a_pymol_selection}")
    except Exception as e:
        logger.error(f"Unexpected error occurred. Exception: {e}")
        return False, ""
    else:
        return True, tmp_chain_color


def color_pymol_selection_atoms_by_element(
        a_pymol_selection: str,
        the_pymol_session_manager: "pymol_session_manager.PyMolSessionManager"
) -> tuple:
    """
    Colors protein chain atoms based on the element for a given protein.

    Args:
        a_pymol_selection (str): The PyMOL selection of atoms for a protein chain.
        the_pymol_session_manager: The PyMOL session manager.

    Returns:
        tuple: A tuple with a single boolean value indicating the success of the operation.
    """
    # <editor-fold desc="Checks">
    if a_pymol_selection is None or a_pymol_selection == "":
        logger.error("a_pymol_selection is either None or an empty string.")
        return (False,)
    if the_pymol_session_manager is None:
        logger.error("the_pymol_session_manager is None.")
        return (False,)

    # </editor-fold>

    try:
        the_pymol_session_manager.color_protein(
            "atomic", f"{a_pymol_selection} and not elem C"
        )
        the_pymol_session_manager.color_protein(
            "grey70", f"{a_pymol_selection} and elem C"
        )
    except Exception as e:
        logger.error(f"Unexpected error occurred. Exception: {e}")
        return (False,)
    else:
        return (True,)


def get_residue_color_config_of_a_given_protein_chain(
        a_protein_name: str,
        a_chain_letter: str,
        the_pymol_session_manager: "pymol_session_manager.PymolSessionManager"
) -> tuple[bool, Optional["residue_color_config.ResidueColorConfig"]]:
    """
    Gets the colors of C-, N-, and O-atoms for the first residue of the given selection.

    Args:
        a_protein_name (str): The name of the protein.
        a_chain_letter (str): The chain letter of the protein.
        the_pymol_session_manager (pymol_session_manager.PymolSessionManager): The PyMOL session manager instance.

    Returns:
        tuple: A tuple with a single boolean value indicating the success of the operation and
        a ResidueColorConfig instance or None if an error occurred.

    """
    # <editor-fold desc="Checks">
    if a_protein_name is None or a_protein_name == "":
        logger.error("a_protein_name is either None or an empty string.")
        return False, None
    if a_chain_letter is None or a_chain_letter == "":
        logger.error("a_chain_letter is either None or an empty string.")
        return False, None
    if a_chain_letter not in constants.chain_dict.values():
        logger.error("a_chain_letter is not part of the chain_dict.")
        return False, None
    if the_pymol_session_manager is None:
        logger.error("the_pymol_session_manager is None.")
        return False, None

    # </editor-fold>

    try:
        tmp_residue_config: "residue_color_config.ResidueColorConfig" = the_pymol_session_manager.get_residue_color_config_of_a_given_selection(
            a_protein_name, a_chain_letter
        )
    except Exception as e:
        logger.error(f"Unexpected error occurred. Exception: {e}")
        return False, None
    else:
        return True, tmp_residue_config


def get_representation_config_of_a_given_protein_chain(
        a_selection_string: str,
        a_chain_letter: str,
        the_pymol_session_manager: "pymol_session_manager.PymolSessionManager"
) -> tuple[bool, Optional[dict]]:
    """
    Returns the representation state of a specific chain in PyMOL.

    Args:
        a_selection_string: A string representing the selection of atoms in the protein chain.
        a_chain_letter: A string representing the letter identifying the protein chain.
        the_pymol_session_manager: An instance of the class PymolSessionManager from the pymol_session_manager module.

    Returns:
        A tuple containing a boolean value indicating the success of the method and an optional dictionary representing the representation state of the protein chain.
    """
    # <editor-fold desc="Checks">
    if a_selection_string is None or a_selection_string == "":
        logger.error("a_selection_string is either None or an empty string.")
        return False, None
    if a_chain_letter is None or a_chain_letter == "":
        logger.error("a_chain_letter is either None or an empty string.")
        return False, None
    if a_chain_letter not in constants.chain_dict.values():
        logger.error("a_chain_letter is not part of the chain_dict.")
        return False, None
    if the_pymol_session_manager is None:
        logger.error("the_pymol_session_manager is None.")
        return False, None

    # </editor-fold>

    try:
        tmp_repr_state = the_pymol_session_manager.get_chain_repr_state(a_selection_string,a_chain_letter)
    except Exception as e:
        logger.error(f"Unexpected error occurred. Exception: {e}")
        return False, None
    else:
        return True, tmp_repr_state
