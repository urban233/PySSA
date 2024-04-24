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
"""Module for all asynchronous functions used in the main presenter."""
import copy
import logging
import os
import pathlib
import shutil
import subprocess
import platform
from typing import TYPE_CHECKING

import pymol
from Bio import SeqRecord, SeqIO
from pymol import cmd

from pyssa.controller import database_manager, interface_manager, pymol_session_manager, watcher
from pyssa.internal.data_structures import project, protein, structure_analysis, structure_prediction
from pyssa.internal.data_structures.data_classes import (
    prediction_protein_info,
    prediction_configuration,
    current_session, database_operation,
)
from pyssa.internal.portal import pymol_io, graphic_operations, protein_operations
from pyssa.internal.thread import database_thread
from pyssa.internal.thread.async_pyssa import custom_signals
from pyssa.io_pyssa import path_util, filesystem_io, bio_data
from pyssa.logging_pyssa import log_handlers
from pyssa.util import analysis_util, exception, exit_codes, constants, enums

if TYPE_CHECKING:
    from pyssa.internal.data_structures import settings
    from pyssa.internal.data_structures import protein_pair

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def clean_protein_new(
    a_protein_name: str,
    a_project: "project.Project",
) -> tuple:
    """Cleans a protein by creating a duplicate and removing all solvent and sugar molecules.

    Args:
        a_protein_name: the name of the protein to clean.
        a_project: the current project.

    Returns:
        a tuple with ("result", a_updates_project_object)
    """
    # TODO: needs checks
    tmp_protein = a_project.search_protein(
        a_protein_name,
    )
    clean_tmp_protein = tmp_protein.clean_protein(new_protein=True)
    constants.PYSSA_LOGGER.info("The protein %s has been cleaned.", clean_tmp_protein.get_molecule_object())
    a_project.add_existing_protein(clean_tmp_protein)
    #a_project.serialize_project(a_project.get_project_xml_path())
    return ("result", a_project)


def clean_protein_update(a_protein: "protein.Protein", the_database_filepath :str) -> tuple:
    """Cleans a protein by removing all solvent and sugar molecules in the current molecule object.

    Args:
        a_protein: the protein object to clean.
        a_project: the current project.

    Returns:
        a tuple with ("result", a_updates_project_object)
    """
    # TODO: needs checks
    # TODO: needs tests!
    a_protein.clean_protein()
    if len(a_protein.get_pdb_data()) == 0:
        logger.error("No PDB data found after cleaning process!")
        raise ValueError("No PDB data found after cleaning process!")

    with database_manager.DatabaseManager(the_database_filepath) as db_manager:
        db_manager.open_project_database()
        # tmp_pdb_atom_data = db_manager.get_pdb_atoms_of_protein(a_protein.get_id())
        # tmp_pdb_atom_dict_1 = [{key.value: value for key, value in zip(enums.PdbAtomEnum, t)} for t in tmp_pdb_atom_data]
        # a_protein.set_pdb_data(tmp_pdb_atom_dict_1)
        db_manager.update_protein_pdb_atom_data(a_protein.get_id(), a_protein.get_pdb_data())
        db_manager.close_project_database()
        a_protein.set_pdb_data([])
    constants.PYSSA_LOGGER.info("The protein %s has been cleaned.", a_protein.get_molecule_object())
    return ("result", a_protein)


def save_selected_protein_sequence_as_fasta_file(
    a_seq_record: "SeqRecord.SeqRecord",
    a_filepath: str,
    the_database_filepath: str,
) -> tuple:
    """Saves a given protein sequence to a fasta file."""
    # with database_manager.DatabaseManager(the_database_filepath) as db_manager:
    #     db_manager.open_project_database()
    #     tmp_pdb_atom_data = db_manager.get_pdb_atoms_of_protein(a_protein.get_id())
    #     tmp_pdb_atom_dict_1 = [{key.value: value for key, value in zip(enums.PdbAtomEnum, t)} for t in tmp_pdb_atom_data]
    #     a_protein.set_pdb_data(tmp_pdb_atom_dict_1)
    #     db_manager.close_project_database()
    try:
        print(a_seq_record.id)
        if a_seq_record.id == "<unknown id>":
            a_seq_record.id = a_seq_record.name
        with open(a_filepath, "w") as file_handler:
            SeqIO.write(a_seq_record, file_handler, "fasta")
    except Exception as e:
        logger.error(f"Saving sequence to fasta file ended with error: {e}")
        return (exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0], exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[1])
    else:
        return (exit_codes.EXIT_CODE_ZERO[0], exit_codes.EXIT_CODE_ZERO[1])


def save_selected_protein_structure_as_pdb_file(
    a_protein: "protein.Protein",
    a_filepath: str,
    the_database_filepath: str,
) -> tuple:
    """Saves a given protein structure to a pdb file."""
    with database_manager.DatabaseManager(the_database_filepath) as db_manager:
        db_manager.open_project_database()
        tmp_pdb_atom_data = db_manager.get_pdb_atoms_of_protein(a_protein.get_id())
        tmp_pdb_atom_dict_1 = [{key.value: value for key, value in zip(enums.PdbAtomEnum, t)} for t in tmp_pdb_atom_data]
        a_protein.set_pdb_data(tmp_pdb_atom_dict_1)
        db_manager.close_project_database()
    try:
        bio_data.build_pdb_file(a_protein.get_pdb_data(), a_filepath)
        # reset pdb data to reduce memory space
        a_protein.set_pdb_data([])
    except Exception as e:
        logger.error(f"Saving protein to pdb file ended with error: {e}")
        return (exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0], exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[1])
    else:
        return (exit_codes.EXIT_CODE_ZERO[0], exit_codes.EXIT_CODE_ZERO[1])


def rename_selected_protein_structure(
    a_protein: "protein.Protein",
    the_new_protein_name: str,
    the_database_filepath: str,
) -> tuple:
    """Renames a certain protein from a project.

    Args:
        a_protein: the protein object to rename.
        the_new_protein_name: the new name for the given protein.
        the_database_filepath: the filepath of the project database.

    Returns:
        a tuple with ("result", an_existing_protein_object)
    """
    tmp_old_protein_name = a_protein.get_molecule_object()

    # Update in memory
    a_protein.set_molecule_object(the_new_protein_name)

    # Update in pymol
    try:
        cmd.set_name(tmp_old_protein_name, a_protein.get_molecule_object())
    except pymol.CmdException:
        logger.error(f"Renaming the protein in PyMOL failed.")
        raise RuntimeError("Renaming the protein in PyMOL failed. Maybe it is not opened in a session?")
    # Update in database
    with database_manager.DatabaseManager(the_database_filepath) as db_manager:
        db_manager.open_project_database()
        db_manager.update_protein_name(
            the_new_protein_name, tmp_old_protein_name, a_protein.get_id()
        )
        db_manager.close_project_database()
    return ("result", a_protein)


def color_protein_pair_by_rmsd_value(a_protein_pair: "protein_pair.ProteinPair", placeholder: int) -> tuple:
    """Colors a given protein pair by their rmsd value.

    Args:
        a_project: a project object containing the protein pair.
        a_results_name: the name of the protein pair.

    Returns:
        a tuple with ("result", an_existing_protein_pair_object)
    """
    graphic_operations.color_protein_pair_by_rmsd(a_protein_pair)
    return ("result", a_protein_pair)


def check_chains_for_analysis(the_protein_1_name: str, the_protein_2_name: str, a_project: "project.Project") -> tuple:
    """Checks if the two proteins have only one chain.

    Args:
        the_protein_1_name: the name of the first protein from the analysis.
        the_protein_2_name: the name of the second protein from the analysis.
        a_project: the current project.

    Returns:
        a tuple with ("result", tmp_is_only_one_chain, tmp_analysis_run_name, tmp_protein_1, tmp_protein_2)
    """
    # TODO: checks needed
    # TODO: tests needed
    tmp_is_only_one_chain: bool = False
    tmp_analysis_run_name: str = ""
    tmp_protein_1 = a_project.search_protein(the_protein_1_name)
    tmp_protein_2 = a_project.search_protein(the_protein_2_name)
    if len(tmp_protein_1.get_protein_sequences()) == 1: # or len(tmp_protein_2.get_protein_sequences()) == 1:
        tmp_is_only_one_chain = True
        tmp_protein_1_first_protein_chain_letter = protein_operations.get_chain_letter_of_first_protein_sequence(
            tmp_protein_1.chains
        )
        tmp_protein_2_first_protein_chain_letter = protein_operations.get_chain_letter_of_first_protein_sequence(
            tmp_protein_2.chains
        )
        tmp_analysis_run_name = (
            f"{tmp_protein_1.get_molecule_object()};{tmp_protein_1_first_protein_chain_letter}"
            f"_vs_{tmp_protein_2.get_molecule_object()};{tmp_protein_2_first_protein_chain_letter}"
        )
    return ("result", tmp_is_only_one_chain, tmp_analysis_run_name, tmp_protein_1, tmp_protein_2)


def check_chains_for_subsequent_analysis(
    the_protein_1_name: str,
    the_protein_2_name: str,
    a_project: "project.Project",
    a_list_with_proteins_to_predict: list[str],
) -> tuple:
    """Checks if the two proteins have only one chain.

    Args:
        the_protein_1_name: the name of the first protein from the analysis.
        the_protein_2_name: the name of the second protein from the analysis.
        a_project: the current project.
        a_list_with_proteins_to_predict: a list of all proteins that should be predicted.
    """
    tmp_protein_1_only_one_chain: bool = False
    tmp_protein_2_only_one_chain: bool = False
    tmp_sub_analysis_name_1: str = ""
    tmp_sub_analysis_name_2: str = ""
    if the_protein_1_name not in a_list_with_proteins_to_predict:
        tmp_protein_1 = a_project.search_protein(the_protein_1_name)
        if len(tmp_protein_1.chains) == 1:
            tmp_protein_1_only_one_chain = True
            tmp_sub_analysis_name_1 = f"{tmp_protein_1.get_molecule_object()};{tmp_protein_1.chains[0].chain_letter}"
    else:
        tmp_protein_1_only_one_chain = True
        tmp_sub_analysis_name_1 = f"{the_protein_1_name};A"

    if the_protein_2_name not in a_list_with_proteins_to_predict:
        tmp_protein_2 = a_project.search_protein(the_protein_2_name)
        if len(tmp_protein_2.chains) == 1:
            tmp_protein_2_only_one_chain = True
            tmp_sub_analysis_name_2 = f"{tmp_protein_2.get_molecule_object()};{tmp_protein_2.chains[0].chain_letter}"
    else:
        tmp_protein_2_only_one_chain = True
        tmp_sub_analysis_name_2 = f"{the_protein_2_name};A"

    if tmp_protein_1_only_one_chain and tmp_protein_2_only_one_chain:
        tmp_analysis_run_name = f"{tmp_sub_analysis_name_1}_vs_{tmp_sub_analysis_name_2}"
    else:
        tmp_analysis_run_name = ""
    return ("result", tmp_analysis_run_name)


def preview_image(a_placeholder_1: int, a_placeholder_2: int) -> tuple:
    # TODO: the renderer should be changeable
    try:
        cmd.ray(2400, 2400, renderer=int(0))
    except pymol.CmdException:
        logger.warning("Unexpected exception.")
    return 0, ""


def create_drawn_image(an_image_filepath: str, the_app_settings: "settings.Settings") -> tuple:
    cmd.bg_color(the_app_settings.image_background_color)
    try:
        cmd.draw(2400, 2400)
        cmd.png(an_image_filepath, dpi=300)
    except pymol.CmdException:
        logger.warning("Unexpected exception.")
    cmd.bg_color("black")
    return 0, ""


def load_protein_pymol_session(a_protein, the_pymol_session_manager, needs_to_be_reinitialized_flag: bool = False) -> tuple:

    #return the_pymol_session_manager, True
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
            db_manager.open_project_database()
            tmp_protein = the_interface_manager.get_current_active_protein_object()
            tmp_protein.save_pymol_session_as_base64_string()
            db_manager.update_pymol_session_of_protein(
                tmp_protein.get_id(),
                tmp_protein.pymol_session
            )
            db_manager.close_project_database()
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
            db_manager.open_project_database()
            tmp_protein_pair = the_interface_manager.get_current_active_protein_pair_object()
            tmp_protein_pair.save_session_of_protein_pair()
            db_manager.update_pymol_session_of_protein_pair(
                tmp_protein_pair.get_id(),
                tmp_protein_pair.pymol_session
            )
            db_manager.close_project_database()
    except Exception as e:
        logger.error(f"Unexpected error occurred. Exception: {e}")
        return 0, False
    else:
        return 0, True
