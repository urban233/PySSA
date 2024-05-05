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
"""Module for all asynchronous functions that are related to the project object."""
import copy
import logging
import pathlib
import platform
import time

from pyssa.controller import database_manager
from pyssa.internal.data_structures import project, protein
from pyssa.internal.data_structures.data_classes import database_operation
from pyssa.internal.thread import database_thread
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants, enums

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def create_new_project(
    the_project_name: str,
    the_workspace_path: pathlib.Path,
    the_watcher,
    the_interface_manager,
) -> tuple:
    """Creates a new project object.

    Args:
        the_project_name: the project name to use for the new project.
        the_workspace_path: the current workspace path.
        an_add_protein_flag: a boolean value if a protein should be added or not.
        protein_source_information: either a filepath or a pdb id

    Returns:
        a tuple with ("results", a_new_project_obj)
    """
    # TODO: checks are needed
    tmp_project = project.Project(the_project_name, the_workspace_path)

    tmp_database_filepath = str(pathlib.Path(f"{the_workspace_path}/{the_project_name}.db"))
    with database_manager.DatabaseManager(tmp_database_filepath) as db_manager:
        tmp_project.set_id(db_manager.insert_new_project(tmp_project.get_project_name(), platform.system()))
    constants.PYSSA_LOGGER.info("Create empty project finished.")
    the_watcher.setup_blacklists(
        tmp_project,
        the_interface_manager.job_manager.get_queue(enums.JobType.PREDICTION),
        the_interface_manager.job_manager.get_queue(enums.JobType.DISTANCE_ANALYSIS),
        the_interface_manager.job_manager.current_prediction_job,
        the_interface_manager.job_manager.current_distance_analysis_job,
    )
    return "result", tmp_project, the_watcher, the_interface_manager


def create_use_project(
    the_project_name: str,
    the_workspace_path: pathlib.Path,
    the_proteins_to_add: list,
    the_watcher,
    the_interface_manager
) -> tuple:
    """Creates a new project object.

    Args:
        the_project_name: the project name to use for the new project.
        the_workspace_path: the current workspace path.
        the_proteins_to_add: a list of protein objects for the new project.

    Returns:
        a tuple with ("results", a_new_project_obj)
    """
    # TODO: checks are needed
    tmp_project = project.Project(the_project_name, the_workspace_path)

    tmp_database_filepath = str(pathlib.Path(f"{the_workspace_path}/{the_project_name}.db"))
    with database_manager.DatabaseManager(tmp_database_filepath) as db_manager:
        tmp_project.set_id(db_manager.insert_new_project(tmp_project.get_project_name(), platform.system()))
        for tmp_protein in the_proteins_to_add:
            tmp_protein_copy = copy.deepcopy(tmp_protein)
            tmp_protein_copy.db_project_id = tmp_project.get_id()
            tmp_project.add_existing_protein(tmp_protein_copy)
            tmp_protein_copy.set_id(db_manager.insert_new_protein(tmp_protein_copy))

    constants.PYSSA_LOGGER.info("Use project finished.")
    the_watcher.setup_blacklists(
        tmp_project,
        the_interface_manager.job_manager.get_queue(enums.JobType.PREDICTION),
        the_interface_manager.job_manager.get_queue(enums.JobType.DISTANCE_ANALYSIS),
        the_interface_manager.job_manager.current_prediction_job,
        the_interface_manager.job_manager.current_distance_analysis_job,
    )
    return "result", tmp_project


def open_project(
        tmp_project_name: str,
        tmp_project_database_filepath: str,
        the_interface_manager: "interface_manager.InterfaceManager",
        the_pymol_session_manager: "pymol_session_manager.PymolSessionManager",
        the_custom_progress_signal: "custom_signals.ProgressSignal",
        the_watcher: "watcher.Watcher"
) -> tuple:
    the_custom_progress_signal.emit_signal("Opening database ...", 10)
    with database_manager.DatabaseManager(tmp_project_database_filepath) as db_manager:
        the_custom_progress_signal.emit_signal("Setting up project ...", 30)
        tmp_project = db_manager.get_project_as_object(
            tmp_project_name,
            the_interface_manager.get_application_settings().workspace_path,
            the_interface_manager.get_application_settings(),
            the_custom_progress_signal
        )
    the_interface_manager.set_new_project(tmp_project)
    the_custom_progress_signal.emit_signal("Reinitializing PyMOL session ...", 96)
    the_pymol_session_manager.reinitialize_session()
    the_watcher.setup_blacklists(
        tmp_project,
        the_interface_manager.job_manager.get_queue(enums.JobType.PREDICTION),
        the_interface_manager.job_manager.get_queue(enums.JobType.DISTANCE_ANALYSIS),
        the_interface_manager.job_manager.current_prediction_job,
        the_interface_manager.job_manager.current_distance_analysis_job,
    )
    # except Exception as e:
    #     logger.error(f"Unexpected error occurred. Exception: {e}")
    #     return 1, 0, 0, 0
    # else:
    return 0, tmp_project, the_interface_manager, the_watcher


def close_project(the_database_thread: "database_thread.DatabaseThread", the_pymol_session_manager) -> tuple:
    try:
        the_database_thread.put_database_operation_into_queue(database_operation.DatabaseOperation(
            enums.SQLQueryType.CLOSE_PROJECT, (0, ""))
        )
        print(the_database_thread.queue_is_running())
        while the_database_thread.queue_is_running():
            time.sleep(0.5)
        the_pymol_session_manager.reinitialize_session()
    except Exception as e:
        logger.error(f"Unknown error occurred while waiting for the database thread to finish: {e}.")
        return False, "Waiting for database thread queue failed!"
    else:
        logger.info("Waiting for database thread queue finished.")
        return True, "Waiting for database thread queue finished."


def search_for_not_matching_proteins(the_main_view_state, the_proteins) -> tuple:
    tmp_proteins = the_main_view_state.get_not_matching_proteins(the_proteins)
    return "result", tmp_proteins


def add_protein_from_pdb_to_project(tmp_protein_name,
                                    the_interface_manager: "interface_manager.InterfaceManager") -> tuple:
    the_main_socket, the_general_purpose_socket = the_interface_manager.job_manager.get_general_purpose_socket_pair()
    with database_manager.DatabaseManager(the_interface_manager.get_current_project().get_database_filepath()) as db_manager:
        tmp_ref_protein = protein.Protein(tmp_protein_name.upper())
        tmp_ref_protein.set_id(db_manager.get_next_id_of_protein_table())
        tmp_ref_protein.db_project_id = the_interface_manager.get_current_project().get_id()
        tmp_ref_protein.add_protein_structure_data_from_pdb_db(tmp_protein_name.upper(), the_main_socket, the_general_purpose_socket)
        tmp_ref_protein.add_id_to_all_chains(db_manager.get_next_id_of_chain_table())

    tmp_ref_protein.create_new_pymol_session(the_main_socket, the_general_purpose_socket)
    the_interface_manager.add_protein_to_proteins_model(tmp_ref_protein)
    return 0, tmp_ref_protein


def add_protein_from_local_filesystem_to_project(tmp_protein_name,
                                                 the_interface_manager):
    pdb_filepath = pathlib.Path(tmp_protein_name)
    tmp_ref_protein = protein.Protein(
        pdb_filepath.name.replace(".pdb", "")
    )
    with database_manager.DatabaseManager(the_interface_manager.get_current_project().get_database_filepath()) as db_manager:
        tmp_ref_protein.set_id(db_manager.get_next_id_of_protein_table())
        tmp_ref_protein.db_project_id = the_interface_manager.get_current_project().get_id()
        the_main_socket, the_general_purpose_socket = the_interface_manager.job_manager.get_general_purpose_socket_pair()
        tmp_ref_protein.add_protein_structure_data_from_local_pdb_file(pdb_filepath, the_main_socket, the_general_purpose_socket)
        tmp_ref_protein.add_id_to_all_chains(db_manager.get_next_id_of_chain_table())
    tmp_ref_protein.create_new_pymol_session(the_main_socket, the_general_purpose_socket)
    the_interface_manager.add_protein_to_proteins_model(tmp_ref_protein)
    return 0, tmp_ref_protein
