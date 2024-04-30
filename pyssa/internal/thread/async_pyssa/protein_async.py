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
"""Module for all asynchronous functions that are related to the protein object."""
import logging
from pyssa.controller import database_manager
from pyssa.io_pyssa import bio_data
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants, enums, exit_codes

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


# def clean_protein_new(
#     a_protein_name: str,
#     a_project: "project.Project",
# ) -> tuple:
#     """Cleans a protein by creating a duplicate and removing all solvent and sugar molecules.
#
#     Args:
#         a_protein_name: the name of the protein to clean.
#         a_project: the current project.
#
#     Returns:
#         a tuple with ("result", a_updates_project_object)
#     """
#     # TODO: needs checks
#     tmp_protein = a_project.search_protein(
#         a_protein_name,
#     )
#     clean_tmp_protein = tmp_protein.clean_protein(new_protein=True)
#     constants.PYSSA_LOGGER.info("The protein %s has been cleaned.", clean_tmp_protein.get_molecule_object())
#     a_project.add_existing_protein(clean_tmp_protein)
#     return ("result", a_project)


def clean_protein_update(a_protein: "protein.Protein",
                         the_database_filepath: str,
                         the_main_socket,
                         the_general_purpose_socket) -> tuple:
    """Cleans a protein by removing all solvent and sugar molecules in the current molecule object.

    Args:
        a_protein: the protein object to clean.
        the_database_filepath: the filepath to the database file of the project.
        the_main_socket: the main socket of the auxiliary pymol.
        the_general_purpose_socket: the general purpose socket of the auxiliary pymol.

    Returns:
        a tuple with ("result", a_updates_project_object)
    """
    # TODO: needs checks
    # TODO: needs tests!
    a_protein.clean_protein(the_main_socket, the_general_purpose_socket)
    if len(a_protein.get_pdb_data()) == 0:
        logger.error("No PDB data found after cleaning process!")
        return ("result", None)

    with database_manager.DatabaseManager(the_database_filepath) as db_manager:
        db_manager.update_protein_pdb_atom_data(a_protein.get_id(), a_protein.get_pdb_data())
        a_protein.set_pdb_data([])
    constants.PYSSA_LOGGER.info("The protein %s has been cleaned.", a_protein.get_molecule_object())
    return ("result", a_protein)


def save_selected_protein_structure_as_pdb_file(
    a_protein: "protein.Protein",
    a_filepath: str,
    the_database_filepath: str,
) -> tuple:
    """Saves a given protein structure to a pdb file."""
    with database_manager.DatabaseManager(the_database_filepath) as db_manager:
        tmp_pdb_atom_data = db_manager.get_pdb_atoms_of_protein(a_protein.get_id())
        tmp_pdb_atom_dict_1 = [{key.value: value for key, value in zip(enums.PdbAtomEnum, t)} for t in tmp_pdb_atom_data]
        a_protein.set_pdb_data(tmp_pdb_atom_dict_1)
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
    # TODO: this needs to be reimplemented using the new PyMOLInterface class (through the pymol_session_manager!)
    # try:
    #     cmd.set_name(tmp_old_protein_name, a_protein.get_molecule_object())
    # except Exception as e:
    #     logger.error(f"Renaming the protein in PyMOL failed. {e}")
    #     raise RuntimeError("Renaming the protein in PyMOL failed. Maybe it is not opened in a session?")
    # Update in database
    with database_manager.DatabaseManager(the_database_filepath) as db_manager:
        db_manager.update_protein_name(
            the_new_protein_name, tmp_old_protein_name, a_protein.get_id()
        )
    return ("result", a_protein)
