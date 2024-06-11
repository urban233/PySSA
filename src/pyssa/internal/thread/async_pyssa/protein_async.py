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
"""Module for all asynchronous functions that are related to the protein object."""
import logging
from typing import Optional

import zmq

from src.pyssa.controller import database_manager
from src.pyssa.internal.data_structures import protein
from src.pyssa.io_pyssa import bio_data
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import constants, enums, exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"

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


def clean_protein_update(
    a_protein: "protein.Protein",
    the_database_filepath: str,
    the_main_socket: zmq.Socket,
    the_general_purpose_socket: zmq.Socket,
) -> tuple[str, Optional[protein.Protein], Optional[bool]]:
  """Cleans a protein by removing all solvent and sugar molecules in the current molecule object.

  Args:
      a_protein (protein.Protein): The protein object to clean.
      the_database_filepath (str): The filepath to the database file of the project.
      the_main_socket (zmq.Socket): The main socket of the auxiliary pymol.
      the_general_purpose_socket (zmq.Socket): The general purpose socket of the auxiliary pymol.

  Returns:
      A tuple with ("result", a_protein or None, tmp_more_than_one_ca_atom or None).
  """
  # <editor-fold desc="Checks">
  if a_protein is None:
    logger.error("a_protein is None.")
    return "", None, None
  if the_database_filepath is None or the_database_filepath == "":
    logger.error("the_database_filepath is either None or an empty string.")
    return "", None, None
  if the_main_socket is None:
    logger.error("the_main_socket is None.")
    return "", None, None
  if the_general_purpose_socket is None:
    logger.error("the_general_purpose_socket is None.")
    return "", None, None

  # </editor-fold>

  try:
    tmp_more_than_one_ca_atom = a_protein.clean_protein(the_main_socket, the_general_purpose_socket)
    if len(a_protein.get_pdb_data()) == 0:
      logger.error("No PDB data found after cleaning process!")
      return "", None, None

    with database_manager.DatabaseManager(the_database_filepath) as db_manager:
      db_manager.update_protein_pdb_atom_data(
        a_protein.get_id(), a_protein.get_pdb_data()
      )
      a_protein.set_pdb_data([])
      for i, tmp_chain in enumerate(a_protein.chains, start=0):
        if tmp_chain.chain_type == enums.ChainTypeEnum.NON_PROTEIN_CHAIN.value:
          db_manager.delete_specific_chain(a_protein.get_id(), tmp_chain.get_id())
    constants.PYSSA_LOGGER.info(
      "The protein %s has been cleaned.", a_protein.get_molecule_object()
    )
  except Exception as e:
    logger.error(e)
    return "", None, None
  else:
    return "result", a_protein, tmp_more_than_one_ca_atom


def save_selected_protein_structure_as_pdb_file(
    a_protein: "protein.Protein",
    a_filepath: str,
    the_database_filepath: str,
) -> tuple[str]:
  """Saves the selected protein structure as a PDB file.

  Args:
      a_protein (protein.Protein): The selected protein object.
      a_filepath (str): The file path where the PDB file should be saved.
      the_database_filepath (str): The file path of the database.

  Returns:
      A tuple containing the result of the operation.
      If successful, the tuple will contain a single string with the value "result". Otherwise, an empty string will be returned.
  """
  # <editor-fold desc="Checks">
  if a_protein is None:
    logger.error("a_protein is None.")
    return ("",)
  if a_filepath is None or a_filepath == "":
    logger.error("a_filepath is either None or an empty string.")
    return ("",)
  if the_database_filepath is None or the_database_filepath == "":
    logger.error("the_database_filepath is either None or an empty string.")
    return ("",)

  # </editor-fold>

  try:
    with database_manager.DatabaseManager(the_database_filepath) as db_manager:
      tmp_pdb_atom_data = db_manager.get_pdb_atoms_of_protein(
          a_protein.get_id()
      )
      tmp_pdb_atom_dict_1 = [
          {key.value: value for key, value in zip(enums.PdbAtomEnum, t)}
          for t in tmp_pdb_atom_data
      ]
      a_protein.set_pdb_data(tmp_pdb_atom_dict_1)
      bio_data.build_pdb_file(a_protein.get_pdb_data(), a_filepath)
      # reset pdb data to reduce memory space
      a_protein.set_pdb_data([])
  except Exception as e:
    logger.error(e)
    return ("",)
  else:
    return ("result",)


def rename_selected_protein_structure(
    a_protein: "protein.Protein",
    the_new_protein_name: str,
    the_database_filepath: str,
) -> tuple:
  """Renames a certain protein from a project.

  Args:
      a_protein (protein.Protein): The protein object to rename.
      the_new_protein_name (str): The new name for the given protein.
      the_database_filepath (str): The filepath of the project database.

  Returns:
      a tuple with ("result" or "", an_existing_protein_object or None)
  """
  # <editor-fold desc="Checks">
  if a_protein is None:
    logger.error("a_protein is None.")
    return "", None
  if the_new_protein_name is None or the_new_protein_name == "":
    logger.error("the_new_protein_name is either None or an empty string.")
    return "", None
  if the_database_filepath is None or the_database_filepath == "":
    logger.error("the_database_filepath is either None or an empty string.")
    return "", None

  # </editor-fold>

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
        the_new_protein_name,
        tmp_old_protein_name,
        a_protein.get_id(),
    )
  return "result", a_protein
