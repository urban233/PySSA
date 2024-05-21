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
"""Module for the project class."""
import collections
import logging
import pathlib
import platform
from typing import TYPE_CHECKING, Optional
from PyQt5 import QtCore
import numpy as np
from Bio import SeqRecord

from pyssa.logging_pyssa import log_handlers
from pyssa.io_pyssa import safeguard
from pyssa.internal.data_structures.data_classes import basic_protein_info
from pyssa.util import exception

if TYPE_CHECKING:
  from pyssa.internal.data_structures import protein
  from pyssa.internal.data_structures import protein_pair

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class Project:
  """Class for the projects used in the plugin."""

  # <editor-fold desc="Class attributes">
  _id: int
  """The id of the project in the database."""

  _project_name: str
  """The name of the project."""

  _workspace: pathlib.Path
  """The absolute path of the current workspace."""

  _operating_system = platform.system()
  """The used OS."""

  proteins: list["protein.Protein"]
  """A list of all protein objects of the project."""

  protein_pairs: list["protein_pair.ProteinPair"]
  """A list of all protein_pair objects of the project."""

  sequences: list[SeqRecord.SeqRecord]
  """A list of all sequence objects of the project."""

  # </editor-fold>

  def __init__(
      self,
      a_project_name: str = "",
      a_workspace_path: pathlib.Path = pathlib.Path(""),
  ) -> None:
    """Constructor.

    Args:
        a_project_name (str): The name of the project.
        a_workspace_path (pathlib.Path): The path of the workspace.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_project_name is None:
      logger.error("a_project_name is None.")
      raise exception.IllegalArgumentError("a_project_name is None.")
    if a_workspace_path is None:
      logger.error("a_workspace_path is None.")
      raise exception.IllegalArgumentError("a_workspace_path is None.")

    # </editor-fold>

    self._project_name: str = a_project_name
    self._workspace: pathlib.Path = a_workspace_path
    self.proteins: list["protein.Protein"] = []
    self.protein_pairs: list["protein_pair.ProteinPair"] = []
    self.sequences: list[SeqRecord.SeqRecord] = []

  def set_id(self, an_id: int) -> None:
    """Sets the ID of the object.

    Args:
        an_id (int): The ID to be set.

    Raises:
        exception.IllegalArgumentError: If an_id is either None or a value less than 0.
    """
    # <editor-fold desc="Checks">
    if an_id is None or an_id < 0:
      logger.error("an_id is None.")
      raise exception.IllegalArgumentError(
          "an_id is either None or a value less than 0."
      )

    # </editor-fold>

    self._id = an_id

  def get_id(self) -> int:
    """Gets the ID of the object.

    Returns:
        The ID of the object.
    """
    return self._id

  # def set_workspace_path(self, a_workspace_path: pathlib.Path) -> None:
  #     """Setter for the workspace path.
  #
  #     Args:
  #         a_workspace_path (pathlib.Path): the new workspace path
  #
  #     Raises:
  #         ValueError: if one of the arguments are None or path does not exist
  #     """
  #     # <editor-fold desc="Checks">
  #     safeguard.Safeguard.check_if_value_is_not_none(a_workspace_path, logger)
  #
  #     if not safeguard.Safeguard.check_filepath(a_workspace_path):
  #         msg = "The given workspace path does not exists!"
  #         logger.error(msg)
  #         raise ValueError(msg)
  #
  #     # </editor-fold>
  #
  #     self._workspace = a_workspace_path

  def add_existing_protein(self, a_protein: "protein.Protein") -> None:
    """Adds an existing protein object to the project.

    Args:
        a_protein (protein.Protein): An existing protein object.

    Raises:
        exception.IllegalArgumentError: If a_protein is None.
    """
    # <editor-fold desc="Checks">
    if a_protein is None:
      logger.error("a_protein is None.")
      raise exception.IllegalArgumentError("a_protein is None.")

    # </editor-fold>

    self.proteins.append(a_protein)

  def add_protein_pair(
      self, a_protein_pair: "protein_pair.ProteinPair"
  ) -> None:
    """Adds an existing protein_pair object to the project.

    Args:
        a_protein_pair (protein_pair.ProteinPair): An existing protein pair object.

    Raises:
        ValueError: If one of the arguments are None
    """
    # <editor-fold desc="Checks">
    if a_protein_pair is None:
      logger.error("a_protein_pair is None.")
      raise exception.IllegalArgumentError("a_protein_pair is None.")

    # </editor-fold>

    self.protein_pairs.append(a_protein_pair)

  def get_number_of_proteins(self) -> int:
    """Returns the number of proteins in the project.

    Returns:
        the number of proteins in the project.
    """
    return len(self.proteins)

  def is_sequence_as_protein_in_project(self, a_sequence_name: str) -> bool:
    """Checks if a given sequence name is present as a protein in the project.

    Args:
        a_sequence_name: A string representing the sequence name to check.

    Returns:
        A boolean value indicating whether the given sequence name is present as a protein in the project. Returns True if the sequence name is found, False otherwise.

    Raises:
        exception.IllegalArgumentError: If a_sequence_name is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_sequence_name is None or a_sequence_name == "":
      logger.error("a_sequence_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_sequence_name is either None or an empty string."
      )

    # </editor-fold>

    for tmp_protein in self.proteins:
      if tmp_protein.get_molecule_object() == a_sequence_name:
        return True
    return False

  def get_project_name(self) -> str:
    """Getter for the project name.

    Returns:
        The project name.
    """
    return self._project_name

  def get_database_filepath(self) -> pathlib.Path:
    """Gets the filepath of the database for the current project.

    Returns:
        The filepath of the database.
    """
    return pathlib.Path(f"{self._workspace}/{self.get_project_name()}.db")

  def search_protein(self, a_protein_name: str) -> Optional["protein.Protein"]:
    """Searches the project for a specific protein name.

    Args:
        a_protein_name (str): The name of the protein to search.

    Returns:
        A protein object or None if the protein was not found.

    Raises:
        exception.IllegalArgumentError: If a_protein_name is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_protein_name is None or a_protein_name == "":
      logger.error("a_protein_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_protein_name is either None or an empty string."
      )

    # </editor-fold>

    for tmp_protein in self.proteins:
      if tmp_protein.get_molecule_object() == a_protein_name:
        return tmp_protein
    print(
        f"No matching protein with the name {a_protein_name} found."
    )  # noqa: RET503
    return None

  def search_protein_pair(
      self, a_protein_pair_name: str
  ) -> Optional["protein_pair.ProteinPair"]:
    """Searches all protein_pairs within the project and returns true if the project contains the pair.

    Args:
        a_protein_pair_name: the name of the protein pair to search

    Raises:
        exception.IllegalArgumentError: If a_protein_pair_name is either None or an empty string.

    Returns:
        A protein pair or None if the protein pair was not found.
    """
    # <editor-fold desc="Checks">
    if a_protein_pair_name is None or a_protein_pair_name == "":
      logger.error("a_protein_pair_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_protein_pair_name is either None or an empty string."
      )

    # </editor-fold>

    for tmp_protein_pair in self.protein_pairs:
      if tmp_protein_pair.name == a_protein_pair_name:
        return tmp_protein_pair
    print(
        f"No matching protein with the name {a_protein_pair_name} found."
    )  # noqa: RET503
    return None

  def search_sequence(self, a_seq_name: str) -> Optional[SeqRecord.SeqRecord]:
    """Searches the project for a specific seq name.

    Args:
        a_seq_name (str): the name of the protein to search

    Returns:
        A seq record or None if the sequence was not found.

    Raises:
        exception.IllegalArgumentError: If a_seq_name is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_seq_name is None or a_seq_name == "":
      logger.error("a_seq_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_seq_name is either None or an empty string."
      )

    # </editor-fold>

    for tmp_seq_record in self.sequences:
      if tmp_seq_record.name == a_seq_name:
        return tmp_seq_record
    print(
        f"No matching sequence with the name {a_seq_name} found."
    )  # noqa: RET503
    return None

  def check_if_protein_is_in_any_protein_pair(
      self, a_protein_name: str
  ) -> bool:
    """Checks if a certain protein is part of an existing protein pair.

    Args:
        a_protein_name (str): the name of the protein to check

    Returns:
        A boolean value indicating whether the protein is part of an existing protein pair.

    Raises:
        exception.IllegalArgumentError: If a_protein_name is None.
    """
    # <editor-fold desc="Checks">
    if a_protein_name is None:
      logger.error("a_protein_name is None.")
      raise exception.IllegalArgumentError("a_protein_name is None.")

    # </editor-fold>

    protein_obj = self.search_protein(a_protein_name)
    tmp_selected_protein_id = protein_obj.get_id()

    for tmp_protein_pair in self.protein_pairs:
      tmp_protein_1_id = tmp_protein_pair.protein_1.get_id()
      tmp_protein_2_id = tmp_protein_pair.protein_2.get_id()
      if (
          tmp_selected_protein_id == tmp_protein_1_id
          or tmp_selected_protein_id == tmp_protein_2_id
      ):
        return True
    return False

  def delete_specific_sequence(self, a_seq_name: str) -> None:
    """Deletes a specific sequence from the list of sequences.

    Args:
        a_seq_name (str): The name of the sequence to be deleted.

    Raises:
        exception.IllegalArgumentError: If a_seq_name is either None or an empty string.
        ValueError: If a_seq_name is not in the list.
    """
    # <editor-fold desc="Checks">
    if a_seq_name is None or a_seq_name == "":
      logger.error("a_seq_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_seq_name is either None or an empty string."
      )

    # </editor-fold>

    tmp_seq_record = self.search_sequence(a_seq_name)
    i = 0
    for tmp_seq in self.sequences:
      if tmp_seq.name == tmp_seq_record.name:
        self.sequences.pop(i)
        return
      i += 1
    else:
      raise ValueError("a_seq_name is not in the list.")

  def delete_specific_protein(self, a_protein_name: str) -> None:
    """Deletes a certain protein from the project based on the protein name.

    Args:
        a_protein_name (str): the name of the protein to delete

    Raises:
        exception.IllegalArgumentError: If a_protein_name is either None or an empty string.
        ValueError: If a_protein_name is not in the list.
    """
    # <editor-fold desc="Checks">
    if a_protein_name is None or a_protein_name == "":
      logger.error("a_protein_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_protein_name is either None or an empty string."
      )

    # </editor-fold>

    protein_obj = self.search_protein(a_protein_name)
    if protein_obj in self.proteins:
      self.proteins.remove(protein_obj)
    else:
      raise ValueError("a_protein_name is not in the list.")

  def delete_specific_protein_pair(self, a_protein_pair_name: str) -> None:
    """Deletes a certain protein from the project based on the protein name.

    Args:
        a_protein_pair_name (str): the name of the protein pair to delete

    Raises:
        exception.IllegalArgumentError: If a_protein_pair_name is either None or an empty string.
        ValueError: If a_protein_pair_name is not in the list.
    """
    # <editor-fold desc="Checks">
    if a_protein_pair_name is None or a_protein_pair_name == "":
      logger.error("a_protein_pair_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_protein_pair_name is either None or an empty string."
      )

    # </editor-fold>

    protein_pair_obj = self.search_protein_pair(a_protein_pair_name)
    if protein_pair_obj in self.protein_pairs:
      self.protein_pairs.remove(protein_pair_obj)
    else:
      raise ValueError("a_protein_pair_name is not in the list.")

  def convert_list_of_proteins_to_list_of_protein_infos(self) -> np.ndarray:
    """Converts a list of proteins to a list of protein infos.

    Returns:
        An numpy array of BasicProteinInfo objects.
    """
    tmp_protein_infos: collections.deque = collections.deque()
    for tmp_protein in self.proteins:
      tmp_protein_infos.append(
          basic_protein_info.BasicProteinInfo(
              tmp_protein.get_molecule_object(),
              tmp_protein.get_id(),
              self._project_name,
          ),
      )
    return np.array(list(tmp_protein_infos))
