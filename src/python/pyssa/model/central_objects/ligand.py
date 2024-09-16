import logging
import pickle

from rdkit import Chem
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)

__docformat__ = "google"


class Ligand:
  """Class for storing a ligand."""

  def __init__(self, a_name, a_structure: Chem.Mol) -> None:
    """Constructor.

    Args:
      a_structure: A Chem.Mol object representing the structure of the ligand

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if a_structure is None:
      default_logging.append_to_log_file(logger, "a_structure is None.", logging.ERROR)
      raise exception.NoneValueError("a_structure is None.")
    # </editor-fold>
    # <editor-fold desc="Instance attributes">
    self.id: int = 0  # Default value is 0 and should be changed if loaded from DB
    """The id of the ligand molecule"""
    self.name: str = a_name
    """The name of the ligand molecule"""
    self.structure: Chem.Mol = a_structure
    """A Chem.Mol object representing the ligand structure information"""
    self.session: bytes = bytes()
    """The session of the ligand."""
    self.cross_references: int = 0
    """The number of cross-references used in complexes"""
    # </editor-fold>

  def get_structure_as_binary(self) -> bytes:
    """Creates a binary representation of the structure."""
    return pickle.dumps(self.structure)
