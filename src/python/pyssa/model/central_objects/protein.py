import logging
import pathlib
import pickle
from typing import Optional
from Bio.PDB import Structure
from Bio.PDB import Polypeptide
from Bio.PDB import PDBIO
from Bio import SeqRecord
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)
__docformat__ = "google"


class Protein:
  """Class for storing a protein."""

  def __init__(self, a_name: str, a_structure: Optional[Structure.Structure] = None):
    """Constructor.

    Args:
      a_name: The name of the protein
      a_structure: The structure of the protein (Default: None)

    Raises:
      exception.NoneValueError: If `a_name` is None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    # </editor-fold>
    # <editor-fold desc="Instance attributes">
    self.structure: Optional[Structure.Structure] = a_structure
    """
    The structure information of the protein
    More info: https://biopython.org/docs/1.75/api/Bio.PDB.Structure.html
    """
    self.name: str = a_name
    """The name of the protein"""
    self.id: int = 0  # Default value is 0 and should be changed if loaded from DB
    """The id of the database entry"""
    self.chain_sequence_map: dict[str, SeqRecord.SeqRecord] = {}
    """A map containing the protein chain with the corresponding sequence"""
    self.selection: str  # Maybe it's better if the session manager does the selection?
    """A selection used in a molecular viewer."""
    self.session: str = "generic_session_string"
    """A base64 encoded string containing the protein session"""
    self.cross_references: int = 0
    """The number of cross-references used in complexes"""
    # </editor-fold>
    if a_structure is not None:
      self._fill_chain_sequence_map()
    else:
      default_logging.append_to_log_file(logger, f"a_structure is None.", logging.WARNING)

  def _fill_chain_sequence_map(self) -> None:
    """Fills the chain sequence map with the chain information of the structure."""
    for tmp_chain in self.structure.get_chains():
      self.chain_sequence_map[tmp_chain.id] = SeqRecord.SeqRecord(name=self.name, seq=Polypeptide.Polypeptide(tmp_chain).get_sequence(), id=tmp_chain.id)

  def dump_to_pdb_file(self, a_filepath: pathlib.Path):
    """Creates a pdb file from the protein structure

    Args:
      a_filepath: The filepath of the pdb file to create
    """
    # <editor-fold desc="Checks">
    if a_filepath is None:
      raise exception.NoneValueError("a_filepath is None.")
    if not a_filepath.parent.exists():
      raise exception.IllegalArgumentError(f"a_filepath does not exist: {a_filepath}")
    # </editor-fold>
    tmp_io = PDBIO()
    tmp_io.set_structure(self.structure)
    tmp_io.save(str(a_filepath))

  def get_structure_as_binary(self) -> bytes:
    """Creates a binary representation of the structure."""
    return pickle.dumps(self.structure)

  def get_session_as_binary(self) -> bytes:
    """Creates a binary representation of the session string."""
    return pickle.dumps(self.session)

  def get_chain_sequence_map_as_binary(self) -> bytes:
    """Creates a binary representation of chain sequence dict."""
    return pickle.dumps(self.chain_sequence_map)
