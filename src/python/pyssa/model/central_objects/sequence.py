import collections
import logging
import pathlib
import pickle
from typing import Optional
from Bio import SeqRecord
from Bio.Seq import Seq
from pyssa.model.preference import central_objects_definitions
from pyssa.model.util import exception
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)

__docformat__ = "google"


class Sequence:
  """Class for storing a protein sequence."""

  def __init__(self, a_name: str, a_list_of_sequences: list[str], a_list_of_chain_ids: Optional[list[str]] = None):
    """Constructor.

    Args:
      a_name: The name of the sequence object
      a_list_of_sequences: A list of sequences
      a_list_of_chain_ids: A list of chain letters (optional, default: None)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If any of the arguments are an empty string or an_id is longer than one char.
    """
    # <editor-fold desc="Checks">
    if a_list_of_sequences is None:
      default_logging.append_to_log_file(logger, "a_list_of_sequences is None.", logging.ERROR)
      raise exception.NoneValueError("a_list_of_sequences is None.")
    if len(a_list_of_sequences) == 0:
      default_logging.append_to_log_file(logger, "a is an empty list.", logging.ERROR)
      raise exception.IllegalArgumentError("a_list_of_sequences is an empty list.")
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if len(a_list_of_chain_ids) == 0:
      default_logging.append_to_log_file(logger, "a is a_list_of_chain_ids empty list.", logging.ERROR)
      raise exception.IllegalArgumentError("a_list_of_chain_ids is an empty list.")
    # </editor-fold>
    # <editor-fold desc="Instance attributes">
    self.name = a_name
    """The name of the sequence object"""
    self.chain_sequence_map: dict[str, SeqRecord.SeqRecord] = {}
    """A map that contains the chain letter and the corresponding sequence as SeqRecord object"""
    self.type: central_objects_definitions.SeqTypesEnum = central_objects_definitions.SeqTypesEnum.PROTEIN  # Set default to PROTEIN as long as PyDD supports no other sequence types
    """A type that defines if the sequence is a protein or non protein sequence"""
    self._fill_chain_sequence_map(a_name, a_list_of_sequences, a_list_of_chain_ids)
    # </editor-fold>

  def _fill_chain_sequence_map(self, a_name: str, a_list_of_sequences: list[str], a_list_of_chain_ids: Optional[list[str]] = None) -> None:
    """Fills the object's chain sequence map.

    Args:
      a_name: The name of the sequence object
      a_list_of_sequences: A list of sequences
      a_list_of_chain_ids: A list of chain letters (optional, default: None)

    Notes:
      Caution: No arguments checks!
      Reason: these arguments are the same as in the constructor and therefore already checked.
    """
    for i, tmp_sequence in enumerate(a_list_of_sequences):
      self.chain_sequence_map[a_list_of_chain_ids[i]] = SeqRecord.SeqRecord(Seq(tmp_sequence), id=a_list_of_chain_ids[i], name=a_name)

  def get_chain_sequence_map_as_binary(self) -> bytes:
    """Creates a binary representation of the chain sequence map."""
    return pickle.dumps(self.chain_sequence_map)

  @staticmethod
  def from_single_seq_record(a_seq_record: SeqRecord.SeqRecord) -> "Sequence":
    """Alternative constructor.

    Args:
      a_seq_record: A SeqRecord object to use for building the sequence object

    Raises:
      exception.NoneValueError: If `a_seq_record` is None.
    """
    # <editor-fold desc="Checks">
    if a_seq_record is None:
      default_logging.append_to_log_file(logger, "a_seq_record is None.", logging.ERROR)
      raise exception.NoneValueError("a_seq_record is None.")
    # </editor-fold>
    tmp_sequence = Sequence(a_seq_record.name, [a_seq_record.seq], ["A"])
    return tmp_sequence

  @staticmethod
  def from_chain_sequence_map(a_name: str, a_chain_sequence_map: dict[str, SeqRecord.SeqRecord]) -> "Sequence":
    """Alternative constructor that uses a chain sequence map.

    Args:
      a_name: The name of the sequence
      a_chain_sequence_map: A map containing the chain and corresponding sequence

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string or a_chain_sequence_map is an empty dict.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if a_chain_sequence_map is None:
      default_logging.append_to_log_file(logger, "a_chain_sequence_map is None.", logging.ERROR)
      raise exception.NoneValueError("a_chain_sequence_map is None.")
    if len(a_chain_sequence_map) == 0:
      default_logging.append_to_log_file(logger, "a_chain_sequence_map is an empty dict.", logging.ERROR)
      raise exception.IllegalArgumentError("a_chain_sequence_map is an empty dict.")
    # </editor-fold>
    tmp_sequences: collections.deque = collections.deque()
    tmp_chain_ids: collections.deque = collections.deque()
    for tmp_key in a_chain_sequence_map:
      tmp_sequences.append(a_chain_sequence_map[tmp_key])
      tmp_chain_ids.append(tmp_key)
    return Sequence(a_name, list(tmp_sequences), list(tmp_chain_ids))
