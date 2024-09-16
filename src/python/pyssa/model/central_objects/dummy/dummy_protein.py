import logging
from typing import Optional

from Bio import SeqRecord

from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)


class DummyProtein:
  """Class for storing a dummy protein which gets predicted."""

  def __init__(self, a_name: str, a_chain_sequence_map: dict) -> None:
    """Constructor.

    Args:
      a_name: The name of the dummy protein.

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
    if a_chain_sequence_map is None:
      default_logging.append_to_log_file(logger, "a_chain_sequence_map is None.", logging.ERROR)
      raise exception.NoneValueError("a_chain_sequence_map is None.")
    # </editor-fold>
    # <editor-fold desc="Instance attributes">
    self.name: str = a_name
    """The name of the dummy protein."""
    self.chain_sequence_map: dict[str, SeqRecord.SeqRecord] = a_chain_sequence_map
    """A map containing the protein chain with the corresponding sequence"""
    # </editor-fold>
