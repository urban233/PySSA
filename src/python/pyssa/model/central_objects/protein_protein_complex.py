import logging
import pickle

from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)


class ProteinProteinComplex:
  """Class for storing a protein-protein complex."""

  def __init__(self,
               a_name: str,
               a_protein_1: "protein.Protein",
               a_protein_2: "protein.Protein") -> None:
    """Constructor.

    Args:
      a_name: The name of the protein-protein complex
      a_protein_1: The first protein of the complex
      a_protein_2: The second protein of the complex

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
    if a_protein_1 is None:
      default_logging.append_to_log_file(logger, "a_protein_1 is None.", logging.ERROR)
      raise exception.NoneValueError("a_protein_1 is None.")
    if a_protein_2 is None:
      default_logging.append_to_log_file(logger, "a_protein_2 is None.", logging.ERROR)
      raise exception.NoneValueError("a_protein_2 is None.")
    # </editor-fold>
    # <editor-fold desc="Instance attributes">
    self.name: str = a_name
    """The name of the protein-protein complex"""
    self.protein_1: "protein.Protein" = a_protein_1
    """The first protein object of the complex"""
    self.protein_2: "protein.Protein" = a_protein_2
    """The second protein object of the complex"""
    self.results = ""  # TODO: What datatype should this be?
    self.session: bytes = bytes()
    """A base64 encoded string containing the protein session"""
    # </editor-fold>

  def get_results_as_binary(self) -> bytes:
    """Creates a binary representation of the results object."""
    return pickle.dumps(self.results)