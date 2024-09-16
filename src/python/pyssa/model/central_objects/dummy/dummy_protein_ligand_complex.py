import logging

from pyssa.model.central_objects import protein
from pyssa.model.central_objects import ligand
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)


class DummyProteinLigandComplex:
  """Class for storing a dummy protein ligand pair."""

  def __init__(self,
               a_name: str,
               a_protein: "protein.Protein",
               a_ligand: "ligand.Ligand") -> None:
    """Constructor.

    Args:
      a_name: The name of the protein-ligand complex
      a_protein: The protein of the complex
      a_ligand: The ligand of the complex

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
    if a_protein is None:
      default_logging.append_to_log_file(logger, "a_protein is None.", logging.ERROR)
      raise exception.NoneValueError("a_protein is None.")
    if a_ligand is None:
      default_logging.append_to_log_file(logger, "a_ligand is None.", logging.ERROR)
      raise exception.NoneValueError("a_ligand is None.")
    # </editor-fold>
    # <editor-fold desc="Instance attributes">
    self.name: str = a_name
    """The name of the protein-ligand dummy complex"""
    self.protein: "protein.Protein" = a_protein
    """The protein object of the complex"""
    self.ligand: "ligand.Ligand" = a_ligand
    """The ligand object of the complex"""
    # </editor-fold>
