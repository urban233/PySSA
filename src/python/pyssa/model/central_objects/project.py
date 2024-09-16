import logging
from typing import TYPE_CHECKING
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)

if TYPE_CHECKING:
  from pyssa.model.central_objects import sequence
  from pyssa.model.central_objects import protein
  from pyssa.model.central_objects import ligand
  from pyssa.model.central_objects import protein_ligand_complex
  from pyssa.model.central_objects import protein_protein_complex
  from pyssa.model.central_objects.dummy import dummy_protein
  from pyssa.model.central_objects.dummy import dummy_protein_ligand_complex
  from pyssa.model.central_objects.dummy import dummy_protein_protein_complex


class Project:
  """Stores all information about a project."""

  def __init__(self, a_project_name: str) -> None:
    """Constructor.

    Args:
      a_project_name: The name of the project.

    Raises:
      exception.NoneValueError: If `a_project_name` is None.
      exception.IllegalArgumentError: If `a_project_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_project_name is None:
      default_logging.append_to_log_file(logger, "a_project_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_project_name is None.")
    if a_project_name == "":
      default_logging.append_to_log_file(logger, "a_project_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_project_name is an empty string.")
    # </editor-fold>
    # <editor-fold desc="Instance attributes">
    self.name = a_project_name
    """A name for the project."""
    self.sequences: list["sequence.Sequence"] = []
    """An array of sequences"""
    self.proteins: list["protein.Protein"] = []
    """An array of proteins"""
    self.ligands: list["ligand.Ligand"] = []
    """An array of ligands"""
    self.protein_ligand_complexes: list["protein_ligand_complex.ProteinLigandComplex"] = []
    """An array of protein ligand complexes"""
    self.protein_protein_complexes: list["protein_protein_complex.ProteinProteinComplex"] = []
    """An array of protein-protein complexes"""
    self.dummy_proteins: list["dummy_protein.DummyProtein"] = []
    """An array of dummy proteins"""
    self.dummy_protein_ligand_complexes: list["dummy_protein_ligand_complex.DummyProteinLigandComplex"] = []
    """An array of dummy protein-ligand complexes"""
    self.dummy_protein_protein_complexes: list["dummy_protein_protein_complex.DummyProteinProteinComplex"] = []
    """An array of dummy protein-protein complexes"""
    # </editor-fold>

  def add_protein(self, a_protein: "protein.Protein") -> None:
    """Adds a protein to the project.

    Args:
      a_protein: The protein instance to add

    Raises:
      exception.NoneValueError: If `a_protein` is None.
    """
    # <editor-fold desc="Checks">
    if a_protein is None:
      default_logging.append_to_log_file(logger, "a_protein is None.", logging.ERROR)
      raise exception.NoneValueError("a_protein is None.")
    # </editor-fold>
    try:
      self.proteins.append(a_protein)
      default_logging.append_to_log_file(logger, f"The protein {a_protein} was added to the project {self.name}.", logging.INFO)
    except Exception as e:
      default_logging.append_to_log_file(
        logger,
        f"An error occurred while adding the protein {a_protein} to the project {self.name}. The error message is {e.__str__()}.",
        logging.ERROR
      )
