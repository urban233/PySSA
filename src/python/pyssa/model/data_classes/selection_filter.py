import logging
import re
from dataclasses import dataclass
from pyssa.model.preference import model_definitions
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)


@dataclass
class SelectionFilter:
  """Class for storing a selection filter string in its parts"""
  # <editor-fold desc="Class attributes">
  # Class attributes are used here because this is a python data class
  full_selection: str
  """Stores the full selection filter string"""
  protein: str
  """Contains the protein part of the selection filter string"""
  chain: str
  """Contains the chain part of the selection filter string"""
  residue_name: str
  """Contains the residue name part of the selection filter string"""
  residue_number: str
  """Contains the residue number part of the selection filter string"""
  atom: str
  """Contains the atom part of the selection filter string"""
  # </editor-fold>

  def __init__(self, a_selection_filter_string: str) -> None:
    """Constructor.

    Args:
      a_selection_filter_string: A string representing the selection filter

    Raises:
      exception.NoneValueError: If `a_selection_filter_string` is None.
      exception.IllegalArgumentError: If `a_selection_filter_string` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_selection_filter_string is None:
      default_logging.append_to_log_file(logger, "a_selection_filter_string is None.", logging.ERROR)
      raise exception.NoneValueError("a_selection_filter_string is None.")
    if a_selection_filter_string == "":
      default_logging.append_to_log_file(logger, "a_selection_filter_string is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_selection_filter_string is an empty string.")
    # </editor-fold>
    # <editor-fold desc="Instance attributes">
    self.full_selection: str = a_selection_filter_string
    self.protein: str = ""
    self.chain: str = ""
    self.residue_name: str = ""
    self.residue_number: str = ""
    self.atom: str = ""
    # </editor-fold>

    tmp_selection_filter_parts = a_selection_filter_string.split(",")
    if not tmp_selection_filter_parts:
      return

    tmp_validators = [
      (r'^[a-zA-Z0-9_-]+$', 'protein'),
      (r'^[A-Z](\+[A-Z])*\+?$', 'chain'),
      (self._validate_residue_name, 'residue_name'),
      (r'^[0-9]+([+-][0-9]+)*[+-]?$', 'residue_number'),
      (r'^[A-Z]{1,3}(\+[A-Z]{1,3})*\+?$', 'atom')
    ]
    if len(tmp_selection_filter_parts) > len(tmp_validators):
      default_logging.append_to_log_file(
        logger,
        f"Invalid selection filter: {a_selection_filter_string}.",
        logging.ERROR
      )
      raise exception.IllegalArgumentError(f"Invalid selection filter: {a_selection_filter_string}.")
    for i, part in enumerate(tmp_selection_filter_parts):
      validator, attr = tmp_validators[i]
      if callable(validator):
        is_valid = validator(part)
      else:
        is_valid = self._validate_input_based_on_regex(part, validator)

      if is_valid:
        setattr(self, attr, part)
      else:
        default_logging.append_to_log_file(
          logger,
          f"Invalid selection filter: {a_selection_filter_string}.",
          logging.ERROR
        )
        raise exception.IllegalArgumentError(f"Invalid selection filter: {a_selection_filter_string}.")

  # <editor-fold desc="Validation methods">
  def _validate_input_based_on_regex(self, an_input: str, a_regex: str) -> bool:
    """Validates a given input string for the selection filter based on a given regex.

    Args:
      an_input: The input to validate
      a_regex: The regex to check against

    Returns:
      True if the protein name is valid, Otherwise: false

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_regex` is an empty string.
    """
    # <editor-fold desc="Checks">
    if an_input is None:
      default_logging.append_to_log_file(logger, "an_input is None.", logging.ERROR)
      raise exception.NoneValueError("an_input is None.")
    if a_regex is None:
      default_logging.append_to_log_file(logger, "a_regex is None.", logging.ERROR)
      raise exception.NoneValueError("a_regex is None.")
    if a_regex == "":
      default_logging.append_to_log_file(logger, "a_regex is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_regex is an empty string.")
    # </editor-fold>
    if re.match(a_regex, an_input) or an_input == "":
      return True
    return False

  def _validate_residue_name(self, a_residue_name: str) -> bool:
    """Validates the residue name for the selection filter string part.

    Args:
      a_residue_name: A residue name that can contain multiple residues separated by +

    Returns:
      True if the residue name is valid, Otherwise: false

    Raises:
      exception.NoneValueError: If `a_residue_name` is None.
    """
    # <editor-fold desc="Checks">
    if a_residue_name is None:
      default_logging.append_to_log_file(logger, "a_residue_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_residue_name is None.")
    # </editor-fold>
    if self.residue_name_part_is_incomplete():
      return True  # TODO: Should be changed to False if try-except is used in protein_model.ProteinModel
    for tmp_residue_name in a_residue_name.split("+"):
      if tmp_residue_name in model_definitions.ModelDefinitions.AMINO_ACID_CODE.keys() or tmp_residue_name == "":
        return True
    return False
  # </editor-fold>

  # <editor-fold desc="Check parts">
  def chain_part_is_empty(self) -> bool:
    """Checks if the chain part of the selection filter string is empty."""
    if self.chain == "":
      return True
    return False

  def residue_name_part_is_empty(self) -> bool:
    """Checks if the residue name part of the selection filter string is empty."""
    if self.residue_name == "":
      return True
    return False

  def residue_name_part_is_incomplete(self) -> bool:
    """Checks if the residue name part is incomplete."""
    tmp_residue_names = self.residue_name.split("+")
    for tmp_residue_name in tmp_residue_names:
      if len(tmp_residue_name) != 3:
        return True
    return False

  def residue_number_part_is_empty(self) -> bool:
    """Checks if the residue number part of the selection filter string is empty."""
    if self.residue_number == "":
      return True
    return False

  def atom_part_is_empty(self) -> bool:
    """Checks if the atom part of the selection filter string is empty."""
    if self.atom == "":
      return True
    return False
  # </editor-fold>

  def multiple_chains_are_selected(self) -> bool:
    """Checks if the selection filter string contains a selection of multiple chains.

    Returns:
      True if the chain part of the filtered selection contains more than one chain that is delimited with a +, Otherwise: false
    """
    if re.match(r'^[A-Z](\+[A-Z])*\+?$', self.chain):
      if self.chain.find("+") != -1:
        return True
    return False

  def get_single_chain_identifiers(self) -> list[str]:
    """Returns a list of all chains of the chain part."""
    if self.chain_part_is_empty():
      return []
    tmp_chains: list[str] = self.chain.split("+")
    if tmp_chains[-1] == "":
      tmp_chains.remove("")
    return tmp_chains

  def get_single_residue_names(self) -> list[str]:
    """Returns a list of all residues of the residue name part."""
    if self.residue_name_part_is_empty():
      return []
    tmp_residue_names: list[str] = self.residue_name.split("+")
    if tmp_residue_names[-1] == "":
      tmp_residue_names.remove("")
    return tmp_residue_names

  def get_single_residue_numbers(self) -> tuple[list[str], list[str]]:
    """Returns a tuple of lists of all residues of the residue name part."""
    if self.residue_number_part_is_empty():
      return [], []
    tmp_residue_no_plus: list[str] = []
    tmp_residue_no_hyphen: list[str] = []
    if self.residue_number.find("+") != -1:
      tmp_residue_no_plus: list[str] = self.residue_number.split("+")
      if tmp_residue_no_plus[-1] == "":
        tmp_residue_no_plus.remove("")
    if self.residue_number.find("-") != -1:
      tmp_residue_no_hyphen: list[str] = self.residue_number.split("-")
      if tmp_residue_no_hyphen[-1] == "":
        tmp_residue_no_hyphen.remove("")
    return tmp_residue_no_plus, tmp_residue_no_hyphen
