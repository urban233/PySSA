import logging
import pathlib
from typing import Optional

from PyQt6 import QtGui
from PyQt6 import QtCore
from Bio import SeqIO

from pyssa.model.central_objects import sequence
from pyssa.model.preference import model_definitions
from pyssa.model.util import exception
from pyssa.model.qmodel import base_tree_model
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)

__docformat__ = "google"


class SequenceModel(base_tree_model.BaseTreeModel):
  """Class for storing sequences of the current opened project."""

  def __init__(self):
    """Constructor."""
    super().__init__()

  @staticmethod
  def from_a_list_of_sequences(a_list_of_sequences: list["sequence.Sequence"]) -> "SequenceModel":
    """Alternative constructor.

    Notes:
      - If the list of sequence is empty, then an empty model with a root node will be created.

    Raises:
      exception.NoneValueError: If `a_list_of_sequences` is None.

    """
    # <editor-fold desc="Checks">
    if a_list_of_sequences is None:
      default_logging.append_to_log_file(logger, "a_list_of_sequences is None.", logging.ERROR)
      raise exception.NoneValueError("a_list_of_sequences is None.")
    # </editor-fold>
    tmp_sequence_model = SequenceModel()
    tmp_sequence_model.create_root_node()
    if len(a_list_of_sequences) == 0:
      default_logging.append_to_log_file(
        logger,
        "Initialized sequence model with an empty list.",
        logging.INFO
      )
      return tmp_sequence_model  # Return empty model if no sequences are given
    try:
      for tmp_sequence in a_list_of_sequences:
        tmp_sequence_model.add_node(
          a_parent_node=tmp_sequence_model.root_node,
          an_item_name=tmp_sequence.name,
          an_item_type_value=model_definitions.TypesEnum.SEQUENCE_TYPE,
          an_item_object_value=tmp_sequence
        )
    except Exception as e:
      default_logging.append_to_log_file(
        logger,
        "A problem occurred while adding the sequence nodes to the model. Following error was raised:\n" + str(e),
        logging.ERROR
      )
    finally:
      # Returning the model ensures at least some of the sequences are in the model
      return tmp_sequence_model

  def add_sequence(self, a_sequence: "sequence.Sequence") -> QtGui.QStandardItem:
    """Adds a sequence to the model.

      Args:
        a_sequence: The sequence to add

      Raises:
        exception.NoneValueError: If `a_sequence` is None.
      """
    # <editor-fold desc="Checks">
    if a_sequence is None:
      default_logging.append_to_log_file(logger, "a_sequence is None.", logging.ERROR)
      raise exception.NoneValueError("a_sequence is None.")
    # </editor-fold>
    tmp_sequence_node = self.add_node(
      a_parent_node=self.root_node,
      an_item_name=a_sequence.name,
      an_item_type_value=model_definitions.TypesEnum.SEQUENCE_TYPE,
      an_item_object_value=a_sequence
    )
    return tmp_sequence_node
