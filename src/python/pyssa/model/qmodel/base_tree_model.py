import logging
from typing import Optional, Any

from PyQt6 import QtGui
from PyQt6 import QtCore

from pyssa.model.preference import model_definitions
from pyssa.model.util import enums
from pyssa.model.util import exception
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)

__docformat__ = "google"


class BaseTreeModel(QtGui.QStandardItemModel):
  """Base class for tree model classes"""

  def __init__(self) -> None:
    """Constructor."""
    super().__init__()
    self.root_node: Optional[QtGui.QStandardItem] = None

  def create_root_node(self) -> None:
    """Creates the root node of the tree model."""
    if self.root_node is not None:
      default_logging.append_to_log_file(logger, "Root node already exists. Nothing to do.", logging.WARNING)
      return
    self.root_node = self.invisibleRootItem()

  def is_empty(self) -> bool:
    """Checks if the model is empty."""
    return True if self.rowCount() == 0 else False

  def add_node(
    self,
    a_parent_node: QtGui.QStandardItem,
    an_item_name: str,
    an_item_type_value: "model_definitions.TypesEnum",
    an_item_object_value: Optional[object] = None
  ) -> QtGui.QStandardItem:
    """Adds a node to the tree model.

    Args:
      a_parent_node: A node where the new node will be added.
      an_item_name: A name for the new node.
      an_item_type_value: A value for the "type" role of the new node.
      an_item_object_value: A value for the "object" role of the new node.

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `an_item_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_parent_node is None:
      default_logging.append_to_log_file(logger, "a_parent_node is None.", logging.ERROR)
      raise exception.NoneValueError("a_parent_node is None.")
    if an_item_name is None:
      default_logging.append_to_log_file(logger, "an_item_name is None.", logging.ERROR)
      raise exception.NoneValueError("an_item_name is None.")
    if an_item_name == "":
      default_logging.append_to_log_file(logger, "an_item_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("an_item_name is an empty string.")
    if an_item_type_value is None:
      default_logging.append_to_log_file(logger, "an_item_type_value is None.", logging.ERROR)
      raise exception.NoneValueError("an_item_type_value is None.")
    if an_item_object_value is None:
      default_logging.append_to_log_file(logger, "an_item_object_value is None.", logging.ERROR)
      raise exception.NoneValueError("an_item_object_value is None.")
    # </editor-fold>
    tmp_item = QtGui.QStandardItem(an_item_name)
    if an_item_object_value is not None:
      tmp_item.setData(an_item_object_value, enums.ModelEnum.OBJECT_ROLE)
    tmp_item.setData(an_item_type_value, enums.ModelEnum.TYPE_ROLE)
    a_parent_node.appendRow(tmp_item)
    return tmp_item

  def remove_node(self, a_model_index: QtCore.QModelIndex) -> None:
    """Removes a node for a given model index.

    Args:
      a_model_index: The index of the item to be removed.

    Raises:
      exception.NoneValueError: If `a_model_index` is None.

    """
    # <editor-fold desc="Checks">
    if a_model_index is None:
      default_logging.append_to_log_file(logger, "a_model_index is None.", logging.ERROR)
      raise exception.NoneValueError("a_model_index is None.")
    # </editor-fold>
    tmp_item = self.itemFromIndex(a_model_index)
    self.removeRow(tmp_item.row())

  def get_index(self, a_row: int, a_parent: Optional[QtCore.QModelIndex] = None) -> QtCore.QModelIndex:
    """Gets an index based on the given row and optionally parent index.

    Args:
      a_row: The row to get the index for
      a_parent: The parent index to set the row in context (Default: None)

    Raises:
      exception.NoneValueError: If `a_row` is None.
      exception.IllegalArgumentError: If `a_row` has a value less than zero.
    """
    # <editor-fold desc="Checks">
    if a_row is None:
      default_logging.append_to_log_file(logger, "a_row is None.", logging.ERROR)
      raise exception.NoneValueError("a_row is None.")
    if a_row < 0:
      default_logging.append_to_log_file(logger, "a_row has a value less than zero.")
      raise exception.IllegalArgumentError("a_row has a value less than zero")
    # </editor-fold>
    if a_parent is None:
      return self.index(a_row, 0)
    return self.index(a_row, 0, a_parent)

  def get_root_node_as_index(self) -> QtCore.QModelIndex:
    """Returns the root node as a QModelIndex."""
    return self.indexFromItem(self.root_node)

  def get_display_data_of_index(self, an_index: QtCore.QModelIndex) -> str:
    """Gets the display data of the index.

    Args:
      an_index: The index to get the type of

    Raises:
      exception.NoneValueError: If `an_index` is None.

    """
    # <editor-fold desc="Checks">
    if an_index is None:
      default_logging.append_to_log_file(logger, "an_index is None.", logging.ERROR)
      raise exception.NoneValueError("an_index is None.")
    # </editor-fold>
    return an_index.data(QtCore.Qt.ItemDataRole.DisplayRole)

  def get_type_data_of_index(self, an_index: QtCore.QModelIndex) -> str:
    """Gets the type of the index.

    Args:
      an_index: The index to get the type of

    Raises:
      exception.NoneValueError: If `an_index` is None.

    """
    # <editor-fold desc="Checks">
    if an_index is None:
      default_logging.append_to_log_file(logger, "an_index is None.", logging.ERROR)
      raise exception.NoneValueError("an_index is None.")
    # </editor-fold>
    return an_index.data(enums.ModelEnum.TYPE_ROLE)

  def get_object_data_of_index(self, an_index: QtCore.QModelIndex) -> Any:
    """Gets the object of the index.

    Args:
      an_index: The index to get the object of

    Raises:
      exception.NoneValueError: If `an_index` is None.

    """
    # <editor-fold desc="Checks">
    if an_index is None:
      default_logging.append_to_log_file(logger, "an_index is None.", logging.ERROR)
      raise exception.NoneValueError("an_index is None.")
    # </editor-fold>
    return an_index.data(enums.ModelEnum.OBJECT_ROLE)

  def create_row_number_iterator(self, an_index: Optional[QtCore.QModelIndex] = None) -> range:
    """Creates a list of row numbers for the given index.

    Args:
      an_index: The index of which the children rows should be used (Default: None)
    """
    if an_index is None:
      return range(self.rowCount())
    return range(self.rowCount(an_index))
