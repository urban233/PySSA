import logging
from typing import Optional, Any

from src.pyssa.gui.qt import QtGui
from src.pyssa.gui.qt import QtCore
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import enums

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)

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
      logger.warning("Root node already exists. Nothing to do.")
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
      an_item_object_value: Optional[object] = None,
  ) -> QtGui.QStandardItem:
    """Adds a new child node to a parent node in a tree model structure.

    This method is designed to simplify the process of adding nodes to a
    QStandardItemModel tree. It validates the input arguments, creates a new
    `QtGui.QStandardItem`, sets its data roles, and appends it as a child to
    the specified parent node.

    Args:
      a_parent_node: The parent node to which the new child node will be added. This must be
        a valid `QStandardItem` instance, and cannot be `None`.

      an_item_name: The name of the new node. This will be displayed in the tree view and
        used as the item's label. The name must be a non-empty string.

      an_item_type_value: A value corresponding to the "type" role (`TYPE_ROLE`) of the new node.
        This is typically used to categorize or identify the type of the item
        within the model.

      an_item_object_value (Optional): An optional value for the "object" role (`OBJECT_ROLE`) of the new node.
        This can be any Python object associated with the node (e.g., metadata,
        references, or user-defined data). If not provided, this role will be
        left unset.

    Returns:
      The newly created and appended `QStandardItem` instance. This can be
      used to modify or retrieve the newly added node after creation.

    Raises:
      NoneValueError: If any of the arguments are None.
      IllegalArgumentError: If `an_item_name` is an empty string.

    Note:
      - This method assumes that the caller has already set up the parent node and
        is managing the tree model externally. It does not handle the creation or
        initialization of the model itself.

      - The `TYPE_ROLE` and `OBJECT_ROLE` are defined in `enums.ModelEnum`. Ensure
        that these roles are correctly implemented in your project to avoid runtime
        errors.

      - To extend the functionality of this method (e.g., adding additional data
        roles), modify the `setData` calls within the method.
    """
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
    tmp_item = self.itemFromIndex(a_model_index)
    self.removeRow(tmp_item.row())

  def get_index(
      self, a_row: int, a_parent: Optional[QtCore.QModelIndex] = None
  ) -> QtCore.QModelIndex:
    """Gets an index based on the given row and optionally parent index.

    Args:
      a_row: The row to get the index for
      a_parent: The parent index to set the row in context (Default: None)

    Raises:
      exception.NoneValueError: If `a_row` is None.
      exception.IllegalArgumentError: If `a_row` has a value less than zero.
    """
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
    return an_index.data(QtCore.Qt.ItemDataRole.DisplayRole)

  def get_type_data_of_index(self, an_index: QtCore.QModelIndex) -> str:
    """Gets the type of the index.

    Args:
      an_index: The index to get the type of

    Raises:
      exception.NoneValueError: If `an_index` is None.

    """
    return an_index.data(enums.ModelEnum.TYPE_ROLE)

  def get_object_data_of_index(self, an_index: QtCore.QModelIndex) -> Any:
    """Gets the object of the index.

    Args:
      an_index: The index to get the object of

    Raises:
      exception.NoneValueError: If `an_index` is None.

    """
    return an_index.data(enums.ModelEnum.OBJECT_ROLE)

  def create_row_number_iterator(
      self, an_index: Optional[QtCore.QModelIndex] = None
  ) -> range:
    """Creates a list of row numbers for the given index.

    Args:
      an_index: The index of which the children rows should be used (Default: None)
    """
    if an_index is None:
      return range(self.rowCount())
    return range(self.rowCount(an_index))
