import threading

from PyQt6 import QtWidgets
from PyQt6 import QtCore
from PyQt6.QtCore import Qt

from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)

__docformat__ = "google"


def is_main_thread() -> bool:
  """Check if the current thread is the main thread.

  Returns:
      The boolean True if the current thread is the main thread, False otherwise.
  """
  if threading.current_thread() == threading.main_thread():
    logger.info("Running in main thread.")
    return True
  logger.info("Running in separate thread.")
  return False


def select_matching_string_in_q_list_view(
        a_string_to_match: str,
        a_q_list_view_to_select: QtWidgets.QListView,
        a_q_line_edit: QtWidgets.QLineEdit,
) -> None:
  """Selects a QListView item from a QListView with a string to match and sets it into a q_line_edit.

  Args:
      a_string_to_match (str): The string to match in the QListView.
      a_q_list_view_to_select (QtWidgets.QListView): The QListView to select in.
      a_q_line_edit (QtWidgets.QLineEdit): The QLineEdit widget associated with the QListView.

  Raises:
      exception.IllegalArgumentError: If either `a_string_to_match`, `a_q_list_view_to_select` or `a_q_line_edit` is None.
      exception.NotMainThreadError: If function is called not from the main thread.
  """
  # <editor-fold desc="Checks">
  if a_string_to_match is None:
    logger.error("a_string_to_match is None.")
    raise exception.IllegalArgumentError("a_string_to_match is None.")
  if a_q_list_view_to_select is None:
    logger.error("a_label a_q_list_view_to_select None.")
    raise exception.IllegalArgumentError("a_q_list_view_to_select is None.")
  if a_q_line_edit is None:
    logger.error("a_q_line_edit is None.")
    raise exception.IllegalArgumentError("a_q_line_edit is None.")
  if not is_main_thread():
    raise exception.NotMainThreadError()

  # </editor-fold>

  if a_string_to_match == "":
    a_q_line_edit.clear()
    a_q_list_view_to_select.selectionModel().clearCurrentIndex()
    a_q_list_view_to_select.selectionModel().clearSelection()
    return

  tmp_exactly_matched_items = a_q_list_view_to_select.model().findItems(
    a_string_to_match,
    Qt.MatchFlag.MatchExactly,
  )
  if tmp_exactly_matched_items == 1:
    tmp_match_items = tmp_exactly_matched_items
  else:
    tmp_match_items = a_q_list_view_to_select.model().findItems(
      a_string_to_match,
      QtCore.Qt.MatchFlag.MatchContains,
    )
  if len(tmp_match_items) == 0:
    return
  # Sets the selection
  a_q_list_view_to_select.selectionModel().setCurrentIndex(
    tmp_match_items[0].index(),
    QtCore.QItemSelectionModel.SelectionFlag.Clear | QtCore.QItemSelectionModel.SelectionFlag.Select,
    )
  a_q_line_edit.setText(tmp_match_items[0].data(Qt.ItemDataRole.DisplayRole))
