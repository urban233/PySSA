#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/urban233/PySSA>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""Module for the sequence list context menu."""
import logging
from typing import Callable

from PyQt5 import QtWidgets

from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class SequenceListContextMenu:
  """A custom context menu wrapper for the sequence list view of the main dialog."""

  def __init__(self) -> None:
    """Constructor."""
    self._context_menu = QtWidgets.QMenu()
    self._rename_sequence_action = QtWidgets.QAction("Rename Selected Protein")
    self._help_action = QtWidgets.QAction("Get Help")

    self._context_menu.addAction(self._rename_sequence_action)
    self._context_menu.addAction(self._help_action)

    self._rename_sequence_action.triggered.connect(self._generic_action)
    self._help_action.triggered.connect(self._generic_action)

  # <editor-fold desc="Connect actions from outside">
  def connect_rename_sequence_action(
      self, the_function_to_connect: Callable
  ) -> None:
    """Connects the `_rename_sequence_action` triggered signal to a given callable.

    Args:
        the_function_to_connect: A callable object that will be connected to the triggered signal of the `_rename_sequence_action` attribute.

    Raises:
        exception.IllegalArgumentError: If the_function_to_connect is None.
    """
    # <editor-fold desc="Checks">
    if the_function_to_connect is None:
      logger.error("the_function_to_connect is None.")
      raise exception.IllegalArgumentError("the_function_to_connect is None.")

    # </editor-fold>

    self._rename_sequence_action.triggered.disconnect()
    self._rename_sequence_action.triggered.connect(the_function_to_connect)

  def connect_help_action(self, the_function_to_connect: Callable) -> None:
    """Connects the `_help_action` triggered signal to a given callable.

    Args:
        the_function_to_connect: A callable object that will be connected to the triggered signal of the `_help_action` attribute.

    Raises:
        exception.IllegalArgumentError: If the_function_to_connect is None.
    """
    # <editor-fold desc="Checks">
    if the_function_to_connect is None:
      logger.error("the_function_to_connect is None.")
      raise exception.IllegalArgumentError("the_function_to_connect is None.")

    # </editor-fold>

    self._help_action.triggered.disconnect()
    self._help_action.triggered.connect(the_function_to_connect)

  # </editor-fold>

  def _generic_action(self) -> None:
    raise NotImplementedError(
        "The called action is not connected to a function!"
    )

  def get_context_menu(self, the_selected_indexes) -> QtWidgets.QMenu:
    """Gets the context menu for the current situation.

    Args:
        the_selected_indexes (list): A list of selected indexes.

    Returns:
        QMenu: The context menu.
    """
    if len(the_selected_indexes) == 1:
      level = 0
      index = the_selected_indexes[0]
      while index.parent().isValid():
        index = index.parent()
        level += 1
    elif len(the_selected_indexes) > 1:
      # More than one sequence is selected
      self._help_action.setVisible(False)
      self._rename_sequence_action.setVisible(False)
      return self._context_menu
    else:
      # No sequences are selected
      self._help_action.setVisible(False)
      self._rename_sequence_action.setVisible(False)
      return self._context_menu
    # One sequence is selected
    self._help_action.setVisible(False)
    self._rename_sequence_action.setVisible(True)

    return self._context_menu
