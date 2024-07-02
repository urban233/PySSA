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
"""Module for the protein tree context menu."""
import logging
from typing import Callable

from PyQt5 import QtWidgets

from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class ProteinTreeContextMenu:
  """A custom context menu wrapper for the protein tree view of the main view window."""

  def __init__(self) -> None:
    """Constructor."""
    self._context_menu = QtWidgets.QMenu()
    self._expand_protein_action = QtWidgets.QAction("Expand Protein")
    self._collapse_protein_action = QtWidgets.QAction("Collapse Protein")
    self._clean_protein_action = QtWidgets.QAction("Clean Selected Protein")
    self._rename_protein_action = QtWidgets.QAction("Rename Selected Protein")
    self._show_sequence_action = QtWidgets.QAction("Show Sequence")
    self._help_action = QtWidgets.QAction("Get Help")

    self._context_menu.addAction(self._expand_protein_action)
    self._context_menu.addAction(self._collapse_protein_action)
    self._context_menu.addSeparator()
    self._context_menu.addAction(self._clean_protein_action)
    self._context_menu.addAction(self._rename_protein_action)
    self._context_menu.addAction(self._show_sequence_action)
    self._context_menu.addAction(self._help_action)

    self._expand_protein_action.triggered.connect(self._generic_action)
    self._collapse_protein_action.triggered.connect(self._generic_action)
    self._clean_protein_action.triggered.connect(self._generic_action)
    self._rename_protein_action.triggered.connect(self._generic_action)
    self._show_sequence_action.triggered.connect(self._generic_action)
    self._help_action.triggered.connect(self._generic_action)

  # <editor-fold desc="Connect actions from outside">
  def connect_expand_protein_action(self, the_function_to_connect: Callable):
    """Connects the `_expand_protein_action` triggered signal to a given callable.

    Args:
        the_function_to_connect: A callable object that will be connected to the triggered signal of the `_expand_protein_action` attribute.

    Raises:
        exception.IllegalArgumentError: If the_function_to_connect is None.
    """
    # <editor-fold desc="Checks">
    if the_function_to_connect is None:
      logger.error("the_function_to_connect is None.")
      raise exception.IllegalArgumentError("the_function_to_connect is None.")

    # </editor-fold>
    self._expand_protein_action.triggered.disconnect()
    self._expand_protein_action.triggered.connect(the_function_to_connect)

  def connect_collapse_protein_action(self, the_function_to_connect: Callable):
    """Connects the `_collapse_protein_action` triggered signal to a given callable.

    Args:
        the_function_to_connect: A callable object that will be connected to the triggered signal of the `_collapse_protein_action` attribute.

    Raises:
        exception.IllegalArgumentError: If the_function_to_connect is None.
    """
    # <editor-fold desc="Checks">
    if the_function_to_connect is None:
      logger.error("the_function_to_connect is None.")
      raise exception.IllegalArgumentError("the_function_to_connect is None.")

    # </editor-fold>

    self._collapse_protein_action.triggered.disconnect()
    self._collapse_protein_action.triggered.connect(the_function_to_connect)

  def connect_clean_protein_action(self, the_function_to_connect: Callable):
    """Connects the `_clean_protein_action` triggered signal to a given callable.

    Args:
        the_function_to_connect: A callable object that will be connected to the triggered signal of the `_clean_protein_action` attribute.

    Raises:
        exception.IllegalArgumentError: If the_function_to_connect is None.
    """
    # <editor-fold desc="Checks">
    if the_function_to_connect is None:
      logger.error("the_function_to_connect is None.")
      raise exception.IllegalArgumentError("the_function_to_connect is None.")

    # </editor-fold>

    self._clean_protein_action.triggered.disconnect()
    self._clean_protein_action.triggered.connect(the_function_to_connect)

  def connect_rename_protein_action(self, the_function_to_connect: Callable):
    """Connects the `_rename_protein_action` triggered signal to a given callable.

    Args:
        the_function_to_connect: A callable object that will be connected to the triggered signal of the `_rename_protein_action` attribute.

    Raises:
        exception.IllegalArgumentError: If the_function_to_connect is None.
    """
    # <editor-fold desc="Checks">
    if the_function_to_connect is None:
      logger.error("the_function_to_connect is None.")
      raise exception.IllegalArgumentError("the_function_to_connect is None.")

    # </editor-fold>

    self._rename_protein_action.triggered.disconnect()
    self._rename_protein_action.triggered.connect(the_function_to_connect)

  def connect_show_sequence_action(self, the_function_to_connect: Callable):
    """Connects the `_show_sequence_action` triggered signal to a given callable.

    Args:
        the_function_to_connect: A callable object that will be connected to the triggered signal of the `_show_sequence_action` attribute.

    Raises:
        exception.IllegalArgumentError: If the_function_to_connect is None.
    """
    # <editor-fold desc="Checks">
    if the_function_to_connect is None:
      logger.error("the_function_to_connect is None.")
      raise exception.IllegalArgumentError("the_function_to_connect is None.")

    # </editor-fold>

    self._show_sequence_action.triggered.disconnect()
    self._show_sequence_action.triggered.connect(the_function_to_connect)

  def connect_help_action(self, the_function_to_connect: Callable):
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

  def get_context_menu(
      self,
      the_selected_indexes,
      the_type: str,
      is_protein_in_any_pair_flag: bool,
      is_protein_in_session_flag: bool,
      is_protein_expanded: bool,
      is_database_thread_running: bool,
  ) -> QtWidgets.QMenu:
    """Gets the context menu for the current situation.

    Args:
        the_selected_indexes (list): List of QModelIndex objects representing the selected indexes.
        the_type (str): Type of the selected object. Can be "chain" or "scene".
        is_protein_in_any_pair_flag (bool): Flag indicating if the selected protein is in any protein pair.
        is_protein_in_session_flag (bool): Flag indicating if the selected protein is in the session.
        is_protein_expanded (bool): Flag indicating if the selected protein is expanded.
        is_database_thread_running (bool): Flag indicating if the database thread is running.

    Returns:
        QMenu: The context menu that should be displayed.

    Raises:
        exception.IllegalArgumentError: If any argument except the_selected_indexes is None.
    """
    # <editor-fold desc="Checks">
    if the_type is None:
      logger.error("the_type is None.")
      raise exception.IllegalArgumentError("the_type is None.")
    if is_protein_in_any_pair_flag is None:
      logger.error("is_protein_in_any_pair_flag is None.")
      raise exception.IllegalArgumentError(
          "is_protein_in_any_pair_flag is None."
      )
    if is_protein_in_session_flag is None:
      logger.error("is_protein_in_session_flag is None.")
      raise exception.IllegalArgumentError(
          "is_protein_in_session_flag is None."
      )
    if is_protein_expanded is None:
      logger.error("is_protein_expanded is None.")
      raise exception.IllegalArgumentError("is_protein_expanded is None.")
    if is_database_thread_running is None:
      logger.error("is_database_thread_running is None.")
      raise exception.IllegalArgumentError("is_database_thread_running is None.")

    if len(the_selected_indexes) > 0:
      level = 0
      index = the_selected_indexes[0]
      while index.parent().isValid():
        index = index.parent()
        level += 1
    else:
      self._help_action.setVisible(True)
      self._expand_protein_action.setVisible(False)
      self._collapse_protein_action.setVisible(False)
      self._clean_protein_action.setVisible(False)
      self._rename_protein_action.setVisible(False)
      self._show_sequence_action.setVisible(False)
      return self._context_menu

    # </editor-fold>

    # Manage of context menu
    self._help_action.setVisible(False)
    if level == 0:
      # A protein is selected
      if is_protein_expanded is True:
        self._expand_protein_action.setVisible(False)
        self._collapse_protein_action.setVisible(True)
      else:
        self._expand_protein_action.setVisible(True)
        self._collapse_protein_action.setVisible(False)
      self._clean_protein_action.setVisible(True)
      if is_database_thread_running is True:
        self._clean_protein_action.setEnabled(False)
      else:
        self._clean_protein_action.setEnabled(True)
      self._rename_protein_action.setVisible(False)
      self._show_sequence_action.setVisible(False)

      # <editor-fold desc="Check if protein is in any protein pair">
      if is_protein_in_any_pair_flag:
        self._rename_protein_action.setEnabled(False)
      else:
        self._rename_protein_action.setEnabled(True)

      if is_protein_in_session_flag:
        self._rename_protein_action.setEnabled(True)
      else:
        self._rename_protein_action.setEnabled(False)

      # </editor-fold>

    elif level == 1:
      # A header is selected (Scenes or Chains)
      self._expand_protein_action.setVisible(False)
      self._collapse_protein_action.setVisible(False)
      self._clean_protein_action.setVisible(False)
      self._rename_protein_action.setVisible(False)
      self._show_sequence_action.setVisible(False)
    elif level == 2:
      # A chain or scene is selected
      if the_type == "chain":
        self._expand_protein_action.setVisible(False)
        self._collapse_protein_action.setVisible(False)
        self._clean_protein_action.setVisible(False)
        self._rename_protein_action.setVisible(False)
        self._show_sequence_action.setVisible(True)
      elif the_type == "scene":
        # the code below could change if other actions are added!
        self._expand_protein_action.setVisible(False)
        self._collapse_protein_action.setVisible(False)
        self._clean_protein_action.setVisible(False)
        self._rename_protein_action.setVisible(False)
        self._show_sequence_action.setVisible(False)
      else:
        self._expand_protein_action.setVisible(False)
        self._collapse_protein_action.setVisible(False)
        self._clean_protein_action.setVisible(False)
        self._rename_protein_action.setVisible(False)
        self._show_sequence_action.setVisible(False)
    return self._context_menu
