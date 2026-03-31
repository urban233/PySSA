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
"""HelpPanelController class for the PySSA frontend application.

Authors: Martin Urban, Hannah Kullik

Version: 1.4.1
"""
import logging

from src.pyssa.gui.ui.views import help_panel
from src.pyssa.gui import main_window
from src.pyssa.gui.qt import QtCore, QtWidgets
from src.pyssa.logging_pyssa import log_levels, log_handlers
from src.pyssa.controller.help_dialog_registry import help_registry
from src.pyssa.controller import help_manager

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class HelpPanelController:
  """Controller class for the HelpPanel."""

  def __init__(
          self,
          the_main_window: "main_window.MainWindow",
          a_help_panel: "help_panel.HelpPanel",
          the_help_manager: "help_manager.HelpManager" = None,
  ):
    """Constructor.

    Args:
        the_main_window: Reference to the main window
        a_help_panel: The help panel view
        the_help_manager: The help manager for opening help dialogs (optional, will create if not provided)
    """
    # <editor-fold desc="Instance attributes">
    self._main_window = the_main_window
    self._panel = a_help_panel
    self._help_manager = the_help_manager if the_help_manager is not None else help_manager.HelpManager()
    self._registry = help_registry
    # </editor-fold>
    self._connect_signals()
    # self._populate_help_panel()  # Disabled - only show hover help

  def _connect_signals(self) -> None:
    """Connects panel signals to their respective handlers."""
    self._panel.panelClosed.connect(self.__slot_close_panel)
    self._panel.helpDialogRequested.connect(self.__slot_open_help_dialog)

  def _populate_help_panel(self) -> None:
    """Populates the help panel with available help topics from the registry."""
    categories = {}
    for dialog in self._registry.get_all_dialogs():
      if dialog.category not in categories:
        categories[dialog.category] = []
      categories[dialog.category].append(dialog)

    self._panel.populate_help_topics(categories)

  def __slot_close_panel(self) -> None:
    """Handler for panel close event."""
    self._main_window.tool_window_layout.set_right_panel_hidden(True)

  def __slot_open_help_dialog(self, dialog_id: str) -> None:
    """Handler for opening a help dialog.

    Args:
        dialog_id: The unique identifier of the help dialog to open
    """
    try:
      self._help_manager.open_help_dialog(dialog_id)
      logger.info(f"Opened help dialog: {dialog_id}")
    except ValueError as e:
      logger.error(f"Failed to open help dialog: {e}")
      QtWidgets.QMessageBox.warning(
        self._main_window,
        "Help Not Available",
        f"The requested help topic could not be found: {dialog_id}"
      )
