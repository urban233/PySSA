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
"""PySSAObjectsPanelController class for the PySSA frontend application.

Authors: Martin Urban, Hannah Kullik

Version: 1.4.0
"""
import logging

from src.pyssa.gui.ui.views import pyssa_objects_panel, help_panel
from src.pyssa.gui import user_pymol, app_state, main_window
from src.pyssa.gui.qt import QtCore, QtWidgets
from src.pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class HelpPanelController:
  """Controller class for the PySSAObjectsPanel."""

  def __init__(
          self,
          the_main_window: "main_window.MainWindow",
          a_help_panel: "help_panel.HelpPanel",
  ):
    """Constructor."""
    # <editor-fold desc="Instance attributes">
    self._main_window = the_main_window
    self._panel = a_help_panel
    # </editor-fold>
    self._connect_signals()

  def _connect_signals(self) -> None:
    # self._panel.import_seq_action.triggered.connect()
    self._panel.panelClosed.connect(self.__slot_close_panel)

  def __slot_close_panel(self):
    self._main_window.tool_window_layout.set_right_panel_hidden(True)
