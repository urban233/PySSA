#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/zielesny/PySSA>
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
"""Module for the help manager."""
import logging
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtCore import QUrl
from src.pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class HelpManager:
  """Class for the HelpManager."""

  def __init__(self) -> None:
    """Constructor."""
    super().__init__()
    self._browser = QWebEngineView()
    self.current_url: QUrl = QUrl("https://www.youtube.com/watch?v=NIwRJjk-LzI")
    self._browser.setUrl(self.current_url)

  def change_url(self, an_url: str) -> None:
    """Changes the URL of the browser window."""
    self.current_url = QUrl(an_url)
    self._browser.setUrl(QUrl(an_url))
    self._browser.show()
