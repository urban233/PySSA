#
# TEA - Task Event-based Async library for Python
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/urban233/TEA>
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
"""Module for the TaskSignal class."""
import pathlib

from PyQt5.QtCore import pyqtSignal
from PyQt5 import QtCore

from tea.util import tea_logging

logger = tea_logging.setup_logger(__file__, pathlib.Path("../../logs/"))
__docformat__ = "google"


class TaskSignals(QtCore.QObject):
  """Signals related to the tasks."""

  started = pyqtSignal()
  """
  A signal indicating that the action has started with its execution.
  """
  finished = pyqtSignal(str)
  """
  A signal indicating that the action is finished whether the task was
  successful or not.
  """
  error = pyqtSignal(tuple)
  """
  A signal indicating that the action caught an exception.
  """
  result = pyqtSignal(tuple)
  """
  A signal containing the results of the action in the form (a_success_flag,
  the_result_object)
  """
