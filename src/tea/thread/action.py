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
"""Module for the Action class."""
import pathlib
import uuid
from typing import Optional, Callable, Any

from PyQt5 import QtCore

from src.tea.thread import task_signals
from src.tea.util import tea_exception
from src.tea.util import tea_logging

logger = tea_logging.setup_logger(__file__, pathlib.Path("../../logs/"))
__docformat__ = "google"


class Action(QtCore.QRunnable):
  """A container for running an asynchronous operation."""

  # <editor-fold desc="Class attributes">
  id: uuid.UUID  # noqa: A003
  """
  An ID for the action instance.
  """
  _signals: "task_signals.TaskSignals"
  """
  The signals used to communicate with the task instance.
  """
  _target: Callable
  """
  The target function to execute.
  """
  _args: tuple
  """
  The arguments for the target function.
  """
  _lock: Optional[QtCore.QMutex]
  """
  A lock used to indicate if the action is executable.
  """
  is_runnable: bool
  """
  A flag used to indicate if the target function can be run.
  """

  # </editor-fold>

  def __init__(self, a_target: Callable, args: tuple = ()) -> None:
    """Constructor.

    Args:
      a_target (Callable): The function that should get executed.
      args (tuple): The arguments for the function.

    Raises:
      tea_exception.IllegalArgumentError: If `a_target` or an argument is None.
    """
    # <editor-fold desc="Checks">
    if a_target is None:
      logger.error("a_target is None")
      raise tea_exception.IllegalArgumentError("a_target is None")
    if args is None:
      logger.error("args is None")
      raise tea_exception.IllegalArgumentError("args is None")

    # </editor-fold>

    super().__init__()
    self.id: str = str(uuid.uuid4())
    self._signals: "task_signals.TaskSignals" = task_signals.TaskSignals()
    self._target: Callable = a_target
    self._args: tuple = args
    self._lock = None
    self.is_runnable: bool = False

  # <editor-fold desc="Public methods">
  def get_finished_signal(
      self,
  ):  # noqa: ANN201; suppressed because type annotation of pyqtsignal is problematic.
    """Gets the `finished` signal of the instance.

    Returns:
        The `finished` signal of the instance.
    """
    return self._signals.finished

  def get_result_signal(
      self,
  ):  # noqa: ANN201; suppressed because type annotation of pyqtsignal is problematic.
    """Gets the `result` signal of the instance.

    Returns:
        The `result` signal of the instance.
    """
    return self._signals.result

  def set_lock(self, a_lock: QtCore.QMutex) -> None:
    """Sets the `_lock` of the instance.

    Args:
      a_lock: The QMutex lock object to be set as the `_lock` attribute of the class.

    Raises:
      tea_exception.IllegalArgumentError: If `a_lock` is None.
    """
    # <editor-fold desc="Checks">
    if a_lock is None:
      logger.error("a_lock is None")
      raise tea_exception.IllegalArgumentError("a_lock is None")

    # </editor-fold>

    self._lock = a_lock

  def connect_started_signal(self, a_function: Callable) -> None:
    """Connects the given function to the started signal.

    Args:
      a_function (Callable): The function to be connected to the 'error' signal.

    Raises:
      tea_exception.IllegalArgumentError: If `a_function` is None.
    """
    # <editor-fold desc="Checks">
    if a_function is None:
      logger.error("a_function is None.")
      raise tea_exception.IllegalArgumentError("a_function is None.")

    # </editor-fold>

    self._signals.started.connect(a_function)

  def connect_error_signal(self, a_function: Callable) -> None:
    """Connects the given function to the `error` signal.

    Args:
      a_function (Callable): The function to be connected to the `error` signal.

    Raises:
      tea_exception.IllegalArgumentError: If `a_function` is None.
    """
    # <editor-fold desc="Checks">
    if a_function is None:
      logger.error("a_function is None.")
      raise tea_exception.IllegalArgumentError("a_function is None.")

    # </editor-fold>

    self._signals.error.connect(a_function)

  def connect_finished_signal(self, a_function: Callable) -> None:
    """Connects a provided function to the 'finished' signal of the object.

    Args:
      a_function (Callable): The function to be connected to the 'finished' signal.

    Raises:
      tea_exception.IllegalArgumentError: If `a_function` is None.
    """
    # <editor-fold desc="Checks">
    if a_function is None:
      logger.error("a_function is None.")
      raise tea_exception.IllegalArgumentError("a_function is None.")

    # </editor-fold>

    self._signals.finished.connect(a_function)

  def connect_result_signal(self, a_function: Callable) -> None:
    """Connects a provided function to the 'result' signal of the object.

    Args:
      a_function (Callable): The function to be connected to the 'result' signal.

    Raises:
      tea_exception.IllegalArgumentError: If `a_function` is None.
    """
    # <editor-fold desc="Checks">
    if a_function is None:
      logger.error("a_function is None.")
      raise tea_exception.IllegalArgumentError("a_function is None.")

    # </editor-fold>

    self._signals.result.connect(a_function)

  def run(self) -> None:
    """Overwrites the run() method of QRunnable and executes the target
    function.
    """
    if self.is_runnable is False:
      raise tea_exception.ActionIsNotRunnableError()

    try:
      with QtCore.QMutexLocker(self._lock):
        if len(self._args) == 0:
          self._signals.started.emit()
          tmp_result: Any = self._target()
        else:
          self._signals.started.emit()
          tmp_result: Any = self._target(*self._args)
    except Exception as e:
      logger.error(e)
      self._signals.finished.emit((str(self.id), False, e))
      #self._signals.error.emit((e,))
    else:
      self._signals.finished.emit((str(self.id), (True, tmp_result)))
      #self._signals.result.emit((True, tmp_result))

  # </editor-fold>
