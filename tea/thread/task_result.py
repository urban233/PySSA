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
"""Module for the TaskResult class."""
import pathlib
import queue
from typing import Optional, Callable, Any

from PyQt5.QtCore import pyqtSignal
from PyQt5 import QtCore

from tea.thread import task, action
from tea.util import tea_logging, tea_exception, tea_enums

logger = tea_logging.setup_logger(__file__, pathlib.Path("../../logs/"))
__docformat__ = "google"


class TaskResult(task.Task):
  """Represents an asynchronous operation that can return a value."""

  # <editor-fold desc="Class attributes">
  t_result: queue.Queue = queue.Queue()
  """Contains the results of a finished task."""
  # </editor-fold>

  def __init__(self):
    """Empty constructor."""
    super().__init__()
    self.t_result: queue.Queue = queue.Queue()
    self._lock = QtCore.QMutex()

  # <editor-fold desc="Private methods">
  def _action_completed_with_result(self, the_action_id: str) -> None:
    """Updates the task status and completion flag when an action is completed.

    Args:
      the_action_id (str): The `id` of an action.

    Raises:
      tea_exception.IllegalArgumentError: If `the_action_id` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if the_action_id is None or the_action_id == "":
      logger.error("the_action_id is either None or an empty string.")
      raise tea_exception.IllegalArgumentError(
          "the_action_id is either None or an empty string."
      )

    # </editor-fold>

    if (
        the_action_id == self.id_of_last_action
        and self.status != tea_enums.TaskStatus.FAILED
    ):
      self.is_completed: bool = True
      self.status: "tea_enums.TaskStatus" = (
          tea_enums.TaskStatus.RAN_TO_COMPLETION
      )
      #self.task_finished.emit((self.id, []))
    elif (
        the_action_id == self.id_of_last_action
        and self.status == tea_enums.TaskStatus.FAILED
    ):
      # Case is handled by method _action_ended_with_error
      return
    else:
      self.is_completed: bool = False
      self.status: "tea_enums.TaskStatus" = tea_enums.TaskStatus.RUNNING

  def _save_single_action_result(
      self,
      an_action_result: tuple[bool, Any],
  ) -> None:
    """Saves a result to a queue and sends a task finished signal.

    Args:
      an_action_result (tuple[bool, Any]): A tuple representing the `action result`. It should be of the form (success: bool, result: Any).

    Raises:
      IllegalArgumentError: If `an_action_result` is either None or an empty tuple.
    """
    # <editor-fold desc="Checks">
    if an_action_result is None or len(an_action_result) == 0:
      logger.error("an_action_result is either None or an empty tuple.")
      raise tea_exception.IllegalArgumentError(
          "an_action_result is either None or an empty tuple.",
      )

    # </editor-fold>

    # TODO: wrap in try-except block
    self.t_result.put(an_action_result)
    if self.is_completed:
      tmp_action_results = []
      while self.t_result.empty() is False:
      #for tmp_action_result in self.t_result.get():
        tmp_action_results.append(self.t_result.get())

      self.task_finished.emit((self.id, tmp_action_results))

  # </editor-fold>

  # <editor-fold desc="Public methods">
  def connect_basic_signals(self, an_action: "action.Action") -> None:
    """Connects the basic signals (started, error and finished) of the given `an_action` object
    with the corresponding member functions of the current object.

    Args:
      an_action: An instance of the `Action` class.

    Raises:
      tea_exception.IllegalArgumentError: If `an_action` is None.
    """
    # <editor-fold desc="Checks">
    if an_action is None:
      logger.error("an_action is None.")
      raise tea_exception.IllegalArgumentError("an_action is None.")

    # </editor-fold>

    an_action.connect_started_signal(self._action_started)
    an_action.connect_error_signal(self._action_ended_with_error)
    an_action.connect_finished_signal(self._action_completed_with_result)

  def connect_result_signal(self, an_action: "action.Action") -> None:
    """Connects the result signal of an action to the internal
    `save_single_action_result` method.

    Args:
      an_action ("action.Action"): The `action` object whose result signal needs to be connected.

    Raises: IllegalArgumentError: If `an_action` is None.
    """
    # <editor-fold desc="Checks">
    if an_action is None:
      logger.error("an_action is None.")
      tea_exception.IllegalArgumentError("an_action is None.")

    # </editor-fold>

    an_action.connect_result_signal(self._save_single_action_result)

  def get_lock(self) -> QtCore.QMutex:
    """Gets the `_lock` of the instance.

    Returns: The `_lock` of the instance.
    """
    return self._lock

  # </editor-fold>

  # <editor-fold desc="Static methods">
  @staticmethod
  def run_action(
      an_action: "action.Action",
      the_threadpool: QtCore.QThreadPool,
      an_await_function: Optional[Callable] = None,
  ) -> "TaskResult":
    """Static constructor for the TaskResult class that takes an action to
    construct the `TaskResult` and starts the task.

    Args:
      an_action (action.Action): The `action` to be executed.
      the_threadpool (QtCore.QThreadPool): The `thread pool` to be used for executing the action.
      an_await_function (Optional[Callable]): The function to be called when a `task` is completed. Defaults to None.

    Returns: A `TaskResult` object representing the running task.

    Raises:
      IllegalArgumentError: If `an_action` or `the_threadpool` is None.
    """
    # <editor-fold desc="Checks">
    if an_action is None:
      logger.error("an_action is None.")
      raise tea_exception.IllegalArgumentError("an_action is None.")
    if the_threadpool is None:
      logger.error("the_threadpool is None.")
      raise tea_exception.IllegalArgumentError("the_threadpool is None.")

    # </editor-fold>

    tmp_task_result: "TaskResult" = TaskResult()
    tmp_task_result.actions.put(an_action)
    tmp_action: "action.Action" = tmp_task_result.actions.get()
    tmp_task_result.connect_basic_signals(tmp_action)
    tmp_task_result.connect_result_signal(tmp_action)
    if an_await_function is not None:
      tmp_task_result.task_finished.connect(an_await_function)
    tmp_action.is_runnable = True
    the_threadpool.start(tmp_action)
    return tmp_task_result

  @staticmethod
  def run_actions(
      the_actions: tuple["action.Action"],
      the_threadpool: QtCore.QThreadPool,
      an_await_function: Optional[Callable] = None,
  ) -> "TaskResult":
    """Static constructor for the TaskResult class that takes actions to
    construct the `TaskResult` and starts the task.

    Args:
      the_actions (tuple["Action"]): `Actions` to be executed.
      the_threadpool (QtCore.QThreadPool): The `thread pool` to be used for executing the actions.
      an_await_function (Optional[Callable]): The function to be called when a `task` is completed. Defaults to None.

    Returns: A `TaskResult` object representing the running task.

    Raises:
      IllegalArgumentError: If `an_action_result` is either None or an empty tuple or `the_threadpool` is None.
    """
    # <editor-fold desc="Checks">
    if the_actions is None or len(the_actions) == 0:
      logger.error("the_actions is either None or an empty tuple.")
      raise tea_exception.IllegalArgumentError(
          "the_actions is either None or an empty tuple.",
      )
    if the_threadpool is None:
      logger.error("the_threadpool is None.")
      raise tea_exception.IllegalArgumentError("the_threadpool is None.")

    # </editor-fold>

    tmp_task_result: "TaskResult" = TaskResult()
    for tmp_action in the_actions:
      tmp_task_result.connect_basic_signals(tmp_action)
      tmp_task_result.connect_result_signal(tmp_action)
      tmp_action.set_lock(tmp_task_result.get_lock())
      tmp_task_result.actions.put(tmp_action)
    if an_await_function is not None:
      tmp_task_result.task_finished.connect(an_await_function)
    while tmp_task_result.actions.empty() is False:
      tmp_action = tmp_task_result.actions.get()
      tmp_action.is_runnable = True
      the_threadpool.start(tmp_action)
    return tmp_task_result

  @staticmethod
  def from_action(an_action: "action.Action", an_await_function: Optional[Callable] = None) -> "TaskResult":
    """Takes a single action to construct the `TaskResult`.

    Args:
      an_action (action.Action): The `action` object to be executed.
      an_await_function (Optional[Callable]): The function that can be called or is None.

    Returns: A `TaskResult` object representing the running task.

    Raises:
      tea_exception.IllegalArgumentError: If `an_action` is None.
    """
    # <editor-fold desc="Checks">
    if an_action is None:
      logger.error("an_action is None.")
      raise tea_exception.IllegalArgumentError("an_action is None.")
    # </editor-fold>

    tmp_task_result = TaskResult()
    tmp_task_result.connect_basic_signals(an_action)
    tmp_task_result.connect_result_signal(an_action)
    if an_await_function is not None:
      tmp_task_result.task_finished.connect(an_await_function)
    an_action.is_runnable = False
    an_action.set_lock(tmp_task_result.lock)
    tmp_task_result.actions.put(an_action)
    return tmp_task_result

  @staticmethod
  def from_actions(the_actions: tuple["action.Action"], an_await_function: Optional[Callable] = None) -> "TaskResult":
    """Takes a tuple of actions to construct the TaskResult.

    Args:
      the_actions (tuple[action.Action]): A tuple of `action` objects to be executed.
      an_await_function (Optional[Callable]): The function that can be called or is None.

    Returns: A `TaskResult` object representing the running task.

    Raises:
      IllegalArgumentError: If `the_actions` is None or an empty tuple or `the_threadpool` is None.
    """
    # <editor-fold desc="Checks">
    if the_actions is None or len(the_actions) == 0:
      logger.error("the_actions is either None or an empty tuple.")
      raise tea_exception.IllegalArgumentError(
          "the_actions is either None or an empty tuple.",
      )
    # </editor-fold>

    tmp_task_result = TaskResult()
    for tmp_action in the_actions:
      tmp_task_result.connect_basic_signals(tmp_action)
      tmp_task_result.connect_result_signal(tmp_action)
      if an_await_function is not None:
        tmp_task_result.task_finished.connect(an_await_function)
      tmp_action.is_runnable = False
      tmp_action.set_lock(tmp_task_result.lock)
      tmp_task_result.actions.put(tmp_action)
    return tmp_task_result

  @staticmethod
  def get_single_action_result(a_complete_result: tuple[str, list[tuple[bool, tuple]]]) -> tuple[bool, tuple]:
    """Gets a single action result.

    Args:
      a_complete_result (tuple[str, list[tuple[bool, tuple]]]): A tuple of `a_complete_result`.

    Returns:
      The single action result.

    Raises:
      tea_exception.IllegalArgumentError: If `a_complete_result` is None or `a_complete_result` has more than one action result.
    """
    # <editor-fold desc="Checks">
    if a_complete_result is None:
      logger.error("a_complete_result is None.")
      raise tea_exception.IllegalArgumentError("a_complete_result is None.")
    if len(a_complete_result[1]) > 1:
      logger.error("There are more than one action result.")
      raise tea_exception.IllegalArgumentError("There are more than one action result.")
    # </editor-fold>

    tmp_results = a_complete_result[1]
    return tmp_results[0]

  # </editor-fold>
