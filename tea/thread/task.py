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
"""Module for the Task class."""
import pathlib
import queue
import uuid
from typing import Optional

from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal

from tea.thread import action
from tea.util import tea_exception, tea_enums
from tea.util import tea_logging

logger = tea_logging.setup_logger(__file__, pathlib.Path("../../logs/"))
__docformat__ = "google"


class Task(QtCore.QObject):
  """Represents an asynchronous operation."""

  # <editor-fold desc="Class attributes">
  id: uuid.UUID  # noqa: A003
  """An ID for the task instance."""

  is_canceled: bool
  """Gets whether this task instance has completed execution due to being canceled."""

  is_completed: bool
  """Gets whether the task ran to completion."""

  task_finished = pyqtSignal(tuple)
  """Signals that the task finished executing."""

  status: Optional["tea_enums.TaskStatus"]
  """Gets the task status of this task."""

  actions: queue.Queue
  """A queue of actions to execute."""

  id_of_last_action: Optional[str]
  """The id of the last action in the queue."""

  _lock: QtCore.QMutex
  """A lock used to indicate if the action is executable."""

  # </editor-fold>

  def __init__(self):
    """Constructor."""
    super().__init__()
    self.id: str = str(uuid.uuid4())
    self.is_canceled: bool = False
    self.is_completed: bool = False
    self._status: Optional["tea_enums.TaskStatus"] = None
    self.actions: queue.Queue = queue.Queue()
    self.id_of_last_action: Optional[str] = None
    self.lock = QtCore.QMutex()

  # <editor-fold desc="Private methods">
  def _action_started(self) -> None:
    """Updates the task status to RUNNING."""
    self.status = tea_enums.TaskStatus.RUNNING

  def _action_completed(self, the_action_id: str) -> None:
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
      self.task_finished.emit((self.id, []))
    elif (
        the_action_id == self.id_of_last_action
        and self.status == tea_enums.TaskStatus.FAILED
    ):
      # Case is handled by method _action_ended_with_error
      return
    else:
      self.is_completed: bool = False
      self.status: "tea_enums.TaskStatus" = tea_enums.TaskStatus.RUNNING

  def _action_ended_with_error(self) -> None:
    """Update the task status and completion status when the action ended with an error."""
    if self.actions.empty():
      self.is_completed: bool = True
      self.status: "tea_enums.TaskStatus" = tea_enums.TaskStatus.FAILED
    else:
      self.is_completed: bool = False
      self.status: "tea_enums.TaskStatus" = tea_enums.TaskStatus.RUNNING

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
    an_action.connect_finished_signal(self._action_completed)

  def run(self, the_threadpool: QtCore.QThreadPool) -> None:
    """Puts the actions from the queue in the `QThreadPool` where they are executed.

    Args:
      the_threadpool (QtCore.QThreadPool): A `QThreadPool` object used for executing the actions.

    Raises:
      tea_exception.IllegalArgumentError: If `the_threadpool` is None.
    """
    # <editor-fold desc="Checks">
    if the_threadpool is None:
      logger.error("the_threadpool is None.")
      raise tea_exception.IllegalArgumentError("the_threadpool is None.")

    # </editor-fold>

    while self.actions.empty() is False:
      tmp_action: "action.Action" = self.actions.get()
      self.connect_basic_signals(tmp_action)
      tmp_action.is_runnable = True
      the_threadpool.start(tmp_action)

  # </editor-fold>

  # <editor-fold desc="Static methods">
  @staticmethod
  def run_action(
      an_action: "action.Action",
      the_threadpool: QtCore.QThreadPool,
  ) -> "Task":
    """Static constructor for the Task class that takes an action to construct the
    `Task` and starts the task.

    Args:
      an_action (action.Action): The `action` to be executed.
      the_threadpool (QtCore.QThreadPool): The `thread pool` to be used for executing the action.

    Returns:
      A `Task` object representing the running task.

    Raises:
      tea_exception.IllegalArgumentError: If `an_action` or the_threadpool is None.
    """
    # <editor-fold desc="Checks">
    if an_action is None:
      logger.error("an_action is None.")
      raise tea_exception.IllegalArgumentError("an_action is None.")
    if the_threadpool is None:
      logger.error("the_threadpool is None.")
      raise tea_exception.IllegalArgumentError("the_threadpool is None.")

    # </editor-fold>

    tmp_task: "Task" = Task()
    try:
      tmp_task.actions.put(an_action)
      tmp_action: "action.Action" = tmp_task.actions.get()
      tmp_task.connect_basic_signals(tmp_action)
      tmp_task.id_of_last_action = tmp_action.id
      tmp_action.is_runnable = True
      the_threadpool.start(tmp_action)
      return tmp_task
    except Exception as e:
      logger.error(e)

  @staticmethod
  def from_action(an_action: "action.Action") -> "Task":
    """Converts a single action into a `Task`.

    Args:
      an_action (action.Action): The `action` object to convert.

    Returns: A new `Task` object representing the provided action.

    Raises:
      tea_exception.IllegalArgumentError: If `an_action` is None.
    """
    # <editor-fold desc="Checks">
    if an_action is None:
      logger.error("an_action is None.")
      raise tea_exception.IllegalArgumentError("an_action is None.")
    # </editor-fold>

    a_task = Task()
    an_action.set_lock(a_task.lock)
    a_task.actions.put(an_action)
    return a_task

  @staticmethod
  def from_actions(the_actions: tuple["action.Action"]) -> "Task":
    """Combines multiple actions into a single Task.

    Args:
      the_actions (tuple[action.Action]): The `action` objects to combine.

    Returns: A new `Task` object representing the combined actions.

    Raises:
      tea_exception.IllegalArgumentError: If `the_actions` is None.
    """
    # <editor-fold desc="Checks">
    if the_actions is None:
      logger.error("the_actions is None.")
      raise tea_exception.IllegalArgumentError("the_actions is None.")
    # </editor-fold>

    a_task = Task()
    for tmp_action in the_actions:
      tmp_action.set_lock(a_task.lock)
      a_task.actions.put(tmp_action)
    return a_task

  # </editor-fold>
