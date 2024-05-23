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
"""Module for the TaskManager class."""
import pathlib
import collections

from tea.thread import task, task_result
from tea.util import tea_logging, tea_exception

logger = tea_logging.setup_logger(__file__, pathlib.Path("../../logs/"))
__docformat__ = "google"


class TaskManager:
  """Manages the task objects."""

  # <editor-fold desc="Class attributes">
  running_tasks: collections.deque
  """Contains the running task objects."""

  running_task_results: collections.deque
  """Contains the running task_result objects."""

  # </editor-fold>

  def __init__(self):
    """Constructor."""
    self.running_tasks: collections.deque = collections.deque()
    self.running_task_results: collections.deque = collections.deque()

  def append_task(self, a_task: "task.Task") -> bool:
    """Appends a `a_task` to the deque.

    Args:
      a_task (task.Task): The `task` to append.

    Returns:
       A boolean indicating if the task was successfully appended or not.

    Raises:
      tea_exception.IllegalArgumentError: If `a_task` is None.
    """
    # <editor-fold desc="Checks">
    if a_task is None:
      logger.error("a_task is None.")
      raise tea_exception.IllegalArgumentError("a_task is None.")
    # </editor-fold>

    try:
      a_task.task_finished.connect(self._remove_task)
      self.running_tasks.append(a_task)
    except Exception as e:
      logger.error(e)
      return False
    else:
      return True

  def append_task_result(self, a_task_result: "task_result.TaskResult") -> bool:
    """Appends a `a_task_result` to the deque.

    Args:
      a_task_result (task_result.TaskResult): The `a_task_result` to append.

    Returns:
       A boolean indicating if the a_task_result was successfully appended or not.

    Raises:
      tea_exception.IllegalArgumentError: If `a_task_result` is None.
    """
    # <editor-fold desc="Checks">
    if a_task_result is None:
      logger.error("a_task_result is None.")
      raise tea_exception.IllegalArgumentError("a_task_result is None.")
    # </editor-fold>

    try:
      a_task_result.task_finished.connect(self._remove_task_result)
      self.running_task_results.append(a_task_result)
    except Exception as e:
      logger.error(e)
      return False
    else:
      return True

  def _remove_task(self, return_value: tuple[id, list]) -> bool:
    """Removes a task by the given ID from the deque.

    Args:
      return_value (tuple[id, list]): A tuple representing the data of the `task_finished` signal.

    Returns:
      A boolean indicating if the task was successfully appended or not.

    Raises:
      tea_exception.IllegalArgumentError: If `return_value` is None.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      raise tea_exception.IllegalArgumentError("return_value is None.")
    # </editor-fold>

    try:
      for tmp_task in self.running_tasks:
        if tmp_task.id == return_value[0]:
          self.running_tasks.remove(tmp_task)
          break
    except Exception as e:
      logger.error(e)
      return False
    else:
      return True

  def _remove_task_result(self, return_value: tuple[id, list]) -> bool:
    """Removes a task result by the given ID from the deque.

    Args:
      return_value (tuple[id, list]): A tuple representing the data of the `task_finished` signal.

    Returns:
      A boolean indicating if the task was successfully appended or not.

    Raises:
      tea_exception.IllegalArgumentError: If `return_value` is None.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      raise tea_exception.IllegalArgumentError("return_value is None.")
    # </editor-fold>

    try:
      for tmp_task_result in self.running_task_results:
        if tmp_task_result.id == return_value[0]:
          self.running_task_results.remove(tmp_task_result)
          break
    except Exception as e:
      logger.error(e)
      return False
    else:
      return True
