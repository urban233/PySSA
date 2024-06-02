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
"""Module for the TaskScheduler class."""
import os
import pathlib
from typing import Union

from PyQt5 import QtCore

from src.tea.thread import task
from src.tea.thread import task_result
from src.tea.util import tea_logging, tea_exception

logger = tea_logging.setup_logger(__file__, pathlib.Path("../../logs/"))
__docformat__ = "google"


class TaskScheduler:
  """Scheduler for Tasks."""

  # <editor-fold desc="Class attributes">
  threadpool: QtCore.QThreadPool
  """The threadpool which is used for running tasks."""
  # </editor-fold>

  def __init__(self) -> None:
    """Constructor."""
    self.threadpool = QtCore.QThreadPool()
    self.threadpool.setMaxThreadCount(os.cpu_count())

  def schedule_task(self, a_task: "task.Task") -> None:
    """Schedule task.

    Args:
      a_task (task.Task): The `Task` object to be scheduled.

    Raises:
      tea_exception.IllegalArgumentError: If `a_task` is None.
    """
    # <editor-fold desc="Checks">
    if a_task is None:
      logger.error("a_task is None.")
      raise tea_exception.IllegalArgumentError("a_task is None.")
    # </editor-fold>

    try:
      while a_task.actions.empty() is False:
        tmp_actions = a_task.actions.get()
        a_task.id_of_last_action = tmp_actions.id
        tmp_actions.is_runnable = True
        self.threadpool.start(tmp_actions)
    except Exception as e:
      logger.error(e)

  def schedule_task_result(
      self, a_task_result: "task_result.TaskResult"
  ) -> None:
    """Schedule task result.

    Args:
      a_task_result (task_result.TaskResult): The `TaskResult` object to be scheduled.

    Raises:
      tea_exception.IllegalArgumentError: If `a_task_result` is None.
    """
    # <editor-fold desc="Checks">
    if a_task_result is None:
      logger.error("a_task_result is None.")
      raise tea_exception.IllegalArgumentError("a_task_result is None.")
    # </editor-fold>

    try:
      while a_task_result.actions.empty() is False:
        tmp_actions = a_task_result.actions.get()
        a_task_result.id_of_last_action = tmp_actions.id
        tmp_actions.is_runnable = True
        self.threadpool.start(tmp_actions)
    except Exception as e:
      logger.error(e)

  def schedule(
      self, a_task_to_schedule: Union["task.Task", "task_result.TaskResult"]
  ) -> None:
    """Schedules a task or a task result for execution.

    Args:
      a_task_to_schedule (Union[task.Task, task_result.TaskResult]): The object to be scheduled, either a `Task` or a `TaskResult`.

    Raises:
      tea_exception.IllegalArgumentError: If the provided item is neither a `Task` nor a `TaskResult`.
    """
    # <editor-fold desc="Checks">
    if a_task_to_schedule is None:
      logger.error("a_task_to_schedule is None.")
      raise tea_exception.IllegalArgumentError("a_task_to_schedule is None.")
    # </editor-fold>

    try:
      while a_task_to_schedule.actions.empty() is False:
        tmp_actions = a_task_to_schedule.actions.get()
        a_task_to_schedule.id_of_last_action = tmp_actions.id
        tmp_actions.is_runnable = True
        self.threadpool.start(tmp_actions)
    except Exception as e:
      logger.error(e)
