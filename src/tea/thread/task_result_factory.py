#
# TEA - Task Event-based Async library for Python
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
"""Module for the TaskResultFactory class."""
import pathlib

from src.tea.thread import task_factory
from src.tea.thread import task_result
from src.tea.thread import task_scheduler
from src.tea.util import tea_logging, tea_exception

logger = tea_logging.setup_logger(__file__, pathlib.Path("../../logs/"))
__docformat__ = "google"


class TaskResultFactory(task_factory.TaskFactory):
  """Factory for creating and running task result objects."""

  @staticmethod
  def run_task_result(
      a_task_result: "task_result.TaskResult",
      a_task_scheduler: "task_scheduler.TaskScheduler",
  ) -> "task_result.TaskResult":
    """Runs the given task.
    Args:
      a_task_result (task_result.TaskResult): The `task` to run.
      a_task_scheduler (task_scheduler.TaskScheduler): The `task_scheduler`.

    Returns:
      task_result.TaskResult: The `task` result.

    Raises:
      tea_exception.IllegalArgumentError: If `a_task_result` or `a_task_scheduler` is None.
    """
    # <editor-fold desc="Checks">
    if a_task_result is None:
      logger.error("Task result is None.")
      raise tea_exception.IllegalArgumentError("Task result is None.")

    if a_task_scheduler is None:
      logger.error("Task scheduler is None.")
      raise tea_exception.IllegalArgumentError("Task scheduler is None.")
    # </editor-fold>

    a_task_scheduler.schedule_task_result(a_task_result)
    return a_task_result
