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
"""Module for all task workers which run in separate threads."""
import logging
import queue
import uuid
from typing import Optional, Callable

from PyQt5.QtCore import pyqtSignal
from PyQt5 import QtCore

from pyssa.logging_pyssa import log_handlers
from pyssa.util import exception, enums

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class _LegacyAction(QtCore.QObject):
    """Actual execution unit."""

    is_finished: bool
    finished = pyqtSignal()
    result = pyqtSignal(tuple)

    def __init__(self, target, args=()) -> None:  # noqa: ANN001
        """Constructor.

        Args:
            target: function that should get executed.
            args: arguments for the function.
        """
        super().__init__()
        self.target = target
        self.args = args
        self.is_finished = False

    def run_action(self) -> None:
        """Runs the function."""
        tmp_result: tuple = self.target(*self.args)
        self.is_finished = True
        self.finished.emit()
        self.result.emit(tmp_result)


class LegacyTask:
    """Container for running an asynchronous operation."""

    def __init__(self, target, args=(), post_func=None) -> None:  # noqa: ANN001
        """Constructor.

        Args:
            target: function that should get executed.
            args: arguments for the function.
            post_func: function that runs after the target finished.

        Note:
            The post function needs a tuple as argument! After the action finished it will return a
            tuple of the results to the post function.

            Be aware that no UI related functions can be used because this would lead to
            an unauthorized memory access violation! Do everything related to the UI in the post function if
            necessary.
        """
        self.action = _LegacyAction(target, args)
        self.thread = QtCore.QThread()

        self.action.moveToThread(self.thread)
        self.action.finished.connect(self.thread.quit)
        self.thread.finished.connect(self.thread.deleteLater)
        if post_func is not None:
            self.action.result.connect(post_func)

    def start(self) -> None:
        """Starts the thread with the action."""
        self.thread.started.connect(self.action.run_action)
        self.thread.start()

    def wait(self) -> None:
        """Waits for the action to finish.

        Note:
            This will block the thread until the thread is finished!
        """
        self.thread.wait()

    def is_finished(self) -> bool:
        """Checks if the action is finished."""
        if self.action.is_finished:
            return True
        return False


class _TaskSignals(QtCore.QObject):
    """Signals related to the tasks."""
    finished = pyqtSignal()
    error = pyqtSignal(tuple)
    result = pyqtSignal(tuple)


class Action(QtCore.QRunnable):
    """A container for running an asynchronous operation."""
    def __init__(self, a_target, args: tuple = ()) -> None:
        """Constructor.

        Args:
            a_target: The function that should get executed.
            args (tuple): The arguments for the function.

        Raises:
             exception.IllegalArgumentError: If an argument is None.
        """
        # <editor-fold desc="Checks">
        if a_target is None:
            logger.error("a_target is None")
            raise exception.IllegalArgumentError("a_target is None")
        if args is None:
            logger.error("args is None")
            raise exception.IllegalArgumentError("args is None")

        # </editor-fold>

        super().__init__()
        self._signals = _TaskSignals()
        self._target = a_target
        self._args: tuple = args
        self.is_runnable = False

    def connect_error_signal(self, a_function: Callable) -> None:
        self._signals.error.connect(a_function)

    def connect_finished_signal(self, a_function: Callable) -> None:
        self._signals.finished.connect(a_function)

    def run(self) -> None:
        """Overwrites the run() method of QRunnable and executes the target function."""
        try:
            if self.is_runnable is False:
                raise exception.ActionIsNotRunnableError()
            tmp_result = self._target(*self._args)
        except Exception as e:
            logger.error(e)
            self._signals.error.emit((e,))
        else:
            self._signals.result.emit((True, tmp_result))
        finally:
            self._signals.finished.emit()


class Task(QtCore.QObject):
    """Represents an asynchronous operation."""
    """
    An ID for the task instance.
    """
    id: uuid.UUID

    """
    Gets whether this task instance has completed execution due to being canceled. 
    """
    is_canceled: bool

    """
    Gets whether the task ran to completion.    
    """
    is_completed: bool

    """
    Gets the task status of this task.
    """
    status: Optional["enums.TaskStatus"]

    """
    A queue of actions to execute.
    """
    actions: queue.Queue

    def __init__(self):
        """Constructor."""
        super().__init__()
        self.id = uuid.uuid4()
        self.is_canceled = False
        self.is_completed = False
        self.status = None
        self.actions = queue.Queue()

    def _action_completed(self) -> None:
        if self.actions.empty() and self.status != enums.TaskStatus.FAILED:
            self.is_completed = True
            self.status = enums.TaskStatus.RAN_TO_COMPLETION
        elif self.actions.empty() and self.status == enums.TaskStatus.FAILED:
            # Case is handled by method _action_ended_with_error
            return
        else:
            self.is_completed = False
            self.status = enums.TaskStatus.RUNNING

    def _action_ended_with_error(self) -> None:
        if self.actions.empty():
            self.is_completed = True
            self.status = enums.TaskStatus.FAILED
        else:
            self.is_completed = False
            self.status = enums.TaskStatus.RUNNING

    def connect_basic_signals(self, an_action: "Action") -> None:
        an_action.connect_error_signal(self._action_ended_with_error)
        an_action.connect_finished_signal(self._action_completed)

    @staticmethod
    def run_action(an_action: "Action", the_threadpool: QtCore.QThreadPool) -> "Task":
        tmp_task: "Task" = Task()
        tmp_task.actions.put(an_action)
        tmp_action: "Action" = tmp_task.actions.get()
        tmp_task.connect_basic_signals(tmp_action)
        the_threadpool.start(tmp_action)
        return tmp_task
