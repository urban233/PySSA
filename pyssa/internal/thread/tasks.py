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
from PyQt5.QtCore import pyqtSignal
from PyQt5 import QtCore


class _Action(QtCore.QObject):
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


class Task:
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
        self.action = _Action(target, args)
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
