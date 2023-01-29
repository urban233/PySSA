#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2022
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
from dataclasses import dataclass
from internal.thread import workers
from PyQt5 import QtCore


@dataclass
class ThreadWorkerPair:
    # first the thread where the worker will be assigned to, then the worker itself
    main_task: str
    thread: QtCore.QThread
    worker: workers.PredictionWorker or workers.AnalysisWorker

    def move_worker_to_thread(self):
        self.worker.moveToThread(self.thread)

    def quit_thread_and_delete_worker(self):
        self.thread.quit()
        self.worker.deleteLater()

    def create_standard_connections(self):
        self.thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.quit_thread_and_delete_worker)

    def define_custom_thread_finished_connection(self, function):
        """

        Args:
            function:
                function which should get executed after the thread has completed its task

        Notes:
            the function needs to have the "self.prediction_thread.deleteLater()" function at the end
        """
        self.thread.finished.connect(function)

    def start_thread(self):
        self.thread.start()

    def setup_and_run_thread(self, function):
        self.move_worker_to_thread()
        self.create_standard_connections()
        self.define_custom_thread_finished_connection(function)
        self.start_thread()
