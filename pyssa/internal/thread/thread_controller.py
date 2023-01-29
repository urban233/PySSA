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
from pyssa.internal.thread import thread_worker_pair


class ThreadController:

    def __init__(self, thread_worker_pairs: dict[str, thread_worker_pair.ThreadWorkerPair]):
        self.thread_worker_pairs: dict[str, thread_worker_pair.ThreadWorkerPair] = thread_worker_pairs

    def delete_specific_thread_later(self, main_task: str):
        for key in self.thread_worker_pairs:
            if self.thread_worker_pairs.get(key).main_task == main_task:
                self.thread_worker_pairs.get(key).thread.deleteLater()

    def get_all_running_threads(self):
        running_threads = []
        for key in self.thread_worker_pairs:
            if self.thread_worker_pairs.get(key).thread.isRunning() is True:
                running_threads.append(self.thread_worker_pairs.get(key))
        return running_threads

    def create_and_add_new_thread_worker_pair(self, main_task, thread, worker):
        for key in self.thread_worker_pairs:
            if key == main_task:
                raise ValueError(f"A thread_worker_pair with the main task {main_task} already exists!")
        self.thread_worker_pairs.update({main_task: thread_worker_pair.ThreadWorkerPair(main_task, thread, worker)})
