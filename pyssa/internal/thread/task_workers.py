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
import pathlib
from PyQt5.QtCore import QObject, pyqtSignal, QThread
from pyssa.util import tools


def setup_worker_for_work(tmp_thread, tmp_worker, return_value_func):
    tmp_worker.moveToThread(tmp_thread)
    tmp_thread.started.connect(tmp_worker.run)
    tmp_worker.finished.connect(tmp_thread.quit)
    tmp_worker.finished.connect(tmp_worker.deleteLater)
    tmp_thread.finished.connect(tmp_thread.deleteLater)
    tmp_worker.return_value.connect(return_value_func)
    return tmp_thread


# Step 1: Create a worker class
class Worker(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal(int)
    return_value = pyqtSignal(tuple)
    workspace_path: str

    def __init__(self, workspace):
        super().__init__()
        self.workspace_path = workspace

    def run(self):
        protein_dict, protein_names = tools.scan_workspace_for_non_duplicate_proteins(pathlib.Path(self.workspace_path))
        self.return_value.emit((protein_dict, protein_names))
        self.finished.emit()
