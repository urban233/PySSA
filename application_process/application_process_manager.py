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
"""Module for the application process manager class."""
import subprocess
import time
from pyssa.util import constants


class ApplicationProcessManager:
    def __init__(self):
        self.pymol_process = None
        self._should_exit = False
        self._is_crashed = False

    def start_pyssa(self):
        self.pymol_process = subprocess.Popen(
            [r"C:\Users\martin\github_repos\PySSA\scripts\batch\start_pyssa.bat"],
            cwd=r"C:\ProgramData\pyssa\win_start\vb_script",
            creationflags=subprocess.CREATE_NO_WINDOW
        )
        if self.pymol_process.poll() is None:
            print("PySSA started correctly.")
        else:
            print("PySSA failed to start.")

    def start_pymol(self):
        self.pymol_process = subprocess.Popen(
            [f"{constants.PLUGIN_PATH}\\scripts\\batch\\start_pymol.bat"],
            creationflags=subprocess.CREATE_NO_WINDOW
        )
        if self.pymol_process.poll() is None:
            print("PyMOL from ApplicationProcessManager class started correctly.")
        else:
            print("PyMOL failed to start.")

    def close_manager(self):
        self._should_exit = True

    def pymol_crashed(self):
        return self._is_crashed

    def check_process(self, placeholder_1, placeholder_2):
        self.start_pymol()
        while True:
            if self.pymol_process.poll() is not None:
                print("PyMOL crashed!")
                self._is_crashed = True
                self.start_pymol()
            time.sleep(2)
            if self._should_exit:
                break
