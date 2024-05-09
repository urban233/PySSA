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
"""Module for the application process manager class."""
import os
import subprocess
import time
from pyssa.gui.ui.custom_dialogs import custom_message_box
from pyssa.util import constants, exception

__docformat__ = "google"


class ApplicationProcessManager:
    """
    Manages the application's processes, including the initiation, checking,
    and closing of the PyMOL and PySSA processes. Instances of this class
    provide an interface for managing these processes.
    """

    def __init__(self, the_reset_pymol_session_func) -> None:
        """
        Constructor.

        Args:
            the_reset_pymol_session_func: The function used to reset the PyMOL session.

        Raises:
            exception.IllegalArgumentError: If the_reset_pymol_session_func is None.

        Returns:
            None
        """
        # <editor-fold desc="Checks">
        if the_reset_pymol_session_func is None:
            raise exception.IllegalArgumentError("the_reset_pymol_session_func is None.")

        # </editor-fold>

        self._reset_pymol_session_func = the_reset_pymol_session_func
        self.pymol_process = None
        self._should_exit = False
        self._is_crashed = False

    def arrange_windows(self):
        """
        Runs a script to arrange the PySSA and PyMOL window on the screen.

        If the script for arranging the windows does not exist, a custom message box is displayed
        with an error message. If the script exists, it is executed using subprocess.Popen(), which
        runs the script in a separate process without displaying a console window.

        Returns:
            None
        """
        if not os.path.exists(constants.ARRANGE_WINDOWS_EXE_FILEPATH):
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "The script for arranging the windows could not be found!",
                "Arrange Windows",
                custom_message_box.CustomMessageBoxIcons.ERROR.value,
            )
            tmp_dialog.exec_()
        else:
            subprocess.Popen([constants.ARRANGE_WINDOWS_EXE_FILEPATH], creationflags=subprocess.CREATE_NO_WINDOW)

    def start_pyssa(self):
        """
        Starts the PySSA process.

        Returns:
            None
        """
        self.pymol_process = subprocess.Popen(
            [r"C:\Users\martin\github_repos\PySSA\scripts\batch\start_pyssa.bat"],
            cwd=r"C:\ProgramData\pyssa\win_start\vb_script",
            creationflags=subprocess.CREATE_NO_WINDOW,
        )
        if self.pymol_process.poll() is None:
            print("PySSA started correctly.")
        else:
            print("PySSA failed to start.")

    def start_pymol(self):
        """
        Starts PyMOL application.

        Returns:
            None
        """
        self.pymol_process = subprocess.Popen(
            [f"{constants.PLUGIN_PATH}\\scripts\\batch\\start_pymol.bat"], creationflags=subprocess.CREATE_NO_WINDOW
        )
        if self.pymol_process.poll() is None:
            print("PyMOL from ApplicationProcessManager class started correctly.")
        else:
            print("PyMOL failed to start.")
            self._is_crashed = True

    def close_manager(self):
        """
        Closes the manager.

        This method sets the internal flag `_should_exit` to `True`, indicating that the manager should exit.

        Returns:
            None
        """
        self._should_exit = True

    def pymol_crashed(self):
        """
        Returns the class attribute _is_crashed if the PyMOL process has crashed.

        Returns:
            bool: True if the PyMOL process has crashed, False otherwise.
        """
        return self._is_crashed

    def check_process(self, placeholder_1, placeholder_2) -> tuple[str, str]:
        """
        Checks the status of a PyMOL process, looped until the process exists or crashes.
        The function also initiates the start of the PyMOL.

        Args:
            placeholder_1 (unknown): Placeholder description for the first parameter.
            placeholder_2 (unknown): Placeholder description for the second parameter.

        Returns:
            tuple: A tuple containing two empty strings needed for the task class.
        """
        self.start_pymol()
        while self._should_exit is False and self._is_crashed is False:
            if self.pymol_process.poll() is not None:
                print("PyMOL crashed!")
                self._is_crashed = True
            else:
                time.sleep(2)
        print("Closing check_process method.")
        return "", ""  # These two empty strings are needed for the task class
