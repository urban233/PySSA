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
import sys
import os
from PyQt5.QtCore import Qt

# appending the plugin path to the python search path
# this is necessary to be able to use import (submodules) statements

# linux, macOS, Windows path
path_list = [
    f"/home/{os.getlogin()}/.local/pyssa/pyssa-mamba-env/lib/python3.10/site-packages/pmg_tk/startup/PySSA"
    f"C:\\ProgramData\\pyssa\\plugin\\Miniconda3\\envs\\pyssa_colab\\Lib\\site-packages\\pymol\\pymol_path\\data\\startup\\PySSA",
]
styles_path_list = [
    pathlib.Path(f"{path_list[0]}/pyssa/gui/ui/styles/styles.css"),
    pathlib.Path(f"{path_list[1]}/pyssa/gui/ui/styles/styles.css"),
]
# appends the os specific python path
if sys.platform.startswith("linux"):
    # Linux path
    sys.path.append(path_list[0])
elif sys.platform.startswith("win32"):
    # Windows path
    sys.path.append(path_list[1])


def __init_plugin__(app=None):
    """This function creates an entry in the PyMOL "Plugin" menu

    Args:
        app (optional):
            None
    """
    from pymol.plugins import addmenuitemqt
    plugin_name = 'PySSA'
    addmenuitemqt(plugin_name, run_plugin_gui)


# global reference to avoid garbage collection of our dialog
mainWindow = None


def run_plugin_gui():
    """This function is the entry point for the plugin to start.

    NOTE:
        Open the custom dialog but WITHOUT using:
        app = Qt.QtWidgets.QApplication(sys.argv)
        app.exec_()
    """
    from .pyssa.main import MainWindow
    # getting the value of the global var mainWindow
    global mainWindow

    if mainWindow is None:
        mainWindow = MainWindow()
        # Open the qss styles file and read in the css-alike styling code
        if sys.platform.startswith("linux"):
            with open(styles_path_list[0], 'r', encoding="utf-8") as file:
                style = file.read()
        elif sys.platform.startswith("win32"):
            with open(styles_path_list[1], 'r', encoding="utf-8") as file:
                style = file.read()
        # Set the stylesheet of the application
        mainWindow.setStyleSheet(style)

    mainWindow.show()
