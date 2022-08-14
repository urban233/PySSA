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

import sys
import os

# appending the plugin path to the python search path
# this is necessary to be able to use import (submodules) statements
# sys.path.append("/home/matt/snap/pymol-oss/107/.pymol/startup/pymol_plugin")
# Linux path
sys.path.append("/home/matt/anaconda3/envs/pymol_plugin/lib/python3.9/site-packages/pmg_tk/startup/pymol_plugin/")
# MacOS path
sys.path.append("/Users/matt/opt/anaconda3/envs/pymol/lib/python3.9/site-packages/pmg_tk/startup/pymol_plugin/")
# Windows path
sys.path.append(f"{os.path.expanduser('~')}\\anaconda3\\envs\\pyssa\\lib\\site-packages\\pmg_tk\\startup\\pyssa\\")


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
    from main import MainWindow
    # getting the value of the global var mainWindow
    global mainWindow

    if mainWindow is None:
        mainWindow = MainWindow()

    mainWindow.show()
