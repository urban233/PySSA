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
"""Module for all message boxes which occur in the settings dialog."""
import PyQt5
from PyQt5.QtWidgets import QMessageBox
from pyssa.gui.utilities import constants


def restart_now_later() -> bool:
    """This is a function which creates a basic QMessageBox with the buttons yes or no.

    Returns:
        True if "yes" is clicked
        False if "no" is clicked
    """
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Information)
    msg.setWindowIcon(PyQt5.QtGui.QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
    msg.setWindowTitle("WSL2 installation")
    msg.setText("Installation is finished! Restart is necessary.")
    btn_later = msg.addButton("Restart later", QMessageBox.ActionRole)
    btn_now = msg.addButton("Restart now", QMessageBox.ActionRole)
    msg.exec_()

    if msg.clickedButton() == btn_now:
        return True
    elif msg.clickedButton() == btn_later:
        msg.close()
        return False
    else:
        print("Unexpected Error.")
        msg.close()
