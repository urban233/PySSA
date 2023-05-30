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

from pyssa.gui.ui.styles import styles
from pyssa.util import constants


def yes_or_no(window_title, text_message, icon) -> bool:
    """This is a function which creates a basic QMessageBox with the buttons yes or no.

    Returns:
        True if "yes" is clicked
        False if "no" is clicked
    """
    msg = QMessageBox()
    msg.setIcon(icon)
    msg.setWindowIcon(PyQt5.QtGui.QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
    styles.set_stylesheet(msg)
    msg.setWindowTitle(window_title)
    msg.setText(text_message)
    btn_yes = msg.addButton("Yes", QMessageBox.ActionRole)
    btn_no = msg.addButton("No", QMessageBox.ActionRole)
    msg.exec_()

    if msg.clickedButton() == btn_yes:
        return True
    elif msg.clickedButton() == btn_no:
        msg.close()
        return False
    else:
        print("Unexpected Error.")
        msg.close()


def ok(window_title, text_message, icon):
    msg = QMessageBox()
    msg.setIcon(icon)
    msg.setWindowTitle(window_title)
    msg.setWindowIcon(PyQt5.QtGui.QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
    styles.set_stylesheet(msg)
    msg.setText(text_message)
    msg.setStandardButtons(QMessageBox.Ok)
    msg.exec_()

    if msg.clickedButton() == msg.setStandardButtons(QMessageBox.Ok):
        return True
    else:
        msg.close()
        return False


def no_buttons(window_title, text_message, icon):
    msg = QMessageBox()
    msg.setIcon(icon)
    msg.setWindowTitle(window_title)
    msg.setWindowIcon(PyQt5.QtGui.QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
    styles.set_stylesheet(msg)
    msg.setText(text_message)
    msg.setStandardButtons(QMessageBox.NoButton)
    return msg
