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
import os
import logging
import subprocess
from pyssa.gui.ui.forms.auto_generated.auto_dialog_settings_global import Ui_Dialog
from pyssa.internal.data_structures import settings
from pyssa.util import constants, gui_utils
from pyssa.gui.ui.styles import styles
from pyssa.internal.thread import workers
from PyQt5.QtWidgets import QMessageBox
import PyQt5
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.gui.ui.messageboxes import settings_boxes

# setup logger
logging.basicConfig(level=logging.DEBUG)


class DialogHelp(QtWidgets.QDialog):
    """This class opens a help dialog."""
    """This variable is for controlling whether the dialog opens or not"""
    ERROR = False

    def __init__(self, html_content: str, parent=None) -> None:
        """Constructor.

        Args:
            args
            kwargs
        """
        QtWidgets.QDialog.__init__(self, parent)
        self.text_browser = QtWidgets.QTextBrowser(self)
        self.text_browser.setGeometry(10, 10, 780, 580)
        self.text_browser.setHtml(html_content)
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
        styles.set_stylesheet(self)
        self.text_browser.setStyleSheet("""background-color: white;""")
        self.setWindowTitle("Help")



