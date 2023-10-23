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
from PyQt5.QtGui import QIcon
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore
from pyssa.gui.ui.forms.auto_generated.auto_dialog_about import Ui_Dialog
from pyssa.gui.ui.styles import styles
from pyssa.util import constants


class DialogAbout(QtWidgets.QDialog):

    def __init__(self, parent=None):
        """Constructor.

        Args:
            args
            kwargs
        """
        QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)

        self.ui.btn_ok.clicked.connect(self.close_dialog)
        original_pixmap = QtGui.QPixmap(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png")
        scaled_pixmap = original_pixmap.scaled(150, 150)
        self.ui.lbl_pyssa_logo.setPixmap(scaled_pixmap)
        self.ui.label_2.setText(f"Version: {constants.VERSION_NUMBER}")
        styles.set_stylesheet(self)
        self.setWindowIcon(QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
        self.setWindowTitle("About")
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)

    # @SLOT
    def close_dialog(self):
        self.close()
