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
"""Module for the add protein pair view."""
from PyQt5.QtCore import Qt
from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from pyssa.gui.ui.forms.auto_generated import auto_add_protein_pair_view

from pyssa.gui.ui.styles import styles
from pyssa.util import constants


class AddProteinPairView(QtWidgets.QDialog):
    """Dialog for adding protein pairs."""
    
    dialogClosed = QtCore.pyqtSignal(tuple)
    """A signal indicating that the dialog is closed."""
    
    def __init__(self, parent=None) -> None:  # noqa: ANN001
        """Constructor."""
        QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = auto_add_protein_pair_view.Ui_Dialog()
        self.ui.setupUi(self)
        styles.color_bottom_frame_button(self.ui.btn_add)
        self.resize(550, 650)

        self.ui.btn_cancel.clicked.connect(self.close)
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        styles.set_stylesheet(self)
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
        self.setWindowTitle("Add Protein Pair")
        self.setModal(True)
