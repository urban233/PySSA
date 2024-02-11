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
"""Module for the add proteins dialog."""
import os
import pymol
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtCore import Qt
from pymol import cmd
from PyQt5 import QtCore
from PyQt5 import QtGui
from PyQt5 import QtWidgets
from pyqtspinner import spinner
from pyssa.gui.ui import icon_resources  # this import is used for the icons! DO NOT DELETE THIS
from pyssa.gui.ui.forms.auto_generated import auto_predict_monomer_view
from pyssa.gui.ui.styles import styles
from pyssa.util import constants, tools, gui_utils


class PredictMonomerView(QtWidgets.QDialog):
    """Class for a dialog to add proteins to a project."""

    """
    A pyqtsignal that is used to hand-over the protein structure information.
    """
    return_value = pyqtSignal(tuple)

    def __init__(self, parent=None) -> None:  # noqa: ANN001
        """Constructor.

        Args:
            parent: The parent.
        """
        QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = auto_predict_monomer_view.Ui_Dialog()
        self.ui.setupUi(self)
        self.wait_spinner = spinner.WaitingSpinner(
            parent=self,
            center_on_parent=True,
            disable_parent_when_spinning=True,
            modality=Qt.ApplicationModal,
            roundness=100.0,
            fade=45.0,
            radius=14,
            lines=8,
            line_length=17,
            line_width=10,
            speed=1.25,
            color=QtGui.QColor(75, 145, 247),
        )
        self.ui.tabWidget.setTabEnabled(1, False)
        self._initialize_ui()
        self.setModal(True)

    def _initialize_ui(self) -> None:
        """Initialize the UI elements."""
        styles.color_bottom_frame_button(self.ui.btn_pred_analysis_mono_start)
        styles.color_bottom_frame_button(self.ui.btn_pred_analysis_mono_go_analysis_setup)
        self.ui.btn_help.setIcon(QtGui.QIcon(":/icons/help_w200.svg"))
        self.ui.btn_help.setIconSize(self.ui.btn_help.icon().actualSize(QtCore.QSize(30, 30)))
        self.ui.btn_help.setText("")
        self.ui.btn_help_2.setIcon(QtGui.QIcon(":/icons/help_w200.svg"))
        self.ui.btn_help_2.setIconSize(self.ui.btn_help.icon().actualSize(QtCore.QSize(30, 30)))
        self.ui.btn_help_2.setText("")
        styles.set_stylesheet(self)
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        self.setWindowTitle("Prediction Monomer")
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
