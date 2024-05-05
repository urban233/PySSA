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
"""Module for the import sequence view."""
from PyQt5.QtCore import pyqtSignal
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5 import QtGui
from pyssa.gui.ui import icon_resources  # this import is used for the icons! DO NOT DELETE THIS
from pyssa.gui.ui.forms.auto_generated import auto_import_sequence_view
from pyssa.gui.ui.styles import styles
from pyssa.util import constants

global_var_add_model = ("", False)


class ImportSequenceView(QtWidgets.QDialog):
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
        self.ui = auto_import_sequence_view.Ui_Dialog()
        self.ui.setupUi(self)
        self.ui.btn_import_sequence.setEnabled(False)
        self.ui.txt_import_sequence.setEnabled(False)
        styles.color_bottom_frame_button(self.ui.btn_import_sequence)
        self.ui.lbl_status.setText("")
        self.ui.btn_choose_fasta_file.setToolTip("Click to add a .fasta file")
        self.setWindowTitle("Import Protein Sequence")
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        styles.set_stylesheet(self)
        self.ui.btn_help.setIcon(QtGui.QIcon(":/icons/help_w200.png"))
        self.ui.btn_help.setIconSize(self.ui.btn_help.icon().actualSize(QtCore.QSize(30, 30)))
        self.ui.btn_help.setText("")
        # fixme: this flag needs to be set if the WhatsThat icon in the window bar should be hidden
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
        self.setModal(True)
        # # check internet connectivity
        # if not tools.check_internet_connectivity():
        #     gui_utils.no_internet_dialog_with_custom_msg(
        #         "You do not have a working internet connection which is "
        #         "necessary for connecting to the PDB!\n"
        #         "However you can add a protein structure from "
        #         "your local filesystem.",
        #     )
        #     self.ui.txt_import_sequence.setEnabled(False)
        #     self.ui.lbl_status.setText("You cannot enter a PDB ID (no internet).")


