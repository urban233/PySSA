#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
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
"""Module for the fasta file import preview view."""
from PyQt5.QtWidgets import QDialog
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from PyQt5 import QtGui

from src.pyssa.gui.ui.custom_delegates import sequence_table_delegate
from src.pyssa.gui.ui.styles import styles
from src.pyssa.util import constants
from src.pyssa.gui.ui import icon_resources  # this import is used for the icons! DO NOT DELETE THIS
from src.pyssa.gui.ui.forms.auto_generated.auto_fasta_file_import_preview_view import Ui_Dialog


class FastaFileImportPreviewView(QDialog):
  """A QDialog that allows users to customize the fasta file import."""

  def __init__(self, parent=None) -> None:
    """Constructor."""
    QtWidgets.QDialog.__init__(self, parent)
    # build ui object
    self.ui = Ui_Dialog()
    self.ui.setupUi(self)

    self.ui.table.setItemDelegate(
        sequence_table_delegate.InputCheckDelegate(self.ui.table)
    )

    self.ui.btn_help.setIcon(QtGui.QIcon(":/icons/help_w200.png"))
    self.ui.btn_help.setIconSize(
        self.ui.btn_help.icon().actualSize(QtCore.QSize(30, 30))
    )
    self.ui.btn_help.setText("")
    self.ui.btn_cancel.clicked.connect(self.close)
    self.setWindowTitle("FASTA File Import Preview")
    self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
    styles.set_stylesheet(self)
    self.resize(900, 600)
    # fixme: this flag needs to be set if the WhatsThat icon in the window bar should be hidden
    self.setWindowFlags(
        self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint
    )
    self.setWindowModality(Qt.WindowModal)
