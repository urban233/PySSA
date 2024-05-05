import sys
from PyQt5.QtWidgets import QDialog, QVBoxLayout, QTableWidget, QPushButton, QHBoxLayout, QLabel
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5 import QtGui

from pyssa.gui.ui.custom_delegates import sequence_table_delegate
from pyssa.gui.ui.styles import styles
from pyssa.util import constants
from pyssa.gui.ui import icon_resources  # this import is used for the icons! DO NOT DELETE THIS
from pyssa.gui.ui.forms.auto_generated.auto_fasta_file_import_preview_view import Ui_Dialog


class FastaFileImportPreviewView(QDialog):
    """A QDialog that allows users to customize the fasta file import."""

    def __init__(self, parent=None):
        QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)

        self.ui.table.setItemDelegate(sequence_table_delegate.InputCheckDelegate(self.ui.table))

        self.ui.btn_help.setIcon(QtGui.QIcon(":/icons/help_w200.png"))
        self.ui.btn_help.setIconSize(self.ui.btn_help.icon().actualSize(QtCore.QSize(30, 30)))
        self.ui.btn_help.setText("")
        self.setWindowTitle("FASTA File Import Preview")
        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        styles.set_stylesheet(self)
        self.resize(900, 600)
        # fixme: this flag needs to be set if the WhatsThat icon in the window bar should be hidden
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
        self.setModal(True)
