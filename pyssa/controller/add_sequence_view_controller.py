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
"""Module for the Open Dialog."""

import glob
import os
from PyQt5 import QtCore
from PyQt5.QtCore import Qt

from pyssa.controller import interface_manager, pymol_session_manager
from pyssa.gui.ui.styles import styles
from pyssa.internal.thread import tasks
from pyssa.internal.thread.async_pyssa import util_async
from pyssa.io_pyssa import safeguard
from pyssa.util import input_validator, constants, tools


class AddSequenceViewController(QtCore.QObject):
    """Class for the Open Project View Controller."""
    return_value = QtCore.pyqtSignal(tuple)

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager") -> None:
        super().__init__()
        self._interface_manager = the_interface_manager
        self._view = the_interface_manager.get_add_sequence_view()
        self._connect_all_ui_elements_to_slot_functions()
        self.restore_default_view()

    def open_help(self, a_page_name: str):
        """Opens the pyssa documentation window if it's not already open.

        Args:
            a_page_name (str): a name of a documentation page to display
        """
        self._interface_manager.update_status_bar("Opening help center ...")
        self._active_task = tasks.Task(
            target=util_async.open_documentation_on_certain_page,
            args=(a_page_name, 0),
            post_func=self.__await_open_help,
        )
        self._active_task.start()

    def __await_open_help(self):
        self._interface_manager.update_status_bar("Opening help center finished.")

    def _open_help_for_dialog(self):
        self.open_help("help/project/open_project/")

    def restore_default_view(self) -> None:
        self._view.ui.le_seq_name.clear()
        self._view.ui.le_protein_seq.clear()

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.le_seq_name.textChanged.connect(self._validate_protein_name)
        self._view.ui.le_protein_seq.textChanged.connect(self._validate_protein_sequence)
        self._view.ui.btn_next.clicked.connect(self._switch_ui_to_sequence_input)
        self._view.ui.btn_back.clicked.connect(self._switch_ui_to_sequence_name_input)
        self._view.ui.btn_add.clicked.connect(self._add_sequence)
        # self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)

    def _validate_protein_name(self) -> None:
        """Validates the input of the protein name in real-time."""
        tools.validate_protein_name(
            self._view.ui.le_seq_name,
            self._view.ui.lbl_status,
            self._view.ui.btn_next,
        )

    def _validate_protein_sequence(self) -> None:
        """Validates the input of the protein sequence in real-time."""
        tools.validate_protein_sequence(
            self._view.ui.le_protein_seq,
            self._view.ui.lbl_status_seq,
            self._view.ui.btn_add,
        )

    def _switch_ui_to_sequence_input(self):
        self._view.ui.le_seq_name.setEnabled(False)
        self._view.ui.btn_next.setEnabled(False)
        self._view.ui.lbl_protein_seq.show()
        self._view.ui.le_protein_seq.show()
        self._view.ui.lbl_status_seq.show()
        self._view.ui.btn_back.show()
        self._activate_add_button()

    def _switch_ui_to_sequence_name_input(self):
        self._view.ui.le_seq_name.setEnabled(True)
        self._view.ui.btn_next.setEnabled(True)
        self._view.ui.lbl_protein_seq.hide()
        self._view.ui.le_protein_seq.hide()
        self._view.ui.lbl_status_seq.hide()
        self._view.ui.btn_back.hide()
        self._view.ui.btn_add.setEnabled(False)

    def _activate_add_button(self) -> None:
        """Activates the open button."""
        if self._view.ui.le_protein_seq.toPlainText() == "":
            self._view.ui.btn_add.setEnabled(False)
        else:
            self._view.ui.btn_add.setEnabled(True)

    def _add_sequence(self) -> None:
        """Adds a protein to the global variable and closes the dialog."""
        self.return_value.emit((self._view.ui.le_seq_name.text(), self._view.ui.le_protein_seq.toPlainText()))
        self._view.close()
