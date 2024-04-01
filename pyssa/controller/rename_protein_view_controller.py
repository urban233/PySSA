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
"""Module for the Rename Protein Dialog."""
import glob
import logging
import os
from PyQt5 import QtCore
from PyQt5.QtCore import Qt

from pyssa.controller import interface_manager
from pyssa.gui.ui.styles import styles
from pyssa.util import input_validator, constants
from pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class RenameProteinViewController(QtCore.QObject):
    """Class for the RenameProteinViewController class"""
    user_input = QtCore.pyqtSignal(tuple)

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager"):
        super().__init__()
        self._interface_manager = the_interface_manager
        self._view = the_interface_manager.get_rename_protein_view()
        self._connect_all_ui_elements_to_slot_functions()

    def restore_ui(self):
        self._view.ui.le_name.clear()
        self._view.ui.lbl_status.setText("")

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.btn_rename.clicked.connect(self._rename_protein)

    def _rename_protein(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Rename' button was clicked.")
        self._view.close()
        self.user_input.emit((self._view.ui.le_name.text(), True))
