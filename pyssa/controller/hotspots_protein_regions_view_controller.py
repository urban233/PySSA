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
"""Module for the Hotspots Dialog."""

import glob
import os
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from pymol import cmd

from pyssa.controller import interface_manager
from pyssa.gui.ui.styles import styles
from pyssa.util import input_validator, constants, session_util


class OpenProjectViewController(QtCore.QObject):
    """Class for the Hotspots Protein Regions View Controller."""
    return_value = QtCore.pyqtSignal(str)

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager") -> None:
        super().__init__()
        self._interface_manager = the_interface_manager
        self._view = the_interface_manager.get_hotspots_protein_regions_view()
        self._connect_all_ui_elements_to_slot_functions()

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.btn_sticks_hide.clicked.connect(self.hide_resi_sticks)
        self._view.ui.btn_sticks_show.clicked.connect(self.show_resi_sticks)
        self._view.ui.btn_position_zoom.clicked.connect(self.zoom_resi_position)
        self._view.ui.cb_disulfide_bonds.clicked.connect(self.show_hide_disulfide_bonds)
        # self._view.ui.btn_info.clicked.connect()

    def hide_resi_sticks(self) -> None:
        """Hides the balls and sticks representation of the pymol selection."""
        session_util.check_if_sele_is_empty()
        cmd.hide(representation="sticks", selection="sele")

    def show_resi_sticks(self) -> None:
        """Shows the pymol selection as sticks."""
        session_util.check_if_sele_is_empty()
        cmd.show(representation="sticks", selection="sele and not hydrogens")
        cmd.select(name="sele", selection="sele and not hydrogens")
        cmd.color(color="atomic", selection="sele and not elem C")
        cmd.set("valence", 0)  # this needs to be better implemented

    def show_hide_disulfide_bonds(self) -> None:
        """Shows and hide all disulfid bonds within the pymol session."""
        # hide disulfid bonds
        if self._view.ui.cb_disulfide_bonds.currentText("hide"):
            tmp_pymol_selection_option: str = "byres (resn CYS and name SG) within 2 of (resn CYS and name SG)"
            cmd.select(
                name="disulfides",
                selection=f"{self._view.ui.box_manage_choose_protein.currentText()} & {tmp_pymol_selection_option}",
            )
            cmd.hide("sticks", "disulfides")
        else:
            # show disulfid bonds
            tmp_pymol_selection_option: str = "byres (resn CYS and name SG) within 2 of (resn CYS and name SG)"
            cmd.select(
                name="disulfides",
                selection=f"{self._view.ui.box_manage_choose_protein.currentText()} & {tmp_pymol_selection_option}",
            )
            cmd.color(color="atomic", selection="disulfides and not elem C")
            cmd.set("valence", 0)  # this needs to be better implemented
            cmd.show("sticks", "disulfides")
            cmd.hide("sticks", "elem H")

    def zoom_resi_position(self) -> None:
        """Zooms to the pymol selection."""
        session_util.check_if_sele_is_empty()
        cmd.zoom(selection="sele", buffer=8.0, state=0, complete=0)
