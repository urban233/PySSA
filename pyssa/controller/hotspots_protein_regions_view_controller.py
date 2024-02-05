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
import typing

from PyQt5 import QtCore
from pymol import cmd

from pyssa.controller import interface_manager
from pyssa.util import session_util


class HotspotsProteinRegionsViewController(QtCore.QObject):
    """Class for the Hotspots Protein Regions View Controller."""
    return_value = QtCore.pyqtSignal(str)

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager") -> None:
        super().__init__()
        self._interface_manager = the_interface_manager
        self._view = the_interface_manager.get_hotspots_protein_regions_view()
        self._connect_all_ui_elements_to_slot_functions()
        self._protein_names: tuple[str, typing.Union[str, int]] = self._get_protein_names_from_tree_view()

    # TODO: Write a function, to intercept the protein from the protein pair, when a protein pair gets loading in the PyMOL session.
    def _get_protein_from_protein_pair(self):
        """Get the protein from the protein pair."""
        if self._interface_manager.get_current_protein_tree_index_type() == "protein_pair":
            return self._interface_manager.get_current_protein_tree_index_object()
        else:
            return self._interface_manager.get_current_protein_tree_index_object()

    def _get_protein_names_from_tree_view(self):
        """Checks if the object is a protein or a chain."""
        if self._interface_manager.get_current_protein_tree_index_type() == "protein":
            return self._interface_manager.get_current_protein_tree_index_object().get_molecule_object()
        elif:

        else:
            return self._interface_manager.get_parent_index_object_of_current_protein_tree_index().get_molecule_object()

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.btn_sticks_show.clicked.connect(self.show_resi_sticks)
        self._view.ui.btn_sticks_hide.clicked.connect(self.hide_resi_sticks)
        self._view.ui.btn_disulfide_bonds_show.clicked.connect(self.show_disulfide_bonds)
        self._view.ui.btn_disulfide_bonds_hide.clicked.connect(self.hide_disulfide_bonds)
        self._view.ui.btn_position_zoom.clicked.connect(self.zoom_resi_position)
        # self._view.ui.btn_help.clicked.connect()

    def show_resi_sticks(self) -> None:
        """Shows the pymol selection as sticks."""
        session_util.check_if_sele_is_empty()
        cmd.show(representation="sticks", selection="sele and not hydrogens")
        cmd.select(name="sele", selection="sele and not hydrogens")
        cmd.color(color="atomic", selection="sele and not elem C")
        cmd.set("valence", 0)  # this needs to be better implemented

    def hide_resi_sticks(self) -> None:
        """Hides the balls and sticks representation of the pymol selection."""
        session_util.check_if_sele_is_empty()
        cmd.hide(representation="sticks", selection="sele")

    def show_disulfide_bonds(self) -> None:
        """Shows all disulfid bonds within the pymol session."""
        tmp_pymol_selection_option: str = "byres (resn CYS and name SG) within 2 of (resn CYS and name SG)"
        cmd.select(
            name="disulfides",
            selection=f"{self._protein_names} & {tmp_pymol_selection_option}",
        )
        cmd.color(color="atomic", selection="disulfides and not elem C")
        cmd.set("valence", 0)  # this needs to be better implemented
        cmd.show("sticks", "disulfides")
        cmd.hide("sticks", "elem H")

    def hide_disulfide_bonds(self) -> None:
        """Hides all disulfid bonds within the pymol session."""
        tmp_pymol_selection_option: str = "byres (resn CYS and name SG) within 2 of (resn CYS and name SG)"
        cmd.select(
            name="disulfides",
            selection=f"{self._protein_names} & {tmp_pymol_selection_option}",
        )
        cmd.hide("sticks", "disulfides")

    def zoom_resi_position(self) -> None:
        """Zooms to the pymol selection."""
        session_util.check_if_sele_is_empty()
        cmd.zoom(selection="sele", buffer=8.0, state=0, complete=0)
