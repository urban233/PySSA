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
"""This module contains helper functions for the pymol session management."""
import pymol
from PyQt5.QtWidgets import QMessageBox
from pymol import cmd
from pyssa.gui.ui.messageboxes import basic_boxes


def check_if_sele_is_empty() -> bool:
    """Checks if a selection is empty."""
    try:
        tmp_selection = cmd.get_model("sele")
    except pymol.CmdException:
        # gets thrown if no sele object exists in pymol
        basic_boxes.ok(
            "PyMOL session",
            "Please select at least one residue from the sequence view.",
            QMessageBox.Information,
        )
        return True
    try:
        tmp_selection.atom[0].resi
    except IndexError:
        # gets thrown if sele object is empty
        basic_boxes.ok(
            "PyMOL session",
            "Please select at least one residue from the sequence view.",
            QMessageBox.Information,
        )
        return True
    return False
