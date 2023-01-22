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
"""Module for graphic operations in pymol"""
import pymol
from pymol import cmd
from pyssa.internal.portal import pymol_safeguard
from pyssa.internal.data_structures import protein

def show_protein_selection_as_balls_and_sticks(protein_obj: protein.Protein):
    """This function show the protein as balls and sticks in representation mode.

    """
    if not pymol_safeguard.PymolSafeguard.check_if_protein_in_session():
        raise pymol.CmdException("No protein is in pymol session.")
    try:
        cmd.show(representation="sticks", selection=protein_obj.selection)
    except pymol.CmdException:
        print("No stick and balls ca be shown in protein.")

# TODO: function for hide and zoom
