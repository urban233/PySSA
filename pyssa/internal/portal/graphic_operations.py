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
from pyssa.util import constants
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pyssa.internal.data_structures import protein


def show_protein_selection_as_balls_and_sticks(selection: str):
    """This function shows the protein as balls and sticks in representation mode.

    Args:
        protein_obj:

    """
    if not pymol_safeguard.PymolSafeguard.check_if_protein_in_session():
        raise pymol.CmdException("No protein is in pymol session.")
    try:
        cmd.show(representation="sticks", selection=selection)
    except pymol.CmdException:
        print("No sticks and balls can be shown in protein.")


def hide_protein_selection_as_balls_and_sticks(selection: str):
    """This function hides the protein from balls and sticks.

    Args:
        protein_obj:

    Returns:

    """
    if not pymol_safeguard.PymolSafeguard.check_if_protein_in_session():
        raise pymol.CmdException("No protein is in pymol session.")
    try:
        cmd.hide(representation="sticks", selection=selection)
    except pymol.CmdException:
        print("No sticks and balls can be hidden in protein.")


def zoom_to_residue_in_protein_position(selection: str):
    """This function zooms to the residue in protein.

    Args:
        protein_obj:

    Returns:

    """
    if not pymol_safeguard.PymolSafeguard.check_if_protein_in_session():
        raise pymol.CmdException("No protein is in pymol session.")
    try:
        cmd.zoom(selection=selection, buffer=8.0, state=0, complete=0)
    except pymol.CmdException:
        print("No residue can be shown in protein.")


def color_protein(pymol_color: str, a_selection_string: str):
    """This function colors a specific protein selection with a given PyMOL color.

    Args:
        pymol_color:
            a color which is available in PyMOL
        a_selection_string:
            a PyMOL conform selection string

    """
    if pymol_color not in constants.PYMOL_COLORS:
        raise ValueError("An illegal color argument.")
    if not pymol_safeguard.PymolSafeguard.check_if_protein_in_session():
        raise pymol.CmdException("No protein is in pymol session.")
    try:
        cmd.color(pymol_color, a_selection_string)
    except pymol.CmdException:
        print("Color process was unsuccessful.")


def setup_default_session_graphic_settings():
    cmd.bg_color(constants.PYMOL_DEFAULT_BACKGROUND_COLOR)
    cmd.set("ray_trace_mode", constants.PYMOL_DEFAULT_RAY_TRACE_MODE)
    cmd.set("antialias", constants.PYMOL_DEFAULT_ANTIALIAS)
    cmd.set("ambient", constants.PYMOL_DEFAULT_AMBIENT)
    cmd.set("cartoon_fancy_helices", constants.PYMOL_DEFAULT_FANCY_HELICES)
    cmd.set("cartoon_discrete_colors", constants.PYMOL_DEFAULT_CARTOON_DISCRETE_COLORS)
    cmd.set("cartoon_sampling", constants.PYMOL_DEFAULT_CARTOON_SAMPLING)
    cmd.set("spec_power", constants.PYMOL_DEFAULT_SPEC_POWER)
    cmd.set("spec_reflect", constants.PYMOL_DEFAULT_SPEC_REFLECT)
    cmd.set("ray_transparency_contrast", constants.PYMOL_DEFAULT_RAY_TRANSPARENCY_CONTRAST)
    cmd.set("ray_transparency_oblique", constants.PYMOL_DEFAULT_RAY_TRANSPARENCY_OBLIQUE)
    cmd.set("ray_transparency_oblique_power", constants.PYMOL_DEFAULT_RAY_OBLIQUE_POWER)
    cmd.set("ray_trace_color", constants.PYMOL_DEFAULT_RAY_COLOR)
    cmd.unset("depth_cue")
