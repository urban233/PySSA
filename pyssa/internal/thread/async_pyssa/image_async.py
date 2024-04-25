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
"""Module for all asynchronous functions that are related to the image creation and preview for PyMOL."""
import logging

from pymol import cmd
from pyssa.logging_pyssa import log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def preview_image(a_placeholder_1: int, a_placeholder_2: int) -> tuple:
    # TODO: the renderer should be changeable
    try:
        cmd.ray(2400, 2400, renderer=int(0))
    except pymol.CmdException:
        logger.warning("Unexpected exception.")
    return 0, ""


def create_drawn_image(an_image_filepath: str, the_app_settings: "settings.Settings") -> tuple:
    cmd.bg_color(the_app_settings.image_background_color)
    try:
        cmd.draw(2400, 2400)
        cmd.png(an_image_filepath, dpi=300)
    except pymol.CmdException:
        logger.warning("Unexpected exception.")
    cmd.bg_color("black")
    return 0, ""
