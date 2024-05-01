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

from pyssa.controller import pymol_session_manager
from pyssa.logging_pyssa import log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def preview_image(the_pymol_session_manager: "pymol_session_manager.PymolSessionManager", a_placeholder_2: int) -> tuple:
    # TODO: the renderer should be changeable
    the_pymol_session_manager.pymol_interface.ray(a_width=2400, a_height=2400, a_renderer=int(0))
    return 0, ""


def create_drawn_image(an_image_filepath: str, the_pymol_session_manager: "pymol_session_manager.PymolSessionManager") -> tuple:
    the_pymol_session_manager.pymol_interface.draw(a_width=2400, a_height=2400, an_antialias_value=2)
    the_pymol_session_manager.pymol_interface.png(an_image_filepath, a_dpi_value=300)
    return 0, ""
