#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
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
"""Module for all icon filepaths."""
import pathlib
from src.pyssa.util import constants


ICONS_PATH = pathlib.Path(constants.PROGRAM_BIN_ROOT_PATH, "assets", "icons")

TREE_RIGHT_ARROW = pathlib.Path(ICONS_PATH, "keyboard_arrow_right_w400.svg")
TREE_DOWN_ARROW = pathlib.Path(ICONS_PATH, "keyboard_arrow_down_w400.svg")
