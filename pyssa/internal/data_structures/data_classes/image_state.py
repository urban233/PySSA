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
from dataclasses import dataclass


@dataclass
class ImageState:
    """Contains the current state of the image page."""

    representation: str
    background_color: str
    renderer: int
    transparent_background: bool
    ray_tracing: bool
    ray_trace_mode: int
    ray_texture: int

    def __init__(self) -> None:
        """Empty constructor."""
        self.representation: str = ""
        self.background_color: str = ""
        self.renderer: int = 0
        self.transparent_background: bool = False
        self.ray_tracing: bool = False
        self.ray_trace_mode: int = 0
        self.ray_texture: int = 0
