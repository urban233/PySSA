#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/zielesny/PySSA>
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
"""Module contains the CurrentSession dataclass."""
from dataclasses import dataclass


@dataclass
class CurrentSession:
    """Class which holds information about the current pymol session.

    This is a dataclass with:
        type (str)
        name (str)
        session (str)

    DO NOT USE THIS CLASS ANYMORE IN PRODUCTION!
    """

    type: str  # noqa: A003
    name: str
    session: str


class CurrentPymolSession:
    """Class which holds information about the current pymol session"""

    session_name: str
    object_type: str

    def __init__(self, a_session_name: str, an_object_type) -> None:
        self.session_name = a_session_name
        self.object_type = an_object_type
