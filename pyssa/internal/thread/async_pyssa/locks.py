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
"""Module for custom lock implementations."""

class PyMOL_LOCK:
    _state: bool

    def __init__(self):
        self._state: bool = False

    def lock(self):
        """Sets the state to True."""
        self._state = True

    def unlock(self):
        """Sets the state to False."""
        self._state = False

    def is_locked(self):
        """Checks if the lock is locked.

        Notes:
            The lock is locked if the state is True.
        """
        return self._state
