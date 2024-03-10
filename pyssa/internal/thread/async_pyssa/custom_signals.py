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
"""Module for all custom signals which can be used as arguments either in the main thread or in separate threads."""
from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal


class ProgressSignal(QtCore.QObject):
    progress = pyqtSignal(tuple)

    def emit_signal(self, a_message, a_value):
        if a_value < 0 or a_value > 100:
            raise ValueError(f"Value must be between 0 and 100 but the given value us {a_value}!")
        self.progress.emit((a_message, a_value))
