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
"""Init file for PyMOL interface."""
import sys
sys.path.append("C:\\ProgramData\\IBCI\\PySSA\\bin\\PySSA")

from src.pyssa_pymol import user_pymol_interface

# This global reference is needed to avoid garbage collection due to reference counting
mainInterface = None


def start_user_pymol_interface():
    """Function to start the PyMOL interface, by instantiating the Interface class."""
    global mainInterface
    mainInterface = user_pymol_interface.UserPyMOLInterface()


# Starting the actual interface
start_user_pymol_interface()
