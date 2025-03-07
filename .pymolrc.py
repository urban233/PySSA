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
"""Init file for PyMOL interface. This exact file is only needed for development purposes!"""
import sys
import pathlib

ROOT_PATH: pathlib.Path = pathlib.Path(__file__).parent
"""The path is .../user_pymol/lib/pymol"""

PYSSA_LIB_PATH: pathlib.Path = pathlib.Path(ROOT_PATH.parent.parent.parent / "lib")
"""
This is the lib path of the main PySSA application and only valid if 
the user_pymol folder is placed under the program root directory!
"""

sys.path.append(str(PYSSA_LIB_PATH))

from src.pyssa_pymol import user_pymol_interface

# This global reference is needed to avoid garbage collection due to reference counting
mainInterface = None


def start_user_pymol_interface() -> None:
    """Function to start the PyMOL interface, by instantiating the Interface class."""
    global mainInterface
    mainInterface = user_pymol_interface.UserPyMOLInterface()


# Starting the actual interface
start_user_pymol_interface()
