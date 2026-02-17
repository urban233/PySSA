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

from pmg_qt.pymol_gl_widget import PyMOLGLWidget

__docformat__ = "google"


class UserPyMOL:
    """Class for manging the interaction of the user with the embedded PyMOL viewer widget."""

    def __init__(
        self, a_pymol_widget: PyMOLGLWidget
    ) -> None:  # Type-annotation does not work in this case!
        """Constructor.

        Args:
          a_pymol_widget: Instance of the PyMOL OpenGL widget used in the main window.

        Note:
          The OpenGL widget is necessary as an argument because this makes the underlying
          cmd module available to the class in a way that the IDE understands.
        """
        # <editor-fold desc="Instance attributes">
        self._embedded_cmd = a_pymol_widget.cmd
        """Embedded instance of the cmd module"""
        # </editor-fold>
        self._development_setup()

    def _development_setup(self) -> None:
        """Sets up the PyMOL viewer widget for development."""
        self._embedded_cmd.set("seq_view", 1)
        # self._embedded_cmd.fetch("1nb1")

    def get_cmd_module(self):
        """Getter for the embedded PyMOL cmd module."""
        return self._embedded_cmd
