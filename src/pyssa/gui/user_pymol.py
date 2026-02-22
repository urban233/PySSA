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
import pathlib
from typing import TYPE_CHECKING

from pmg_qt.pymol_gl_widget import PyMOLGLWidget
from src.pyssa.io_pyssa import binary_data
from src.pyssa.util import constants

if TYPE_CHECKING:
    from src.pyssa.internal.data_structures import protein_pair, protein

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
        self._current_object: "protein.Protein | protein_pair.ProteinPair | None" = None
        # </editor-fold>
        self._development_setup()

    def _development_setup(self) -> None:
        """Sets up the PyMOL viewer widget for development."""
        self._embedded_cmd.set("seq_view", 1)
        # self._embedded_cmd.fetch("1nb1")

    def get_cmd_module(self):
        """Getter for the embedded PyMOL cmd module."""
        return self._embedded_cmd

    def get_currently_loaded_object(self) -> "protein.Protein | protein_pair.ProteinPair | None":
        return self._current_object

    def load_session(self, a_pymol_session: str, a_current_object: "protein.Protein | protein_pair.ProteinPair") -> None:
        tmp_cache_dir = pathlib.Path(constants.CACHE_PYMOL_SESSION_DIR)
        tmp_cache_dir.mkdir(parents=True, exist_ok=True)
        tmp_session_path = tmp_cache_dir / f"load_session_cache.pse"
        binary_data.write_binary_file_from_base64_string(
            tmp_session_path, a_pymol_session
        )
        self._embedded_cmd.load(tmp_session_path)
        self._current_object = a_current_object

    def save_session(self) -> str:
        """Saves the current PyMOL session as a Base64-encoded string."""
        tmp_cache_dir = pathlib.Path(constants.CACHE_PYMOL_SESSION_DIR)
        tmp_cache_dir.mkdir(parents=True, exist_ok=True)
        tmp_session_path = tmp_cache_dir / f"save_session_cache.pse"
        self._embedded_cmd.save(tmp_session_path)
        return binary_data.create_base64_string_from_file(tmp_session_path)

    def reinitialize_session(self):
        self._embedded_cmd.reinitialize()
        self._current_object = None
