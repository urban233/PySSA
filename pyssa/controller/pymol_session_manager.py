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
"""Module for the PyMOL session manager class."""
import os.path
import pathlib

from pymol import cmd

from pyssa.controller import interface_manager
from pyssa.gui.ui.views import main_view
from pyssa.internal.data_structures import protein, protein_pair
from pyssa.internal.data_structures.data_classes import database_operation
from pyssa.internal.portal import pymol_io
from pyssa.io_pyssa import binary_data
from pyssa.util import enums, constants, exception


class PymolSessionManager:

    session_name: str
    session_object_type: str
    session_objects: list
    current_scene_name: str
    all_scenes: list[str]

    frozen_scene_name: str
    frozen_session_base64_string_filepath: pathlib.Path
    frozen_protein_object: "protein.Protein"
    frozen_protein_pair_object: "protein_pair.ProteinPair"

    def __init__(self) -> None:
        self.session_name = ""
        self.session_object_type = ""
        self.session_objects = []
        self.current_scene_name: str = ""
        self.all_scenes: list[str] = []

        self.frozen_protein_object = None
        self.frozen_protein_pair_object = None

    def reinitialize_session(self) -> None:
        """Reinitialize the pymol session and class attributes."""
        # reset class attributes
        self.session_name = ""
        self.session_object_type = ""
        self.session_objects: list = []
        # reset actual pymol session
        cmd.reinitialize()

    def freeze_current_protein_pymol_session(self,
                                             a_protein: "protein.Protein") -> "database_operation.DatabaseOperation":
        if not self.is_the_current_session_empty() and self.session_object_type == "protein":
            self.frozen_scene_name = self.current_scene_name
            self.frozen_protein_object = a_protein

            tmp_database_operation = database_operation.DatabaseOperation(
                enums.SQLQueryType.UPDATE_PYMOL_SESSION_PROTEIN,
                (0, self.frozen_protein_object)
            )
            return tmp_database_operation

    def unfreeze_current_protein_pymol_session(self):
        if self.frozen_protein_object is not None:
            self.load_protein_session(self.frozen_protein_object)
            self.frozen_protein_object = None
            self.frozen_scene_name = ""

    def freeze_current_protein_pair_pymol_session(
            self, a_protein_pair: "protein_pair.ProteinPair"
    ) -> "database_operation.DatabaseOperation":
        if not self.is_the_current_session_empty() and self.session_object_type == "protein_pair":
            self.frozen_scene_name = self.current_scene_name
            self.frozen_protein_pair_object = a_protein_pair
            tmp_database_operation = database_operation.DatabaseOperation(
                enums.SQLQueryType.UPDATE_PYMOL_SESSION_PROTEIN_PAIR,
                (0, self.frozen_protein_pair_object)
            )
            return tmp_database_operation

    def unfreeze_current_protein_pair_pymol_session(self):
        if self.frozen_protein_pair_object is not None:
            self.load_protein_pair_session(self.frozen_protein_pair_object)
            self.frozen_protein_pair_object = None
            self.frozen_scene_name = ""

    def load_protein_session(self, a_protein: "protein.Protein"):
        """Loads a pymol session of a single protein."""
        self.session_name = a_protein.get_molecule_object()
        self.session_object_type = "protein"
        self.session_objects = [a_protein]
        a_protein.load_protein_pymol_session()

        # <editor-fold desc="Integrity check">
        if not self._check_session_integrity(self.session_name):
            raise exception.ProteinNotFoundInPyMOLSession(f"Loading the PyMOL session failed, because the protein {self.session_name} can not be found in the PyMOL object list.")

        # </editor-fold>

    def load_protein_pair_session(self, a_protein_pair: "protein_pair.ProteinPair"):
        """Loads a pymol session of a protein pair."""
        self.session_name = a_protein_pair.name
        self.session_object_type = "protein_pair"
        self.session_objects = [a_protein_pair]
        a_protein_pair.load_pymol_session()

        # <editor-fold desc="Integrity check">
        if a_protein_pair.protein_1.get_molecule_object() == a_protein_pair.protein_2.get_molecule_object():
            if not self._check_session_integrity(f"{a_protein_pair.protein_1.get_molecule_object()}_1"):
                raise exception.ProteinNotFoundInPyMOLSession(
                    f"Loading the PyMOL session failed, because the protein {a_protein_pair.protein_1.get_molecule_object()}_1 can not be found in the PyMOL object list.")
            if not self._check_session_integrity(f"{a_protein_pair.protein_2.get_molecule_object()}_2"):
                raise exception.ProteinNotFoundInPyMOLSession(
                    f"Loading the PyMOL session failed, because the protein {a_protein_pair.protein_1.get_molecule_object()}_1 can not be found in the PyMOL object list.")
        else:
            if not self._check_session_integrity(a_protein_pair.protein_1.get_molecule_object()):
                raise exception.ProteinNotFoundInPyMOLSession(f"Loading the PyMOL session failed, because the protein {a_protein_pair.protein_1.get_molecule_object()} can not be found in the PyMOL object list.")
            if not self._check_session_integrity(a_protein_pair.protein_2.get_molecule_object()):
                raise exception.ProteinNotFoundInPyMOLSession(f"Loading the PyMOL session failed, because the protein {a_protein_pair.protein_2.get_molecule_object()} can not be found in the PyMOL object list.")

        # </editor-fold>

    def load_scene(self, a_scene_name):
        cmd.scene(a_scene_name, "recall")

    def load_current_scene(self):
        cmd.scene(self.current_scene_name, "recall")

    def save_current_pymol_session_as_pse_cache_file(self):
        tmp_session_path = pathlib.Path(
            f"{constants.CACHE_PYMOL_SESSION_DIR}/{self.session_name}.pse",
        )
        tmp_session_info = pymol_io.convert_pymol_session_to_base64_string(self.session_name)
        binary_data.write_binary_file_from_base64_string(tmp_session_path, tmp_session_info)
        return tmp_session_path

    def save_current_pymol_session_as_base64_string(self) -> str:
        return pymol_io.convert_pymol_session_to_base64_string("temp_active")

    def is_the_current_session_empty(self) -> bool:
        """Checks if the manager is in an empty session state."""
        if self.session_name == "" and self.session_object_type == "" and self.session_objects == []:
            return True
        else:
            return False

    def is_the_current_protein_in_session(self, the_name_of_the_selected_protein) -> bool:
        """Checks if the current protein is in the session."""
        if self.session_object_type == "protein" and self.session_name == the_name_of_the_selected_protein:
            return True
        else:
            return False

    def is_the_current_protein_pair_in_session(self, the_name_of_the_selected_protein_pair) -> bool:
        """Checks if the current protein pair is in the session."""
        if self.session_object_type == "protein_pair" and self.session_name == the_name_of_the_selected_protein_pair:
            return True
        else:
            return False

    def get_all_scenes_in_current_session(self):
        self.all_scenes.clear()
        self.all_scenes = cmd.get_scene_list()

    def set_all_scenes_for_current_session(self, all_scenes):
        self.all_scenes.clear()
        self.all_scenes = all_scenes

    @staticmethod
    def _check_session_integrity(a_protein_name) -> bool:
        """Checks if the current session is consistent with the manager."""
        tmp_pymol_objects = cmd.get_names()
        if a_protein_name in tmp_pymol_objects:
            return True
        else:
            return False

    def show_sequence_view(self) -> None:
        cmd.set("seq_view", 1)

    def hide_sequence_view(self) -> None:
        cmd.set("seq_view", 0)
