import os.path
import pathlib

from pymol import cmd

from pyssa.controller import interface_manager
from pyssa.gui.ui.views import main_view
from pyssa.internal.data_structures import protein, protein_pair
from pyssa.internal.data_structures.data_classes import database_operation
from pyssa.internal.portal import pymol_io
from pyssa.util import enums, constants


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

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager") -> None:
        self.session_name = ""
        self.session_object_type = ""
        self.session_objects = []
        self.current_scene_name: str = ""
        self.all_scenes: list[str] = []

        self.frozen_protein_object = None
        self.frozen_protein_pair_object = None

        self._interface_manager: "interface_manager.InterfaceManager" = the_interface_manager

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
            raise RuntimeError(f"Loading the PyMOL session failed, because the protein {self.session_name} can not be found in the PyMOL object list.")

        # </editor-fold>

    def load_protein_pair_session(self, a_protein_pair: "protein_pair.ProteinPair"):
        """Loads a pymol session of a protein pair."""
        self.session_name = a_protein_pair.name
        self.session_object_type = "protein_pair"
        self.session_objects = [a_protein_pair]
        a_protein_pair.load_pymol_session()

        # <editor-fold desc="Integrity check">
        if not self._check_session_integrity(a_protein_pair.protein_1.get_molecule_object()):
            raise RuntimeError(f"Loading the PyMOL session failed, because the protein {a_protein_pair.protein_1.get_molecule_object()} can not be found in the PyMOL object list.")
        if not self._check_session_integrity(a_protein_pair.protein_2.get_molecule_object()):
            raise RuntimeError(f"Loading the PyMOL session failed, because the protein {a_protein_pair.protein_2.get_molecule_object()} can not be found in the PyMOL object list.")

        # </editor-fold>

    def load_scene(self, a_scene_name):
        cmd.scene(a_scene_name, "recall")

    def save_current_pymol_session_as_base64_string(self) -> str:
        return pymol_io.convert_pymol_session_to_base64_string("temp_active")

    def is_the_current_session_empty(self) -> bool:
        """Checks if the manager is in an empty session state."""
        if self.session_name == "" and self.session_object_type == "" and self.session_objects == []:
            return True
        else:
            return False

    def is_the_current_protein_in_session(self) -> bool:
        """Checks if the current protein is in the session."""
        # if self._interface_manager.get_current_protein_tree_index_type() == "protein":
        #     tmp_protein_name: str = self._interface_manager.get_current_active_protein_object().get_molecule_object()
        # elif self._interface_manager.get_current_protein_tree_index_type() == "chain":
        #     tmp_protein_name: str = self._interface_manager.get_parent_index_object_of_current_protein_tree_index().get_molecule_object()
        # else:
        #     raise ValueError("Unknown type!")

        if self.session_object_type == "protein" and self.session_name == self._interface_manager.get_current_active_protein_object().get_molecule_object():
            return True
        else:
            return False

    def is_the_current_protein_pair_in_session(self) -> bool:
        """Checks if the current protein pair is in the session."""
        if self.session_object_type == "protein_pair" and self.session_name == self._interface_manager.get_current_active_protein_pair_object().name:
            return True
        else:
            return False

    def get_all_scenes_in_current_session(self):
        self.all_scenes.clear()
        self.all_scenes = cmd.get_scene_list()

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
