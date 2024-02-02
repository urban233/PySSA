from pymol import cmd

from pyssa.controller import interface_manager
from pyssa.gui.ui.views import main_view
from pyssa.internal.data_structures import protein, protein_pair


class PymolSessionManager:

    session_name: str
    session_object_type: str
    session_objects: list

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager") -> None:
        self.session_name = ""
        self.session_object_type = ""
        self.session_objects = []
        self._interface_manager: "interface_manager.InterfaceManager" = the_interface_manager

    def reinitialize_session(self) -> None:
        """Reinitialize the pymol session and class attributes."""
        # reset class attributes
        self.session_name = ""
        self.session_object_type = ""
        self.session_objects: list = []
        # reset actual pymol session
        cmd.reinitialize()

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

    def is_the_current_session_empty(self) -> bool:
        """Checks if the manager is in an empty session state."""
        if self.session_name == "" and self.session_object_type == "" and self.session_objects == []:
            return True
        else:
            return False

    def is_the_current_protein_in_session(self) -> bool:
        """Checks if the current protein is in the session."""
        if self._interface_manager.get_current_protein_tree_index_type() == "protein":
            tmp_protein_name: str = self._interface_manager.get_current_protein_tree_index_object().get_molecule_object()
        elif self._interface_manager.get_current_protein_tree_index_type() == "chain":
            tmp_protein_name: str = self._interface_manager.get_parent_index_object_of_current_protein_tree_index().get_molecule_object()
        else:
            raise ValueError("Unknown type!")

        if self.session_object_type == "protein" and self.session_name == tmp_protein_name:
            return True
        else:
            return False

    def is_the_current_protein_pair_in_session(self) -> bool:
        """Checks if the current protein pair is in the session."""
        if self._interface_manager.get_current_protein_pair_tree_index_object() == "protein":
            tmp_protein_pair_name: str = self._interface_manager.get_current_protein_pair_tree_index().get_molecule_object()
        elif self._interface_manager.get_current_protein_pair_tree_index_type() == "chain":
            tmp_protein_pair_name: str = self._interface_manager.get_parent_index_object_of_current_protein_pair_tree_index().get_molecule_object()
        else:
            raise ValueError("Unknown type!")

        if self.session_object_type == "protein" and self.session_name == tmp_protein_pair_name:
            return True
        else:
            return False

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
