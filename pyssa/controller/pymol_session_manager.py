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

from application_process import application_process_manager
from pyssa.gui.ui.custom_dialogs import custom_message_box
from pyssa.internal.data_structures import protein, protein_pair
from pyssa.io_pyssa import binary_data, path_util
from pyssa.util import constants, exception, protein_pair_util
from pyssa_pymol import pymol_interface


class PymolSessionManager:

    session_name: str
    session_object_type: str
    session_objects: list
    current_scene_name: str
    all_scenes: list[str]

    def __init__(self, the_app_process_manager: "application_process_manager.ApplicationProcessManager") -> None:
        self.pymol_interface: "pymol_interface.PyMOLInterface" = pymol_interface.PyMOLInterface()
        self.session_name = ""
        self.session_object_type = ""
        self.session_objects = []
        self.current_scene_name: str = ""
        self.all_scenes: list[str] = []
        self._app_process_manager = the_app_process_manager

    # <editor-fold desc="Private methods">
    def _check_session_integrity(self, a_protein_name) -> bool:
        """Checks if the current session is consistent with the manager."""
        tmp_result = self.pymol_interface.get_all_object_names()
        if tmp_result["success"]:
            if a_protein_name in tmp_result["data"]:
                return True
            else:
                return False
        else:
            print(tmp_result["message"])
            return False

    def _load_pymol_session(self, a_pymol_session: str) -> None:
        """Loads a pymol session based on the given base64 data."""
        if self.session_name != "":
            tmp_session_path = pathlib.Path(
                f"{constants.CACHE_PYMOL_SESSION_DIR}/session_of_{self.session_name}.pse",
            )
        else:
            tmp_session_path = pathlib.Path(
                f"{constants.CACHE_PYMOL_SESSION_DIR}/temp_session_file.pse",
            )
        binary_data.write_binary_file_from_base64_string(tmp_session_path, a_pymol_session)
        if self._app_process_manager.pymol_crashed() is True:
            pass
        else:
            tmp_result = self.pymol_interface.load_pymol_session(tmp_session_path)
            if tmp_result["success"]:
                print("Session loaded.")
            else:
                print(tmp_result["message"])
                print("Session loaded failed!")

    def _convert_pymol_session_to_base64_string(self, pymol_molecule_object: str) -> str:
        """This function converts a pymol session file into a base64 string.

        Args:
            pymol_molecule_object (str): PyMOL molecule object to be converted.
        """
        session_filepath = pathlib.Path(f"{constants.SCRATCH_DIR}/{pymol_molecule_object}_session.pse")
        tmp_result = self.pymol_interface.save_pymol_session(session_filepath)
        if tmp_result["success"]:
            base64_string = binary_data.create_base64_string_from_file(session_filepath)
            os.remove(session_filepath)
            return base64_string
        else:
            print(tmp_result["message"])
            return ""
    # </editor-fold>

    # <editor-fold desc="Public methods">
    # <editor-fold desc="Non-cmd methods">
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
    # </editor-fold>

    # <editor-fold desc="Cmd-depended methods">
    def reinitialize_session(self) -> None:
        """Reinitialize the pymol session and class attributes."""
        # reset class attributes
        self.session_name = ""
        self.session_object_type = ""
        self.session_objects: list = []
        # reset actual pymol session
        self.pymol_interface.reinitialize_session()

    def load_protein_session(self, a_protein: "protein.Protein"):
        """Loads a pymol session of a single protein."""
        self.session_name = a_protein.get_molecule_object()
        self.session_object_type = "protein"
        self.session_objects = [a_protein]
        self._load_pymol_session(a_protein.pymol_session)

        # <editor-fold desc="Integrity check">
        if not self._check_session_integrity(a_protein.get_molecule_object()):
            raise exception.ProteinNotFoundInPyMOLSession(f"Loading the PyMOL session failed, because the protein {self.session_name} can not be found in the PyMOL object list.")

        # </editor-fold>

    def load_protein_pair_session(self, a_protein_pair: "protein_pair.ProteinPair"):
        """Loads a pymol session of a protein pair."""
        self.session_name = a_protein_pair.name
        self.session_object_type = "protein_pair"
        self.session_objects = [a_protein_pair]
        self._load_pymol_session(a_protein_pair.pymol_session)

        # <editor-fold desc="Integrity check">
        if a_protein_pair.protein_1.get_molecule_object() == a_protein_pair.protein_2.get_molecule_object():
            if not self._check_session_integrity(f"{a_protein_pair.protein_1.get_molecule_object()}_1"):
                raise exception.ProteinNotFoundInPyMOLSession(
                    f"Loading the PyMOL session failed, because the protein {a_protein_pair.protein_1.get_molecule_object()}_1 can not be found in the PyMOL object list.")
            if not self._check_session_integrity(f"{a_protein_pair.protein_2.get_molecule_object()}_2"):
                raise exception.ProteinNotFoundInPyMOLSession(
                    f"Loading the PyMOL session failed, because the protein {a_protein_pair.protein_2.get_molecule_object()}_2 can not be found in the PyMOL object list.")
        else:
            if not self._check_session_integrity(a_protein_pair.protein_1.get_molecule_object()):
                raise exception.ProteinNotFoundInPyMOLSession(f"Loading the PyMOL session failed, because the protein {a_protein_pair.protein_1.get_molecule_object()} can not be found in the PyMOL object list.")
            if not self._check_session_integrity(a_protein_pair.protein_2.get_molecule_object()):
                raise exception.ProteinNotFoundInPyMOLSession(f"Loading the PyMOL session failed, because the protein {a_protein_pair.protein_2.get_molecule_object()} can not be found in the PyMOL object list.")

        # </editor-fold>

    def load_scene(self, a_scene_name):
        tmp_result = self.pymol_interface.load_scene(a_scene_name)
        return tmp_result["success"]

    def load_current_scene(self):
        tmp_result = self.pymol_interface.load_scene(self.current_scene_name)
        return tmp_result["success"]

    def save_current_pymol_session_as_pse_cache_file(self):
        tmp_session_path = pathlib.Path(
            f"{constants.CACHE_PYMOL_SESSION_DIR}/{self.session_name}.pse",
        )
        tmp_session_base64 = self._convert_pymol_session_to_base64_string(self.session_name)
        binary_data.write_binary_file_from_base64_string(tmp_session_path, tmp_session_base64)
        return tmp_session_path

    def save_current_session_as_base64(self):
        tmp_session_filepath = pathlib.Path(f"{constants.CACHE_PYMOL_SESSION_DIR}/session_of_{self.session_name}.pse")
        tmp_session_base64 = self._convert_pymol_session_to_base64_string(self.session_name)
        binary_data.write_binary_file_from_base64_string(tmp_session_filepath, tmp_session_base64)
        return binary_data.create_base64_string_from_file(tmp_session_filepath)
        # tmp_result = self.pymol_interface.save_pymol_session(str(session_filepath))
        # if tmp_result["success"]:
        #     base64_string = binary_data.create_base64_string_from_file(path_util.FilePath(session_filepath))
        #     #os.remove(session_filepath)
        #     return base64_string
        # else:
        #     print(tmp_result["message"])
        #     return ""

    def get_all_scenes_in_current_session(self):
        self.all_scenes.clear()
        tmp_result = self.pymol_interface.get_scene_list()
        if tmp_result["success"]:
            self.all_scenes = tmp_result["data"]
        else:
            self.all_scenes = []

    def set_all_scenes_for_current_session(self, all_scenes):
        self.all_scenes.clear()
        self.all_scenes = all_scenes

    def show_sequence_view(self) -> None:
        self.pymol_interface.set_custom_setting("seq_view", 1)

    def hide_sequence_view(self) -> None:
        self.pymol_interface.set_custom_setting("seq_view", 0)
    # </editor-fold>

    def show_specific_representation(self, a_representation, a_selection_string):
        self.pymol_interface.show_custom_representation(a_representation, a_selection_string)

    def hide_specific_representation(self, a_representation, a_selection_string):
        self.pymol_interface.hide_custom_representation(a_representation, a_selection_string)

    def get_residue_colors(self, a_selection_string: str):
        tmp_result = self.pymol_interface.get_residue_colors(a_selection_string)
        if tmp_result["success"]:
            return tmp_result["data"]
        return None

    def get_chain_color(self, a_selection_string: str, chain_letter: str):
        tmp_result = self.pymol_interface.get_chain_color(a_selection_string, chain_letter)
        print(tmp_result)
        if tmp_result["success"]:
            return tmp_result["data"]
        return None

    def get_chain_repr_state(self, a_selection_string: str, chain_letter: str) -> dict:
        tmp_result = self.pymol_interface.get_chain_repr_state(a_selection_string, chain_letter)
        if tmp_result["success"]:
            return tmp_result["data"]
        return None

    def show_protein_selection_as_balls_and_sticks(self, selection: str) -> None:
        self.pymol_interface.show_custom_representation("sticks", selection)

    def hide_protein_selection_as_balls_and_sticks(self, selection: str) -> None:
        self.pymol_interface.hide_custom_representation("sticks", selection)

    def zoom_to_residue_in_protein_position(self, selection: str) -> None:
        self.pymol_interface.zoom_with_custom_parameters(selection)

    def color_protein(self, pymol_color: str, a_selection_string: str) -> None:
        """Colors a specific protein selection with a given PyMOL color.

        Args:
            pymol_color: a color which is available in PyMOL
            a_selection_string: a PyMOL conform selection string

        """
        if pymol_color == "":
            return
        if pymol_color not in constants.PYMOL_COLORS:
            raise ValueError(f"An illegal color argument. {pymol_color}")
        self.pymol_interface.color_selection(pymol_color, a_selection_string)

    def color_protein_pair_by_rmsd(self, a_protein_pair: "protein_pair.ProteinPair") -> None:
        """Colors a specific protein pair based on their rmsd value.

        Args:
            a_protein_pair: the protein pair to color.
        """
        cutoff_1 = 0.5
        cutoff_2 = 1.0
        cutoff_3 = 2
        cutoff_4 = 4
        cutoff_5 = 6

        color_1 = "br0"
        color_2 = "br2"
        color_3 = "br4"
        color_4 = "br6"
        color_5 = "br8"
        color_6 = "red"

        self.pymol_interface.color_selection("hydrogen", a_protein_pair.protein_2.get_molecule_object())

        i: int = 0
        for distance_value in a_protein_pair.distance_analysis.analysis_results.distance_data.get("distance"):
            if distance_value <= cutoff_1:
                atom_info = protein_pair_util.get_chain_and_position(
                    a_protein_pair.distance_analysis.analysis_results.distance_data,
                    i,
                )
                # create two atoms for the get_distance command
                atom1: str = (
                    f"/{a_protein_pair.protein_1.get_molecule_object()}//" f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                )
                atom2: str = (
                    f"/{a_protein_pair.protein_2.get_molecule_object()}//" f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                )
                # coloring
                self.pymol_interface.color_selection(color_1, atom1)
                self.pymol_interface.color_selection(color_1, atom2)
                i += 1

            elif distance_value <= cutoff_2:
                atom_info = protein_pair_util.get_chain_and_position(
                    a_protein_pair.distance_analysis.analysis_results.distance_data,
                    i,
                )
                # create two atoms for the get_distance command
                atom1: str = (
                    f"/{a_protein_pair.protein_1.get_molecule_object()}//" f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                )
                atom2: str = (
                    f"/{a_protein_pair.protein_2.get_molecule_object()}//" f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                )
                # coloring
                self.pymol_interface.color_selection(color_2, atom1)
                self.pymol_interface.color_selection(color_2, atom2)
                i += 1

            elif distance_value <= cutoff_3:
                atom_info = protein_pair_util.get_chain_and_position(
                    a_protein_pair.distance_analysis.analysis_results.distance_data,
                    i,
                )
                # create two atoms for the get_distance command
                atom1: str = (
                    f"/{a_protein_pair.protein_1.get_molecule_object()}//"
                    f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}/CA"
                )
                atom2: str = (
                    f"/{a_protein_pair.protein_2.get_molecule_object()}//"
                    f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}/CA"
                )
                # coloring
                self.pymol_interface.color_selection(color_3, atom1)
                self.pymol_interface.color_selection(color_3, atom2)
                i += 1

            elif distance_value <= cutoff_4:
                atom_info = protein_pair_util.get_chain_and_position(
                    a_protein_pair.distance_analysis.analysis_results.distance_data,
                    i,
                )
                # create two atoms for the get_distance command
                atom1: str = (
                    f"/{a_protein_pair.protein_1.get_molecule_object()}//" f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                )
                atom2: str = (
                    f"/{a_protein_pair.protein_2.get_molecule_object()}//" f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                )
                # coloring
                self.pymol_interface.color_selection(color_4, atom1)
                self.pymol_interface.color_selection(color_4, atom2)
                i += 1

            elif distance_value <= cutoff_5:
                atom_info = protein_pair_util.get_chain_and_position(
                    a_protein_pair.distance_analysis.analysis_results.distance_data,
                    i,
                )
                # create two atoms for the get_distance command
                atom1: str = (
                    f"/{a_protein_pair.protein_1.get_molecule_object()}//" f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                )
                atom2: str = (
                    f"/{a_protein_pair.protein_2.get_molecule_object()}//" f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                )
                # coloring
                self.pymol_interface.color_selection(color_5, atom1)
                self.pymol_interface.color_selection(color_5, atom2)
                i += 1

            elif distance_value > cutoff_5:
                atom_info = protein_pair_util.get_chain_and_position(
                    a_protein_pair.distance_analysis.analysis_results.distance_data,
                    i,
                )
                # create two atoms for the get_distance command
                atom1: str = (
                    f"/{a_protein_pair.protein_1.get_molecule_object()}//" f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                )
                atom2: str = (
                    f"/{a_protein_pair.protein_2.get_molecule_object()}//" f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                )
                # coloring
                self.pymol_interface.color_selection(color_6, f"({atom1})")
                self.pymol_interface.color_selection(color_6, f"({atom2})")
                i += 1

    def setup_default_session_graphic_settings(self) -> None:
        """This functions modifies the pymol session to look fancy."""
        self.pymol_interface.set_background_color(constants.PYMOL_DEFAULT_BACKGROUND_COLOR)
        self.pymol_interface.set_default_graphic_settings()

    def setup_default_image_graphic_settings(self, ray_shadows: bool, opaque_background: int = 0) -> None:
        """Sets up the default image graphic settings for PyMOL.

        Args:
            ray_shadows (bool): false if no shadows, true if shadows should be displayed.
            opaque_background (int, optional): 0 for a transparent background and 1 for a white background.
        """
        if not ray_shadows:
            opt_ray_shadows: str = "off"
        else:
            opt_ray_shadows: str = "on"
        self.pymol_interface.set_background_color(constants.PYMOL_DEFAULT_BACKGROUND_COLOR)
        self.pymol_interface.set_custom_setting("ray_trace_mode", constants.PYMOL_DEFAULT_RAY_TRACE_MODE)
        self.pymol_interface.set_custom_setting("antialias", constants.PYMOL_DEFAULT_ANTIALIAS)
        self.pymol_interface.set_custom_setting("ray_shadows", opt_ray_shadows)
        self.pymol_interface.set_custom_setting("ray_opaque_background", opaque_background)

    def setup_default_graphic_settings_for_interesting_regions(self) -> None:
        """Sets up the default graphic settings for interesting regions."""
        self.pymol_interface.set_background_color(constants.PYMOL_DEFAULT_BACKGROUND_COLOR)
        self.pymol_interface.set_custom_setting("label_size", 14)
        self.pymol_interface.set_custom_setting("label_font_id", 13)
        self.pymol_interface.set_custom_setting("label_color", "hotpink")
        self.pymol_interface.set_custom_setting("depth_cue", 0)
        # interacts directly with molecule objects in the session
        self.pymol_interface.hide_custom_representation("cartoon", "all")
        self.pymol_interface.show_custom_representation("ribbon", "all")

    def check_if_sele_is_empty(self) -> bool:
        """Checks if a selection is empty."""
        tmp_result = self.pymol_interface.get_model("sele")
        if tmp_result["success"]:
            tmp_selection = tmp_result["data"]
        else:
            # gets thrown if no sele object exists in pymol
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Please select at least one residue from the sequence view.",
                "PyMOL Selection",
                custom_message_box.CustomMessageBoxIcons.INFORMATION.value
            )
            tmp_dialog.exec_()
            return True
        try:
            tmp_selection.atom[0].resi
        except IndexError:
            # gets thrown if sele object is empty
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "Please select at least one residue from the sequence view.",
                "PyMOL Selection",
                custom_message_box.CustomMessageBoxIcons.INFORMATION.value
            )
            tmp_dialog.exec_()
            return True
        return False
    # </editor-fold>
