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
from src.pyssa.util import constants, pyssa_keys

import numpy as np

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

    def color_protein_pair_by_rmsd(
            self, a_protein_pair: "protein_pair.ProteinPair"
    ) -> None:
        """Colors a specific protein pair based on their rmsd value.
    
        Args:
            a_protein_pair (protein_pair.ProteinPair): The protein pair to color.
    
        Raises:
            exception.IllegalArgumentError: If `a_protein_pair` is None.
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
        
        self._embedded_cmd.color("hydrogen", a_protein_pair.protein_2.get_molecule_object())
            
        i: int = 0
        for (
                distance_value
        ) in a_protein_pair.distance_analysis.analysis_results.distance_data.get(
            "distance"
        ):
            if distance_value <= cutoff_1:
                atom_info = _get_chain_and_position(
                    a_protein_pair.distance_analysis.analysis_results.distance_data,
                    i,
                )
                # create two atoms for the get_distance command
                atom1: str = (
                    f"/{a_protein_pair.protein_1.get_molecule_object()}//"
                    f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                )
                atom2: str = (
                    f"/{a_protein_pair.protein_2.get_molecule_object()}//"
                    f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                )
                # coloring
                self._embedded_cmd.color(color_1, atom1)
                self._embedded_cmd.color(color_1, atom2)
                i += 1
    
            elif distance_value <= cutoff_2:
                atom_info = _get_chain_and_position(
                    a_protein_pair.distance_analysis.analysis_results.distance_data,
                    i,
                )
                # create two atoms for the get_distance command
                atom1: str = (
                    f"/{a_protein_pair.protein_1.get_molecule_object()}//"
                    f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                )
                atom2: str = (
                    f"/{a_protein_pair.protein_2.get_molecule_object()}//"
                    f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                )
                # coloring
                self._embedded_cmd.color(color_2, atom1)
                self._embedded_cmd.color(color_2, atom2)
                i += 1
    
            elif distance_value <= cutoff_3:
                atom_info = _get_chain_and_position(
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
                self._embedded_cmd.color(color_3, atom1)
                self._embedded_cmd.color(color_3, atom2)
                i += 1
    
            elif distance_value <= cutoff_4:
                atom_info = _get_chain_and_position(
                    a_protein_pair.distance_analysis.analysis_results.distance_data,
                    i,
                )
                # create two atoms for the get_distance command
                atom1: str = (
                    f"/{a_protein_pair.protein_1.get_molecule_object()}//"
                    f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                )
                atom2: str = (
                    f"/{a_protein_pair.protein_2.get_molecule_object()}//"
                    f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                )
                # coloring
                self._embedded_cmd.color(color_4, atom1)
                self._embedded_cmd.color(color_4, atom2)
                i += 1
    
            elif distance_value <= cutoff_5:
                atom_info = _get_chain_and_position(
                    a_protein_pair.distance_analysis.analysis_results.distance_data,
                    i,
                )
                # create two atoms for the get_distance command
                atom1: str = (
                    f"/{a_protein_pair.protein_1.get_molecule_object()}//"
                    f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                )
                atom2: str = (
                    f"/{a_protein_pair.protein_2.get_molecule_object()}//"
                    f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                )
                # coloring
                self._embedded_cmd.color(color_5, atom1)
                self._embedded_cmd.color(color_5, atom2)
                i += 1
    
            elif distance_value > cutoff_5:
                atom_info = _get_chain_and_position(
                    a_protein_pair.distance_analysis.analysis_results.distance_data,
                    i,
                )
                # create two atoms for the get_distance command
                atom1: str = (
                    f"/{a_protein_pair.protein_1.get_molecule_object()}//"
                    f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
                )
                atom2: str = (
                    f"/{a_protein_pair.protein_2.get_molecule_object()}//"
                    f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
                )
                # coloring
                self._embedded_cmd.color(color_6, f"({atom1})")
                self._embedded_cmd.color(color_6, f"({atom2})")
                i += 1


def _get_chain_and_position(
        results_hashtable: dict[str, np.ndarray], index: int
) -> tuple:
    """This function gets the chain and the residue postion from the results hash table.
  
    Args:
        results_hashtable (dict[str, np.ndarray]): A hash table which contains all information from the distance calculation.
        index (int): An interator for the results hash table index.
  
    Returns:
        A tuple of ref_chain, ref_pos, ref_resi, model_chain, model_pos, model_resi for a specific index.
    """
    ref_chain_array = results_hashtable.get(
        pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN
    )
    ref_chain = ref_chain_array[index]
    ref_pos_array = results_hashtable.get(
        pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION
    )
    ref_pos = ref_pos_array[index]
    ref_resi_array = results_hashtable.get(pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI)
    ref_resi = ref_resi_array[index]
    model_chain_array = results_hashtable.get(
        pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN
    )
    model_chain = model_chain_array[index]
    model_pos_array = results_hashtable.get(
        pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION
    )
    model_pos = model_pos_array[index]
    model_resi_array = results_hashtable.get(
        pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI
    )
    model_resi = model_resi_array[index]
    return ref_chain, ref_pos, ref_resi, model_chain, model_pos, model_resi
