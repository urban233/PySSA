#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
# 
# Source code is available at <https://github.com/zielesny/PySSA>
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
"""Module contains the MainViewState dataclass."""
from dataclasses import dataclass
from PyQt5 import QtCore

from pyssa.internal.data_structures import protein, protein_pair


@dataclass
class MainViewState:
    """Contains the current state of the main view."""
    selected_sequence: QtCore.QModelIndex
    selected_chain_proteins: QtCore.QModelIndex
    current_proteins: list[str]
    current_protein_pairs: list[str]
    selected_chain_protein_pairs: QtCore.QModelIndex

    def __init__(self,
                 the_sequence_list_view,
                 _show_sequence_information,
                 the_proteins_tree_view,
                 __slot_get_information_about_selected_object_in_protein_branch,
                 the_protein_pairs_tree_view,
                 __slot_get_information_about_selected_object_in_protein_pair_branch
                 ):
        self.selected_sequence: QtCore.QModelIndex = None
        self.selected_chain_proteins: QtCore.QModelIndex = None
        self.selected_chain_protein_pairs: QtCore.QModelIndex = None

        self.sequence_list_view = the_sequence_list_view
        self._show_sequence_information = _show_sequence_information
        self.proteins_tree_view = the_proteins_tree_view
        self.__slot_get_information_about_selected_object_in_protein_branch = __slot_get_information_about_selected_object_in_protein_branch
        self.protein_pairs_tree_view = the_protein_pairs_tree_view
        self.__slot_get_information_about_selected_object_in_protein_pair_branch = __slot_get_information_about_selected_object_in_protein_pair_branch

        self.current_proteins: list[str] = []
        self.current_protein_pairs: list[str] = []

    def restore_main_view_state(self):
        """Restores the main view state."""
        # Restore sequences tab
        if self.selected_sequence is not None:
            self.sequence_list_view.setCurrentIndex(self.selected_sequence)
            self._show_sequence_information()
            self.selected_sequence = None
        # Restore proteins tab
        if self.selected_chain_proteins is not None:
            self.proteins_tree_view.setCurrentIndex(
                self.selected_chain_proteins
            )
            self.__slot_get_information_about_selected_object_in_protein_branch()
            self.selected_chain_proteins = None
        # Restore protein pairs tab
        if self.selected_chain_protein_pairs is not None:
            self.protein_pairs_tree_view.setCurrentIndex(
                self.selected_chain_protein_pairs
            )
            self.__slot_get_information_about_selected_object_in_protein_pair_branch()
            self.selected_chain_protein_pairs = None

    def set_proteins_list(self, a_list_of_proteins: list[protein.Protein]):
        for tmp_protein in a_list_of_proteins:
            self.current_proteins.append(tmp_protein.get_molecule_object())

    def get_not_matching_proteins(self, proteins_to_check: list["protein.Protein"]) -> list["protein.Protein"]:
        """Returns all proteins that are not in the current_proteins list of this class."""
        not_matching_proteins = []
        for tmp_protein in proteins_to_check:
            if tmp_protein.get_molecule_object() not in self.current_proteins:
                not_matching_proteins.append(tmp_protein)
        return not_matching_proteins

    def set_protein_pairs_list(self, a_list_of_protein_pairs: list["protein_pair.ProteinPair"]):
        for tmp_protein_pair in a_list_of_protein_pairs:
            self.current_protein_pairs.append(tmp_protein_pair.name)

    def get_not_matching_protein_pairs(self, protein_pairs_to_check: list["protein_pair.ProteinPair"]) -> list["protein_pair.ProteinPair"]:
        """Returns all proteins that are not in the current_proteins list of this class."""
        not_matching_protein_pairs = []
        for tmp_protein_pair in protein_pairs_to_check:
            if tmp_protein_pair.name not in self.current_protein_pairs:
                not_matching_protein_pairs.append(tmp_protein_pair)
        return not_matching_protein_pairs
