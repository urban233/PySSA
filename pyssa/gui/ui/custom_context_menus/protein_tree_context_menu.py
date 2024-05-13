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
"""Module for the protein tree context menu."""
from PyQt5 import QtWidgets


class ProteinTreeContextMenu:
    """A custom context menu wrapper for the protein tree view of the main dialog."""

    def __init__(self):
        self._context_menu = QtWidgets.QMenu()
        self._expand_protein_action = QtWidgets.QAction("Expand protein")
        self._clean_protein_action = QtWidgets.QAction("Clean Selected Protein")
        self._rename_protein_action = QtWidgets.QAction("Rename Selected Protein")
        self._show_sequence_action = QtWidgets.QAction("Show Sequence")
        self._help_action = QtWidgets.QAction("Get Help")

        self._context_menu.addAction(self._expand_protein_action)
        self._context_menu.addSeparator()
        self._context_menu.addAction(self._clean_protein_action)
        self._context_menu.addAction(self._rename_protein_action)
        self._context_menu.addAction(self._show_sequence_action)
        self._context_menu.addAction(self._help_action)

        self._expand_protein_action.triggered.connect(self._generic_action)
        self._clean_protein_action.triggered.connect(self._generic_action)
        self._rename_protein_action.triggered.connect(self._generic_action)
        self._show_sequence_action.triggered.connect(self._generic_action)
        self._help_action.triggered.connect(self._generic_action)

    # <editor-fold desc="Connect actions from outside">
    def connect_expand_protein_action(self, the_function_to_connect):
        self._expand_protein_action.triggered.disconnect()
        self._expand_protein_action.triggered.connect(the_function_to_connect)

    def connect_clean_protein_action(self, the_function_to_connect):
        self._clean_protein_action.triggered.disconnect()
        self._clean_protein_action.triggered.connect(the_function_to_connect)
    
    def connect_rename_protein_action(self, the_function_to_connect):
        self._rename_protein_action.triggered.disconnect()
        self._rename_protein_action.triggered.connect(the_function_to_connect)
        
    def connect_show_sequence_action(self, the_function_to_connect):
        self._show_sequence_action.triggered.disconnect()
        self._show_sequence_action.triggered.connect(the_function_to_connect)

    def connect_help_action(self, the_function_to_connect):
        self._help_action.triggered.disconnect()
        self._help_action.triggered.connect(the_function_to_connect)
    # </editor-fold>
    
    def _generic_action(self):
        raise NotImplementedError("The called action is not connected to a function!")
    
    def get_context_menu(self,
                         the_selected_indexes,
                         the_type: str,
                         is_protein_in_any_pair_flag: bool,
                         is_protein_in_session_flag: bool):
        # <editor-fold desc="Checks">
        if len(the_selected_indexes) > 0:
            level = 0
            index = the_selected_indexes[0]
            while index.parent().isValid():
                index = index.parent()
                level += 1
        else:
            self._help_action.setVisible(True)
            self._expand_protein_action.setVisible(False)
            self._clean_protein_action.setVisible(False)
            self._rename_protein_action.setVisible(False)
            self._show_sequence_action.setVisible(False)
            return self._context_menu
        # </editor-fold>

        # Manage of context menu
        self._help_action.setVisible(False)
        if level == 0:
            # A protein is selected
            self._clean_protein_action.setVisible(True)
            self._rename_protein_action.setVisible(False)
            self._show_sequence_action.setVisible(False)

            # <editor-fold desc="Check if protein is in any protein pair">
            if is_protein_in_any_pair_flag:
                self._rename_protein_action.setEnabled(False)
            else:
                self._rename_protein_action.setEnabled(True)

            if is_protein_in_session_flag:
                self._rename_protein_action.setEnabled(True)
            else:
                self._rename_protein_action.setEnabled(False)

            # </editor-fold>

        elif level == 1:
            # A header is selected (Scenes or Chains)
            self._clean_protein_action.setVisible(False)
            self._rename_protein_action.setVisible(False)
            self._show_sequence_action.setVisible(False)
        elif level == 2:
            # A chain or scene is selected
            if the_type == "chain":
                self._clean_protein_action.setVisible(False)
                self._rename_protein_action.setVisible(False)
                self._show_sequence_action.setVisible(True)
            elif the_type == "scene":
                # the code below could change if other actions are added!
                self._clean_protein_action.setVisible(False)
                self._rename_protein_action.setVisible(False)
                self._show_sequence_action.setVisible(False)
            else:
                self._clean_protein_action.setVisible(False)
                self._rename_protein_action.setVisible(False)
                self._show_sequence_action.setVisible(False)
        return self._context_menu
