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
"""Module for the protein pair tree context menu."""
import logging
from typing import Callable

from PyQt5 import QtWidgets

from pyssa.logging_pyssa import log_handlers
from pyssa.util import exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class ProteinPairTreeContextMenu:
    """A custom context menu wrapper for the protein pair tree view of the main dialog."""

    def __init__(self) -> None:
        """Constructor."""
        self._context_menu = QtWidgets.QMenu()
        self._expand_protein_pair_action = QtWidgets.QAction("Expand Protein Pair")
        self._collapse_protein_pair_action = QtWidgets.QAction("Collapse Protein Pair")
        self._open_results_summary_action = QtWidgets.QAction("Open Results Summary")
        self._color_based_on_rmsd_action = QtWidgets.QAction("Color Based On RMSD")
        self._help_action = QtWidgets.QAction("Get Help")

        self._context_menu.addAction(self._expand_protein_pair_action)
        self._context_menu.addAction(self._collapse_protein_pair_action)
        self._context_menu.addSeparator()
        self._context_menu.addAction(self._open_results_summary_action)
        self._context_menu.addAction(self._color_based_on_rmsd_action)
        self._context_menu.addAction(self._help_action)

        self._expand_protein_pair_action.triggered.connect(self._generic_action)
        self._collapse_protein_pair_action.triggered.connect(self._generic_action)
        self._open_results_summary_action.triggered.connect(self._generic_action)
        self._color_based_on_rmsd_action.triggered.connect(self._generic_action)
        self._help_action.triggered.connect(self._generic_action)

    # <editor-fold desc="Connect actions from outside">
    def connect_expand_protein_pair_action(self, the_function_to_connect: Callable) -> None:
        """Connects the `_expand_protein_pair_action` triggered signal to a given callable.

        Args:
            the_function_to_connect: A callable object that will be connected to the triggered signal of the `_expand_protein_pair_action` attribute.
        
        Raises:
            exception.IllegalArgumentError: If the_function_to_connect is None.
        """
        # <editor-fold desc="Checks">
        if the_function_to_connect is None:
            logger.error("the_function_to_connect is None.")
            raise exception.IllegalArgumentError("the_function_to_connect is None.")
        
        # </editor-fold>
        
        self._expand_protein_pair_action.triggered.disconnect()
        self._expand_protein_pair_action.triggered.connect(the_function_to_connect)

    def connect_collapse_protein_pair_action(self, the_function_to_connect: Callable) -> None:
        """Connects the `_expand_protein_pair_action` triggered signal to a given callable.

        Args:
            the_function_to_connect: A callable object that will be connected to the triggered signal of the `_expand_protein_pair_action` attribute.

        Raises:
            exception.IllegalArgumentError: If the_function_to_connect is None.
        """
        # <editor-fold desc="Checks">
        if the_function_to_connect is None:
            logger.error("the_function_to_connect is None.")
            raise exception.IllegalArgumentError("the_function_to_connect is None.")

        # </editor-fold>
        
        self._collapse_protein_pair_action.triggered.disconnect()
        self._collapse_protein_pair_action.triggered.connect(the_function_to_connect)

    def connect_open_results_summary_action(self, the_function_to_connect: Callable) -> None:
        """Connects the `_open_results_summary_action` triggered signal to a given callable.

        Args:
            the_function_to_connect: A callable object that will be connected to the triggered signal of the `_open_results_summary_action` attribute.

        Raises:
            exception.IllegalArgumentError: If the_function_to_connect is None.
        """
        # <editor-fold desc="Checks">
        if the_function_to_connect is None:
            logger.error("the_function_to_connect is None.")
            raise exception.IllegalArgumentError("the_function_to_connect is None.")

        # </editor-fold>
        
        self._open_results_summary_action.triggered.disconnect()
        self._open_results_summary_action.triggered.connect(the_function_to_connect)
    
    def connect_color_based_on_rmsd_action(self, the_function_to_connect: Callable) -> None:
        """Connects the `_color_based_on_rmsd_action` triggered signal to a given callable.

        Args:
            the_function_to_connect: A callable object that will be connected to the triggered signal of the `_color_based_on_rmsd_action` attribute.

        Raises:
            exception.IllegalArgumentError: If the_function_to_connect is None.
        """
        # <editor-fold desc="Checks">
        if the_function_to_connect is None:
            logger.error("the_function_to_connect is None.")
            raise exception.IllegalArgumentError("the_function_to_connect is None.")

        # </editor-fold>
        
        self._color_based_on_rmsd_action.triggered.disconnect()
        self._color_based_on_rmsd_action.triggered.connect(the_function_to_connect)

    def connect_help_action(self, the_function_to_connect: Callable) -> None:
        """Connects the `_help_action` triggered signal to a given callable.

        Args:
            the_function_to_connect: A callable object that will be connected to the triggered signal of the `_help_action` attribute.

        Raises:
            exception.IllegalArgumentError: If the_function_to_connect is None.
        """
        # <editor-fold desc="Checks">
        if the_function_to_connect is None:
            logger.error("the_function_to_connect is None.")
            raise exception.IllegalArgumentError("the_function_to_connect is None.")

        # </editor-fold>
        
        self._help_action.triggered.disconnect()
        self._help_action.triggered.connect(the_function_to_connect)
    
    # </editor-fold>
    
    def _generic_action(self) -> None:
        raise NotImplementedError("The called action is not connected to a function!")
    
    def get_context_menu(self,
                         the_selected_indexes,
                         is_protein_pair_in_current_session_flag: bool,
                         is_protein_pair_expanded: bool) -> QtWidgets.QMenu:
        """Gets the context menu for the current situation.

        Args:
            the_selected_indexes (list): A list of selected indexes.
            is_protein_pair_in_current_session_flag (bool): Flag indicating if the protein pair is in the current session.
            is_protein_pair_expanded (bool): Flag indicating if the protein pair is expanded.

        Returns:
            QMenu: The context menu.
        
        Raises:
            exception.IllegalArgumentError: If either `is_protein_pair_in_current_session_flag` or `is_protein_pair_expanded` is None.
        """
        # <editor-fold desc="Checks">
        if is_protein_pair_in_current_session_flag is None:
            logger.error("is_protein_pair_in_current_session_flag is None.")
            raise exception.IllegalArgumentError("is_protein_pair_in_current_session_flag is None.")
        if is_protein_pair_expanded is None:
            logger.error("is_protein_pair_expanded is None.")
            raise exception.IllegalArgumentError("is_protein_pair_expanded is None.")
        
        if len(the_selected_indexes) > 0:
            level = 0
            index = the_selected_indexes[0]
            while index.parent().isValid():
                index = index.parent()
                level += 1
        else:
            # No protein pair is selected
            self._expand_protein_pair_action.setVisible(False)
            self._collapse_protein_pair_action.setVisible(False)
            self._help_action.setVisible(True)
            self._open_results_summary_action.setVisible(False)
            self._color_based_on_rmsd_action.setVisible(False)
            return self._context_menu
        # </editor-fold>

        self._help_action.setVisible(False)
        if level == 0:
            # protein pair level
            self._open_results_summary_action.setVisible(True)
            self._color_based_on_rmsd_action.setVisible(True)
            self._expand_protein_pair_action.setVisible(True)

            if is_protein_pair_expanded is True:
                self._expand_protein_pair_action.setVisible(False)
                self._collapse_protein_pair_action.setVisible(True)
            else:
                self._expand_protein_pair_action.setVisible(True)
                self._collapse_protein_pair_action.setVisible(False)

            if is_protein_pair_in_current_session_flag:
                self._color_based_on_rmsd_action.setEnabled(True)
                self._expand_protein_pair_action.setVisible(True)
            else:
                self._color_based_on_rmsd_action.setEnabled(False)
        else:
            self._open_results_summary_action.setVisible(False)
            self._color_based_on_rmsd_action.setVisible(False)
            self._expand_protein_pair_action.setVisible(False)
            self._collapse_protein_pair_action.setVisible(False)
        return self._context_menu
