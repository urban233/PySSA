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
"""Module for the add protein pair view controller."""
import copy
import glob
import logging
import os
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5.QtCore import Qt
from pyssa.controller import interface_manager, watcher
from pyssa.internal.data_structures import chain, protein
from pyssa.logging_pyssa import log_levels, log_handlers
from pyssa.model import proteins_model
from pyssa.util import constants, enums

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class AddProteinPairViewController(QtCore.QObject):
    """Class for the AddProteinPairViewController class"""
    user_input = QtCore.pyqtSignal(tuple)

    def __init__(self,
                 the_interface_manager: "interface_manager.InterfaceManager",
                 the_watcher: "watcher.Watcher",
                 the_existing_analysis_runs: list[str],
                 the_protein_pairs: list[str],
                 a_list_of_extra_proteins: list[protein.Protein] = None):
        super().__init__()
        self._interface_manager = the_interface_manager
        self._watcher = the_watcher
        self._view = the_interface_manager.get_add_protein_pair_view()
        self._temporary_model = proteins_model.TemporaryProteinsModel()  # is needed to not mess up the main model!
        self._temporary_model.build_model_from_scratch(self._interface_manager.get_current_project().proteins)
        self._existing_analysis_runs = the_existing_analysis_runs
        self._existing_protein_pairs = the_protein_pairs
        self._number_of_prot_1_selected_chains: int = 1
        self.restore_ui()

        self._view.ui.tree_prot_1.setModel(self._temporary_model)
        self._view.ui.tree_prot_2.setModel(self._temporary_model)
        if a_list_of_extra_proteins is not None:
            self._add_additional_proteins_to_model(a_list_of_extra_proteins)
        self._hide_scenes_nodes()
        self._hide_non_protein_chains()
        self._connect_all_ui_elements_to_slot_functions()

    def restore_ui(self):
        """Restores the UI."""
        self._view.ui.lbl_prot_2.hide()
        self._view.ui.tree_prot_2.hide()
        self._view.ui.btn_back.hide()
        self._view.ui.btn_next.setEnabled(False)
        self._view.ui.btn_add.setEnabled(False)
        self.__slot_show_tree_prot_1()
        self._view.ui.tree_prot_1.setModel(None)
        self._view.ui.tree_prot_2.setModel(None)
        self._view.ui.tree_prot_1.setHeaderHidden(True)
        self._view.ui.tree_prot_2.setHeaderHidden(True)
        self._view.ui.tree_prot_1.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self._view.ui.tree_prot_2.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self._view.ui.tree_prot_1.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self._view.ui.tree_prot_2.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)

    def _add_additional_proteins_to_model(self, a_list_of_extra_proteins):
        for tmp_protein in a_list_of_extra_proteins:
            self._view.ui.tree_prot_1.model().add_temporary_protein(tmp_protein)

    def _hide_scenes_nodes(self):
        for tmp_row in range(self._view.ui.tree_prot_1.model().rowCount()):
            self._view.ui.tree_prot_1.setRowHidden(0, self._view.ui.tree_prot_1.model().index(tmp_row, 0), True)
            self._view.ui.tree_prot_2.setRowHidden(0, self._view.ui.tree_prot_2.model().index(tmp_row, 0), True)

    def _get_chain_indexes_from_tree_model(self, a_model: QtGui.QStandardItemModel) -> list:
        tmp_chain_indexes = []
        for row in range(a_model.rowCount()):
            tmp_model_index = a_model.index(row, 0)
            for sub_row in range(a_model.rowCount(tmp_model_index)):
                sub_model_index = a_model.index(sub_row, 0, tmp_model_index)
                if sub_model_index.data(Qt.DisplayRole) == "Chains":
                    for tmp_chain_row in range(a_model.rowCount(sub_model_index)):
                        tmp_chain_indexes.append(a_model.index(tmp_chain_row, 0, sub_model_index))
        return tmp_chain_indexes

    def _hide_non_protein_chains(self):
        for tmp_chain_index in self._get_chain_indexes_from_tree_model(self._view.ui.tree_prot_1.model()):
            tmp_chain: "chain.Chain" = tmp_chain_index.data(enums.ModelEnum.OBJECT_ROLE)
            if tmp_chain.chain_type == "non_protein_chain":
                self._view.ui.tree_prot_1.setRowHidden(tmp_chain_index.row(), tmp_chain_index.parent(), True)
                self._view.ui.tree_prot_2.setRowHidden(tmp_chain_index.row(), tmp_chain_index.parent(), True)

    def _get_first_protein_chain(self, an_index):
        i = 0
        tmp_loop_flag = True
        while tmp_loop_flag:
            tmp_chain = an_index.child(1, 0).child(i, 0).data(enums.ModelEnum.OBJECT_ROLE)
            if tmp_chain.chain_type == "protein_chain":
                return tmp_chain.chain_letter
            elif tmp_chain.chain_type == "non_protein_chain":
                i += 1
            else:
                tmp_loop_flag = False

    def _get_protein_name_and_chains(self, the_selected_indexes) -> tuple[str, list]:
        """Gets the name and the selected chains from a QTreeView selection."""
        tmp_protein_name = ""
        tmp_protein_chains = []
        for tmp_index in the_selected_indexes:
            if tmp_index.data(enums.ModelEnum.TYPE_ROLE) == "protein":
                tmp_protein_name = tmp_index.data(Qt.DisplayRole)
                tmp_protein_chains.append(self._get_first_protein_chain(tmp_index))
            elif tmp_index.data(enums.ModelEnum.TYPE_ROLE) == "header":
                tmp_protein_name = tmp_index.parent().data(Qt.DisplayRole)
                tmp_protein_chains.append(self._get_first_protein_chain(tmp_index.parent()))
            elif tmp_index.data(enums.ModelEnum.TYPE_ROLE) == "chain":
                tmp_protein_name = tmp_index.parent().parent().data(Qt.DisplayRole)
                tmp_protein_chains.append(tmp_index.data(Qt.DisplayRole))

        return tmp_protein_name, tmp_protein_chains

    def _create_analysis_run_name(self) -> str:
        """Creates the name of the analysis run."""
        prot_1_name, prot_1_chains = self._get_protein_name_and_chains(self._view.ui.tree_prot_1.selectedIndexes())
        prot_2_name, prot_2_chains = self._get_protein_name_and_chains(self._view.ui.tree_prot_2.selectedIndexes())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        # if prot_1_name == prot_2_name:
        #     # identical proteins getting compared and need therefore be indexed with _1 and _2
        #     tmp_analysis_run_name = f"{prot_1_name}_1;{prot_1_chains}_vs_{prot_2_name}_2;{prot_2_chains}"
        # else:
        tmp_analysis_run_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        return tmp_analysis_run_name

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.tree_prot_1.expanded.connect(self.__slot_collapse_all_tree_prot_1)
        self._view.ui.tree_prot_2.expanded.connect(self.__slot_collapse_all_tree_prot_2)
        self._view.ui.tree_prot_1.selectionModel().selectionChanged.connect(
            self.__slot_handle_selection_change_for_tree_prot_1
        )
        self._view.ui.tree_prot_2.selectionModel().selectionChanged.connect(
            self.__slot_handle_selection_change_for_tree_prot_2
        )
        self._view.ui.btn_next.clicked.connect(self.__slot_show_tree_prot_2)
        self._view.ui.btn_back.clicked.connect(self.__slot_show_tree_prot_1)
        self._view.ui.btn_add.clicked.connect(self.__slot_add_protein_pair)

    # <editor-fold desc="Slot methods">
    def __slot_handle_selection_change_for_tree_prot_1(self, selected, deselected):
        if selected.isEmpty():
            self._view.ui.btn_next.setEnabled(False)
            return

        tmp_selection_model = self._view.ui.tree_prot_1.selectionModel()
        tmp_selection = tmp_selection_model.selection()

        # <editor-fold desc="Checks for selection of multiple proteins">
        i = 0
        invalid = QtCore.QItemSelection()
        for index in tmp_selection.indexes():
            if index.data(enums.ModelEnum.TYPE_ROLE) == "protein":
                if i > 0:
                    invalid.select(index, index)
                i += 1
        if i > 1:
            tmp_selection_model.select(invalid, QtCore.QItemSelectionModel.Deselect)
            self._view.ui.btn_next.setEnabled(True)
            return
        # </editor-fold>

        # <editor-fold desc="Checks for selection of chains and proteins">
        parent = self._view.ui.tree_prot_1.currentIndex().parent()
        invalid = QtCore.QItemSelection()
        for index in tmp_selection.indexes():
            if index.parent() == parent:
                continue
            invalid.select(index, index)
        tmp_selection_model.select(invalid, QtCore.QItemSelectionModel.Deselect)
        self._view.ui.btn_next.setEnabled(True)
        # </editor-fold>

    def __slot_handle_selection_change_for_tree_prot_2(self, selected, deselected):
        if selected.isEmpty():
            self._view.ui.btn_add.setEnabled(False)
            return

        tmp_selection_model = self._view.ui.tree_prot_2.selectionModel()
        tmp_selection = tmp_selection_model.selection()

        # <editor-fold desc="Checks for selection of multiple proteins">
        i = 0
        invalid = QtCore.QItemSelection()
        for index in tmp_selection.indexes():
            if index.data(enums.ModelEnum.TYPE_ROLE) == "protein":
                if i > 0:
                    invalid.select(index, index)
                i += 1
        if i > 1:
            tmp_selection_model.select(invalid, QtCore.QItemSelectionModel.Deselect)
            tmp_analysis_run_name = self._create_analysis_run_name()
            tmp_analysis_run_name_without_semicolon = tmp_analysis_run_name.replace(";", "_")
            tmp_analysis_run_name_without_semicolon_and_comma = tmp_analysis_run_name_without_semicolon.replace(",", "_")
            if len(tmp_selection.indexes()) != self._number_of_prot_1_selected_chains:
                self._view.ui.btn_add.setEnabled(False)
            elif tmp_analysis_run_name in self._existing_analysis_runs:
                self._view.ui.btn_add.setEnabled(False)
            elif tmp_analysis_run_name_without_semicolon_and_comma in self._existing_protein_pairs:
                self._view.ui.btn_add.setEnabled(False)
            elif tmp_analysis_run_name_without_semicolon_and_comma in self._watcher.protein_pair_names_blacklist:
                self._view.ui.btn_add.setEnabled(False)
            else:
                self._view.ui.btn_add.setEnabled(True)
            return
        # </editor-fold>

        # <editor-fold desc="Checks for selection of chains and proteins">
        parent = self._view.ui.tree_prot_2.currentIndex().parent()
        invalid = QtCore.QItemSelection()
        for index in tmp_selection.indexes():
            if index.parent() == parent:
                continue
            invalid.select(index, index)
        tmp_selection_model.select(invalid, QtCore.QItemSelectionModel.Deselect)
        tmp_analysis_run_name = self._create_analysis_run_name()
        tmp_analysis_run_name_without_semicolon = tmp_analysis_run_name.replace(";", "_")
        tmp_analysis_run_name_without_semicolon_and_comma = tmp_analysis_run_name_without_semicolon.replace(",", "_")
        if len(tmp_selection.indexes()) != self._number_of_prot_1_selected_chains:
            self._view.ui.btn_add.setEnabled(False)
        elif tmp_analysis_run_name in self._existing_analysis_runs:
            self._view.ui.btn_add.setEnabled(False)
        elif tmp_analysis_run_name_without_semicolon_and_comma in self._existing_protein_pairs:
            self._view.ui.btn_add.setEnabled(False)
        elif tmp_analysis_run_name_without_semicolon_and_comma in self._watcher.protein_pair_names_blacklist:
            self._view.ui.btn_add.setEnabled(False)
        else:
            self._view.ui.btn_add.setEnabled(True)
        # </editor-fold>

    def __slot_collapse_all_tree_prot_1(self, the_selected_index: QtCore.QModelIndex):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "An object of the tree view 1 was expanded.")
        tmp_type = the_selected_index.data(enums.ModelEnum.TYPE_ROLE)
        if tmp_type == "protein":
            tmp_index_to_check = the_selected_index
            self._view.ui.tree_prot_1.setExpanded(tmp_index_to_check.child(1, 0), True)
        elif tmp_type == "header":
            tmp_index_to_check = the_selected_index.parent()
        elif tmp_type == "chain":
            tmp_index_to_check = the_selected_index.parent().parent()
        else:
            tmp_index_to_check = None

        for tmp_row in range(self._view.ui.tree_prot_1.model().rowCount()):
            tmp_index = self._view.ui.tree_prot_1.model().index(tmp_row, 0)
            if tmp_index != tmp_index_to_check and tmp_index.data(enums.ModelEnum.TYPE_ROLE) == "protein":
                self._view.ui.tree_prot_1.collapse(self._view.ui.tree_prot_1.model().index(tmp_row, 0))

    def __slot_collapse_all_tree_prot_2(self, the_selected_index: QtCore.QModelIndex):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "An object of the tree view 1 was expanded.")
        tmp_type = the_selected_index.data(enums.ModelEnum.TYPE_ROLE)
        if tmp_type == "protein":
            tmp_index_to_check = the_selected_index
            self._view.ui.tree_prot_2.setExpanded(tmp_index_to_check.child(1, 0), True)
        elif tmp_type == "header":
            tmp_index_to_check = the_selected_index.parent()
        elif tmp_type == "chain":
            tmp_index_to_check = the_selected_index.parent().parent()
        else:
            tmp_index_to_check = None

        for tmp_row in range(self._view.ui.tree_prot_2.model().rowCount()):
            tmp_index = self._view.ui.tree_prot_2.model().index(tmp_row, 0)
            if tmp_index != tmp_index_to_check and tmp_index.data(enums.ModelEnum.TYPE_ROLE) == "protein":
                self._view.ui.tree_prot_2.collapse(self._view.ui.tree_prot_2.model().index(tmp_row, 0))

    def __slot_show_tree_prot_2(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Next' button was clicked.")
        # UI changes
        self._view.ui.lbl_prot_1.setEnabled(False)
        self._view.ui.tree_prot_1.setEnabled(False)
        self._view.ui.btn_next.hide()
        self._view.ui.lbl_prot_2.show()
        self._view.ui.tree_prot_2.show()
        self._view.ui.btn_back.show()
        self._number_of_prot_1_selected_chains = len(self._view.ui.tree_prot_1.selectedIndexes())
        self._view.ui.lbl_prot_2.setText(
            f"Select second protein structure with {self._number_of_prot_1_selected_chains} chains"
        )

    def __slot_show_tree_prot_1(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Back' button was clicked.")
        # UI changes
        self._view.ui.lbl_prot_1.setEnabled(True)
        self._view.ui.tree_prot_1.setEnabled(True)
        self._view.ui.btn_next.show()
        self._view.ui.lbl_prot_2.hide()
        self._view.ui.tree_prot_2.hide()
        self._view.ui.btn_back.hide()
        self._view.ui.btn_add.setEnabled(False)

    def __slot_add_protein_pair(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Add' button was clicked.")
        tmp_item = QtWidgets.QListWidgetItem(self._create_analysis_run_name())
        self._view.close()
        self.user_input.emit((tmp_item, True))
    # </editor-fold>
