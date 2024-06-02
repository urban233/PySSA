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
import logging
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5.QtCore import Qt
from src.pyssa.controller import interface_manager, watcher
from src.pyssa.internal.data_structures import chain, protein
from src.pyssa.logging_pyssa import log_levels, log_handlers
from src.pyssa.model import proteins_model
from src.pyssa.util import enums, exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class AddProteinPairViewController(QtCore.QObject):
  """Class for the AddProteinPairViewController."""

  user_input = QtCore.pyqtSignal(tuple)
  """Singal used to transfer data back to the previous window."""

  def __init__(
      self,
      the_interface_manager: "interface_manager.InterfaceManager",
      the_watcher: "watcher.Watcher",
      the_existing_analysis_runs: list[str],
      the_protein_pairs: list[str],
      a_list_of_extra_proteins: list[protein.Protein] = None,
  ) -> None:
    """Constructor.

    Args:
        the_interface_manager (interface_manager.InterfaceManager): An instance of the InterfaceManager class.
        the_watcher (watcher.Watcher): An instance of the Watcher class.
        the_existing_analysis_runs (list[str]): A list of strings representing existing analysis runs.
        the_protein_pairs (list[str]): A list of strings representing protein pairs.
        a_list_of_extra_proteins (list[protein.Protein]): (optional) A list of Protein objects.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_interface_manager is None:
      logger.error("the_interface_manager is None.")
      raise exception.IllegalArgumentError("the_interface_manager is None.")
    if the_watcher is None:
      logger.error("the_watcher is None.")
      raise exception.IllegalArgumentError("the_watcher is None.")
    if the_existing_analysis_runs is None:
      logger.error("the_existing_analysis_runs is None.")
      raise exception.IllegalArgumentError(
          "the_existing_analysis_runs is None."
      )
    if the_protein_pairs is None:
      logger.error("the_protein_pairs is None.")
      raise exception.IllegalArgumentError("the_protein_pairs is None.")

    # </editor-fold>

    super().__init__()
    self._interface_manager = the_interface_manager
    self._watcher = the_watcher
    self._view = the_interface_manager.get_add_protein_pair_view()
    self._temporary_model = (
        proteins_model.TemporaryProteinsModel()
    )  # is needed to not mess up the main model!
    self._temporary_model.build_model_from_scratch(
        self._interface_manager.get_current_project().proteins
    )
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

  def restore_ui(self) -> None:
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
    self._view.ui.tree_prot_1.setEditTriggers(
        QtWidgets.QAbstractItemView.NoEditTriggers
    )
    self._view.ui.tree_prot_2.setEditTriggers(
        QtWidgets.QAbstractItemView.NoEditTriggers
    )
    self._view.ui.tree_prot_1.setSelectionMode(
        QtWidgets.QAbstractItemView.ExtendedSelection
    )
    self._view.ui.tree_prot_2.setSelectionMode(
        QtWidgets.QAbstractItemView.ExtendedSelection
    )

  def _add_additional_proteins_to_model(
      self, a_list_of_extra_proteins: list
  ) -> None:
    """Adds additional proteins to the temporary model.

    Args:
        a_list_of_extra_proteins (list): A list of additional proteins to be added to the model.

    Raises:
        exception.IllegalArgumentError: If `a_list_of_extra_proteins` is None.
    """
    # <editor-fold desc="Checks">
    if a_list_of_extra_proteins is None:
      logger.error("a_list_of_extra_proteins is None.")
      raise exception.IllegalArgumentError("a_list_of_extra_proteins is None.")

    # </editor-fold>

    for tmp_protein in a_list_of_extra_proteins:
      self._view.ui.tree_prot_1.model().add_temporary_protein(tmp_protein)

  def _hide_scenes_nodes(self) -> None:
    """Hides the nodes in the tree_prot_1 and tree_prot_2 views.

    This method iterates through the rows in the model of tree_prot_1 and tree_prot_2 views,
    and sets the row at index (0, tmp_row) to be hidden.
    """
    for tmp_row in range(self._view.ui.tree_prot_1.model().rowCount()):
      self._view.ui.tree_prot_1.setRowHidden(
          0, self._view.ui.tree_prot_1.model().index(tmp_row, 0), True
      )
      self._view.ui.tree_prot_2.setRowHidden(
          0, self._view.ui.tree_prot_2.model().index(tmp_row, 0), True
      )

  def _get_chain_indexes_from_tree_model(
      self, a_model: QtGui.QStandardItemModel
  ) -> list:
    """Gets a list of chain indexes from the given model.

    Args:
        a_model (QtGui.QStandardItemModel): The tree model from which to retrieve the chain indexes.

    Returns:
        list: A list of QModelIndex objects representing the chain indexes in the tree model.

    Raises:
        exception.IllegalArgumentError: If `a_model` is None.
    """
    # <editor-fold desc="Checks">
    if a_model is None:
      logger.error("a_model is None.")
      raise exception.IllegalArgumentError("a_model is None.")

    # </editor-fold>

    tmp_chain_indexes = []
    for row in range(a_model.rowCount()):
      tmp_model_index = a_model.index(row, 0)
      for sub_row in range(a_model.rowCount(tmp_model_index)):
        sub_model_index = a_model.index(sub_row, 0, tmp_model_index)
        if sub_model_index.data(Qt.DisplayRole) == "Chains":
          for tmp_chain_row in range(a_model.rowCount(sub_model_index)):
            tmp_chain_indexes.append(
                a_model.index(tmp_chain_row, 0, sub_model_index)
            )
    return tmp_chain_indexes

  def _hide_non_protein_chains(self) -> None:
    """Hide the rows in both tree_prot_1 and tree_prot_2 that correspond to non-protein chains."""
    for tmp_chain_index in self._get_chain_indexes_from_tree_model(
        self._view.ui.tree_prot_1.model()
    ):
      tmp_chain: "chain.Chain" = tmp_chain_index.data(
          enums.ModelEnum.OBJECT_ROLE
      )
      if tmp_chain.chain_type == "non_protein_chain":
        self._view.ui.tree_prot_1.setRowHidden(
            tmp_chain_index.row(), tmp_chain_index.parent(), True
        )
        self._view.ui.tree_prot_2.setRowHidden(
            tmp_chain_index.row(), tmp_chain_index.parent(), True
        )

  def _get_first_protein_chain(self, an_index: QtCore.QModelIndex) -> str:
    """Searches for the first protein chain within the given index object.

    It iterates over the child models of the index object until it finds a protein chain or reaches the end of the child models.

    Args:
        an_index (QtCore.QModelIndex): The index object to search for the first protein chain.

    Returns:
         A chain letter or an empty string if no matching chain is found.

    Raises:
        exception.IllegalArgumentError: If `an_index` is None.
    """
    # <editor-fold desc="Checks">
    if an_index is None:
      logger.error("an_index is None.")
      raise exception.IllegalArgumentError("an_index is None.")

    # </editor-fold>

    i = 0
    tmp_loop_flag = True
    while tmp_loop_flag is True:
      tmp_chain = (
          an_index.child(1, 0).child(i, 0).data(enums.ModelEnum.OBJECT_ROLE)
      )
      if tmp_chain.chain_type == "protein_chain":
        return tmp_chain.chain_letter
      elif tmp_chain.chain_type == "non_protein_chain":
        i += 1
      else:
        tmp_loop_flag = False
    return ""

  def _get_protein_name_and_chains(
      self, the_selected_indexes: list
  ) -> tuple[str, list]:
    """Gets the name and the selected chains from a QTreeView selection.

    Args:
        the_selected_indexes (list): A list of selected indexes.

    Returns:
        A tuple containing the protein name (str) and a list of protein chains (list[str]).

    Raises:
        exception.IllegalArgumentError: If `the_selected_indexes` is None.
    """
    # <editor-fold desc="Checks">
    if the_selected_indexes is None:
      logger.error("the_selected_indexes is None.")
      raise exception.IllegalArgumentError("the_selected_indexes is None.")

    # </editor-fold>

    tmp_protein_name = ""
    tmp_protein_chains = []
    for tmp_index in the_selected_indexes:
      if tmp_index.data(enums.ModelEnum.TYPE_ROLE) == "protein":
        tmp_protein_name = tmp_index.data(Qt.DisplayRole)
        tmp_protein_chains.append(self._get_first_protein_chain(tmp_index))
      elif tmp_index.data(enums.ModelEnum.TYPE_ROLE) == "header":
        tmp_protein_name = tmp_index.parent().data(Qt.DisplayRole)
        tmp_protein_chains.append(
            self._get_first_protein_chain(tmp_index.parent())
        )
      elif tmp_index.data(enums.ModelEnum.TYPE_ROLE) == "chain":
        tmp_protein_name = tmp_index.parent().parent().data(Qt.DisplayRole)
        tmp_protein_chains.append(tmp_index.data(Qt.DisplayRole))

    return tmp_protein_name, tmp_protein_chains

  def _create_analysis_run_name(self) -> str:
    """Creates the name of the analysis run."""
    prot_1_name, prot_1_chains = self._get_protein_name_and_chains(
        self._view.ui.tree_prot_1.selectedIndexes()
    )
    prot_2_name, prot_2_chains = self._get_protein_name_and_chains(
        self._view.ui.tree_prot_2.selectedIndexes()
    )
    prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
    prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
    # if prot_1_name == prot_2_name:
    #     # identical proteins getting compared and need therefore be indexed with _1 and _2
    #     tmp_analysis_run_name = f"{prot_1_name}_1;{prot_1_chains}_vs_{prot_2_name}_2;{prot_2_chains}"
    # else:
    tmp_analysis_run_name = (
        f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
    )
    return tmp_analysis_run_name

  def _connect_all_ui_elements_to_slot_functions(self) -> None:
    """Connects all UI elements to their corresponding slot functions in the class."""
    self._view.ui.tree_prot_1.expanded.connect(
        self.__slot_collapse_all_tree_prot_1
    )
    self._view.ui.tree_prot_2.expanded.connect(
        self.__slot_collapse_all_tree_prot_2
    )
    self._view.ui.tree_prot_1.selectionModel().selectionChanged.connect(
        self.__slot_handle_selection_change_for_tree_prot_1,
    )
    self._view.ui.tree_prot_2.selectionModel().selectionChanged.connect(
        self.__slot_handle_selection_change_for_tree_prot_2,
    )
    self._view.ui.btn_next.clicked.connect(self.__slot_show_tree_prot_2)
    self._view.ui.btn_back.clicked.connect(self.__slot_show_tree_prot_1)
    self._view.ui.btn_add.clicked.connect(self.__slot_add_protein_pair)

  # <editor-fold desc="Slot methods">
  def __slot_handle_selection_change_for_tree_prot_1(
      self, selected: QtCore.QItemSelection, deselected: QtCore.QItemSelection
  ) -> None:
    """Handles the selection change event for the 'tree_prot_1' tree view.

    Args:
        selected (QtCore.QItemSelection): The selected items in the tree view.
        deselected (QtCore.QItemSelection): The deselected items in the tree view.

    Raises:
        exception.IllegalArgumentError: If `selected` is None.
    """
    # <editor-fold desc="Checks">
    if selected is None:
      logger.error("selected is None.")
      raise exception.IllegalArgumentError("selected is None.")

    # </editor-fold>

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

  def __slot_handle_selection_change_for_tree_prot_2(
      self, selected: QtCore.QItemSelection, deselected: QtCore.QItemSelection
  ) -> None:
    """Handles the selection change event for the 'tree_prot_2' tree view.

    Args:
        selected (QtCore.QItemSelection): The selected items in the tree view.
        deselected (QtCore.QItemSelection): The deselected items in the tree view.

    Raises:
        exception.IllegalArgumentError: If `selected` is None.
    """
    # <editor-fold desc="Checks">
    if selected is None:
      logger.error("selected is None.")
      raise exception.IllegalArgumentError("selected is None.")

    # </editor-fold>

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
      tmp_analysis_run_name_without_semicolon = tmp_analysis_run_name.replace(
          ";", "_"
      )
      tmp_analysis_run_name_without_semicolon_and_comma = (
          tmp_analysis_run_name_without_semicolon.replace(",", "_")
      )
      if len(tmp_selection.indexes()) != self._number_of_prot_1_selected_chains:
        self._view.ui.btn_add.setEnabled(False)
      elif tmp_analysis_run_name in self._existing_analysis_runs:
        self._view.ui.btn_add.setEnabled(False)
      elif (
          tmp_analysis_run_name_without_semicolon_and_comma
          in self._existing_protein_pairs
      ):
        self._view.ui.btn_add.setEnabled(False)
      elif (
          tmp_analysis_run_name_without_semicolon_and_comma
          in self._watcher.protein_pair_names_blacklist
      ):
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
    tmp_analysis_run_name_without_semicolon = tmp_analysis_run_name.replace(
        ";", "_"
    )
    tmp_analysis_run_name_without_semicolon_and_comma = (
        tmp_analysis_run_name_without_semicolon.replace(",", "_")
    )
    if len(tmp_selection.indexes()) != self._number_of_prot_1_selected_chains:
      self._view.ui.btn_add.setEnabled(False)
    elif tmp_analysis_run_name in self._existing_analysis_runs:
      self._view.ui.btn_add.setEnabled(False)
    elif (
        tmp_analysis_run_name_without_semicolon_and_comma
        in self._existing_protein_pairs
    ):
      self._view.ui.btn_add.setEnabled(False)
    elif (
        tmp_analysis_run_name_without_semicolon_and_comma
        in self._watcher.protein_pair_names_blacklist
    ):
      self._view.ui.btn_add.setEnabled(False)
    else:
      self._view.ui.btn_add.setEnabled(True)
    # </editor-fold>

  def __slot_collapse_all_tree_prot_1(
      self, the_selected_index: QtCore.QModelIndex
  ) -> None:
    """Collapses all items in the tree view 1 except for the selected item.

    Args:
        the_selected_index (QtCore.QModelIndex): The index of the selected item in the tree view.

    Raises:
        exception.IllegalArgumentError: If `the_selected_index` is None.
    """
    # <editor-fold desc="Checks">
    if the_selected_index is None:
      logger.error("the_selected_index is None.")
      raise exception.IllegalArgumentError("the_selected_index is None.")
    # </editor-fold>

    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "An object of the tree view 1 was expanded.",
    )
    tmp_type = the_selected_index.data(enums.ModelEnum.TYPE_ROLE)
    if tmp_type == "protein":
      tmp_index_to_check = the_selected_index
      self._view.ui.tree_prot_1.setExpanded(
          tmp_index_to_check.child(1, 0), True
      )
    elif tmp_type == "header":
      tmp_index_to_check = the_selected_index.parent()
    elif tmp_type == "chain":
      tmp_index_to_check = the_selected_index.parent().parent()
    else:
      tmp_index_to_check = None

    for tmp_row in range(self._view.ui.tree_prot_1.model().rowCount()):
      tmp_index = self._view.ui.tree_prot_1.model().index(tmp_row, 0)
      if (
          tmp_index != tmp_index_to_check
          and tmp_index.data(enums.ModelEnum.TYPE_ROLE) == "protein"
      ):
        self._view.ui.tree_prot_1.collapse(
            self._view.ui.tree_prot_1.model().index(tmp_row, 0)
        )

  def __slot_collapse_all_tree_prot_2(
      self, the_selected_index: QtCore.QModelIndex
  ) -> None:
    """Collapses all items in the tree view 2 except for the selected item.

    Args:
        the_selected_index (QtCore.QModelIndex): The index of the selected item in the tree view.

    Raises:
        exception.IllegalArgumentError: If `the_selected_index` is None.
    """
    # <editor-fold desc="Checks">
    if the_selected_index is None:
      logger.error("the_selected_index is None.")
      raise exception.IllegalArgumentError("the_selected_index is None.")
    # </editor-fold>

    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "An object of the tree view 1 was expanded.",
    )
    tmp_type = the_selected_index.data(enums.ModelEnum.TYPE_ROLE)
    if tmp_type == "protein":
      tmp_index_to_check = the_selected_index
      self._view.ui.tree_prot_2.setExpanded(
          tmp_index_to_check.child(1, 0), True
      )
    elif tmp_type == "header":
      tmp_index_to_check = the_selected_index.parent()
    elif tmp_type == "chain":
      tmp_index_to_check = the_selected_index.parent().parent()
    else:
      tmp_index_to_check = None

    for tmp_row in range(self._view.ui.tree_prot_2.model().rowCount()):
      tmp_index = self._view.ui.tree_prot_2.model().index(tmp_row, 0)
      if (
          tmp_index != tmp_index_to_check
          and tmp_index.data(enums.ModelEnum.TYPE_ROLE) == "protein"
      ):
        self._view.ui.tree_prot_2.collapse(
            self._view.ui.tree_prot_2.model().index(tmp_row, 0)
        )

  def __slot_show_tree_prot_2(self) -> None:
    """Performs the necessary UI changes and updates the label and tree widget for the second protein structure selection."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Next' button was clicked."
    )
    # UI changes
    self._view.ui.lbl_prot_1.setEnabled(False)
    self._view.ui.tree_prot_1.setEnabled(False)
    self._view.ui.btn_next.hide()
    self._view.ui.lbl_prot_2.show()
    self._view.ui.tree_prot_2.show()
    self._view.ui.btn_back.show()
    self._number_of_prot_1_selected_chains = len(
        self._view.ui.tree_prot_1.selectedIndexes()
    )
    self._view.ui.lbl_prot_2.setText(
        f"Select second protein structure with {self._number_of_prot_1_selected_chains} chains",
    )

  def __slot_show_tree_prot_1(self) -> None:
    """Performs the necessary UI changes and updates the label and tree widget for the first protein structure selection."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Back' button was clicked."
    )
    # UI changes
    self._view.ui.lbl_prot_1.setEnabled(True)
    self._view.ui.tree_prot_1.setEnabled(True)
    self._view.ui.btn_next.show()
    self._view.ui.lbl_prot_2.hide()
    self._view.ui.tree_prot_2.hide()
    self._view.ui.btn_back.hide()
    self._view.ui.btn_add.setEnabled(False)

  def __slot_add_protein_pair(self) -> None:
    """Emits a signal with the protein pair information and closes the dialog window."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Add' button was clicked."
    )
    tmp_item = QtWidgets.QListWidgetItem(self._create_analysis_run_name())
    self._view.close()
    self.user_input.emit((tmp_item, True))

  # </editor-fold>
