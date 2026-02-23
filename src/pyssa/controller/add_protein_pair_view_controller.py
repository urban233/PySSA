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
"""Module for the add protein pair view controller."""
import logging
from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.qt import QtWidgets
from src.pyssa.gui.qt import QtGui
from src.pyssa.gui.qt import Qt
from typing import TYPE_CHECKING
if TYPE_CHECKING:
  from src.pyssa.gui import app_state
from src.pyssa.gui import name_registry as name_registry_module
from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.internal.data_structures import chain, protein
from src.pyssa.logging_pyssa import log_levels, log_handlers
from src.pyssa.model import psa_objects_model
from src.pyssa.util import enums, exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class AddProteinPairViewController(QtCore.QObject):
  """Class for the AddProteinPairViewController."""

  def __init__(
      self,
      the_app_state: "app_state.AppState",
      a_list_of_used_run_names: list[str],
      a_list_of_used_protein_pair_names: list[str],
      a_list_of_extra_proteins: list[protein.Protein] = None,
      on_add_callback = None,
      a_parent=None
  ) -> None:
    """Constructor.

    Args:
        the_app_state (app_state.AppState): An instance of the AppState class.
        a_list_of_used_run_names (list[str]): A list of strings representing existing analysis runs.
        a_list_of_used_protein_pair_names (list[str]): A list of strings representing protein pairs.
        a_list_of_extra_proteins (list[protein.Protein]): (optional) A list of Protein objects.
        a_parent: Parent widget to pass to the view.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_app_state is None:
      logger.error("the_app_state is None.")
      raise exception.IllegalArgumentError("the_app_state is None.")
    if a_list_of_used_run_names is None:
      logger.error("a_list_of_used_run_names is None.")
      raise exception.IllegalArgumentError(
          "a_list_of_used_run_names is None."
      )
    if a_list_of_used_protein_pair_names is None:
      logger.error("a_list_of_used_protein_pair_names is None.")
      raise exception.IllegalArgumentError("a_list_of_used_protein_pair_names is None.")

    # </editor-fold>

    super().__init__()
    self._app_state = the_app_state
    self._on_add_callback = on_add_callback
    from src.pyssa.gui.ui.views import add_protein_pair_view
    self._view = add_protein_pair_view.AddProteinPairView(a_parent)
    self._local_model = psa_objects_model.PSAObjectsModel()
    if self._app_state.has_open_project():
      main_model = self._app_state.pyssa_objects_model
      main_proteins_index = main_model.get_proteins_section_index()
      main_proteins_item = main_model.itemFromIndex(main_proteins_index)
      
      local_proteins_index = self._local_model.get_proteins_section_index()
      local_proteins_item = self._local_model.itemFromIndex(local_proteins_index)

      if main_proteins_item is not None and local_proteins_item is not None:
        for row in range(main_proteins_item.rowCount()):
          protein_item = main_proteins_item.child(row, 0)
          if protein_item is not None:
            tmp_protein = protein_item.data(enums.ModelEnum.OBJECT_ROLE)
            if self._check_if_protein_has_protein_chains(tmp_protein):
              local_proteins_item.appendRow(self._deep_copy_item(protein_item))
    self._existing_analysis_runs = a_list_of_used_run_names
    self._existing_protein_pairs = a_list_of_used_protein_pair_names
    self._number_of_prot_1_selected_chains: int = 1
    self.restore_ui()
    self.temporary_model_is_valid = False

    self._view.ui.tree_prot_1.setModel(self._local_model)
    self._view.ui.tree_prot_2.setModel(self._local_model)
    proteins_section_index = self._local_model.get_proteins_section_index()
    self._view.ui.tree_prot_1.setRootIndex(proteins_section_index)
    self._view.ui.tree_prot_2.setRootIndex(proteins_section_index)
    
    if a_list_of_extra_proteins is not None:
      self._add_additional_proteins_to_model(a_list_of_extra_proteins)
    if self._local_model.rowCount(proteins_section_index) == 0:
      tmp_dialog = custom_message_box.CustomMessageBoxOk(
        "All proteins in the project have only non-protein chains!",
        "No Protein Chains",
        custom_message_box.CustomMessageBoxIcons.ERROR.value,
      )
      tmp_dialog.exec()
    else:
      self.temporary_model_is_valid = True
      self._hide_scenes_nodes()
      self._hide_non_protein_chains()
      self._connect_all_ui_elements_to_slot_functions()

  def get_view(self):
    return self._view

  def _check_if_protein_has_protein_chains(self, a_protein: "protein.Protein") -> bool:
    """Checks if at least one chain in the protein is an actual protein chain."""
    if a_protein is None:
      return False
    for tmp_chain_in_protein in a_protein.chains:
      if tmp_chain_in_protein.chain_type == "protein_chain":
        return True
    return False

  def _deep_copy_item(self, item: QtGui.QStandardItem) -> QtGui.QStandardItem:
    """Creates a deep copy of a QStandardItem and its children."""
    new_item = item.clone()
    for row in range(item.rowCount()):
        for col in range(item.columnCount()):
            child = item.child(row, col)
            if child is not None:
                new_item.setChild(row, col, self._deep_copy_item(child))
    return new_item

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
      if self._check_if_protein_has_protein_chains(tmp_protein):
        self._local_model.add_temporary_protein(tmp_protein)

  def _hide_scenes_nodes(self) -> None:
    """Hides the nodes in the tree_prot_1 and tree_prot_2 views."""
    tmp_parent = self._local_model.get_proteins_section_index()
    for tmp_row in range(self._local_model.rowCount(tmp_parent)):
      tmp_protein_index = self._local_model.index(tmp_row, 0, tmp_parent)
      self._view.ui.tree_prot_1.setRowHidden(
          0, tmp_protein_index, True
      )
      self._view.ui.tree_prot_2.setRowHidden(
          0, tmp_protein_index, True
      )

  def _get_chain_indexes_from_tree_model(
      self, a_model: QtGui.QStandardItemModel, a_parent: QtCore.QModelIndex = None
  ) -> list:
    """Gets a list of chain indexes from the given model.

    Args:
        a_model (QtGui.QStandardItemModel): The tree model from which to retrieve the chain indexes.
        a_parent (QtCore.QModelIndex): An optional parent index for nested searches.

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
    tmp_parent = a_parent if a_parent is not None else QtCore.QModelIndex()
    for row in range(a_model.rowCount(tmp_parent)):
      tmp_model_index = a_model.index(row, 0, tmp_parent)
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
    tmp_parent = self._local_model.get_proteins_section_index()
    for tmp_chain_index in self._get_chain_indexes_from_tree_model(
        self._local_model, tmp_parent
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

    Iterates over the chain children of the protein node until a protein chain
    is found.  Uses `model.index()` for child navigation because
    `QModelIndex.child()` was removed in PyQt6.

    Args:
        an_index (QtCore.QModelIndex): The protein-level index to search.

    Returns:
        The chain letter of the first protein chain, or an empty string when
        no protein chain exists under the given index.

    Raises:
        exception.IllegalArgumentError: If `an_index` is None.
    """
    # <editor-fold desc="Checks">
    if an_index is None:
      logger.error("an_index is None.")
      raise exception.IllegalArgumentError("an_index is None.")

    # </editor-fold>

    tmp_model = an_index.model()
    # Row 1 under a protein node is always the "Chains" header.
    chains_header_index = tmp_model.index(1, 0, an_index)
    if not chains_header_index.isValid():
      return ""
    chain_count = tmp_model.rowCount(chains_header_index)
    for i in range(chain_count):
      chain_index = tmp_model.index(i, 0, chains_header_index)
      tmp_chain = chain_index.data(enums.ModelEnum.OBJECT_ROLE)
      if tmp_chain is not None and tmp_chain.chain_type == "protein_chain":
        return tmp_chain.chain_letter
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
      elif tmp_index.data(enums.ModelEnum.TYPE_ROLE) == "residue":
        tmp_protein_name = tmp_index.parent().parent().parent().data(Qt.DisplayRole)
        tmp_protein_chains.append(tmp_index.parent().data(Qt.DisplayRole))
      elif tmp_index.data(enums.ModelEnum.TYPE_ROLE) == "atom":
        tmp_protein_name = tmp_index.parent().parent().parent().parent().data(Qt.DisplayRole)
        tmp_protein_chains.append(tmp_index.parent().parent().data(Qt.DisplayRole))

    return tmp_protein_name, list(set(tmp_protein_chains))

  def _create_analysis_run_name(self) -> str:
    """Creates the name of the analysis run based on the current tree selections.

    Returns:
        The analysis run name as a formatted string of the form
        ``prot1;chain(s)_vs_prot2;chain(s)``.
    """
    prot_1_name, prot_1_chains = self._get_protein_name_and_chains(
        self._view.ui.tree_prot_1.selectedIndexes()
    )
    prot_2_name, prot_2_chains = self._get_protein_name_and_chains(
        self._view.ui.tree_prot_2.selectedIndexes()
    )
    prot_1_chains_str = ",".join([str(elem) for elem in prot_1_chains])
    prot_2_chains_str = ",".join([str(elem) for elem in prot_2_chains])
    return f"{prot_1_name};{prot_1_chains_str}_vs_{prot_2_name};{prot_2_chains_str}"

  def _evaluate_and_update_add_button(self, the_selection: QtCore.QItemSelection) -> None:
    """Evaluates the current tree_prot_2 selection and enables or disables the Add button.

    The Add button is enabled only when the number of selected chains matches
    the number of chains selected in tree_prot_1, and the resulting analysis
    run name does not conflict with any existing analysis run, protein pair, or
    reserved name.

    Args:
        the_selection (QtCore.QItemSelection): The current selection in tree_prot_2.
    """
    tmp_analysis_run_name = self._create_analysis_run_name()
    tmp_run_name_normalised = tmp_analysis_run_name.replace(";", "_").replace(",", "_")
    
    _, prot_2_chains = self._get_protein_name_and_chains(
        self._view.ui.tree_prot_2.selectedIndexes()
    )

    if len(prot_2_chains) != self._number_of_prot_1_selected_chains:
      self._view.ui.btn_add.setEnabled(False)
    elif tmp_analysis_run_name in self._existing_analysis_runs:
      self._view.ui.btn_add.setEnabled(False)
    elif tmp_run_name_normalised in self._existing_protein_pairs:
      self._view.ui.btn_add.setEnabled(False)
    elif self._app_state.name_registry.is_reserved(
        name_registry_module.PROTEIN_PAIR, tmp_run_name_normalised
    ):
      self._view.ui.btn_add.setEnabled(False)
    else:
      self._view.ui.btn_add.setEnabled(True)

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

    # <editor-fold desc="Deselect additional protein nodes when multiple are chosen">
    protein_count = 0
    invalid = QtCore.QItemSelection()
    for index in tmp_selection.indexes():
      if index.data(enums.ModelEnum.TYPE_ROLE) == "protein":
        if protein_count > 0:
          invalid.select(index, index)
        protein_count += 1
    if protein_count > 1:
      tmp_selection_model.select(invalid, QtCore.QItemSelectionModel.Deselect)
      self._evaluate_and_update_add_button(tmp_selection)
      return
    # </editor-fold>

    # <editor-fold desc="Restrict selection to a single parent level">
    parent = self._view.ui.tree_prot_2.currentIndex().parent()
    invalid = QtCore.QItemSelection()
    for index in tmp_selection.indexes():
      if index.parent() != parent:
        invalid.select(index, index)
    tmp_selection_model.select(invalid, QtCore.QItemSelectionModel.Deselect)
    # </editor-fold>

    self._evaluate_and_update_add_button(tmp_selection)

  def __slot_collapse_all_tree_prot_1(
      self, the_selected_index: QtCore.QModelIndex
  ) -> None:
    """Collapses all protein nodes in tree_prot_1 except for the expanded one.

    When a protein node is expanded the Chains header (row 1) is also expanded
    automatically.  Uses `model.index()` for child navigation because
    `QModelIndex.child()` was removed in PyQt6.

    Args:
        the_selected_index (QtCore.QModelIndex): The index of the item that was
            just expanded.

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
    tmp_model = self._view.ui.tree_prot_1.model()
    tmp_type = the_selected_index.data(enums.ModelEnum.TYPE_ROLE)
    if tmp_type == "protein":
      tmp_index_to_check = the_selected_index
      # Row 1 under a protein node is the "Chains" header.
      chains_header_index = tmp_model.index(1, 0, tmp_index_to_check)
      self._view.ui.tree_prot_1.setExpanded(chains_header_index, True)
    elif tmp_type == "header":
      tmp_index_to_check = the_selected_index.parent()
    elif tmp_type == "chain":
      tmp_index_to_check = the_selected_index.parent().parent()
    else:
      tmp_index_to_check = None

    tmp_parent = self._local_model.get_proteins_section_index()
    for tmp_row in range(tmp_model.rowCount(tmp_parent)):
      tmp_index = tmp_model.index(tmp_row, 0, tmp_parent)
      if (
          tmp_index != tmp_index_to_check
          and tmp_index.data(enums.ModelEnum.TYPE_ROLE) == "protein"
      ):
        self._view.ui.tree_prot_1.collapse(tmp_index)

  def __slot_collapse_all_tree_prot_2(
      self, the_selected_index: QtCore.QModelIndex
  ) -> None:
    """Collapses all protein nodes in tree_prot_2 except for the expanded one.

    When a protein node is expanded the Chains header (row 1) is also expanded
    automatically.  Uses `model.index()` for child navigation because
    `QModelIndex.child()` was removed in PyQt6.

    Args:
        the_selected_index (QtCore.QModelIndex): The index of the item that was
            just expanded.

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
        "An object of the tree view 2 was expanded.",
    )
    tmp_model = self._view.ui.tree_prot_2.model()
    tmp_type = the_selected_index.data(enums.ModelEnum.TYPE_ROLE)
    if tmp_type == "protein":
      tmp_index_to_check = the_selected_index
      # Row 1 under a protein node is the "Chains" header.
      chains_header_index = tmp_model.index(1, 0, tmp_index_to_check)
      self._view.ui.tree_prot_2.setExpanded(chains_header_index, True)
    elif tmp_type == "header":
      tmp_index_to_check = the_selected_index.parent()
    elif tmp_type == "chain":
      tmp_index_to_check = the_selected_index.parent().parent()
    else:
      tmp_index_to_check = None

    tmp_parent = self._local_model.get_proteins_section_index()
    for tmp_row in range(tmp_model.rowCount(tmp_parent)):
      tmp_index = tmp_model.index(tmp_row, 0, tmp_parent)
      if (
          tmp_index != tmp_index_to_check
          and tmp_index.data(enums.ModelEnum.TYPE_ROLE) == "protein"
      ):
        self._view.ui.tree_prot_2.collapse(tmp_index)

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
    
    _, prot_1_chains = self._get_protein_name_and_chains(
        self._view.ui.tree_prot_1.selectedIndexes()
    )
    self._number_of_prot_1_selected_chains = len(prot_1_chains)
    
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
    if self._on_add_callback:
        self._on_add_callback((tmp_item, True))

  # </editor-fold>
