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
"""Module for the use project view controller."""
import logging
import os
import pathlib
import subprocess

import pygetwindow
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from PyQt5 import QtWidgets
from src.pyssa.controller import database_manager
from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.internal.thread import tasks
from src.pyssa.internal.thread.async_pyssa import util_async
from src.pyssa.util import input_validator, constants, enums, exception
from src.pyssa.util import gui_utils
from src.pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class UseProjectViewController(QtCore.QObject):
  """Class for the UseProjectViewController."""

  user_input = QtCore.pyqtSignal(tuple)
  """Singal used to transfer data back to the previous window."""

  def __init__(
      self, the_interface_manager: "interface_manager.InterfaceManager"
  ) -> None:
    """Constructor.

    Args:
        the_interface_manager (interface_manager.InterfaceManager): The InterfaceManager object.

    Raises:
        exception.IllegalArgumentError: If `the_interface_manager` is None.
    """
    # <editor-fold desc="Checks">
    if the_interface_manager is None:
      logger.error("the_interface_manager is None.")
      raise exception.IllegalArgumentError("the_interface_manager is None.")

    # </editor-fold>

    super().__init__()
    self._interface_manager = the_interface_manager
    self._view = the_interface_manager.get_use_project_view()
    self._initialize_ui()
    self._fill_projects_list_view()
    self._fill_projects_combobox()
    self._project_names: set = self._convert_model_into_set()
    self._connect_all_ui_elements_to_slot_functions()

  def _open_help_for_dialog(self) -> None:
    """Opens the help dialog for the corresponding dialog."""
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Help' button was clicked."
    )
    self._interface_manager.help_manager.open_use_project_page()

  def _convert_model_into_set(self) -> set:
    """Converts the model data into a set of project names.

    Returns:
        set: A set containing the project names from the model.
    """
    tmp_project_names = []
    for tmp_row in range(
        self._view.ui.list_use_existing_projects.model().rowCount()
    ):
      tmp_project_names.append(
          self._view.ui.list_use_existing_projects.model()
          .index(tmp_row, 0)
          .data(Qt.DisplayRole),
      )
    return set(tmp_project_names)

  def _initialize_ui(self) -> None:
    """Initializes the user interface for the application."""
    gui_elements_to_show = [
        self._view.ui.btn_use_next,
        self._view.ui.list_use_existing_projects,
        self._view.ui.label,
        self._view.ui.lbl_use_project_name,
    ]
    gui_utils.show_gui_elements(gui_elements_to_show)
    self._view.ui.txt_use_project_name.setEnabled(True)

    gui_elements_to_hide = [
        self._view.ui.lbl_use_search,
        self._view.ui.lbl_use_status_search,
        self._view.ui.txt_use_search,
        self._view.ui.lbl_choose_project,
        self._view.ui.cb_choose_project,
        self._view.ui.btn_use_add_available_protein_structures,
        self._view.ui.lbl_use_available_protein_structures,
        self._view.ui.list_use_available_protein_structures,
        self._view.ui.btn_use_remove_selected_protein_structures,
        self._view.ui.lbl_use_selected_protein_structures,
        self._view.ui.list_use_selected_protein_structures,
        self._view.ui.btn_use_back,
        self._view.ui.btn_use_create_new_project,
    ]
    gui_utils.hide_gui_elements(gui_elements_to_hide)
    gui_utils.enable_text_box(
        self._view.ui.txt_use_project_name, self._view.ui.lbl_use_project_name
    )
    self._view.ui.txt_use_project_name.clear()
    self._view.ui.txt_use_project_name.setStyleSheet(
        """QLineEdit {color: #000000; border-color: #DCDBE3;}""",
    )
    self._view.ui.lbl_use_status_project_name.setText("")
    self._view.ui.lbl_use_status_project_name.setStyleSheet(
        "color: #ba1a1a; font-size: 11px;"
    )
    self._view.ui.txt_use_search.clear()
    self._view.ui.list_use_selected_protein_structures.clear()
    self._view.ui.lbl_use_status_search.setText("")
    self._view.ui.btn_use_next.setEnabled(False)
    self._temporary_redesign()

  def _temporary_redesign(self) -> None:
    """Changes some parts of the UI temporarily."""
    self._view.ui.lbl_use_search.hide()
    self._view.ui.lbl_use_status_search.hide()
    self._view.ui.txt_use_search.hide()

  def _fill_projects_list_view(self) -> None:
    """Lists all projects."""
    self._view.ui.list_use_existing_projects.setModel(
        self._interface_manager.get_workspace_projects()
    )

  def _connect_all_ui_elements_to_slot_functions(self) -> None:
    """Connects all UI elements to their corresponding slot functions in the class."""
    self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)
    self._view.ui.txt_use_project_name.textChanged.connect(
        self.validate_use_project_name
    )
    self._view.ui.btn_use_next.clicked.connect(
        self.show_protein_selection_for_use
    )
    self._view.ui.btn_use_add_available_protein_structures.clicked.connect(
        self.add_protein_structure_to_new_project,
    )
    self._view.ui.cb_choose_project.currentIndexChanged.connect(
        self._list_all_proteins_of_selected_project
    )
    self._view.ui.list_use_available_protein_structures.itemClicked.connect(
        self.use_enable_add
    )
    self._view.ui.btn_use_remove_selected_protein_structures.clicked.connect(
        self.remove_protein_structure_to_new_project,
    )
    self._view.ui.list_use_selected_protein_structures.itemClicked.connect(
        self.use_enable_remove
    )
    self._view.ui.btn_use_back.clicked.connect(
        self.hide_protein_selection_for_use
    )
    self._view.ui.btn_use_create_new_project.clicked.connect(
        self.create_use_project
    )

  def validate_use_project_name(self, the_entered_text: str) -> None:
    """Validates the entered text for project name.

    Args:
        the_entered_text (str): The text entered for project name.

    Raises:
        exception.IllegalArgumentError: If `the_entered_text` is None.
    """
    # <editor-fold desc="Checks">
    if the_entered_text is None:
      logger.error("the_entered_text is None.")
      raise exception.IllegalArgumentError("the_entered_text is None.")

    # </editor-fold>

    logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "A text was entered.")
    projects_list_view = self._view.ui.list_use_existing_projects
    # Deselect any current item in the list view
    if projects_list_view.currentIndex().isValid():
      projects_list_view.selectionModel().clearSelection()

    tmp_validate_flag, tmp_stylesheet_string, tmp_message = (
        input_validator.validate_input_for_project_name(
            the_entered_text,
            self._project_names,
        )
    )
    self._view.ui.txt_use_project_name.setStyleSheet(tmp_stylesheet_string)

    if tmp_validate_flag:
      self._view.ui.lbl_use_status_project_name.setText("")
      self._view.ui.btn_use_next.setEnabled(True)
    else:
      self._view.ui.lbl_use_status_project_name.setText(tmp_message)
      self._view.ui.btn_use_next.setEnabled(False)

  def validate_use_search(self) -> None:
    """Validates the input of the protein name in real-time."""
    message = "Protein structure does not exists."
    input_validator.InputValidator.validate_search_input(
        self._view.ui.list_use_available_protein_structures,
        self._view.ui.txt_use_search,
        self._view.ui.lbl_use_status_search,
        status_message=message,
    )

  def add_protein_structure_to_new_project(self) -> None:
    """Adds the selected protein to the list which is used to create the new project."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Add' protein to new project button was clicked.",
    )
    prot_to_add = QtWidgets.QListWidgetItem(
        self._view.ui.list_use_available_protein_structures.currentItem().text()
    )
    prot_to_add.setData(
        enums.ModelEnum.OBJECT_ROLE,
        self._view.ui.list_use_available_protein_structures.currentItem().data(
            enums.ModelEnum.OBJECT_ROLE
        ),
    )
    self._view.ui.list_use_selected_protein_structures.addItem(prot_to_add)
    self._view.ui.list_use_available_protein_structures.takeItem(
        self._view.ui.list_use_available_protein_structures.currentRow(),
    )
    self._view.ui.btn_use_add_available_protein_structures.setEnabled(False)
    if self._view.ui.list_use_available_protein_structures.count() > 0:
      try:
        self._view.ui.list_use_available_protein_structures.currentItem().setSelected(
            False
        )
      except AttributeError:
        constants.PYSSA_LOGGER.debug(
            "No selection in use available proteins list on Use page."
        )

    self._view.ui.btn_use_create_new_project.setEnabled(True)

  def remove_protein_structure_to_new_project(self) -> None:
    """Removes the selected protein from the list which is used to create the new project."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Remove' button was clicked."
    )
    prot_to_remove = (
        self._view.ui.list_use_selected_protein_structures.currentItem()
    )
    self._view.ui.list_use_selected_protein_structures.takeItem(
        self._view.ui.list_use_selected_protein_structures.currentRow(),
    )
    self._view.ui.list_use_available_protein_structures.addItem(prot_to_remove)
    self._view.ui.btn_use_remove_selected_protein_structures.setEnabled(False)
    if self._view.ui.list_use_selected_protein_structures.count() > 0:
      try:
        self._view.ui.list_use_selected_protein_structures.currentItem().setSelected(
            False
        )
      except AttributeError:
        constants.PYSSA_LOGGER.debug(
            "No selection in use selected proteins list on Use page."
        )

    if self._view.ui.list_use_selected_protein_structures.count() == 0:
      self._view.ui.btn_use_create_new_project.setEnabled(False)
    else:
      self._view.ui.btn_use_create_new_project.setEnabled(True)

  def _are_duplicate_proteins_in_selected_protein_list(self) -> bool:
    """Checks if there are duplicate proteins in the selected protein list.

    Returns:
        True if there are duplicate proteins, False otherwise.
    """
    try:
      target_string = (
          self._view.ui.list_use_available_protein_structures.currentItem().text()
      )
    except AttributeError:
      return False

    tmp_occurrences = 0
    for tmp_row in range(
        self._view.ui.list_use_selected_protein_structures.count()
    ):
      if (
          self._view.ui.list_use_selected_protein_structures.item(
              tmp_row
          ).text()
          == target_string
      ):
        tmp_occurrences += 1
    if tmp_occurrences > 0:
      return True
    return False

  def show_protein_selection_for_use(self) -> None:
    """Shows the two lists for the protein selection."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Next' button was clicked."
    )
    gui_elements_to_show = [
        self._view.ui.lbl_use_search,
        self._view.ui.lbl_use_status_search,
        self._view.ui.txt_use_search,
        self._view.ui.lbl_choose_project,
        self._view.ui.cb_choose_project,
        self._view.ui.btn_use_add_available_protein_structures,
        self._view.ui.lbl_use_available_protein_structures,
        self._view.ui.list_use_available_protein_structures,
        self._view.ui.btn_use_remove_selected_protein_structures,
        self._view.ui.lbl_use_selected_protein_structures,
        self._view.ui.list_use_selected_protein_structures,
        self._view.ui.btn_use_back,
        self._view.ui.btn_use_create_new_project,
        self._view.ui.lbl_use_project_name,
    ]
    gui_utils.show_gui_elements(gui_elements_to_show)
    self._view.ui.txt_use_project_name.setEnabled(False)
    gui_elements_to_hide = [
        self._view.ui.btn_use_next,
        self._view.ui.list_use_existing_projects,
        self._view.ui.label,
    ]
    gui_utils.hide_gui_elements(gui_elements_to_hide)
    gui_utils.disable_text_box(
        self._view.ui.txt_use_project_name, self._view.ui.lbl_use_project_name
    )
    self._view.ui.btn_use_add_available_protein_structures.setEnabled(False)
    self._view.ui.btn_use_remove_selected_protein_structures.setEnabled(False)
    self._view.ui.btn_use_create_new_project.setEnabled(False)
    self._temporary_redesign()
    self._list_all_proteins_of_selected_project()

    if self._are_duplicate_proteins_in_selected_protein_list():
      self._view.ui.btn_use_create_new_project.setEnabled(False)
    else:
      self._view.ui.btn_use_create_new_project.setEnabled(True)

  def hide_protein_selection_for_use(self) -> None:
    """Hides the two lists for the protein selection."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Back' button was clicked."
    )
    gui_elements_to_show = [
        self._view.ui.btn_use_next,
        self._view.ui.list_use_existing_projects,
        self._view.ui.label,
        self._view.ui.lbl_use_project_name,
    ]
    gui_utils.show_gui_elements(gui_elements_to_show)
    self._view.ui.txt_use_project_name.setEnabled(True)

    gui_elements_to_hide = [
        self._view.ui.lbl_use_search,
        self._view.ui.lbl_use_status_search,
        self._view.ui.txt_use_search,
        self._view.ui.lbl_choose_project,
        self._view.ui.cb_choose_project,
        self._view.ui.btn_use_add_available_protein_structures,
        self._view.ui.lbl_use_available_protein_structures,
        self._view.ui.list_use_available_protein_structures,
        self._view.ui.btn_use_remove_selected_protein_structures,
        self._view.ui.lbl_use_selected_protein_structures,
        self._view.ui.list_use_selected_protein_structures,
        self._view.ui.btn_use_back,
        self._view.ui.btn_use_create_new_project,
    ]
    gui_utils.hide_gui_elements(gui_elements_to_hide)
    gui_utils.enable_text_box(
        self._view.ui.txt_use_project_name, self._view.ui.lbl_use_project_name
    )
    self._temporary_redesign()

  def _fill_projects_combobox(self) -> None:
    """Fills the combo box with the available projects from the workspace."""
    gui_utils.fill_combo_box(
        self._view.ui.cb_choose_project,
        self._interface_manager.get_workspace_projects_as_list(),
    )
    self._view.ui.cb_choose_project.setCurrentIndex(
        self._view.ui.cb_choose_project.findText(
            self._interface_manager.get_current_project().get_project_name()
        ),
    )

  def _list_all_proteins_of_selected_project(self) -> None:
    """Lists all proteins of the selected project."""
    if self._view.ui.cb_choose_project.currentText() == "":
      return
    tmp_database_filepath: str = str(
        pathlib.Path(
            f"{self._interface_manager.get_application_settings().get_workspace_path()}/{self._view.ui.cb_choose_project.currentText()}.db",
        ),
    )
    with database_manager.DatabaseManager(tmp_database_filepath) as db_manager:
      tmp_project = db_manager.get_project_as_object(
          self._view.ui.cb_choose_project.currentText(),
          self._interface_manager.get_application_settings().get_workspace_path(),
          self._interface_manager.get_application_settings(),
      )
      for tmp_protein in tmp_project.proteins:
        tmp_pdb_atom_db_data = db_manager.get_pdb_atoms_of_protein(
            tmp_protein.get_id()
        )
        pdb_atom_dict = [
            {key.value: value for key, value in zip(enums.PdbAtomEnum, t)}
            for t in tmp_pdb_atom_db_data
        ]
        tmp_protein.set_pdb_data(pdb_atom_dict)
    self._view.ui.list_use_available_protein_structures.clear()
    for tmp_protein in tmp_project.proteins:
      item = QtWidgets.QListWidgetItem(tmp_protein.get_molecule_object())
      item.setData(enums.ModelEnum.OBJECT_ROLE, tmp_protein)
      self._view.ui.list_use_available_protein_structures.addItem(item)

    if self._are_duplicate_proteins_in_selected_protein_list():
      self._view.ui.btn_use_create_new_project.setEnabled(False)
    else:
      self._view.ui.btn_use_create_new_project.setEnabled(True)

  def use_enable_add(self) -> None:
    """Enables the add button."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "A protein from the list of available proteins was clicked.",
    )
    if self._are_duplicate_proteins_in_selected_protein_list():
      self._view.ui.btn_use_add_available_protein_structures.setEnabled(False)
    else:
      self._view.ui.btn_use_add_available_protein_structures.setEnabled(True)

  def use_enable_remove(self) -> None:
    """Enables the remove button."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "A protein from the list of proteins to add to the new project was clicked.",
    )
    self._view.ui.btn_use_remove_selected_protein_structures.setEnabled(True)

  def create_use_project(self) -> None:
    """Uses the project by sending the `user_input` signal and closing the dialog."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Create' button was clicked."
    )
    tmp_proteins: list = []
    for tmp_row in range(
        self._view.ui.list_use_selected_protein_structures.count()
    ):
      tmp_proteins.append(
          self._view.ui.list_use_selected_protein_structures.item(tmp_row).data(
              enums.ModelEnum.OBJECT_ROLE
          )
      )
    self.user_input.emit(
        (self._view.ui.txt_use_project_name.text(), tmp_proteins)
    )
    self._view.close()
