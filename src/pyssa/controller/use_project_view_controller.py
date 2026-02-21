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
import pathlib

from src.pyssa.gui import app_state
from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.qt import Qt
from src.pyssa.gui.qt import QtWidgets
from src.pyssa.gui.ui.views import use_project_view

from src.pyssa.io_pyssa.db_pyssa import ProjectDatabase
from src.pyssa.util import input_validator, constants, enums, exception
from src.pyssa.util import gui_utils
from src.pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class UseProjectViewController(QtCore.QObject):
  """Class for the UseProjectViewController."""

  def __init__(
      self, the_app_state: "app_state.AppState"
  ) -> None:
    """Constructor.

    Args:
        the_app_state (app_state.AppState): The AppState object.

    Raises:
        exception.IllegalArgumentError: If `the_app_state` is None.
    """
    # <editor-fold desc="Checks">
    if the_app_state is None:
      logger.error("the_app_state is None.")
      raise exception.IllegalArgumentError("the_app_state is None.")

    # </editor-fold>

    super().__init__()
    self._app_state = the_app_state
    self._view = use_project_view.UseProjectView()
    self._initialize_ui()
    self._fill_projects_list_view()
    self._fill_projects_combobox()
    self._project_names: set = self._convert_model_into_set()
    self._connect_all_ui_elements_to_slot_functions()

  def get_view(self):
    return self._view

  def restore_default_view(self) -> None:
    self._initialize_ui()
    self._fill_projects_list_view()
    self._fill_projects_combobox()

  def _open_help_for_dialog(self) -> None:
    """Opens the help dialog for the corresponding dialog."""
    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Help' button was clicked."
    )

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
        self._app_state.workspace.get_model()
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
        self._app_state.workspace.get_projects_as_string_list(),
    )
    current_project_name = ""
    if self._app_state.has_open_project():
        current_project_name = self._app_state.project.get_project_name()

    self._view.ui.cb_choose_project.setCurrentIndex(
        self._view.ui.cb_choose_project.findText(
            current_project_name
        ),
    )

  def _list_all_proteins_of_selected_project(self) -> None:
    """Lists all proteins of the selected project."""
    if self._view.ui.cb_choose_project.currentText() == "":
      return
    
    selected_project_name = self._view.ui.cb_choose_project.currentText()
    tmp_database_filepath: str = str(
       self._app_state.workspace.construct_project_db_path(selected_project_name)
    )

    source_db = ProjectDatabase(db_path=tmp_database_filepath, project_id=selected_project_name)
    
    try:
      tmp_project = source_db.service.load_project(
          selected_project_name, 
          self._app_state.get_settings().workspace_path,
          self._app_state.get_settings()
      )
    except Exception as e:
       logger.error(f"Failed to load source project in _list_all_proteins: {e}")
       source_db.close()
       return
    
    source_db.close()

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

    project_name = self._view.ui.txt_use_project_name.text()
    db_path = str(self._app_state.workspace.construct_project_db_path(project_name))

    self._set_ui_loading(True)

    from src.pyssa.io_pyssa.db_pyssa import ProjectDatabase
    from src.pyssa.model import psa_objects_model
    from src.pyssa.internal.thread.thread_api import thread_runtime
    import platform
    import copy

    def create_project(progress_callback, is_cancelled):
      tmp_db = ProjectDatabase(db_path=db_path, project_id=project_name)
      if is_cancelled():
        tmp_db.close()
        raise InterruptedError("Cancelled before creating project database.")
      
      tmp_db.initialise_schema()
      project_id = tmp_db.insert_project(name=project_name, os=platform.system())

      from src.pyssa.internal.data_structures import project
      import pathlib
      tmp_project = project.Project(project_name, pathlib.Path(self._app_state.get_settings().workspace_path))
      tmp_project.set_id(project_id)

      if is_cancelled():
        tmp_db.close()
        raise InterruptedError("Cancelled after creating project.")

      # Insert the copied proteins into the db, and associate them with the new project
      for idx, tmp_protein in enumerate(tmp_proteins):
        if progress_callback:
           progress_callback(f"Adding protein {idx+1}/{len(tmp_proteins)}...", int(idx/len(tmp_proteins) * 100))

        tmp_protein_copy = copy.deepcopy(tmp_protein)
        tmp_protein_copy.db_project_id = project_id
        tmp_project.add_existing_protein(tmp_protein_copy)
        tmp_protein_copy.set_id(tmp_db.insert_protein_full(tmp_protein_copy))

      tmp_pyssa_objects_model = psa_objects_model.PSAObjectsModel()
      tmp_pyssa_objects_model.build_model(tmp_project)
      return tmp_project, tmp_db, tmp_pyssa_objects_model

    def on_success(result):
      tmp_project, tmp_db, tmp_pyssa_objects_model = result
      self._app_state.pyssa_objects_model = tmp_pyssa_objects_model
      self._app_state.open_project(tmp_project, tmp_db)

    def on_error(exc):
      logger.exception("Failed to create use project.", exc_info=exc)
      self._set_ui_loading(False)
      QtWidgets.QMessageBox.critical(
        self._view,
        "Failed to create project",
        f"Could not create the project:\n{exc}",
      )

    (
      thread_runtime.get_singleton_thread_runtime()
      .run(create_project)
      .on_success(on_success)
      .on_error(on_error)
    )
    self._view.close()

  def _set_ui_loading(self, loading: bool) -> None:
     self._view.ui.btn_use_create_new_project.setEnabled(not loading)
     self._view.ui.btn_use_back.setEnabled(not loading)
     self._view.ui.list_use_selected_protein_structures.setEnabled(not loading)
     self._view.ui.list_use_available_protein_structures.setEnabled(not loading)

