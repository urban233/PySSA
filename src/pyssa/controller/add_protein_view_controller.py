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
"""Module for the add protein view controller."""
import logging
import os
from typing import TYPE_CHECKING

from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.qt import Qt
from src.pyssa.gui.qt import QtWidgets
from src.pyssa.gui.ui.views import add_protein_view
from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.internal.thread import tasks
from src.pyssa.internal.thread.async_pyssa import validate_async
from src.pyssa.logging_pyssa import log_levels, log_handlers
from src.pyssa.util import constants, tools, exception

if TYPE_CHECKING:
  from src.pyssa.gui import app_state

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class AddProteinViewController(QtCore.QObject):
  """Class for the AddProteinViewController."""

  def __init__(
      self, the_app_state: "app_state.AppState", a_parent=None
  ) -> None:
    """Constructor.

    Args:
        the_app_state (app_state.AppState): The AppState object.
        a_parent: Parent widget to pass to the view.

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
    self._view = add_protein_view.AddProteinView(a_parent)
    self._active_task = None
    self._connect_all_ui_elements_to_slot_functions()

  def get_view(self):
    # # check internet connectivity
    # if not tools.check_internet_connectivity():
    #   tmp_dialog = custom_message_box.CustomMessageBoxOk(
    #       "You do not have a working internet connection which is "
    #       "necessary for connecting to the PDB!\n"
    #       "However you can add a protein structure from "
    #       "your local filesystem.",
    #       "Internet Connection",
    #       custom_message_box.CustomMessageBoxIcons.ERROR.value,
    #   )
    #   tmp_dialog.exec()
    #   self._view.ui.txt_add_protein.setEnabled(False)
    #   self._view.ui.lbl_status.setText(
    #       "You cannot enter a PDB ID (no working internet connection)."
    #   )
    return self._view

  def _connect_all_ui_elements_to_slot_functions(self) -> None:
    """Connects all UI elements to their corresponding slot functions in the class."""
    self._view.ui.btn_choose_protein.clicked.connect(
        self.__slot_load_protein_from_filesystem
    )
    self._view.ui.btn_add_protein.clicked.connect(self.__slot_add_protein)
    self._view.ui.txt_add_protein.textChanged.connect(
        self.__slot_validate_input
    )
    self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)

  def restore_ui(self) -> None:
    """Restores the UI."""
    self._view.ui.txt_add_protein.clear()
    self._view.ui.txt_add_protein.setStyleSheet(
        """QLineEdit {color: #000000; border-color: #DCDBE3;}""",
    )
    self._view.ui.lbl_status.setText("")
    self._view.ui.btn_add_protein.setEnabled(False)
    self._view.setMinimumWidth(500)

  def _open_help_for_dialog(self) -> None:
    """Opens the help dialog for the corresponding dialog."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Help' button was clicked."
    )
    # self._interface_manager.help_manager.open_protein_import_page()

  # @SLOT
  def __slot_validate_input(self, the_entered_text: str) -> None:
    """Checks if the entered reference protein is valid or not.

    Args:
        the_entered_text (str): The text entered by the user.

    Raises:
        exception.IllegalArgumentError: If `the_entered_text` is None.
    """
    if the_entered_text is None:
      logger.error("the_entered_text is None.")
      raise exception.IllegalArgumentError("the_entered_text is None.")

    logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "A text was entered.")
    self._view.ui.lbl_status.setStyleSheet(
        """QLabel {color: #ba1a1a;}""",
    )
    if len(the_entered_text) == 0:
      # empty line edit
      self._view.ui.txt_add_protein.setStyleSheet(
          """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}""",
      )
      self._view.ui.lbl_status.setText(
          "Please enter a PDB id or choose an existing .pdb file from your filesystem!"
      )
      self._view.ui.btn_add_protein.setEnabled(False)
    elif len(the_entered_text) < 4:
      # length of text is too small
      self._view.ui.txt_add_protein.setStyleSheet(
          """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}""",
      )
      self._view.ui.btn_add_protein.setEnabled(False)
      self._view.ui.lbl_status.setText("Please enter more characters.")
    elif len(the_entered_text) > 4 and not os.path.exists(the_entered_text):
      self._view.ui.txt_add_protein.setText(
          the_entered_text[: len(the_entered_text) - 1]
      )
    # checks if a pdb id was entered
    else:
      # QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
      self._view.ui.txt_add_protein.setStyleSheet(
        """QLineEdit {color: #000000; border-color: #DCDBE3;}""",
      )
      self._view.ui.lbl_status.setStyleSheet(
        """QLabel {color: #367AF6;}""",
      )
      self._view.ui.lbl_status.setText("")
      self._view.ui.btn_add_protein.setEnabled(True)


  def __await__slot_validate_input(self, return_value: tuple) -> None:
    """Validates the input entered by the user.

    Args:
        return_value (tuple): A tuple containing information about the validation result.
            The tuple should have the following structure:
                - return_value[0] (int): The type of input (1 for pdb id entered, 2 for filepath entered).
                - return_value[1] (bool): Flag indicating whether the input is valid or not.
                - return_value[2] (str): The name entered by the user.

    Raises:
        exception.IllegalArgumentError: If `return_value` is None.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      raise exception.IllegalArgumentError("return_value is None.")

    # </editor-fold>

    if return_value[0] == -1:
      QtWidgets.QApplication.restoreOverrideCursor()
      return

    tmp_type: int = return_value[0]
    tmp_is_valid: bool = return_value[1]
    tmp_name: str = return_value[2]
    self._view.ui.lbl_status.setStyleSheet(
        """QLabel {color: #ba1a1a;}""",
    )
    if tmp_is_valid:
      self._view.ui.txt_add_protein.setStyleSheet(
          """QLineEdit {color: #000000; border-color: #DCDBE3;}""",
      )
      self._view.ui.lbl_status.setText("")
      self._view.ui.btn_add_protein.setEnabled(True)
    elif not tmp_is_valid and tmp_type == 1:  # pdb id entered
      self._view.ui.txt_add_protein.setStyleSheet(
          """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}""",
      )
      # self._view.ui.lbl_status.setText(
      #   "There is an invalid PDB ID or NO internet connection is available!"
      # )
      self._view.ui.btn_add_protein.setEnabled(False)
    elif not tmp_is_valid and tmp_type == 2:  # filepath entered
      self._view.ui.txt_add_protein.setStyleSheet(
          """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}""",
      )
      self._view.ui.lbl_status.setText("Invalid filepath!")
      self._view.ui.btn_add_protein.setEnabled(False)
    else:
      constants.PYSSA_LOGGER.error(
          "There is an unknown case, while validating the add protein view user input!"
      )
    if tmp_name in [p.get_molecule_object() for p in self._app_state.project.proteins]:
      self._view.ui.txt_add_protein.setStyleSheet(
          """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}""",
      )
      self._view.ui.lbl_status.setText(
          "Protein already exists in current project!"
      )
      self._view.ui.btn_add_protein.setEnabled(False)
    QtWidgets.QApplication.restoreOverrideCursor()
    self._view.ui.txt_add_protein.setFocus()

  def __slot_load_protein_from_filesystem(self) -> None:
    """Loads a protein from the filesystem into the textbox."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Load protein from filesystem' button was clicked.",
    )
    try:
      # open file dialog
      file_name = QtWidgets.QFileDialog.getOpenFileName(
          self._view,
          "Open existing protein",
          QtCore.QDir.homePath(),
          "PDB Files (*.pdb)",
      )
      if file_name == ("", ""):
        self._view.ui.lbl_status.setText("No file has been selected.")
      else:
        # display path in text box
        self._view.ui.txt_add_protein.setText(str(file_name[0]))
        self._view.ui.btn_add_protein.setEnabled(True)
    except FileNotFoundError:
      self._view.ui.lbl_status.setText("Loading the protein structure failed!")

  def __slot_add_protein(self) -> None:
    """Adds a protein to the project and closes the dialog."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Add' button was clicked."
    )
    from src.pyssa.internal.thread.thread_api import thread_runtime
    from src.pyssa.internal.data_structures import protein
    import pathlib

    tmp_protein_name = self._view.ui.txt_add_protein.text()
    tmp_name_len = len(tmp_protein_name)

    if not tmp_protein_name:
        logger.error("No protein name to import.")
        QtWidgets.QMessageBox.critical(self._view, "Import Error", "No data received!")
        return

    self._view.ui.btn_add_protein.setEnabled(False)

    def import_task(progress_callback, is_cancelled):
        # We can load directly through user_pymol which AppState doesn't hold directly,
        # but PyMOL runs in the same process instance namespace basically via the global cmd module.
        # It's cleaner to access PyMOL directly here or via an established interface.
        if tmp_name_len == 4:
            pdb_name = tmp_protein_name.upper()
            tmp_ref_protein = protein.Protein(pdb_name)
            tmp_ref_protein.set_id(0) # Let the model map it later
            tmp_ref_protein.db_project_id = self._app_state.project.get_id()
            if not tools.check_internet_connectivity():
              raise RuntimeError("No internet connection is available!")
            try:
              tmp_ref_protein.add_protein_structure_data_from_pdb_db(pdb_name)
            except Exception as tmp_exception:
              logger.error(tmp_exception)
              raise RuntimeError("PDB ID is invalid!")
        else:
            pdb_filepath = pathlib.Path(tmp_protein_name)
            pdb_name = pdb_filepath.name.replace(".pdb", "")
            tmp_ref_protein = protein.Protein(pdb_name)
            tmp_ref_protein.set_id(0)
            tmp_ref_protein.db_project_id = self._app_state.project.get_id()
            tmp_ref_protein.add_protein_structure_data_from_local_pdb_file(pdb_filepath)

        tmp_ref_protein.create_new_pymol_session()
        tmp_protein_id = self._app_state.hot_db.insert_protein_full(tmp_ref_protein)
        tmp_ref_protein.set_id(tmp_protein_id)
        return tmp_ref_protein
        
    def on_success(result: "protein.Protein"):
        self._app_state.project.add_existing_protein(result)
        # Use incremental update instead of full rebuild to preserve tree state
        self._app_state.pyssa_objects_model.add_protein(result)
        self._app_state.hot_db.insert_protein_full(result)
        self._app_state.status_bar_manager.show_permanent_message("", False)
        self._app_state.status_bar_manager.show_temporary_message("Protein imported.")
        
    def on_error(an_exception):
        logger.exception("Failed to insert protein.", exc_info=an_exception)
        QtWidgets.QMessageBox.critical(
          self._view,
          "Import Error",
          f"An error occurred: {an_exception}"
        )
        self._view.ui.btn_add_protein.setEnabled(True)
        self._app_state.status_bar_manager.show_error_message(
          "Import protein failed!",
          True
        )

    thread_runtime.get_singleton_thread_runtime().run(import_task).on_success(on_success).on_error(on_error)
    self._app_state.status_bar_manager.show_permanent_message(
      "Importing protein ...", True
    )
    self._view.close()
