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
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from PyQt5 import QtWidgets
from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.internal.thread import tasks
from src.pyssa.internal.thread.async_pyssa import validate_async
from src.pyssa.logging_pyssa import log_levels, log_handlers
from src.pyssa.util import constants, tools, exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class AddProteinViewController(QtCore.QObject):
  """Class for the AddProteinViewController."""

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
    self._view = the_interface_manager.get_add_protein_view()
    self._connect_all_ui_elements_to_slot_functions()
    # check internet connectivity
    if not tools.check_internet_connectivity():
      tmp_dialog = custom_message_box.CustomMessageBoxOk(
          "You do not have a working internet connection which is "
          "necessary for connecting to the PDB!\n"
          "However you can add a protein structure from "
          "your local filesystem.",
          "Internet Connection",
          custom_message_box.CustomMessageBoxIcons.ERROR.value,
      )
      tmp_dialog.exec_()
      self._view.ui.txt_add_protein.setEnabled(False)
      self._view.ui.lbl_status.setText(
          "You cannot enter a PDB ID (no working internet connection)."
      )

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
    self._interface_manager.help_manager.open_protein_import_page()

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
      self._active_task = tasks.LegacyTask(
          target=validate_async.validate_add_protein_view_input,
          args=(
              the_entered_text,
              0,
          ),
          post_func=self.__await__slot_validate_input,
      )
      self._active_task.start()
      QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
      self._view.ui.txt_add_protein.setStyleSheet(
          """QLineEdit {color: #000000; border-color: #DCDBE3;}""",
      )
      self._view.ui.lbl_status.setStyleSheet(
          """QLabel {color: #367AF6;}""",
      )
      self._view.ui.lbl_status.setText("Checking input ...")

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
      self._view.ui.lbl_status.setText("Invalid PDB id!")
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
    if tmp_name in self._interface_manager.watcher.protein_names_blacklist:
      self._view.ui.txt_add_protein.setStyleSheet(
          """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}""",
      )
      self._view.ui.lbl_status.setText(
          "Protein already exists in current project!"
      )
      self._view.ui.btn_add_protein.setEnabled(False)
      self._view.ui.lbl_status.setText("")
      self._view.ui.btn_add_protein.setEnabled(True)
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
    """Adds a protein to the global variable and closes the dialog."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Add' button was clicked."
    )
    self._view.close()
    self.user_input.emit(
        (
            self._view.ui.txt_add_protein.text(),
            len(self._view.ui.txt_add_protein.text()),
        )
    )

  def _validate_scene_name(self, text: str) -> None:
    """Validates the given scene name and updates the UI elements accordingly.

    Args:
        text (str): The scene name to be validated.

    Raises:
        exception.IllegalArgumentError: If `text` is None.
    """
    # <editor-fold desc="Checks">
    if text is None:
      logger.error("text is None.")
      raise exception.IllegalArgumentError("text is None.")

    # </editor-fold>

    new_text = "".join(char for char in text)
    self._view.line_edit_scene_name.setText(new_text)
    if new_text in self._all_current_scenes:
      self._view.btn_add_scene.setEnabled(False)
      self._view.line_edit_scene_name.setToolTip(
          "This scene name already exists. Please enter another name."
      )
      self._view.line_edit_scene_name.setStyleSheet(
          """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}""",
      )
    else:
      self._view.btn_add_scene.setEnabled(True)
      self._view.line_edit_scene_name.setToolTip("")
      self._view.line_edit_scene_name.setStyleSheet(
          """QLineEdit {color: #000000; border-color: #DCDBE3;}""",
      )
