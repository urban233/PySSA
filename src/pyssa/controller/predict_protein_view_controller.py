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
"""Module for the predict protein view controller."""
import logging
import os
import subprocess

import pygetwindow
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtCore import Qt
from src.pyssa.controller import add_protein_pair_view_controller, advanced_prediction_configurations_view_controller
from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.internal.data_structures import protein, chain
from src.pyssa.internal.data_structures.data_classes import prediction_protein_info, prediction_configuration
from src.pyssa.internal.thread import tasks
from src.pyssa.internal.thread.async_pyssa import util_async
from src.pyssa.util import tools, constants, prediction_util, enums, exception
from src.pyssa.util import gui_utils
from src.pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class PredictProteinViewController(QtCore.QObject):
  """Class for the PredictProteinViewController."""

  job_input = pyqtSignal(tuple)
  """Singal used to transfer data back to the previous window."""

  def __init__(
      self,
      the_interface_manager: "interface_manager.InterfaceManager",
      the_watcher: "watcher.Watcher",
      the_selected_indexes: list,
      a_prediction_type: str,
  ) -> None:
    """Constructor.

    Args:
        the_interface_manager (interface_manager.InterfaceManager): An instance of the interface_manager.InterfaceManager class.
        the_watcher (watcher.Watcher): An instance of the watcher.Watcher class.
        the_selected_indexes (list): A list of selected indexes.
        a_prediction_type (str): A string representing the prediction type.

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
    if the_selected_indexes is None:
      logger.error("the_selected_indexes is None.")
      raise exception.IllegalArgumentError("the_selected_indexes is None.")
    if a_prediction_type is None:
      logger.error("a_prediction_type is None.")
      raise exception.IllegalArgumentError("a_prediction_type is None.")

    # </editor-fold>

    super().__init__()
    self._interface_manager = the_interface_manager
    self._watcher = the_watcher
    self._view: "predict_protein_view.PredictProteinView" = (
        the_interface_manager.get_predict_protein_view()
    )
    self.prediction_configuration = (
        prediction_configuration.PredictionConfiguration(True, "pdb70")
    )
    self.temporary_protein_objs = []
    self.has_internet_connection = True  # flag to determine from the main_view_controller.py if there is a working internet connection
    try:
      self.restore_ui_defaults()
    except exception.NoInternetConnectionError:
      logger.warning("No internet connection")
      self.has_internet_connection = False
    self._fill_protein_to_predict_table_with_sequences(
        the_selected_indexes, a_prediction_type, the_watcher
    )
    self._connect_all_ui_elements_to_slot_functions()

  # <editor-fold desc="Util methods">
  def open_help(self, a_page_name: str) -> None:
    """Opens the pyssa documentation window if it's not already open.

    Args:
        a_page_name (str): a name of a documentation page to display

    Raises:
        exception.IllegalArgumentError: If `a_page_name` is None.
    """
    # <editor-fold desc="Checks">
    if a_page_name is None:
      logger.error("a_page_name is None.")
      raise exception.IllegalArgumentError("a_page_name is None.")

    # </editor-fold>

    try:
      self._interface_manager.status_bar_manager.show_temporary_message(
          "Opening help center ..."
      )
      if (
          len(
              pygetwindow.getWindowsWithTitle(
                  constants.WINDOW_TITLE_OF_HELP_CENTER
              )
          )
          != 1
      ):
        self._interface_manager.documentation_window = None
      self._active_task = tasks.LegacyTask(
          target=util_async.open_documentation_on_certain_page,
          args=(a_page_name, self._interface_manager.documentation_window),
          post_func=self.__await_open_help,
      )
    except Exception as e:
      logger.error(f"Error while opening help center: {e}")
    else:
      self._active_task.start()

  def __await_open_help(self, return_value: tuple) -> None:
    """Opens the help center and performs necessary actions based on the return value.

    Args:
        return_value (tuple): The return value from opening the help center.
    """
    # <editor-fold desc="Checks">
    if return_value[0] == "":
      self._interface_manager.status_bar_manager.show_error_message(
          "Opening help center failed!"
      )
      return

    # </editor-fold>

    try:
      self._interface_manager.documentation_window = return_value[2]
      if not os.path.exists(constants.HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH):
        tmp_dialog = custom_message_box.CustomMessageBoxOk(
            "The script for bringing the documentation window in front could not be found!",
            "Documentation",
            custom_message_box.CustomMessageBoxIcons.ERROR.value,
        )
        tmp_dialog.exec_()
      else:
        self._interface_manager.documentation_window.restore()
        subprocess.run([constants.HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH])
        self._interface_manager.status_bar_manager.show_temporary_message(
            "Opening help center finished."
        )
    except Exception as e:
      logger.error(f"Error while opening help center: {e}")
      self._interface_manager.status_bar_manager.show_error_message(
          "Opening help center failed!"
      )

  def _open_help_for_dialog(self) -> None:
    """Opens help dialog."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Help' button was clicked."
    )
    self.open_help("help/protein_structure_prediction/colabfold/")

  def _show_prediction_configuration(self) -> None:
    """Opens the prediction configuration dialog window."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Edit' advanced configurations button was clicked.",
    )
    self._external_controller = advanced_prediction_configurations_view_controller.AdvancedPredictionConfigurationsViewController(
        self._interface_manager,
        self.prediction_configuration,
    )
    self._external_controller.user_input.connect(
        self._post_show_prediction_configuration
    )
    self._interface_manager.get_advanced_prediction_configurations_view().show()

  def _post_show_prediction_configuration(self, user_input: tuple) -> None:
    """Updates the prediction configuration with the new values.

    Args:
        user_input (tuple): A tuple containing two elements. The first element is ignored. The second element is a temporary configuration object.

    Raises:
        exception.IllegalArgumentError: If `user_input` is None.
    """
    # <editor-fold desc="Checks">
    if user_input is None:
      logger.error("user_input is None.")
      raise exception.IllegalArgumentError("user_input is None.")

    # </editor-fold>

    _, tmp_config = user_input
    self.prediction_configuration.amber_force_field = (
        tmp_config.amber_force_field
    )
    self.prediction_configuration.templates = tmp_config.templates

  def restore_ui_defaults(self) -> None:
    """Restores the default UI."""
    self._view.ui.btn_go_to_analysis_setup.setText("Predict")
    self._view.ui.lbl_go_to_analysis_setup.setText("Protein Structure(s)")
    self._view.ui.checkbox_add_analysis.setChecked(False)
    # checks internet connection
    if not tools.check_internet_connectivity():
      tmp_dialog = custom_message_box.CustomMessageBoxOk(
          "You do not have a working internet connection\nbut that is necessary for this operation!",
          "Internet Connection",
          custom_message_box.CustomMessageBoxIcons.ERROR.value,
      )
      tmp_dialog.exec_()
      raise exception.NoInternetConnectionError(
          "No working internet connection."
      )

    self._view.ui.list_analysis_overview.clear()
    self._view.ui.btn_analysis_remove.hide()
    self._view.ui.table_proteins_to_predict.clear()
    self._view.ui.table_proteins_to_predict.setRowCount(0)
    self._view.ui.table_proteins_to_predict.setHorizontalHeaderItem(
        0,
        QtWidgets.QTableWidgetItem("Chain"),
    )
    self._view.ui.table_proteins_to_predict.setHorizontalHeaderItem(
        1,
        QtWidgets.QTableWidgetItem("Sequence"),
    )
    self._view.ui.table_proteins_to_predict.resizeColumnsToContents()
    gui_elements_to_show = [
        self._view.ui.lbl_proteins_to_predict,
        self._view.ui.table_proteins_to_predict,
        self._view.ui.checkbox_add_analysis,
        self._view.ui.btn_prediction_remove,
        self._view.ui.lbl_advanced_config,
        self._view.ui.btn_edit_advanced_config,
        self._view.ui.btn_go_to_analysis_setup,
        self._view.ui.lbl_go_to_analysis_setup,
    ]
    gui_utils.show_gui_elements(gui_elements_to_show)

    if self._view.ui.tab_widget.currentIndex() == 1:
      self._view.ui.tab_widget.setCurrentIndex(0)
    self._view.ui.tab_widget.setTabEnabled(1, False)
    self._view.ui.tab_widget.setTabEnabled(0, True)
    self._view.ui.table_proteins_to_predict.setEnabled(True)
    self._view.ui.btn_go_to_analysis_setup.setEnabled(True)
    self._view.resize(700, 800)

  def _check_if_prediction_and_analysis_should_be_done(self) -> None:
    """Checks if the 'Add analysis' checkbox was clicked and updates the GUI accordingly."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Add analysis' checkbox was clicked.",
    )
    if self._view.ui.checkbox_add_analysis.isChecked():
      self._view.ui.btn_go_to_analysis_setup.setText("Go")
      self._view.ui.lbl_go_to_analysis_setup.setText("To Analysis Setup")
    else:
      self._view.ui.btn_go_to_analysis_setup.setText("Predict")
      self._view.ui.lbl_go_to_analysis_setup.setText("Protein Structure(s)")

  # </editor-fold>

  def _connect_all_ui_elements_to_slot_functions(self) -> None:
    """Connects all UI elements to their corresponding slot functions in the class."""
    self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)
    self._view.ui.btn_help_2.clicked.connect(self._open_help_for_dialog)
    self._view.ui.checkbox_add_analysis.clicked.connect(
        self._check_if_prediction_and_analysis_should_be_done
    )

    # Prediction tab
    self._view.ui.btn_prediction_remove.clicked.connect(
        self._remove_protein_to_predict
    )
    self._view.ui.btn_edit_advanced_config.clicked.connect(
        self._show_prediction_configuration
    )
    self._view.ui.btn_go_to_analysis_setup.clicked.connect(self._switch_tab)
    self._view.ui.table_proteins_to_predict.itemSelectionChanged.connect(
        self._enable_remove_button_of_prediction_table
    )

    # Analysis tab
    self._view.ui.btn_analysis_add.clicked.connect(self._add_protein_pair)
    self._view.ui.btn_analysis_remove.clicked.connect(
        self._remove_analysis_run_from_list
    )
    self._view.ui.btn_analysis_back_4.clicked.connect(self._switch_tab)
    self._view.ui.btn_start_prediction_analysis.clicked.connect(
        self._start_prediction_analysis
    )

  # <editor-fold desc="Prediction + analysis">

  # <editor-fold desc="Prediction section">
  def _fill_protein_to_predict_table_with_sequences(
      self,
      tmp_selected_indices: list,
      a_prediction_type: str,
      the_watcher: "watcher.Watcher",
  ) -> None:
    """Fills the prediction table with the selected sequences.

    Args:
        tmp_selected_indices (list): List of QModelIndex objects representing the selected indices in a table.
        a_prediction_type (str): String representing the type of prediction ('monomer' or 'multimer').
        the_watcher (watcher.Watcher): An instance of the `watcher.Watcher` class.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None or if `a_prediction_type` is an empty string.
        ValueError: If prediction type is unknown.
    """
    # <editor-fold desc="Checks">
    if tmp_selected_indices is None:
      logger.error("tmp_selected_indices is None.")
      raise exception.IllegalArgumentError("tmp_selected_indices is None.")
    if a_prediction_type is None or a_prediction_type == "":
      logger.error("a_prediction_type is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_prediction_type is either None or an empty string."
      )
    if the_watcher is None:
      logger.error("the_watcher is None.")
      raise exception.IllegalArgumentError("the_watcher is None.")

    # </editor-fold>

    tmp_sequences_to_predict_monomer: list = []
    tmp_sequences_to_predict_multimer: list = []
    for tmp_model_index in tmp_selected_indices:
      if (
          tmp_model_index.data(enums.ModelEnum.TYPE_ROLE)
          == enums.ModelTypeEnum.MONOMER_SEQ
      ):
        tmp_sequences_to_predict_monomer.append(
            tmp_model_index.data(enums.ModelEnum.OBJECT_ROLE)
        )
      elif (
          tmp_model_index.data(enums.ModelEnum.TYPE_ROLE)
          == enums.ModelTypeEnum.MULTIMER_SEQ
      ):
        tmp_sequences_to_predict_multimer.append(
            tmp_model_index.data(enums.ModelEnum.OBJECT_ROLE)
        )
    self._view.ui.table_proteins_to_predict.setColumnCount(2)

    if a_prediction_type == "monomer":
      tmp_sequences_to_predict = tmp_sequences_to_predict_monomer
    elif a_prediction_type == "multimer":
      tmp_sequences_to_predict = tmp_sequences_to_predict_multimer
    else:
      raise ValueError(f"Unknown prediction type: {a_prediction_type}")

    tmp_row_no = 0
    self.temporary_protein_objs.clear()
    for tmp_seq_record in tmp_sequences_to_predict:
      if self._interface_manager.get_current_project().is_sequence_as_protein_in_project(
          tmp_seq_record.name
      ):
        # Continues if a protein with the name of the given sequence already exists
        continue
      if the_watcher.is_protein_name_on_blacklist(tmp_seq_record.name):
        continue
      tmp_seqs = tmp_seq_record.seq.split(",")
      tmp_chain_no = 0
      for tmp_seq in tmp_seqs:
        self._view.ui.table_proteins_to_predict.insertRow(tmp_row_no)
        tmp_seq_name_item = QtWidgets.QTableWidgetItem(tmp_seq_record.name)
        tmp_chain_letter_item = QtWidgets.QTableWidgetItem(
            constants.chain_dict.get(tmp_chain_no)
        )
        tmp_seq_item = QtWidgets.QTableWidgetItem(tmp_seq)
        tmp_seq_name_item.setFlags(
            tmp_seq_name_item.flags() & ~Qt.ItemIsEditable
        )
        tmp_chain_letter_item.setFlags(
            tmp_chain_letter_item.flags() & ~Qt.ItemIsEditable
        )
        tmp_seq_item.setFlags(tmp_seq_item.flags() & ~Qt.ItemIsEditable)

        self._view.ui.table_proteins_to_predict.setVerticalHeaderItem(
            tmp_row_no, tmp_seq_name_item
        )
        self._view.ui.table_proteins_to_predict.setItem(
            tmp_row_no, 0, tmp_chain_letter_item
        )
        self._view.ui.table_proteins_to_predict.setItem(
            tmp_row_no, 1, tmp_seq_item
        )
        tmp_chain_no += 1
        tmp_row_no += 1
      tmp_protein = protein.Protein(tmp_seq_record.name)
      tmp_chains = []
      for tmp_chain_number in range(tmp_chain_no):
        tmp_chain = chain.Chain(
            constants.chain_dict.get(tmp_chain_number),
            tmp_seqs[tmp_chain_number],
            "protein_chain",
        )
        tmp_chains.append(tmp_chain)
      tmp_protein.chains = tmp_chains
      self.temporary_protein_objs.append(tmp_protein)

    self._view.ui.table_proteins_to_predict.resizeColumnsToContents()
    self._check_if_proteins_to_predict_table_is_empty()

  def _check_if_proteins_to_predict_table_is_empty(self) -> None:
    """Checks if the list of proteins to predict is empty."""
    if self._view.ui.table_proteins_to_predict.rowCount() == 0:
      gui_elements_to_show = [
          self._view.ui.lbl_proteins_to_predict,
          self._view.ui.table_proteins_to_predict,
      ]
      gui_elements_to_hide = [
          self._view.ui.btn_prediction_remove,
          self._view.ui.lbl_advanced_config,
          self._view.ui.btn_edit_advanced_config,
          self._view.ui.btn_go_to_analysis_setup,
          self._view.ui.lbl_go_to_analysis_setup,
          self._view.ui.checkbox_add_analysis,
      ]
      gui_utils.show_gui_elements(gui_elements_to_show)
      gui_utils.hide_gui_elements(gui_elements_to_hide)
      self._view.ui.btn_go_to_analysis_setup.setEnabled(False)
      self._view.ui.btn_prediction_remove.setEnabled(False)
    else:
      self._view.ui.btn_start_prediction_analysis.setEnabled(True)
      gui_elements_to_show = [
          self._view.ui.lbl_proteins_to_predict,
          self._view.ui.table_proteins_to_predict,
          self._view.ui.btn_prediction_remove,
          self._view.ui.lbl_advanced_config,
          self._view.ui.btn_edit_advanced_config,
          self._view.ui.btn_go_to_analysis_setup,
          self._view.ui.lbl_go_to_analysis_setup,
          self._view.ui.checkbox_add_analysis,
      ]
      gui_elements_to_hide = []
      gui_utils.show_gui_elements(gui_elements_to_show)
      gui_utils.hide_gui_elements(gui_elements_to_hide)
      self._view.ui.btn_go_to_analysis_setup.setEnabled(True)
      self._view.ui.btn_prediction_remove.setEnabled(False)

  def _remove_protein_to_predict(self) -> None:
    """Removes the selected protein from the list of proteins to predict."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Remove' button on the 'Prediction Tab' was clicked.",
    )
    tmp_table = self._view.ui.table_proteins_to_predict
    for tmp_protein in self.temporary_protein_objs:
      if (
          tmp_protein.get_molecule_object()
          == tmp_table.verticalHeaderItem(tmp_table.currentRow()).text()
          and len(tmp_protein.chains) > 1
      ):
        tmp_protein.chains.pop(
            constants.chain_dict_reverse.get(
                tmp_table.item(tmp_table.currentRow(), 0).text()
            )
        )
        i = 0
        for tmp_chain in tmp_protein.chains:
          tmp_chain.chain_letter = constants.chain_dict.get(i)
          i += 1
      elif (
          tmp_protein.get_molecule_object()
          == tmp_table.verticalHeaderItem(tmp_table.currentRow()).text()
          and len(tmp_protein.chains) == 1
      ):
        self.temporary_protein_objs.remove(tmp_protein)

    if tmp_table.rowCount() == 1:
      tmp_table.removeRow(0)
    else:
      tmp_table.removeRow(
          tmp_table.currentRow(),
      )
      prot_name = tmp_table.verticalHeaderItem(
          tmp_table.currentRow(),
      ).text()
      for i in range(tmp_table.rowCount()):
        if tmp_table.verticalHeaderItem(i).text() == prot_name:
          tmp_chain_letter_item = QtWidgets.QTableWidgetItem(
              constants.chain_dict.get(i)
          )
          tmp_chain_letter_item.setFlags(
              tmp_chain_letter_item.flags() & ~Qt.ItemIsEditable
          )
          tmp_table.setItem(
              i,
              0,
              tmp_chain_letter_item,
          )
    self._check_if_proteins_to_predict_table_is_empty()
    self._view.ui.btn_prediction_remove.setEnabled(False)

  def _enable_remove_button_of_prediction_table(self) -> None:
    """Enables the remove button of the list of proteins to predict."""
    self._view.ui.btn_prediction_remove.setEnabled(True)

  # </editor-fold>

  def _switch_tab(self) -> None:
    """Switches the tabs from prediction to analysis and vice versa."""
    if self._view.ui.btn_go_to_analysis_setup.text() == "Predict":
      logger.log(
          log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Predict' button was clicked."
      )
      tmp_prediction_runs: list[
          prediction_protein_info.PredictionProteinInfo
      ] = prediction_util.get_prediction_name_and_seq_from_table(
          self._view.ui.table_proteins_to_predict
      )
      self._view.close()
      self.job_input.emit(
          (
              "job_input",
              tmp_prediction_runs,
              self.prediction_configuration,
              False,
          )
      )
    else:
      logger.log(
          log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Go' button was clicked."
      )
      if self._view.ui.tab_widget.currentIndex() == 0:
        # goes from prediction to analysis
        self._view.ui.tab_widget.setCurrentIndex(1)
        gui_elements_to_show = [
            self._view.ui.lbl_analysis_overview,
            self._view.ui.list_analysis_overview,
            self._view.ui.btn_analysis_add,
            self._view.ui.btn_analysis_back_4,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_analysis_remove,
            # self._view.ui.lbl_analysis_protein_struct_1,
            # self._view.ui.lbl_analysis_protein_struct_2,
            # self._view.ui.lbl_analysis_batch_vs_3,
            # self._view.ui.lbl_analysis_protein_1_chains,
            # self._view.ui.list_analysis_protein_1_chains,
            # self._view.ui.btn_analysis_back_2,
            # self._view.ui.btn_analysis_next_2,
            # self._view.ui.box_analysis_protein_struct_1,
            # self._view.ui.box_analysis_protein_struct_2,
            # self._view.ui.btn_analysis_back,
            # self._view.ui.btn_analysis_next,
            # self._view.ui.lbl_analysis_protein_2_chains,
            # self._view.ui.list_analysis_protein_2_chains,
            # self._view.ui.btn_analysis_back_3,
            # self._view.ui.btn_analysis_next_3,
            self._view.ui.btn_start_prediction_analysis,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.tab_widget.setTabEnabled(1, True)
        self._view.ui.tab_widget.setTabEnabled(0, False)
        if self._view.ui.list_analysis_overview.count() > 0:
          self._view.ui.btn_analysis_remove.show()
          self._view.ui.btn_analysis_remove.setEnabled(False)
          self._view.ui.btn_start_prediction_analysis.show()
      else:
        # goes from analysis to prediction
        self._view.ui.tab_widget.setCurrentIndex(0)
        if self._view.ui.list_analysis_overview.count() > 0:
          gui_elements_to_show = [
              self._view.ui.lbl_proteins_to_predict,
              self._view.ui.table_proteins_to_predict,
              self._view.ui.btn_go_to_analysis_setup,
              self._view.ui.lbl_go_to_analysis_setup,
          ]
          gui_elements_to_hide = [
              self._view.ui.btn_prediction_remove,
              self._view.ui.lbl_advanced_config,
              self._view.ui.btn_edit_advanced_config,
          ]
          gui_utils.show_gui_elements(gui_elements_to_show)
          gui_utils.hide_gui_elements(gui_elements_to_hide)
        else:
          gui_elements_to_show = [
              self._view.ui.lbl_proteins_to_predict,
              self._view.ui.table_proteins_to_predict,
              self._view.ui.btn_prediction_remove,
              self._view.ui.lbl_advanced_config,
              self._view.ui.btn_edit_advanced_config,
              self._view.ui.btn_go_to_analysis_setup,
              self._view.ui.lbl_go_to_analysis_setup,
          ]
          gui_elements_to_hide = []
          gui_utils.show_gui_elements(gui_elements_to_show)
          gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.tab_widget.setTabEnabled(0, True)
        self._view.ui.tab_widget.setTabEnabled(1, False)

  # <editor-fold desc="Analysis section">
  def _get_all_current_analysis_runs(self) -> list[str]:
    """Retrieves a list of all current analysis runs.

    Returns:
        A list of strings representing the current analysis runs.
    """
    tmp_analysis_runs = []
    for tmp_row in range(self._view.ui.list_analysis_overview.count()):
      tmp_analysis_runs.append(
          self._view.ui.list_analysis_overview.item(tmp_row).text()
      )
    return tmp_analysis_runs

  def _get_all_current_protein_pair_names(self) -> list[str]:
    """Retrieves the names of all current protein pairs.

    Returns:
        A list of strings representing the names of all current protein pairs.
    """
    tmp_protein_pair_names = []
    for (
        tmp_protein_pair
    ) in self._interface_manager.get_current_project().protein_pairs:
      tmp_protein_pair_names.append(tmp_protein_pair.name)
    return tmp_protein_pair_names

  def _add_protein_pair(self) -> None:
    """Instantiates the AddProteinPairViewController class and shows the 'AddProteinPair' view."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Add' button was clicked."
    )
    self._external_controller = (
        add_protein_pair_view_controller.AddProteinPairViewController(
            self._interface_manager,
            self._watcher,
            self._get_all_current_analysis_runs(),
            self._get_all_current_protein_pair_names(),
            a_list_of_extra_proteins=self.temporary_protein_objs,
        )
    )

    self._external_controller.user_input.connect(self._post_add_protein_pair)
    self._interface_manager.get_add_protein_pair_view().show()

  def _post_add_protein_pair(self, return_value: tuple) -> None:
    """Adds the protein pair to the list widget and updates the UI.

    Args:
        return_value (tuple): The return value from the method.

    Raises:
        exception.IllegalArgumentError: If `return_value` is None.
    """
    # <editor-fold desc="Checks">
    if return_value is None:
      logger.error("return_value is None.")
      raise exception.IllegalArgumentError("return_value is None.")

    # </editor-fold>

    tmp_item, _ = return_value
    self._view.ui.list_analysis_overview.addItem(tmp_item)
    self._view.ui.btn_analysis_remove.show()
    self._view.ui.btn_analysis_remove.setEnabled(False)
    self._view.ui.btn_start_prediction_analysis.show()
    self._view.ui.btn_start_prediction_analysis.setEnabled(True)

  def _enable_remove_button_for_analysis_list(self) -> None:
    """Enables the remove button."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "An analysis run from the list was clicked.",
    )
    self._view.ui.btn_analysis_remove.setEnabled(True)

  def _remove_analysis_run_from_list(self) -> None:
    """Removes the selected protein pair from the list of protein pairs to analyze."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Remove' button on the 'Analysis Tab' was clicked.",
    )
    self._view.ui.list_analysis_overview.takeItem(
        self._view.ui.list_analysis_overview.currentRow(),
    )
    gui_elements_to_show = [
        self._view.ui.lbl_analysis_overview,
        self._view.ui.list_analysis_overview,
        self._view.ui.btn_analysis_add,
        self._view.ui.btn_analysis_back_4,
    ]
    gui_elements_to_hide = [
        self._view.ui.btn_analysis_remove,
        self._view.ui.lbl_analysis_protein_struct_1,
        self._view.ui.lbl_analysis_protein_struct_2,
        self._view.ui.lbl_analysis_batch_vs_3,
        self._view.ui.lbl_analysis_protein_1_chains,
        self._view.ui.list_analysis_protein_1_chains,
        self._view.ui.btn_analysis_back_2,
        self._view.ui.btn_analysis_next_2,
        self._view.ui.box_analysis_protein_struct_1,
        self._view.ui.box_analysis_protein_struct_2,
        self._view.ui.lbl_analysis_protein_2_chains,
        self._view.ui.list_analysis_protein_2_chains,
        self._view.ui.btn_analysis_back_3,
        self._view.ui.btn_analysis_next_3,
        self._view.ui.btn_start_prediction_analysis,
    ]
    gui_utils.show_gui_elements(gui_elements_to_show)
    gui_utils.hide_gui_elements(gui_elements_to_hide)
    if self._view.ui.list_analysis_overview.count() > 0:
      self._view.ui.btn_analysis_remove.show()
      self._view.ui.btn_analysis_remove.setEnabled(False)
      self._view.ui.btn_start_prediction_analysis.show()

  def _start_prediction_analysis(self) -> None:
    """Starts the process batch based on selected items in the prediction and distance analysis overview."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Start' button was clicked."
    )
    tmp_prediction_runs: list[prediction_protein_info.PredictionProteinInfo] = (
        prediction_util.get_prediction_name_and_seq_from_table(
            self._view.ui.table_proteins_to_predict
        )
    )
    self._view.close()
    self.job_input.emit(
        ("job_input", tmp_prediction_runs, self.prediction_configuration, True),
    )

  # </editor-fold>
  # </editor-fold>
