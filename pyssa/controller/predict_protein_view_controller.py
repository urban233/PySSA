import os
import subprocess

from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal

from pyssa.controller import interface_manager
from pyssa.gui.ui.custom_dialogs import custom_message_box
from pyssa.gui.ui.dialogs import dialog_advanced_prediction_configurations
from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.gui.ui.styles import styles
from pyssa.gui.ui.views import distance_analysis_view, predict_monomer_view, predict_multimer_view, predict_protein_view
from pyssa.internal.data_structures.data_classes import prediction_protein_info, prediction_configuration
from pyssa.internal.thread import tasks
from pyssa.internal.thread.async_pyssa import util_async
from pyssa.io_pyssa import safeguard
from pyssa.presenter import main_presenter_async
from pyssa.util import gui_utils, tools, constants, exit_codes, prediction_util


class PredictProteinViewController(QtCore.QObject):
    job_input = pyqtSignal(tuple)

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager"):
        super().__init__()
        self._interface_manager = the_interface_manager
        self._view: "predict_protein_view.PredictProteinView" = the_interface_manager.get_predict_protein_view()
        self.prediction_configuration = prediction_configuration.PredictionConfiguration(True, "pdb70")
        self.restore_ui_defaults()
        self._connect_all_ui_elements_to_slot_functions()

    # <editor-fold desc="Util methods">
    def open_help(self, a_page_name: str):
        """Opens the pyssa documentation window if it's not already open.

        Args:
            a_page_name (str): a name of a documentation page to display
        """
        self._interface_manager.update_status_bar("Opening help center ...")
        self._active_task = tasks.Task(
            target=util_async.open_documentation_on_certain_page,
            args=(a_page_name, self._interface_manager.documentation_window),
            post_func=self.__await_open_help,
        )
        self._active_task.start()

    def __await_open_help(self, return_value):
        self._interface_manager.documentation_window = return_value[2]
        if not os.path.exists(constants.HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH):
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "The script for bringing the documentation window in front could not be found!", "Documentation",
                custom_message_box.CustomMessageBoxIcons.ERROR.value
            )
            tmp_dialog.exec_()
        else:
            self._interface_manager.documentation_window.restore()
            subprocess.run([constants.HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH])
            self._interface_manager.update_status_bar("Opening help center finished.")

    def _open_help_for_dialog(self):
        self.open_help("help/protein_structure_prediction/colabfold_multimer/")

    def _show_prediction_configuration(self) -> None:
        """Opens the prediction configuration dialog window."""
        config = dialog_advanced_prediction_configurations.DialogAdvancedPredictionConfigurations(
            self.prediction_configuration,
        )
        config.exec_()
        self.prediction_configuration.amber_force_field = config.prediction_config.amber_force_field
        self.prediction_configuration.templates = config.prediction_config.templates

    def restore_ui_defaults(self) -> None:
        """Displays the multimer prediction + analysis page."""
        self._view.ui.btn_go_to_analysis_setup.setText("Predict")
        self._view.ui.lbl_go_to_analysis_setup.setText("Protein Structure(s)")
        # checks internet connection
        if not tools.check_internet_connectivity():
            gui_utils.no_internet_dialog()
            return

        self._view.ui.list_analysis_overview.clear()
        self._view.ui.btn_analysis_remove.hide()
        self._view.ui.table_proteins_to_predict.clear()
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
        ]
        gui_elements_to_hide = [
            self._view.ui.checkbox_add_analysis,
            self._view.ui.btn_prediction_remove,
            self._view.ui.lbl_advanced_config,
            self._view.ui.btn_edit_advanced_config,
            self._view.ui.btn_go_to_analysis_setup,
            self._view.ui.lbl_go_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self._view.ui.tab_widget.currentIndex() == 1:
            self._view.ui.tab_widget.setCurrentIndex(0)
        self._view.ui.tab_widget.setTabEnabled(1, False)
        self._view.ui.tab_widget.setTabEnabled(0, True)
        self._view.ui.table_proteins_to_predict.setEnabled(True)

    def _check_if_prediction_and_analysis_should_be_done(self):
        if self._view.ui.checkbox_add_analysis.isChecked():
            self._view.ui.btn_go_to_analysis_setup.setText("Go")
            self._view.ui.lbl_go_to_analysis_setup.setText("To Analysis Setup")
        else:
            self._view.ui.btn_go_to_analysis_setup.setText("Predict")
            self._view.ui.lbl_go_to_analysis_setup.setText("Protein Structure(s)")
    # </editor-fold>

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)
        self._view.ui.btn_help_2.clicked.connect(self._open_help_for_dialog)
        self._view.ui.checkbox_add_analysis.clicked.connect(self._check_if_prediction_and_analysis_should_be_done)

        # Prediction tab
        self._view.ui.btn_prediction_remove.clicked.connect(self._remove_protein_to_predict)
        self._view.ui.btn_edit_advanced_config.clicked.connect(self._show_prediction_configuration)
        self._view.ui.btn_go_to_analysis_setup.clicked.connect(self._switch_tab)
        self._view.ui.table_proteins_to_predict.itemSelectionChanged.connect(self._enable_remove_button_of_prediction_table)

        # Analysis tab
        self._view.ui.btn_analysis_add.clicked.connect(self._enable_selection_of_proteins_for_analysis)
        self._view.ui.btn_analysis_remove.clicked.connect(self._remove_analysis_run_from_list)
        self._view.ui.btn_analysis_back.clicked.connect(self._disable_selection_of_proteins_for_analysis)
        self._view.ui.btn_analysis_next.clicked.connect(self._show_chains_of_protein_structure_1)
        self._view.ui.btn_analysis_back_2.clicked.connect(self._hide_chains_of_protein_structure_1)
        self._view.ui.btn_analysis_next_2.clicked.connect(self._show_chains_of_protein_structure_2)
        self._view.ui.btn_analysis_back_3.clicked.connect(self._hide_chains_of_protein_structure_2)
        self._view.ui.btn_analysis_next_3.clicked.connect(self._add_analysis_run_to_list_widget)
        self._view.ui.box_analysis_protein_struct_1.currentIndexChanged.connect(
            self._check_if_protein_structs_are_filled,
        )
        self._view.ui.box_analysis_protein_struct_2.currentIndexChanged.connect(
            self._check_if_protein_structs_are_filled,
        )
        self._view.ui.list_analysis_protein_1_chains.itemSelectionChanged.connect(
            self._count_selected_chains_for_protein_struct_1,
        )
        self._view.ui.list_analysis_protein_2_chains.itemSelectionChanged.connect(
            self._check_if_same_number_of_chains_selected,
        )
        self._view.ui.btn_analysis_back_4.clicked.connect(self._switch_tab)
        self._view.ui.list_analysis_overview.clicked.connect(
            self._enable_remove_button_for_analysis_list,
        )
        self._view.ui.btn_start_prediction_analysis.clicked.connect(self._start_prediction_analysis)

    # <editor-fold desc="Prediction + analysis">

    # <editor-fold desc="Prediction section">
    def _check_if_proteins_to_predict_table_is_empty(self) -> None:
        """Checks if the list of proteins to predict is empty."""
        if self._view.ui.table_proteins_to_predict.rowCount() == 0:
            gui_elements_to_show = [
                self._view.ui.lbl_proteins_to_predict,
                self._view.ui.table_proteins_to_predict,
            ]
            gui_elements_to_hide = [
                self._view.ui.btn_prediction_remove,
                self._view.ui.lbl_pred_analysis_multi_prot_name,
                self._view.ui.txt_pred_analysis_multi_prot_name,
                self._view.ui.lbl_pred_analysis_multi_prot_name_status,
                self._view.ui.btn_pred_analysis_multi_back,
                self._view.ui.btn_pred_analysis_multi_next,
                self._view.ui.lbl_pred_analysis_multi_prot_seq,
                self._view.ui.txt_pred_analysis_multi_prot_seq,
                self._view.ui.lbl_pred_analysis_multi_prot_seq_status,
                self._view.ui.lbl_pred_multi_prot_seq_add_2,
                self._view.ui.btn_pred_analysis_multi_prot_seq_add,
                self._view.ui.lbl_pred_analysis_multi_prot_seq_overview,
                self._view.ui.list_pred_analysis_multi_prot_seq_overview,
                self._view.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
                self._view.ui.lbl_pred_analysis_multi_prot_to_predict_2,
                self._view.ui.btn_pred_analysis_multi_back_2,
                self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
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
            ]
            gui_elements_to_hide = [
                self._view.ui.lbl_pred_analysis_multi_prot_name,
                self._view.ui.txt_pred_analysis_multi_prot_name,
                self._view.ui.lbl_pred_analysis_multi_prot_name_status,
                self._view.ui.btn_pred_analysis_multi_back,
                self._view.ui.btn_pred_analysis_multi_next,
                self._view.ui.lbl_pred_analysis_multi_prot_seq,
                self._view.ui.txt_pred_analysis_multi_prot_seq,
                self._view.ui.lbl_pred_analysis_multi_prot_seq_status,
                self._view.ui.lbl_pred_multi_prot_seq_add_2,
                self._view.ui.btn_pred_analysis_multi_prot_seq_add,
                self._view.ui.lbl_pred_analysis_multi_prot_seq_overview,
                self._view.ui.list_pred_analysis_multi_prot_seq_overview,
                self._view.ui.btn_pred_analysis_multi_prot_seq_overview_remove,
                self._view.ui.lbl_pred_analysis_multi_prot_to_predict_2,
                self._view.ui.btn_pred_analysis_multi_back_2,
                self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.btn_prediction_remove.setEnabled(False)

    def _remove_protein_to_predict(self) -> None:
        """Removes the selected protein from the list of proteins to predict."""
        if self._view.ui.table_proteins_to_predict.rowCount() == 1:
            self._view.ui.table_proteins_to_predict.removeRow(0)
        else:
            self._view.ui.table_proteins_to_predict.removeRow(
                self._view.ui.table_proteins_to_predict.currentRow(),
            )
            prot_name = self._view.ui.table_proteins_to_predict.verticalHeaderItem(
                self._view.ui.table_proteins_to_predict.currentRow(),
            ).text()
            for i in range(self._view.ui.table_proteins_to_predict.rowCount()):
                if self._view.ui.table_proteins_to_predict.verticalHeaderItem(
                        i).text() == prot_name:
                    self._view.ui.table_proteins_to_predict.setItem(
                        i,
                        0,
                        QtWidgets.QTableWidgetItem(constants.chain_dict.get(i)),
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
            tmp_prediction_runs: list[
                prediction_protein_info.PredictionProteinInfo
            ] = prediction_util.get_prediction_name_and_seq_from_table(
                self._view.ui.table_proteins_to_predict)
            self._view.close()
            self.job_input.emit(("job_input", tmp_prediction_runs, self.prediction_configuration, False))
        else:
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
                    self._view.ui.lbl_analysis_protein_struct_1,
                    self._view.ui.lbl_analysis_protein_struct_2,
                    self._view.ui.lbl_analysis_batch_vs_3,
                    self._view.ui.lbl_analysis_protein_1_chains,
                    self._view.ui.list_analysis_protein_1_chains,
                    self._view.ui.btn_analysis_back_2,
                    self._view.ui.btn_analysis_next_2,
                    self._view.ui.box_analysis_protein_struct_1,
                    self._view.ui.box_analysis_protein_struct_2,
                    self._view.ui.btn_analysis_back,
                    self._view.ui.btn_analysis_next,
                    self._view.ui.lbl_analysis_protein_2_chains,
                    self._view.ui.list_analysis_protein_2_chains,
                    self._view.ui.btn_analysis_back_3,
                    self._view.ui.btn_analysis_next_3,
                    
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
                    self._view.ui.cb_pred_analysis_multi_images.hide()
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
                    gui_elements_to_hide = [

                    ]
                    gui_utils.show_gui_elements(gui_elements_to_show)
                    gui_utils.hide_gui_elements(gui_elements_to_hide)
                self._view.ui.tab_widget.setTabEnabled(0, True)
                self._view.ui.tab_widget.setTabEnabled(1, False)

    # <editor-fold desc="Analysis section">
    def _enable_selection_of_proteins_for_analysis(self) -> None:
        """Shows the gui elements to choose the two proteins."""
        gui_elements_to_show = [
            self._view.ui.lbl_analysis_overview,
            self._view.ui.list_analysis_overview,
            self._view.ui.lbl_analysis_protein_struct_1,
            self._view.ui.box_analysis_protein_struct_1,
            self._view.ui.lbl_analysis_batch_vs_3,
            self._view.ui.lbl_analysis_protein_struct_2,
            self._view.ui.box_analysis_protein_struct_2,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_analysis_remove,
            self._view.ui.btn_analysis_add,
            self._view.ui.lbl_analysis_protein_1_chains,
            self._view.ui.list_analysis_protein_1_chains,
            self._view.ui.btn_analysis_back_2,
            self._view.ui.btn_analysis_next_2,
            self._view.ui.lbl_analysis_protein_2_chains,
            self._view.ui.list_analysis_protein_2_chains,
            self._view.ui.btn_analysis_back_3,
            self._view.ui.btn_analysis_next_3,
            self._view.ui.btn_start_prediction_analysis,
            self._view.ui.btn_analysis_back_4,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.lbl_analysis_protein_struct_1.clear()
        self._view.ui.lbl_analysis_protein_struct_2.clear()
        self._view.ui.lbl_analysis_protein_struct_1.setText("Protein structure 1")
        self._view.ui.lbl_analysis_protein_struct_2.setText("Protein structure 2")
        self._fill_protein_struct_combo_boxes()
        if self._view.ui.list_analysis_overview.count() > 0:
            try:
                self._view.ui.list_analysis_overview.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection in struction analysis overview.")

    def _disable_selection_of_proteins_for_analysis(self) -> None:
        """Hides the gui elements to choose the two proteins."""
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
            self._view.ui.cb_pred_analysis_multi_images.hide()

    def _show_chains_of_protein_structure_1(self) -> None:
        """Shows the gui elements to select the chains in protein 1."""
        self._view.wait_spinner.start()
        tmp_proteins_to_predict: list[str] = []
        for i in range(self._view.ui.table_proteins_to_predict.rowCount()):
            tmp_proteins_to_predict.append(
                self._view.ui.table_proteins_to_predict.verticalHeaderItem(i).text(),
            )
        self._active_task = tasks.Task(
            target=main_presenter_async.check_chains_for_subsequent_analysis,
            args=(
                self._view.ui.box_analysis_protein_struct_1.currentText(),
                self._view.ui.box_analysis_protein_struct_2.currentText(),
                self._interface_manager.get_current_project(),
                tmp_proteins_to_predict,
            ),
            post_func=self.__await_show_chains_of_protein_structure_1,
        )
        self._active_task.start()

    def __await_show_chains_of_protein_structure_1(self, result: tuple) -> None:
        # fixme: something is not quite right with the detection of the chains!
        _, tmp_analysis_name = result
        if tmp_analysis_name != "":
            gui_elements_to_show = [
                self._view.ui.btn_analysis_remove,
                self._view.ui.btn_analysis_add,
                self._view.ui.lbl_analysis_overview,
                self._view.ui.list_analysis_overview,
                self._view.ui.btn_start_prediction_analysis,
                self._view.ui.btn_analysis_back_4,
            ]
            gui_elements_to_hide = [
                self._view.ui.box_analysis_protein_struct_1,
                self._view.ui.box_analysis_protein_struct_2,
                self._view.ui.lbl_analysis_protein_struct_1,
                self._view.ui.lbl_analysis_protein_struct_2,
                self._view.ui.lbl_analysis_batch_vs_3,
                self._view.ui.lbl_analysis_protein_1_chains,
                self._view.ui.list_analysis_protein_1_chains,
                self._view.ui.lbl_analysis_protein_2_chains,
                self._view.ui.list_analysis_protein_2_chains,
                self._view.ui.btn_analysis_back,
                self._view.ui.btn_analysis_next,
                self._view.ui.btn_analysis_back_2,
                self._view.ui.btn_analysis_next_2,
                self._view.ui.btn_analysis_back_3,
                self._view.ui.btn_analysis_next_3,

            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            item = QtWidgets.QListWidgetItem(tmp_analysis_name)
            self._view.ui.list_analysis_overview.addItem(item)
            self._view.ui.btn_analysis_remove.setEnabled(False)
        else:
            gui_elements_to_show = [
                self._view.ui.lbl_analysis_overview,
                self._view.ui.list_analysis_overview,
                self._view.ui.lbl_analysis_protein_struct_1,
                self._view.ui.lbl_analysis_protein_struct_2,
                self._view.ui.lbl_analysis_batch_vs_3,
                self._view.ui.lbl_analysis_protein_1_chains,
                self._view.ui.list_analysis_protein_1_chains,
                self._view.ui.btn_analysis_back_2,
                self._view.ui.btn_analysis_next_2,
            ]
            gui_elements_to_hide = [
                self._view.ui.btn_analysis_remove,
                self._view.ui.btn_analysis_add,
                self._view.ui.box_analysis_protein_struct_1,
                self._view.ui.box_analysis_protein_struct_2,
                self._view.ui.btn_analysis_back,
                self._view.ui.btn_analysis_next,
                self._view.ui.lbl_analysis_protein_2_chains,
                self._view.ui.list_analysis_protein_2_chains,
                self._view.ui.btn_analysis_back_3,
                self._view.ui.btn_analysis_next_3,

                self._view.ui.btn_start_prediction_analysis,
                self._view.ui.btn_analysis_back_4,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.lbl_analysis_protein_struct_1.setText(
                self._view.ui.box_analysis_protein_struct_1.currentText(),
            )
            self._view.ui.lbl_analysis_protein_struct_2.setText(
                self._view.ui.box_analysis_protein_struct_2.currentText(),
            )
            self._view.ui.list_analysis_protein_1_chains.clear()
            self._view.ui.btn_analysis_next_2.setEnabled(False)
            self._view.ui.list_analysis_protein_1_chains.setEnabled(True)

            for i in range(self._view.ui.table_proteins_to_predict.rowCount()):
                if (
                        self._view.ui.table_proteins_to_predict.verticalHeaderItem(i).text()
                        == self._view.ui.box_analysis_protein_struct_1.currentText()
                ):
                    self._view.ui.list_analysis_protein_1_chains.addItem(
                        self._view.ui.table_proteins_to_predict.item(i, 0).text(),
                    )
            if self._view.ui.list_analysis_protein_1_chains.count() == 0:
                tmp_protein = self._interface_manager.get_current_project().search_protein(
                    self._view.ui.box_analysis_protein_struct_1.currentText(),
                )
                for tmp_chain in tmp_protein.chains:
                    if tmp_chain.chain_type == "protein_chain":
                        self._view.ui.list_analysis_protein_1_chains.addItem(tmp_chain.chain_letter)
            if self._view.ui.list_analysis_protein_1_chains.count() == 1:
                self._view.ui.lbl_analysis_protein_1_chains.setText(
                    f"Select chain in protein structure {self._view.ui.lbl_analysis_protein_struct_1.text()}.",
                )
            else:
                self._view.ui.lbl_analysis_protein_1_chains.setText(
                    f"Select chains in protein structure {self._view.ui.lbl_analysis_protein_struct_1.text()}.",
                )
        self._view.wait_spinner.stop()

    def _show_chains_of_protein_structure_2(self) -> None:
        """Shows the gui elements to select the chains in protein 2."""
        gui_elements_to_show = [
            self._view.ui.lbl_analysis_overview,
            self._view.ui.list_analysis_overview,
            self._view.ui.lbl_analysis_protein_struct_1,
            self._view.ui.lbl_analysis_protein_struct_2,
            self._view.ui.lbl_analysis_batch_vs_3,
            self._view.ui.lbl_analysis_protein_1_chains,
            self._view.ui.list_analysis_protein_1_chains,
            self._view.ui.lbl_analysis_protein_2_chains,
            self._view.ui.list_analysis_protein_2_chains,
            self._view.ui.btn_analysis_back_3,
            self._view.ui.btn_analysis_next_3,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_analysis_remove,
            self._view.ui.btn_analysis_add,
            self._view.ui.box_analysis_protein_struct_1,
            self._view.ui.box_analysis_protein_struct_2,
            
            self._view.ui.btn_analysis_back_2,
            self._view.ui.btn_analysis_next_2,
            
            self._view.ui.btn_start_prediction_analysis,
            self._view.ui.btn_analysis_back_4,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.list_analysis_protein_2_chains.clear()
        self._view.ui.list_analysis_protein_1_chains.setEnabled(False)
        self._view.ui.btn_analysis_next_3.setEnabled(False)

        for i in range(self._view.ui.table_proteins_to_predict.rowCount()):
            if (
                    self._view.ui.table_proteins_to_predict.verticalHeaderItem(i).text()
                    == self._view.ui.box_analysis_protein_struct_2.currentText()
            ):
                self._view.ui.list_analysis_protein_2_chains.addItem(
                    self._view.ui.table_proteins_to_predict.item(i, 0).text(),
                )
        if self._view.ui.list_analysis_protein_2_chains.count() == 0:
            tmp_protein = self._interface_manager.get_current_project().search_protein(
                self._view.ui.box_analysis_protein_struct_2.currentText(),
            )
            for tmp_chain in tmp_protein.chains:
                if tmp_chain.chain_type == "protein_chain":
                    self._view.ui.list_analysis_protein_2_chains.addItem(tmp_chain.chain_letter)
        if len(self._view.ui.list_analysis_protein_1_chains.selectedItems()) == 1:
            self._view.ui.lbl_analysis_protein_2_chains.setText(
                f"Select 1 chain in protein structure {self._view.ui.lbl_analysis_protein_struct_2.text()}.",
            )
        else:
            self._view.ui.lbl_analysis_protein_2_chains.setText(
                f"Select {len(self._view.ui.list_analysis_protein_1_chains.selectedItems())} chains in "
                f"protein structure {self._view.ui.lbl_analysis_protein_struct_2.text()}.",
            )

    def _hide_chains_of_protein_structure_1(self) -> None:
        """Hides the gui elements to select the chains in protein 1."""
        gui_elements_to_show = [
            self._view.ui.lbl_analysis_overview,
            self._view.ui.list_analysis_overview,
            self._view.ui.lbl_analysis_protein_struct_1,
            self._view.ui.box_analysis_protein_struct_1,
            self._view.ui.lbl_analysis_batch_vs_3,
            self._view.ui.lbl_analysis_protein_struct_2,
            self._view.ui.box_analysis_protein_struct_2,
            self._view.ui.btn_analysis_back_4,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_analysis_remove,
            self._view.ui.btn_analysis_add,
            self._view.ui.lbl_analysis_protein_1_chains,
            self._view.ui.list_analysis_protein_1_chains,
            self._view.ui.btn_analysis_back_2,
            self._view.ui.btn_analysis_next_2,
            self._view.ui.lbl_analysis_protein_2_chains,
            self._view.ui.list_analysis_protein_2_chains,
            self._view.ui.btn_analysis_back_3,
            self._view.ui.btn_analysis_next_3,
            self._view.ui.btn_start_prediction_analysis,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.lbl_analysis_protein_struct_1.setText("Protein structure 1")
        self._view.ui.lbl_analysis_protein_struct_2.setText("Protein structure 2")

    def _add_analysis_run_to_list_widget(self) -> None:
        """Adds the protein pair to the list of protein pairs to analyze."""
        gui_elements_to_show = [
            self._view.ui.btn_analysis_remove,
            self._view.ui.btn_analysis_add,
            self._view.ui.lbl_analysis_overview,
            self._view.ui.list_analysis_overview,
            self._view.ui.btn_start_prediction_analysis,
            self._view.ui.btn_analysis_back_4,
        ]
        gui_elements_to_hide = [
            self._view.ui.box_analysis_protein_struct_1,
            self._view.ui.box_analysis_protein_struct_2,
            self._view.ui.lbl_analysis_protein_struct_1,
            self._view.ui.lbl_analysis_protein_struct_2,
            self._view.ui.lbl_analysis_batch_vs_3,
            self._view.ui.lbl_analysis_protein_1_chains,
            self._view.ui.list_analysis_protein_1_chains,
            self._view.ui.lbl_analysis_protein_2_chains,
            self._view.ui.list_analysis_protein_2_chains,
            self._view.ui.btn_analysis_back_2,
            self._view.ui.btn_analysis_next_2,
            self._view.ui.btn_analysis_back_3,
            self._view.ui.btn_analysis_next_3,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        prot_1_name = self._view.ui.lbl_analysis_protein_struct_1.text()
        prot_1_chains = []
        for chain in self._view.ui.list_analysis_protein_1_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self._view.ui.lbl_analysis_protein_struct_2.text()
        prot_2_chains = []
        for chain in self._view.ui.list_analysis_protein_2_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        item = QtWidgets.QListWidgetItem(analysis_name)
        self._view.ui.list_analysis_overview.addItem(item)
        self._view.ui.btn_analysis_remove.setEnabled(False)

    def _hide_chains_of_protein_structure_2(self) -> None:
        """Hides the gui elements to select the chains in protein 2."""
        gui_elements_to_show = [
            self._view.ui.lbl_analysis_overview,
            self._view.ui.list_analysis_overview,
            self._view.ui.lbl_analysis_protein_struct_1,
            self._view.ui.lbl_analysis_protein_struct_2,
            self._view.ui.lbl_analysis_batch_vs_3,
            self._view.ui.lbl_analysis_protein_1_chains,
            self._view.ui.list_analysis_protein_1_chains,
            self._view.ui.btn_analysis_back_2,
            self._view.ui.btn_analysis_next_2,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_analysis_remove,
            self._view.ui.btn_analysis_add,
            self._view.ui.box_analysis_protein_struct_1,
            self._view.ui.box_analysis_protein_struct_2,
            
            self._view.ui.btn_analysis_back_3,
            self._view.ui.btn_analysis_next_3,
            
            self._view.ui.btn_start_prediction_analysis,
            self._view.ui.lbl_analysis_protein_2_chains,
            self._view.ui.list_analysis_protein_2_chains,
            self._view.ui.btn_analysis_back_4,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.list_analysis_protein_1_chains.setEnabled(True)

        # tmp_protein = self._interface_manager.get_current_project().search_protein(self._view.ui.box_analysis_protein_struct_2.currentText())
        # for tmp_chain in tmp_protein.chains:
        #     if tmp_chain.chain_type == "protein_chain":
        #         self._view.ui.list_analysis_protein_1_chains.addItem(tmp_chain.chain_letter)

    def _enable_remove_button_for_analysis_list(self) -> None:
        """Enables the remove button."""
        self._view.ui.btn_analysis_remove.setEnabled(True)

    def _fill_protein_struct_combo_boxes(self) -> None:
        """Fills the combo boxes with the protein names."""
        protein_names = []
        for i in range(self._view.ui.table_proteins_to_predict.rowCount()):
            protein_names.append(
                self._view.ui.table_proteins_to_predict.verticalHeaderItem(i).text())
        for tmp_protein in self._interface_manager.get_current_project().proteins:
            protein_names.append(tmp_protein.get_molecule_object())
        protein_names.insert(0, "")
        protein_names = list(set(protein_names))
        self._view.ui.box_analysis_protein_struct_1.clear()
        self._view.ui.box_analysis_protein_struct_2.clear()
        gui_utils.fill_combo_box(self._view.ui.box_analysis_protein_struct_1, protein_names)
        gui_utils.fill_combo_box(self._view.ui.box_analysis_protein_struct_2, protein_names)

    def _remove_analysis_run_from_list(self) -> None:
        """Removes the selected protein pair from the list of protein pairs to analyze."""
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

    def _check_if_same_number_of_chains_selected(self) -> None:
        """Checks if the same number of chains were selected."""
        self._view.ui.btn_analysis_next_3.setEnabled(False)
        if self.no_of_selected_chains == len(self._view.ui.list_analysis_protein_2_chains.selectedItems()):
            self._view.ui.btn_analysis_next_3.setEnabled(True)

        prot_1_name = self._view.ui.lbl_analysis_protein_struct_1.text()
        prot_1_chains = []
        for chain in self._view.ui.list_analysis_protein_1_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self._view.ui.lbl_analysis_protein_struct_2.text()
        prot_2_chains = []
        for chain in self._view.ui.list_analysis_protein_2_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        for tmp_row in range(self._view.ui.list_analysis_overview.count()):
            if analysis_name == self._view.ui.list_analysis_overview.item(tmp_row).text():
                self._view.ui.btn_analysis_next_3.setEnabled(False)
                return

    def _check_if_protein_structs_are_filled(self) -> None:
        """Checks if two proteins were selected."""
        prot_1 = self._view.ui.box_analysis_protein_struct_1.itemText(
            self._view.ui.box_analysis_protein_struct_1.currentIndex(),
        )
        prot_2 = self._view.ui.box_analysis_protein_struct_2.itemText(
            self._view.ui.box_analysis_protein_struct_2.currentIndex(),
        )
        if prot_1 != "" and prot_2 != "":
            self._view.ui.btn_analysis_next.setEnabled(True)
        else:
            self._view.ui.btn_analysis_next.setEnabled(False)

    def _count_selected_chains_for_protein_struct_1(self) -> None:
        """Counts the number of chains in protein 1."""
        self.no_of_selected_chains = len(self._view.ui.list_analysis_protein_1_chains.selectedItems())
        if self.no_of_selected_chains > 0:
            self._view.ui.btn_analysis_next_2.setEnabled(True)
        else:
            self._view.ui.btn_analysis_next_2.setEnabled(False)

    def _start_prediction_analysis(self):
        tmp_prediction_runs: list[
            prediction_protein_info.PredictionProteinInfo
        ] = prediction_util.get_prediction_name_and_seq_from_table(
            self._view.ui.table_proteins_to_predict)
        self._view.close()
        self.job_input.emit(
            ("job_input",
             tmp_prediction_runs,
             self.prediction_configuration,
             True)
        )

    # </editor-fold>
    # </editor-fold>
