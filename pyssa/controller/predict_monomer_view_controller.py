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
from pyssa.gui.ui.views import distance_analysis_view, predict_monomer_view
from pyssa.internal.data_structures.data_classes import prediction_protein_info, prediction_configuration
from pyssa.internal.thread import tasks
from pyssa.internal.thread.async_pyssa import util_async
from pyssa.io_pyssa import safeguard
from pyssa.presenter import main_presenter_async
from pyssa.util import gui_utils, tools, constants, exit_codes, prediction_util


class PredictMonomerViewController(QtCore.QObject):

    job_input = pyqtSignal(tuple)

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager"):
        super().__init__()
        self._interface_manager = the_interface_manager
        self._view: "predict_monomer_view.PredictMonomerView" = the_interface_manager.get_predict_monomer_view()
        self.prediction_configuration = prediction_configuration.PredictionConfiguration(True, "pdb70")
        self.display_monomer_pred_analysis()
        self._connect_all_ui_elements_to_slot_functions()

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
        self.open_help("help/protein_structure_prediction/colabfold_monomer/")

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)
        self._view.ui.btn_help_2.clicked.connect(self._open_help_for_dialog)
        self._view.ui.checkbox_add_analysis.clicked.connect(self.check_if_prediction_and_analysis_should_be_done)
        # <editor-fold desc="Monomer Prediction + Analysis page">
        # <editor-fold desc="Prediction section">
        self._view.ui.btn_pred_analysis_mono_seq_to_predict.clicked.connect(self.mono_pred_analysis_add_seq_to_predict)
        self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove.clicked.connect(
            self.mono_pred_analysis_remove_protein_to_predict,
        )
        self._view.ui.btn_pred_analysis_mono_next.clicked.connect(self.mono_pred_analysis_next)
        self._view.ui.btn_pred_analysis_mono_back.clicked.connect(self.mono_pred_analysis_back)
        self._view.ui.btn_pred_analysis_mono_add_protein.clicked.connect(self.mono_pred_analysis_add_protein)
        self._view.ui.btn_pred_analysis_mono_back_2.clicked.connect(self.mono_pred_analysis_back_2)
        self._view.ui.txt_pred_analysis_mono_prot_name.textChanged.connect(
            self.mono_pred_analysis_validate_protein_name,
        )
        self._view.ui.txt_pred_analysis_mono_seq_name.textChanged.connect(
            self.mono_pred_analysis_validate_protein_sequence,
        )
        self._view.ui.btn_pred_mono_advanced_config_2.clicked.connect(self.show_prediction_configuration)
        self._view.ui.btn_pred_analysis_mono_go_analysis_setup.clicked.connect(self.switch_monomer_pred_analysis_tab)
        self._view.ui.table_pred_analysis_mono_prot_to_predict.itemSelectionChanged.connect(
            self.mono_pred_analysis_prediction_overview_item_clicked,
        )

        # </editor-fold>

        # <editor-fold desc="Analysis section">
        self._view.ui.btn_pred_analysis_mono_add.clicked.connect(self.mono_pred_analysis_structure_analysis_add)
        self._view.ui.btn_pred_analysis_mono_remove.clicked.connect(self.remove_mono_pred_analysis_analysis_run)
        self._view.ui.btn_pred_analysis_mono_back_3.clicked.connect(self.mono_pred_analysis_structure_analysis_back_3)
        self._view.ui.btn_pred_analysis_mono_next_2.clicked.connect(self.mono_pred_analysis_structure_analysis_next_2)
        self._view.ui.btn_pred_analysis_mono_back_4.clicked.connect(self.mono_pred_analysis_structure_analysis_back_4)
        self._view.ui.btn_pred_analysis_mono_next_3.clicked.connect(self.mono_pred_analysis_structure_analysis_next_3)
        self._view.ui.btn_pred_analysis_mono_back_5.clicked.connect(self.mono_pred_analysis_structure_analysis_back_5)
        self._view.ui.btn_pred_analysis_mono_next_4.clicked.connect(self.mono_pred_analysis_structure_analysis_next_4)
        self._view.ui.box_pred_analysis_mono_prot_struct_1.currentIndexChanged.connect(
            self.check_mono_pred_analysis_if_prot_structs_are_filled,
        )
        self._view.ui.box_pred_analysis_mono_prot_struct_2.currentIndexChanged.connect(
            self.check_mono_pred_analysis_if_prot_structs_are_filled,
        )
        self._view.ui.list_pred_analysis_mono_ref_chains.itemSelectionChanged.connect(
            self.count_mono_pred_analysis_selected_chains_for_prot_struct_1,
        )
        self._view.ui.list_pred_analysis_mono_model_chains.itemSelectionChanged.connect(
            self.check_mono_pred_analysis_if_same_no_of_chains_selected,
        )
        self._view.ui.btn_pred_analysis_mono_back_pred_setup.clicked.connect(self.switch_monomer_pred_analysis_tab)
        self._view.ui.btn_pred_analysis_mono_start.clicked.connect(self.start_monomer_prediction_analysis)
        self._view.ui.list_pred_analysis_mono_overview.clicked.connect(
            self.mono_pred_analysis_structure_analysis_overview_clicked,
        )

        # </editor-fold>

        # </editor-fold>

    def mono_pred_analysis_structure_analysis_next_2(self) -> None:
        """Shows the gui elements for the chain selection of protein 1."""
        self._view.wait_spinner.start()
        tmp_proteins_to_predict: list[str] = []
        for i in range(self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount()):
            tmp_proteins_to_predict.append(
                self._view.ui.table_pred_analysis_mono_prot_to_predict.verticalHeaderItem(i).text(),
            )
        self._active_task = tasks.Task(
            target=main_presenter_async.check_chains_for_subsequent_analysis,
            args=(
                self._view.ui.box_pred_analysis_mono_prot_struct_1.currentText(),
                self._view.ui.box_pred_analysis_mono_prot_struct_2.currentText(),
                self._interface_manager.get_current_project(),
                tmp_proteins_to_predict,
            ),
            post_func=self.__await_mono_pred_analysis_structure_analysis_next_2,
        )
        self._active_task.start()

    def __await_mono_pred_analysis_structure_analysis_next_2(self, result: tuple) -> None:
        _, tmp_analysis_name = result
        if tmp_analysis_name != "":
            gui_elements_to_show = [
                self._view.ui.btn_pred_analysis_mono_remove,
                self._view.ui.btn_pred_analysis_mono_add,
                self._view.ui.lbl_pred_analysis_mono_overview,
                self._view.ui.list_pred_analysis_mono_overview,
                
                self._view.ui.btn_pred_analysis_mono_start,
                self._view.ui.btn_pred_analysis_mono_back_pred_setup,
            ]
            gui_elements_to_hide = [
                self._view.ui.box_pred_analysis_mono_prot_struct_1,
                self._view.ui.box_pred_analysis_mono_prot_struct_2,
                self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
                self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
                self._view.ui.lbl_analysis_batch_vs_2,
                self._view.ui.lbl_pred_analysis_mono_ref_chains,
                self._view.ui.list_pred_analysis_mono_ref_chains,
                self._view.ui.lbl_pred_analysis_mono_model_chains,
                self._view.ui.list_pred_analysis_mono_model_chains,
                self._view.ui.btn_pred_analysis_mono_back_3,
                self._view.ui.btn_pred_analysis_mono_next_2,
                self._view.ui.btn_pred_analysis_mono_back_4,
                self._view.ui.btn_pred_analysis_mono_next_3,
                self._view.ui.btn_pred_analysis_mono_back_5,
                self._view.ui.btn_pred_analysis_mono_next_4,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            item = QtWidgets.QListWidgetItem(tmp_analysis_name)
            self._view.ui.list_pred_analysis_mono_overview.addItem(item)
            self._view.ui.btn_pred_analysis_mono_remove.setEnabled(False)
        else:
            gui_elements_to_show = [
                self._view.ui.lbl_pred_analysis_mono_overview,
                self._view.ui.list_pred_analysis_mono_overview,
                self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
                self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
                self._view.ui.lbl_analysis_batch_vs_2,
                self._view.ui.lbl_pred_analysis_mono_ref_chains,
                self._view.ui.list_pred_analysis_mono_ref_chains,
                self._view.ui.btn_pred_analysis_mono_back_4,
                self._view.ui.btn_pred_analysis_mono_next_3,
            ]
            gui_elements_to_hide = [
                self._view.ui.btn_pred_analysis_mono_remove,
                self._view.ui.btn_pred_analysis_mono_add,
                self._view.ui.box_pred_analysis_mono_prot_struct_1,
                self._view.ui.box_pred_analysis_mono_prot_struct_2,
                self._view.ui.btn_pred_analysis_mono_back_3,
                self._view.ui.btn_pred_analysis_mono_next_2,
                self._view.ui.lbl_pred_analysis_mono_model_chains,
                self._view.ui.list_pred_analysis_mono_model_chains,
                self._view.ui.btn_pred_analysis_mono_back_5,
                self._view.ui.btn_pred_analysis_mono_next_4,

                
                self._view.ui.btn_pred_analysis_mono_start,
                self._view.ui.btn_pred_analysis_mono_back_pred_setup,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.lbl_pred_analysis_mono_prot_struct_1.setText(
                self._view.ui.box_pred_analysis_mono_prot_struct_1.currentText(),
            )
            self._view.ui.lbl_pred_analysis_mono_prot_struct_2.setText(
                self._view.ui.box_pred_analysis_mono_prot_struct_2.currentText(),
            )
            self._view.ui.list_pred_analysis_mono_ref_chains.clear()
            self._view.ui.btn_pred_analysis_mono_next_3.setEnabled(False)
            self._view.ui.list_pred_analysis_mono_ref_chains.setEnabled(True)

            for i in range(self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount()):
                if (
                    self._view.ui.table_pred_analysis_mono_prot_to_predict.verticalHeaderItem(i).text()
                    == self._view.ui.box_pred_analysis_mono_prot_struct_1.currentText()
                ):
                    self._view.ui.list_pred_analysis_mono_ref_chains.addItem(
                        self._view.ui.table_pred_analysis_mono_prot_to_predict.item(i, 0).text(),
                    )
            if self._view.ui.list_pred_analysis_mono_ref_chains.count() == 0:
                tmp_protein = self._interface_manager.get_current_project().search_protein(
                    self._view.ui.box_pred_analysis_mono_prot_struct_1.currentText(),
                )
                for tmp_chain in tmp_protein.chains:
                    if tmp_chain.chain_type == "protein_chain":
                        self._view.ui.list_pred_analysis_mono_ref_chains.addItem(tmp_chain.chain_letter)
            if self._view.ui.list_pred_analysis_mono_ref_chains.count() == 1:
                self._view.ui.lbl_pred_analysis_mono_ref_chains.setText(
                    f"Select chain in protein structure {self._view.ui.lbl_pred_analysis_mono_prot_struct_1.text()}.",
                )
            else:
                self._view.ui.lbl_pred_analysis_mono_ref_chains.setText(
                    f"Select chains in protein structure {self._view.ui.lbl_pred_analysis_mono_prot_struct_1.text()}.",
                )
        self._view.wait_spinner.stop()

    def start_monomer_prediction_analysis(self):
        tmp_prediction_runs: list[
            prediction_protein_info.PredictionProteinInfo
        ] = prediction_util.get_prediction_name_and_seq_from_table(
            self._view.ui.table_pred_analysis_mono_prot_to_predict)
        self._view.close()
        self.job_input.emit(
            ("job_input",
             tmp_prediction_runs,
             self.prediction_configuration,
             True)
        )

    def check_if_prediction_and_analysis_should_be_done(self):
        if self._view.ui.checkbox_add_analysis.isChecked():
            self._view.ui.btn_pred_analysis_mono_go_analysis_setup.setText("Go")
            self._view.ui.lbl_pred_analysis_mono_to_analysis_setup.setText("To Analysis Setup")
        else:
            self._view.ui.btn_pred_analysis_mono_go_analysis_setup.setText("Predict")
            self._view.ui.lbl_pred_analysis_mono_to_analysis_setup.setText("Protein Structure(s)")

    def switch_monomer_pred_analysis_tab(self) -> None:
        """Switches the tabs from prediction to analysis and vice versa."""
        if self._view.ui.btn_pred_analysis_mono_go_analysis_setup.text() == "Predict":
            tmp_prediction_runs: list[
                prediction_protein_info.PredictionProteinInfo
            ] = prediction_util.get_prediction_name_and_seq_from_table(self._view.ui.table_pred_analysis_mono_prot_to_predict)
            self._view.close()
            self.job_input.emit(("job_input", tmp_prediction_runs, self.prediction_configuration, False))
        else:
            if self._view.ui.tabWidget.currentIndex() == 0:
                # goes from prediction to analysis
                self._view.ui.tabWidget.setTabEnabled(1, True)
                self._view.ui.tabWidget.setTabEnabled(0, False)
                self._view.ui.tabWidget.setCurrentIndex(1)
                gui_elements_to_show = [
                    self._view.ui.lbl_pred_analysis_mono_overview,
                    self._view.ui.list_pred_analysis_mono_overview,
                    self._view.ui.btn_pred_analysis_mono_add,
                    self._view.ui.btn_pred_analysis_mono_back_pred_setup,
                ]
                gui_elements_to_hide = [
                    self._view.ui.btn_pred_analysis_mono_remove,
                    self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
                    self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
                    self._view.ui.lbl_analysis_batch_vs_2,
                    self._view.ui.lbl_pred_analysis_mono_ref_chains,
                    self._view.ui.list_pred_analysis_mono_ref_chains,
                    self._view.ui.btn_pred_analysis_mono_back_4,
                    self._view.ui.btn_pred_analysis_mono_next_3,
                    self._view.ui.box_pred_analysis_mono_prot_struct_1,
                    self._view.ui.box_pred_analysis_mono_prot_struct_2,
                    self._view.ui.btn_pred_analysis_mono_back_3,
                    self._view.ui.btn_pred_analysis_mono_next_2,
                    self._view.ui.lbl_pred_analysis_mono_model_chains,
                    self._view.ui.list_pred_analysis_mono_model_chains,
                    self._view.ui.btn_pred_analysis_mono_back_5,
                    self._view.ui.btn_pred_analysis_mono_next_4,

                    
                    self._view.ui.btn_pred_analysis_mono_start,
                ]
                gui_utils.show_gui_elements(gui_elements_to_show)
                gui_utils.hide_gui_elements(gui_elements_to_hide)
                if self._view.ui.list_pred_analysis_mono_overview.count() > 0:
                    self._view.ui.btn_pred_analysis_mono_remove.show()
                    self._view.ui.btn_pred_analysis_mono_remove.setEnabled(False)
                    self._view.ui.btn_pred_analysis_mono_start.show()
                    
                    
            else:
                # goes from analysis to prediction
                if self._view.ui.list_pred_analysis_mono_overview.count() > 0:
                    gui_elements_to_show = [
                        self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
                        self._view.ui.table_pred_analysis_mono_prot_to_predict,
                        self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
                        self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
                    ]
                    gui_elements_to_hide = [
                        self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
                        self._view.ui.btn_pred_analysis_mono_seq_to_predict,
                        self._view.ui.lbl_pred_analysis_mono_prot_name,
                        self._view.ui.txt_pred_analysis_mono_prot_name,
                        self._view.ui.lbl_pred_analysis_mono_prot_name_status,
                        self._view.ui.btn_pred_analysis_mono_back,
                        self._view.ui.btn_pred_analysis_mono_next,
                        self._view.ui.lbl_pred_analysis_mono_seq_name,
                        self._view.ui.txt_pred_analysis_mono_seq_name,
                        self._view.ui.lbl_pred_analysis_mono_seq_name_status,
                        self._view.ui.btn_pred_analysis_mono_back_2,
                        self._view.ui.btn_pred_analysis_mono_add_protein,
                        self._view.ui.lbl_pred_mono_advanced_config_2,
                        self._view.ui.btn_pred_mono_advanced_config_2,
                    ]
                    gui_utils.show_gui_elements(gui_elements_to_show)
                    gui_utils.hide_gui_elements(gui_elements_to_hide)
                else:
                    gui_elements_to_show = [
                        self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
                        self._view.ui.table_pred_analysis_mono_prot_to_predict,
                        self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
                        self._view.ui.btn_pred_analysis_mono_seq_to_predict,
                        self._view.ui.lbl_pred_mono_advanced_config_2,
                        self._view.ui.btn_pred_mono_advanced_config_2,
                        self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
                        self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
                    ]
                    gui_elements_to_hide = [
                        self._view.ui.lbl_pred_analysis_mono_prot_name,
                        self._view.ui.txt_pred_analysis_mono_prot_name,
                        self._view.ui.lbl_pred_analysis_mono_prot_name_status,
                        self._view.ui.btn_pred_analysis_mono_back,
                        self._view.ui.btn_pred_analysis_mono_next,
                        self._view.ui.lbl_pred_analysis_mono_seq_name,
                        self._view.ui.txt_pred_analysis_mono_seq_name,
                        self._view.ui.lbl_pred_analysis_mono_seq_name_status,
                        self._view.ui.btn_pred_analysis_mono_back_2,
                        self._view.ui.btn_pred_analysis_mono_add_protein,
                    ]
                    gui_utils.show_gui_elements(gui_elements_to_show)
                    gui_utils.hide_gui_elements(gui_elements_to_hide)
                    self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove.setEnabled(False)
                self._view.ui.tabWidget.setTabEnabled(0, True)
                self._view.ui.tabWidget.setTabEnabled(1, False)
                self._view.ui.tabWidget.setCurrentIndex(0)

    def show_prediction_configuration(self) -> None:
        """Opens the prediction configuration dialog window."""
        config = dialog_advanced_prediction_configurations.DialogAdvancedPredictionConfigurations(
            self.prediction_configuration,
        )
        config.exec_()
        self.prediction_configuration.amber_force_field = config.prediction_config.amber_force_field
        self.prediction_configuration.templates = config.prediction_config.templates

    # <editor-fold desc="Monomer prediction + analysis">
    def _init_mono_pred_analysis_page(self) -> None:
        """Clears all text boxes and sets default values for the gui elements."""
        # <editor-fold desc="Prediction section">
        self._view.ui.txt_pred_analysis_mono_prot_name.clear()
        self._view.ui.txt_pred_analysis_mono_seq_name.clear()
        for i in range(self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount()):
            self._view.ui.table_pred_analysis_mono_prot_to_predict.removeRow(i)
        # sets up defaults: Prediction
        self._view.ui.btn_pred_analysis_mono_next.setEnabled(False)
        self._view.ui.btn_pred_analysis_mono_add_protein.setEnabled(False)
        self._view.ui.lbl_pred_analysis_mono_prot_name_status.setText("")
        self._view.ui.lbl_pred_analysis_mono_seq_name_status.setText("")

        # </editor-fold>

        # <editor-fold desc="Analysis section">
        self._view.ui.list_pred_analysis_mono_overview.clear()
        self._view.ui.btn_pred_analysis_mono_remove.hide()

        # </editor-fold>

    def display_monomer_pred_analysis(self) -> None:
        """Displays the monomer prediction + analysis page."""
        self._view.ui.btn_pred_analysis_mono_go_analysis_setup.setText("Predict")
        self._view.ui.lbl_pred_analysis_mono_to_analysis_setup.setText("Protein Structure(s)")
        # checks internet connection
        if not tools.check_internet_connectivity():
            gui_utils.no_internet_dialog()
            return

        self._init_mono_pred_analysis_page()
        self._view.ui.table_pred_analysis_mono_prot_to_predict.clear()
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setHorizontalHeaderItem(
            0,
            QtWidgets.QTableWidgetItem("Chain"),
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setHorizontalHeaderItem(
            1,
            QtWidgets.QTableWidgetItem("Sequence"),
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.resizeColumnsToContents()
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
            self._view.ui.table_pred_analysis_mono_prot_to_predict,
            self._view.ui.btn_pred_analysis_mono_seq_to_predict,
        ]
        gui_elements_to_hide = [
            self._view.ui.checkbox_add_analysis,
            self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            self._view.ui.lbl_pred_analysis_mono_prot_name,
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_prot_name_status,
            self._view.ui.btn_pred_analysis_mono_back,
            self._view.ui.btn_pred_analysis_mono_next,
            self._view.ui.lbl_pred_analysis_mono_seq_name,
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name_status,
            self._view.ui.btn_pred_analysis_mono_back_2,
            self._view.ui.btn_pred_analysis_mono_add_protein,
            self._view.ui.lbl_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self._view.ui.tabWidget.currentIndex() == 1:
            self._view.ui.tabWidget.setCurrentIndex(0)
        self._view.ui.tabWidget.setTabEnabled(1, False)
        self._view.ui.tabWidget.setTabEnabled(0, True)
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setEnabled(True)

    # <editor-fold desc="Sections">

    # <editor-fold desc="Prediction section">
    def mono_pred_analysis_validate_protein_name(self) -> None:
        """Validates the input of the protein name in real-time."""
        if safeguard.Safeguard.check_if_value_is_in_table_v_header(
                self._view.ui.txt_pred_analysis_mono_prot_name.text(),
                self._view.ui.table_pred_analysis_mono_prot_to_predict,
        ):
            self._view.ui.lbl_pred_analysis_mono_prot_name_status.setText("Protein name already used.")
            self._view.ui.btn_pred_analysis_mono_next.setEnabled(False)
        else:
            self._view.ui.btn_pred_analysis_mono_next.setEnabled(True)
            tools.validate_protein_name(
                self._view.ui.txt_pred_analysis_mono_prot_name,
                self._view.ui.lbl_pred_analysis_mono_prot_name_status,
                self._view.ui.btn_pred_analysis_mono_next,
            )

    def mono_pred_analysis_validate_protein_sequence(self) -> None:
        """Validates the input of the protein sequence in real-time."""
        tools.validate_protein_sequence(
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name_status,
            self._view.ui.btn_pred_analysis_mono_add_protein,
        )

    def mono_pred_analysis_check_if_table_is_empty(self) -> None:
        """Checks if the table proteins to predict is empty."""
        if self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() == 0:
            self._view.ui.btn_pred_analysis_mono_go_analysis_setup.setEnabled(False)
            gui_elements_to_show = [
                self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
                self._view.ui.table_pred_analysis_mono_prot_to_predict,
                self._view.ui.btn_pred_analysis_mono_seq_to_predict,
            ]
            gui_elements_to_hide = [
                self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
                self._view.ui.lbl_pred_analysis_mono_prot_name,
                self._view.ui.txt_pred_analysis_mono_prot_name,
                self._view.ui.lbl_pred_analysis_mono_prot_name_status,
                self._view.ui.btn_pred_analysis_mono_back,
                self._view.ui.btn_pred_analysis_mono_next,
                self._view.ui.lbl_pred_analysis_mono_seq_name,
                self._view.ui.txt_pred_analysis_mono_seq_name,
                self._view.ui.lbl_pred_analysis_mono_seq_name_status,
                self._view.ui.btn_pred_analysis_mono_back_2,
                self._view.ui.btn_pred_analysis_mono_add_protein,
                self._view.ui.lbl_pred_mono_advanced_config_2,
                self._view.ui.btn_pred_mono_advanced_config_2,
                self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
                self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
                self._view.ui.checkbox_add_analysis,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
        else:
            gui_elements_to_show = [
                self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
                self._view.ui.table_pred_analysis_mono_prot_to_predict,
                self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
                self._view.ui.btn_pred_analysis_mono_seq_to_predict,
                self._view.ui.lbl_pred_mono_advanced_config_2,
                self._view.ui.btn_pred_mono_advanced_config_2,
                self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
                self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
            ]
            gui_elements_to_hide = [
                self._view.ui.lbl_pred_analysis_mono_prot_name,
                self._view.ui.txt_pred_analysis_mono_prot_name,
                self._view.ui.lbl_pred_analysis_mono_prot_name_status,
                self._view.ui.btn_pred_analysis_mono_back,
                self._view.ui.btn_pred_analysis_mono_next,
                self._view.ui.lbl_pred_analysis_mono_seq_name,
                self._view.ui.txt_pred_analysis_mono_seq_name,
                self._view.ui.lbl_pred_analysis_mono_seq_name_status,
                self._view.ui.btn_pred_analysis_mono_back_2,
                self._view.ui.btn_pred_analysis_mono_add_protein,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.btn_pred_analysis_mono_go_analysis_setup.setEnabled(True)

    def setup_defaults_monomer_prediction_analysis(self) -> None:
        """Sets up default values for the prediction tab."""
        # clears everything
        self._view.ui.txt_pred_analysis_mono_prot_name.clear()
        self._view.ui.txt_pred_analysis_mono_seq_name.clear()
        # sets up defaults: Prediction
        self._view.ui.btn_pred_analysis_mono_next.setEnabled(False)
        self._view.ui.btn_pred_analysis_mono_add_protein.setEnabled(False)
        self._view.ui.lbl_pred_analysis_mono_prot_name_status.setText("")
        self._view.ui.lbl_pred_analysis_mono_seq_name_status.setText("")

    def mono_pred_analysis_add_seq_to_predict(self) -> None:
        """Shows the gui elements for the protein name."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
            self._view.ui.table_pred_analysis_mono_prot_to_predict,
            self._view.ui.lbl_pred_analysis_mono_prot_name,
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_prot_name_status,
            self._view.ui.btn_pred_analysis_mono_back,
            self._view.ui.btn_pred_analysis_mono_next,
        ]
        gui_utils.enable_text_box(
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_prot_name,
        )
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            self._view.ui.btn_pred_analysis_mono_seq_to_predict,
            self._view.ui.lbl_pred_analysis_mono_seq_name,
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name_status,
            self._view.ui.btn_pred_analysis_mono_back_2,
            self._view.ui.btn_pred_analysis_mono_add_protein,
            self._view.ui.lbl_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
        ]
        gui_utils.disable_text_box(
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name,
        )
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.btn_pred_analysis_mono_next.setEnabled(False)
        self._view.ui.txt_pred_analysis_mono_prot_name.clear()
        if self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() > 0:
            try:
                self._view.ui.table_pred_analysis_mono_prot_to_predict.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection on Local Monomer Prediction in overview table.")

    def mono_pred_analysis_back(self) -> None:
        """Hides the gui elements for the protein name."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
            self._view.ui.table_pred_analysis_mono_prot_to_predict,
            self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            self._view.ui.btn_pred_analysis_mono_seq_to_predict,
        ]
        gui_elements_to_hide = [
            self._view.ui.lbl_pred_analysis_mono_prot_name,
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_prot_name_status,
            self._view.ui.btn_pred_analysis_mono_back,
            self._view.ui.btn_pred_analysis_mono_next,
            self._view.ui.lbl_pred_analysis_mono_seq_name,
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name_status,
            self._view.ui.btn_pred_analysis_mono_back_2,
            self._view.ui.btn_pred_analysis_mono_add_protein,
            self._view.ui.lbl_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.mono_pred_analysis_check_if_table_is_empty()
        self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove.setEnabled(False)

    def mono_pred_analysis_next(self) -> None:
        """Shows the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
            self._view.ui.table_pred_analysis_mono_prot_to_predict,
            self._view.ui.lbl_pred_analysis_mono_prot_name,
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name,
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name_status,
            self._view.ui.btn_pred_analysis_mono_back_2,
            self._view.ui.btn_pred_analysis_mono_add_protein,
        ]
        gui_utils.enable_text_box(
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name,
        )
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            self._view.ui.btn_pred_analysis_mono_seq_to_predict,
            self._view.ui.lbl_pred_analysis_mono_prot_name_status,
            self._view.ui.btn_pred_analysis_mono_back,
            self._view.ui.btn_pred_analysis_mono_next,
            self._view.ui.lbl_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
        ]
        gui_utils.disable_text_box(
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_prot_name,
        )
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.txt_pred_analysis_mono_seq_name.clear()

    def mono_pred_analysis_back_2(self) -> None:
        """Hides the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
            self._view.ui.table_pred_analysis_mono_prot_to_predict,
            self._view.ui.lbl_pred_analysis_mono_prot_name,
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_prot_name_status,
            self._view.ui.btn_pred_analysis_mono_back,
            self._view.ui.btn_pred_analysis_mono_next,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            self._view.ui.btn_pred_analysis_mono_seq_to_predict,
            self._view.ui.lbl_pred_analysis_mono_seq_name,
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name_status,
            self._view.ui.btn_pred_analysis_mono_back_2,
            self._view.ui.btn_pred_analysis_mono_add_protein,
            self._view.ui.lbl_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
        ]
        gui_utils.enable_text_box(
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_prot_name,
        )
        gui_utils.disable_text_box(
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name,
        )
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)

    def mono_pred_analysis_add_protein(self) -> None:
        """Adds the protein to the list of proteins to predict."""
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setRowCount(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() + 1,
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.insertRow(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() + 1,
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setItem(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() - 1,
            0,
            QtWidgets.QTableWidgetItem("A"),
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setItem(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() - 1,
            1,
            QtWidgets.QTableWidgetItem(self._view.ui.txt_pred_analysis_mono_seq_name.toPlainText()),
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setVerticalHeaderItem(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() - 1,
            QtWidgets.QTableWidgetItem(self._view.ui.txt_pred_analysis_mono_prot_name.text()),
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.resizeColumnsToContents()
        self.mono_pred_analysis_check_if_table_is_empty()
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_prot_to_predict,
            self._view.ui.table_pred_analysis_mono_prot_to_predict,
            self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove,
            self._view.ui.btn_pred_analysis_mono_seq_to_predict,
            self._view.ui.lbl_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_mono_advanced_config_2,
            self._view.ui.btn_pred_analysis_mono_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_mono_to_analysis_setup,
            self._view.ui.checkbox_add_analysis,
        ]
        gui_utils.enable_text_box(
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_prot_name,
        )
        gui_elements_to_hide = [
            self._view.ui.lbl_pred_analysis_mono_prot_name,
            self._view.ui.txt_pred_analysis_mono_prot_name,
            self._view.ui.lbl_pred_analysis_mono_prot_name_status,
            self._view.ui.btn_pred_analysis_mono_back,
            self._view.ui.btn_pred_analysis_mono_next,
            self._view.ui.lbl_pred_analysis_mono_seq_name,
            self._view.ui.txt_pred_analysis_mono_seq_name,
            self._view.ui.lbl_pred_analysis_mono_seq_name_status,
            self._view.ui.btn_pred_analysis_mono_back_2,
            self._view.ui.btn_pred_analysis_mono_add_protein,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.btn_pred_analysis_mono_go_analysis_setup.setEnabled(True)
        self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove.setEnabled(False)
        self.setup_defaults_monomer_prediction_analysis()

    def mono_pred_analysis_prediction_overview_item_clicked(self) -> None:
        """Enables the remove button."""
        self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove.setEnabled(True)

    def mono_pred_analysis_add_protein_to_predict(self) -> None:
        """Needs to be removed."""
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setRowCount(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() + 1,
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.insertRow(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() + 1,
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setItem(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() - 1,
            0,
            QtWidgets.QTableWidgetItem("A"),
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setItem(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() - 1,
            1,
            QtWidgets.QTableWidgetItem(self._view.ui.txt_pred_analysis_mono_seq_name.toPlainText()),
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.setVerticalHeaderItem(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount() - 1,
            QtWidgets.QTableWidgetItem(self._view.ui.txt_pred_analysis_mono_prot_name.text()),
        )
        self._view.ui.table_pred_analysis_mono_prot_to_predict.resizeColumnsToContents()
        self.mono_pred_analysis_check_if_table_is_empty()
        self.setup_defaults_monomer_prediction_analysis()

    def mono_pred_analysis_remove_protein_to_predict(self) -> None:
        """Removes the selected protein from the list of proteins to predict."""
        self._view.ui.table_pred_analysis_mono_prot_to_predict.removeRow(
            self._view.ui.table_pred_analysis_mono_prot_to_predict.currentRow(),
        )
        self.mono_pred_analysis_check_if_table_is_empty()
        self._view.ui.btn_pred_analysis_mono_seq_to_predict_remove.setEnabled(False)

    # </editor-fold>

    # <editor-fold desc="Analysis section">
    def mono_pred_analysis_structure_analysis_add(self) -> None:
        """Shows the gui elements for the selection of the two proteins."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_overview,
            self._view.ui.list_pred_analysis_mono_overview,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
            self._view.ui.box_pred_analysis_mono_prot_struct_1,
            self._view.ui.lbl_analysis_batch_vs_2,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
            self._view.ui.box_pred_analysis_mono_prot_struct_2,
            self._view.ui.btn_pred_analysis_mono_back_3,
            self._view.ui.btn_pred_analysis_mono_next_2,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_mono_remove,
            self._view.ui.btn_pred_analysis_mono_add,
            self._view.ui.lbl_pred_analysis_mono_ref_chains,
            self._view.ui.list_pred_analysis_mono_ref_chains,
            self._view.ui.btn_pred_analysis_mono_back_4,
            self._view.ui.btn_pred_analysis_mono_next_3,
            self._view.ui.lbl_pred_analysis_mono_model_chains,
            self._view.ui.list_pred_analysis_mono_model_chains,
            self._view.ui.btn_pred_analysis_mono_back_5,
            self._view.ui.btn_pred_analysis_mono_next_4,
            
            self._view.ui.btn_pred_analysis_mono_start,
            self._view.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.lbl_pred_analysis_mono_prot_struct_1.clear()
        self._view.ui.lbl_pred_analysis_mono_prot_struct_2.clear()
        self._view.ui.lbl_pred_analysis_mono_prot_struct_1.setText("Protein structure 1")
        self._view.ui.lbl_pred_analysis_mono_prot_struct_2.setText("Protein structure 2")
        self.fill_mono_pred_analysis_protein_boxes()
        if self._view.ui.list_pred_analysis_mono_overview.count() > 0:
            try:
                self._view.ui.list_pred_analysis_mono_overview.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection in struction analysis overview.")

    def mono_pred_analysis_structure_analysis_back_3(self) -> None:
        """Hides the gui elements to select the two proteins."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_overview,
            self._view.ui.list_pred_analysis_mono_overview,
            self._view.ui.btn_pred_analysis_mono_add,
            self._view.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_mono_remove,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_2,
            self._view.ui.lbl_pred_analysis_mono_ref_chains,
            self._view.ui.list_pred_analysis_mono_ref_chains,
            self._view.ui.btn_pred_analysis_mono_back_4,
            self._view.ui.btn_pred_analysis_mono_next_3,
            self._view.ui.box_pred_analysis_mono_prot_struct_1,
            self._view.ui.box_pred_analysis_mono_prot_struct_2,
            self._view.ui.btn_pred_analysis_mono_back_3,
            self._view.ui.btn_pred_analysis_mono_next_2,
            self._view.ui.lbl_pred_analysis_mono_model_chains,
            self._view.ui.list_pred_analysis_mono_model_chains,
            self._view.ui.btn_pred_analysis_mono_back_5,
            self._view.ui.btn_pred_analysis_mono_next_4,
            
            self._view.ui.btn_pred_analysis_mono_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self._view.ui.list_pred_analysis_mono_overview.count() > 0:
            self._view.ui.btn_pred_analysis_mono_remove.show()
            self._view.ui.btn_pred_analysis_mono_remove.setEnabled(False)
            self._view.ui.btn_pred_analysis_mono_start.show()
            
            

    def mono_pred_analysis_structure_analysis_next_3(self) -> None:
        """Shows the gui elements to select the chains of protein 2."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_overview,
            self._view.ui.list_pred_analysis_mono_overview,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_2,
            self._view.ui.lbl_pred_analysis_mono_ref_chains,
            self._view.ui.list_pred_analysis_mono_ref_chains,
            self._view.ui.lbl_pred_analysis_mono_model_chains,
            self._view.ui.list_pred_analysis_mono_model_chains,
            self._view.ui.btn_pred_analysis_mono_back_5,
            self._view.ui.btn_pred_analysis_mono_next_4,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_mono_remove,
            self._view.ui.btn_pred_analysis_mono_add,
            self._view.ui.box_pred_analysis_mono_prot_struct_1,
            self._view.ui.box_pred_analysis_mono_prot_struct_2,
            self._view.ui.btn_pred_analysis_mono_back_3,
            self._view.ui.btn_pred_analysis_mono_next_2,
            self._view.ui.btn_pred_analysis_mono_back_4,
            self._view.ui.btn_pred_analysis_mono_next_3,
            
            self._view.ui.btn_pred_analysis_mono_start,
            self._view.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.list_pred_analysis_mono_model_chains.clear()
        self._view.ui.list_pred_analysis_mono_ref_chains.setEnabled(False)
        self._view.ui.btn_pred_analysis_mono_next_4.setEnabled(False)

        for i in range(self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount()):
            if (
                    self._view.ui.table_pred_analysis_mono_prot_to_predict.verticalHeaderItem(i).text()
                    == self._view.ui.box_pred_analysis_mono_prot_struct_2.currentText()
            ):
                self._view.ui.list_pred_analysis_mono_model_chains.addItem(
                    self._view.ui.table_pred_analysis_mono_prot_to_predict.item(i, 0).text(),
                )
        if self._view.ui.list_pred_analysis_mono_model_chains.count() == 0:
            tmp_protein = self._interface_manager.get_current_project().search_protein(
                self._view.ui.box_pred_analysis_mono_prot_struct_2.currentText(),
            )
            for tmp_chain in tmp_protein.chains:
                if tmp_chain.chain_type == "protein_chain":
                    self._view.ui.list_pred_analysis_mono_model_chains.addItem(tmp_chain.chain_letter)
        if self._view.ui.list_pred_analysis_mono_model_chains.count() == 1:
            self._view.ui.lbl_pred_analysis_mono_model_chains.setText(
                f"Select chain in protein structure {self._view.ui.lbl_pred_analysis_mono_prot_struct_2.text()}.",
            )
        else:
            self._view.ui.lbl_pred_analysis_mono_model_chains.setText(
                f"Select {len(self._view.ui.list_pred_analysis_mono_model_chains.selectedItems())} chains "
                f"in protein structure {self._view.ui.lbl_pred_analysis_mono_prot_struct_2.text()}.",
            )

    def mono_pred_analysis_structure_analysis_back_4(self) -> None:
        """Hides the gui elements to select the chains in protein 1."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_overview,
            self._view.ui.list_pred_analysis_mono_overview,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
            self._view.ui.box_pred_analysis_mono_prot_struct_1,
            self._view.ui.lbl_analysis_batch_vs_2,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
            self._view.ui.box_pred_analysis_mono_prot_struct_2,
            self._view.ui.btn_pred_analysis_mono_back_3,
            self._view.ui.btn_pred_analysis_mono_next_2,
            self._view.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_mono_remove,
            self._view.ui.btn_pred_analysis_mono_add,
            self._view.ui.lbl_pred_analysis_mono_ref_chains,
            self._view.ui.list_pred_analysis_mono_ref_chains,
            self._view.ui.btn_pred_analysis_mono_back_4,
            self._view.ui.btn_pred_analysis_mono_next_3,
            self._view.ui.lbl_pred_analysis_mono_model_chains,
            self._view.ui.list_pred_analysis_mono_model_chains,
            self._view.ui.btn_pred_analysis_mono_back_5,
            self._view.ui.btn_pred_analysis_mono_next_4,
            
            self._view.ui.btn_pred_analysis_mono_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.lbl_pred_analysis_mono_prot_struct_1.setText("Protein structure 1")
        self._view.ui.lbl_pred_analysis_mono_prot_struct_2.setText("Protein structure 2")

    def mono_pred_analysis_structure_analysis_next_4(self) -> None:
        """Adds the protein pair to the list of protein pairs to analyze."""
        gui_elements_to_show = [
            self._view.ui.btn_pred_analysis_mono_remove,
            self._view.ui.btn_pred_analysis_mono_add,
            self._view.ui.lbl_pred_analysis_mono_overview,
            self._view.ui.list_pred_analysis_mono_overview,
            
            self._view.ui.btn_pred_analysis_mono_start,
            self._view.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self._view.ui.box_pred_analysis_mono_prot_struct_1,
            self._view.ui.box_pred_analysis_mono_prot_struct_2,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_2,
            self._view.ui.lbl_pred_analysis_mono_ref_chains,
            self._view.ui.list_pred_analysis_mono_ref_chains,
            self._view.ui.lbl_pred_analysis_mono_model_chains,
            self._view.ui.list_pred_analysis_mono_model_chains,
            self._view.ui.btn_pred_analysis_mono_back_3,
            self._view.ui.btn_pred_analysis_mono_next_2,
            self._view.ui.btn_pred_analysis_mono_back_4,
            self._view.ui.btn_pred_analysis_mono_next_3,
            self._view.ui.btn_pred_analysis_mono_back_5,
            self._view.ui.btn_pred_analysis_mono_next_4,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        prot_1_name = self._view.ui.lbl_pred_analysis_mono_prot_struct_1.text()
        prot_1_chains = []
        for chain in self._view.ui.list_pred_analysis_mono_ref_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self._view.ui.lbl_pred_analysis_mono_prot_struct_2.text()
        prot_2_chains = []
        for chain in self._view.ui.list_pred_analysis_mono_model_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        item = QtWidgets.QListWidgetItem(analysis_name)
        self._view.ui.list_pred_analysis_mono_overview.addItem(item)
        self._view.ui.btn_pred_analysis_mono_remove.setEnabled(False)

    def mono_pred_analysis_structure_analysis_back_5(self) -> None:
        """Hides the gui elements to select the chains in protein 2."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_overview,
            self._view.ui.list_pred_analysis_mono_overview,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_2,
            self._view.ui.lbl_pred_analysis_mono_ref_chains,
            self._view.ui.list_pred_analysis_mono_ref_chains,
            self._view.ui.btn_pred_analysis_mono_back_4,
            self._view.ui.btn_pred_analysis_mono_next_3,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_mono_remove,
            self._view.ui.btn_pred_analysis_mono_add,
            self._view.ui.box_pred_analysis_mono_prot_struct_1,
            self._view.ui.box_pred_analysis_mono_prot_struct_2,
            self._view.ui.btn_pred_analysis_mono_back_3,
            self._view.ui.btn_pred_analysis_mono_next_2,
            self._view.ui.btn_pred_analysis_mono_back_5,
            self._view.ui.btn_pred_analysis_mono_next_4,
            
            self._view.ui.btn_pred_analysis_mono_start,
            self._view.ui.lbl_pred_analysis_mono_model_chains,
            self._view.ui.list_pred_analysis_mono_model_chains,
            self._view.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.list_pred_analysis_mono_ref_chains.setEnabled(True)

        # tmp_protein = self._interface_manager.get_current_project().search_protein(self._view.ui.box_pred_analysis_mono_prot_struct_2.currentText())
        # for tmp_chain in tmp_protein.chains:
        #     if tmp_chain.chain_type == "protein_chain":
        #         self._view.ui.list_pred_analysis_mono_ref_chains.addItem(tmp_chain.chain_letter)

    def mono_pred_analysis_structure_analysis_overview_clicked(self) -> None:
        """Enables the remove button."""
        self._view.ui.btn_pred_analysis_mono_remove.setEnabled(True)

    def fill_mono_pred_analysis_protein_boxes(self) -> None:
        """Fills the combo box of the protein structures."""
        protein_names = []
        for i in range(self._view.ui.table_pred_analysis_mono_prot_to_predict.rowCount()):
            protein_names.append(
                self._view.ui.table_pred_analysis_mono_prot_to_predict.verticalHeaderItem(i).text())
        for tmp_protein in self._interface_manager.get_current_project().proteins:
            protein_names.append(tmp_protein.get_molecule_object())
        protein_names.insert(0, "")
        self._view.ui.box_pred_analysis_mono_prot_struct_1.clear()
        self._view.ui.box_pred_analysis_mono_prot_struct_2.clear()
        gui_utils.fill_combo_box(self._view.ui.box_pred_analysis_mono_prot_struct_1, protein_names)
        gui_utils.fill_combo_box(self._view.ui.box_pred_analysis_mono_prot_struct_2, protein_names)

    def remove_mono_pred_analysis_analysis_run(self) -> None:
        """Removes a selected protein pair form the list of protein pairs to analyze."""
        self._view.ui.list_pred_analysis_mono_overview.takeItem(
            self._view.ui.list_pred_analysis_mono_overview.currentRow(),
        )
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_mono_overview,
            self._view.ui.list_pred_analysis_mono_overview,
            self._view.ui.btn_pred_analysis_mono_add,
            self._view.ui.btn_pred_analysis_mono_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_mono_remove,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_1,
            self._view.ui.lbl_pred_analysis_mono_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_2,
            self._view.ui.lbl_pred_analysis_mono_ref_chains,
            self._view.ui.list_pred_analysis_mono_ref_chains,
            self._view.ui.btn_pred_analysis_mono_back_4,
            self._view.ui.btn_pred_analysis_mono_next_3,
            self._view.ui.box_pred_analysis_mono_prot_struct_1,
            self._view.ui.box_pred_analysis_mono_prot_struct_2,
            self._view.ui.btn_pred_analysis_mono_back_3,
            self._view.ui.btn_pred_analysis_mono_next_2,
            self._view.ui.lbl_pred_analysis_mono_model_chains,
            self._view.ui.list_pred_analysis_mono_model_chains,
            self._view.ui.btn_pred_analysis_mono_back_5,
            self._view.ui.btn_pred_analysis_mono_next_4,
            
            self._view.ui.btn_pred_analysis_mono_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self._view.ui.list_pred_analysis_mono_overview.count() > 0:
            self._view.ui.btn_pred_analysis_mono_remove.show()
            self._view.ui.btn_pred_analysis_mono_remove.setEnabled(False)
            self._view.ui.btn_pred_analysis_mono_start.show()
            
            
        # if self._view.ui.list_pred_analysis_mono_overview.count() == 0:
        #
        #     self._view.ui.btn_pred_analysis_mono_back_pred_setup.show()
        #     self._view.ui.btn_pred_analysis_mono_remove.hide()

    def check_mono_pred_analysis_if_same_no_of_chains_selected(self) -> None:
        """Checks if the same number of chains were selected."""
        self._view.ui.btn_pred_analysis_mono_next_4.setEnabled(False)
        if self.no_of_selected_chains == len(self._view.ui.list_pred_analysis_mono_model_chains.selectedItems()):
            self._view.ui.btn_pred_analysis_mono_next_4.setEnabled(True)

        prot_1_name = self._view.ui.lbl_pred_analysis_mono_prot_struct_1.text()
        prot_1_chains = []
        for chain in self._view.ui.list_pred_analysis_mono_ref_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self._view.ui.lbl_pred_analysis_mono_prot_struct_2.text()
        prot_2_chains = []
        for chain in self._view.ui.list_pred_analysis_mono_model_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        for tmp_row in range(self._view.ui.list_pred_analysis_mono_overview.count()):
            if analysis_name == self._view.ui.list_pred_analysis_mono_overview.item(tmp_row).text():
                self._view.ui.btn_pred_analysis_mono_next_4.setEnabled(False)
                return

    def check_mono_pred_analysis_if_prot_structs_are_filled(self) -> None:
        """Checks if two proteins were selected."""
        prot_1 = self._view.ui.box_pred_analysis_mono_prot_struct_1.itemText(
            self._view.ui.box_pred_analysis_mono_prot_struct_1.currentIndex(),
        )
        prot_2 = self._view.ui.box_pred_analysis_mono_prot_struct_2.itemText(
            self._view.ui.box_pred_analysis_mono_prot_struct_2.currentIndex(),
        )
        if prot_1 != "" and prot_2 != "":
            self._view.ui.btn_pred_analysis_mono_next_2.setEnabled(True)
        else:
            self._view.ui.btn_pred_analysis_mono_next_2.setEnabled(False)

    def count_mono_pred_analysis_selected_chains_for_prot_struct_1(self) -> None:
        """Counts the number of chains selected in protein 1."""
        self.no_of_selected_chains = len(self._view.ui.list_pred_analysis_mono_ref_chains.selectedItems())
        if self.no_of_selected_chains > 0:
            self._view.ui.btn_pred_analysis_mono_next_3.setEnabled(True)
        else:
            self._view.ui.btn_pred_analysis_mono_next_3.setEnabled(False)

    # </editor-fold>

    # </editor-fold>

    # </editor-fold>
