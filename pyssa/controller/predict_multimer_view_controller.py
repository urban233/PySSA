import os
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal

from pyssa.controller import interface_manager
from pyssa.gui.ui.dialogs import dialog_advanced_prediction_configurations
from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.gui.ui.styles import styles
from pyssa.gui.ui.views import distance_analysis_view, predict_monomer_view, predict_multimer_view
from pyssa.internal.data_structures.data_classes import prediction_protein_info, prediction_configuration
from pyssa.internal.thread import tasks
from pyssa.io_pyssa import safeguard
from pyssa.presenter import main_presenter_async
from pyssa.util import gui_utils, tools, constants, exit_codes, prediction_util


class PredictMultimerViewController(QtCore.QObject):

    job_input = pyqtSignal(tuple)

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager"):
        super().__init__()
        self._interface_manager = the_interface_manager
        self._view: predict_multimer_view.PredictMultimerView = the_interface_manager.get_predict_multimer_view()
        self.prediction_configuration = prediction_configuration.PredictionConfiguration(True, "pdb70")
        self.display_multimer_pred_analysis()
        self._connect_all_ui_elements_to_slot_functions()

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.checkbox_add_analysis.clicked.connect(self.check_if_prediction_and_analysis_should_be_done)

        # <editor-fold desc="Multimer Prediction + Analysis page">
        # <editor-fold desc="Prediction section">
        self._view.ui.btn_pred_analysis_multi_prot_to_predict_add.clicked.connect(self.multi_pred_analysis_add)
        self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove.clicked.connect(
            self.multi_pred_analysis_remove_protein_to_predict,
        )
        self._view.ui.btn_pred_analysis_multi_next.clicked.connect(self.multi_pred_analysis_next)
        self._view.ui.btn_pred_analysis_multi_back.clicked.connect(self.multi_pred_analysis_back)
        self._view.ui.btn_pred_analysis_multi_prot_seq_add.clicked.connect(
            self.multi_pred_analysis_add_sequence_to_list,
        )
        self._view.ui.btn_pred_analysis_multi_prot_seq_overview_remove.clicked.connect(
            self.multi_pred_analysis_remove_sequence_to_list,
        )
        self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2.clicked.connect(
            self.multi_pred_analysis_add_protein_to_predict,
        )
        self._view.ui.btn_pred_analysis_multi_back_2.clicked.connect(self.multi_pred_analysis_back_2)

        self._view.ui.txt_pred_analysis_multi_prot_name.textChanged.connect(
            self.multi_pred_analysis_validate_protein_name,
        )
        self._view.ui.txt_pred_analysis_multi_prot_seq.textChanged.connect(
            self.multi_pred_analysis_validate_protein_sequence,
        )
        self._view.ui.btn_pred_analysis_multi_advanced_config.clicked.connect(self.show_prediction_configuration)
        self._view.ui.btn_pred_analysis_multi_go_analysis_setup.clicked.connect(self.switch_multimer_pred_analysis_tab)

        self._view.ui.list_pred_analysis_multi_prot_seq_overview.clicked.connect(
            self.multi_pred_analysis_prot_seq_overview_item_changed,
        )
        self._view.ui.table_pred_analysis_multi_prot_to_predict.itemSelectionChanged.connect(
            self.multi_pred_analysis_prot_to_predict_item_changed,
        )

        # </editor-fold>

        # <editor-fold desc="Analysis section">
        self._view.ui.btn_pred_analysis_multi_add.clicked.connect(self.multi_pred_analysis_structure_analysis_add)
        self._view.ui.btn_pred_analysis_multi_remove.clicked.connect(self.remove_multi_pred_analysis_analysis_run)
        self._view.ui.btn_pred_analysis_multi_back_3.clicked.connect(self.multi_pred_analysis_structure_analysis_back_3)
        self._view.ui.btn_pred_analysis_multi_next_2.clicked.connect(self.multi_pred_analysis_structure_analysis_next_2)
        self._view.ui.btn_pred_analysis_multi_back_4.clicked.connect(self.multi_pred_analysis_structure_analysis_back_4)
        self._view.ui.btn_pred_analysis_multi_next_3.clicked.connect(self.multi_pred_analysis_structure_analysis_next_3)
        self._view.ui.btn_pred_analysis_multi_back_5.clicked.connect(self.multi_pred_analysis_structure_analysis_back_5)
        self._view.ui.btn_pred_analysis_multi_next_4.clicked.connect(self.multi_pred_analysis_structure_analysis_next_4)
        self._view.ui.box_pred_analysis_multi_prot_struct_1.currentIndexChanged.connect(
            self.check_multi_pred_analysis_if_prot_structs_are_filled,
        )
        self._view.ui.box_pred_analysis_multi_prot_struct_2.currentIndexChanged.connect(
            self.check_multi_pred_analysis_if_prot_structs_are_filled,
        )
        self._view.ui.list_pred_analysis_multi_ref_chains.itemSelectionChanged.connect(
            self.count_multi_pred_analysis_selected_chains_for_prot_struct_1,
        )
        self._view.ui.list_pred_analysis_multi_model_chains.itemSelectionChanged.connect(
            self.check_multi_pred_analysis_if_same_no_of_chains_selected,
        )
        self._view.ui.btn_pred_analysis_multi_back_pred_setup.clicked.connect(self.switch_multimer_pred_analysis_tab)
        #self._view.ui.btn_pred_analysis_multi_start.clicked.connect(self.start_multimer_prediction_analysis)
        self._view.ui.list_pred_analysis_multi_overview.clicked.connect(
            self.multi_pred_analysis_structure_analysis_overview_clicked,
        )
        # </editor-fold>

        # </editor-fold>

    def multi_pred_analysis_structure_analysis_next_2(self) -> None:
        """Shows the gui elements to select the chains in protein 1."""
        self._view.wait_spinner.start()
        tmp_proteins_to_predict: list[str] = []
        for i in range(self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount()):
            tmp_proteins_to_predict.append(
                self._view.ui.table_pred_analysis_multi_prot_to_predict.verticalHeaderItem(i).text(),
            )
        self._active_task = tasks.Task(
            target=main_presenter_async.check_chains_for_subsequent_analysis,
            args=(
                self._view.ui.box_pred_analysis_multi_prot_struct_1.currentText(),
                self._view.ui.box_pred_analysis_multi_prot_struct_2.currentText(),
                self._interface_manager.get_current_project(),
                tmp_proteins_to_predict,
            ),
            post_func=self.__await_multi_pred_analysis_structure_analysis_next_2,
        )
        self._active_task.start()

    def __await_multi_pred_analysis_structure_analysis_next_2(self, result: tuple) -> None:
        # fixme: something is not quite right with the detection of the chains!
        _, tmp_analysis_name = result
        if tmp_analysis_name != "":
            gui_elements_to_show = [
                self._view.ui.btn_pred_analysis_multi_remove,
                self._view.ui.btn_pred_analysis_multi_add,
                self._view.ui.lbl_pred_analysis_multi_overview,
                self._view.ui.list_pred_analysis_multi_overview,
                self._view.ui.btn_pred_analysis_multi_start,
                self._view.ui.btn_pred_analysis_multi_back_pred_setup,
            ]
            gui_elements_to_hide = [
                self._view.ui.box_pred_analysis_multi_prot_struct_1,
                self._view.ui.box_pred_analysis_multi_prot_struct_2,
                self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
                self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
                self._view.ui.lbl_analysis_batch_vs_3,
                self._view.ui.lbl_pred_analysis_multi_ref_chains,
                self._view.ui.list_pred_analysis_multi_ref_chains,
                self._view.ui.lbl_pred_analysis_multi_model_chains,
                self._view.ui.list_pred_analysis_multi_model_chains,
                self._view.ui.btn_pred_analysis_multi_back_3,
                self._view.ui.btn_pred_analysis_multi_next_2,
                self._view.ui.btn_pred_analysis_multi_back_4,
                self._view.ui.btn_pred_analysis_multi_next_3,
                self._view.ui.btn_pred_analysis_multi_back_5,
                self._view.ui.btn_pred_analysis_multi_next_4,
                self._view.ui.cb_pred_analysis_multi_images,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            item = QtWidgets.QListWidgetItem(tmp_analysis_name)
            self._view.ui.list_pred_analysis_multi_overview.addItem(item)
            self._view.ui.btn_pred_analysis_multi_remove.setEnabled(False)
            styles.color_button_ready(self._view.ui.btn_pred_analysis_multi_start)
        else:
            gui_elements_to_show = [
                self._view.ui.lbl_pred_analysis_multi_overview,
                self._view.ui.list_pred_analysis_multi_overview,
                self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
                self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
                self._view.ui.lbl_analysis_batch_vs_3,
                self._view.ui.lbl_pred_analysis_multi_ref_chains,
                self._view.ui.list_pred_analysis_multi_ref_chains,
                self._view.ui.btn_pred_analysis_multi_back_4,
                self._view.ui.btn_pred_analysis_multi_next_3,
            ]
            gui_elements_to_hide = [
                self._view.ui.btn_pred_analysis_multi_remove,
                self._view.ui.btn_pred_analysis_multi_add,
                self._view.ui.box_pred_analysis_multi_prot_struct_1,
                self._view.ui.box_pred_analysis_multi_prot_struct_2,
                self._view.ui.btn_pred_analysis_multi_back_3,
                self._view.ui.btn_pred_analysis_multi_next_2,
                self._view.ui.lbl_pred_analysis_multi_model_chains,
                self._view.ui.list_pred_analysis_multi_model_chains,
                self._view.ui.btn_pred_analysis_multi_back_5,
                self._view.ui.btn_pred_analysis_multi_next_4,
                self._view.ui.cb_pred_analysis_multi_images,
                self._view.ui.btn_pred_analysis_multi_start,
                self._view.ui.btn_pred_analysis_multi_back_pred_setup,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.lbl_pred_analysis_multi_prot_struct_1.setText(
                self._view.ui.box_pred_analysis_multi_prot_struct_1.currentText(),
            )
            self._view.ui.lbl_pred_analysis_multi_prot_struct_2.setText(
                self._view.ui.box_pred_analysis_multi_prot_struct_2.currentText(),
            )
            self._view.ui.list_pred_analysis_multi_ref_chains.clear()
            self._view.ui.btn_pred_analysis_multi_next_3.setEnabled(False)
            self._view.ui.list_pred_analysis_multi_ref_chains.setEnabled(True)

            for i in range(self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount()):
                if (
                    self._view.ui.table_pred_analysis_multi_prot_to_predict.verticalHeaderItem(i).text()
                    == self._view.ui.box_pred_analysis_multi_prot_struct_1.currentText()
                ):
                    self._view.ui.list_pred_analysis_multi_ref_chains.addItem(
                        self._view.ui.table_pred_analysis_multi_prot_to_predict.item(i, 0).text(),
                    )
            if self._view.ui.list_pred_analysis_multi_ref_chains.count() == 0:
                tmp_protein = self._interface_manager.get_current_project().search_protein(
                    self._view.ui.box_pred_analysis_multi_prot_struct_1.currentText(),
                )
                for tmp_chain in tmp_protein.chains:
                    if tmp_chain.chain_type == "protein_chain":
                        self._view.ui.list_pred_analysis_multi_ref_chains.addItem(tmp_chain.chain_letter)
            if self._view.ui.list_pred_analysis_multi_ref_chains.count() == 1:
                self._view.ui.lbl_pred_analysis_multi_ref_chains.setText(
                    f"Select chain in protein structure {self._view.ui.lbl_pred_analysis_multi_prot_struct_1.text()}.",
                )
            else:
                self._view.ui.lbl_pred_analysis_multi_ref_chains.setText(
                    f"Select chains in protein structure {self._view.ui.lbl_pred_analysis_multi_prot_struct_1.text()}.",
                )
        self._view.wait_spinner.stop()

    # <editor-fold desc="Multimer prediction + analysis">
    def _init_multi_pred_analysis_page(self) -> None:
        """Clears the text boxes and sets the default values for the gui elements."""
        # <editor-fold desc="Prediction section">
        # clears everything
        self._view.ui.txt_pred_analysis_multi_prot_name.clear()
        self._view.ui.txt_pred_analysis_multi_prot_seq.clear()
        self._view.ui.list_pred_analysis_multi_prot_seq_overview.clear()
        for i in range(self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount() - 1, -1, -1):
            self._view.ui.table_pred_analysis_multi_prot_to_predict.removeRow(i)

        # sets up defaults: Prediction
        self._view.ui.btn_pred_analysis_multi_next.setEnabled(False)
        self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2.setEnabled(False)
        self._view.ui.lbl_pred_analysis_multi_prot_name_status.setText("")
        self._view.ui.lbl_pred_analysis_multi_prot_seq_status.setText("")

        # </editor-fold>

        # <editor-fold desc="Analysis section">
        self._view.ui.list_pred_analysis_multi_overview.clear()
        self._view.ui.btn_pred_analysis_multi_remove.hide()

        # </editor-fold>

        # self.multi_pred_analysis_show_protein_overview()

    def display_multimer_pred_analysis(self) -> None:
        """Displays the multimer prediction + analysis page."""
        self._view.ui.btn_pred_analysis_multi_go_analysis_setup.setText("Predict")
        self._view.ui.lbl_pred_analysis_multi_to_analysis_setup.setText("Protein Structure(s)")
        # checks internet connection
        if not tools.check_internet_connectivity():
            gui_utils.no_internet_dialog()
            return

        self._init_multi_pred_analysis_page()
        self._view.ui.table_pred_analysis_multi_prot_to_predict.clear()
        self._view.ui.table_pred_analysis_multi_prot_to_predict.setHorizontalHeaderItem(
            0,
            QtWidgets.QTableWidgetItem("Chain"),
        )
        self._view.ui.table_pred_analysis_multi_prot_to_predict.setHorizontalHeaderItem(
            1,
            QtWidgets.QTableWidgetItem("Sequence"),
        )
        self._view.ui.table_pred_analysis_multi_prot_to_predict.resizeColumnsToContents()
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_prot_to_predict,
            self._view.ui.table_pred_analysis_multi_prot_to_predict,
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_add,
        ]
        gui_elements_to_hide = [
            self._view.ui.checkbox_add_analysis,
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove,
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
            self._view.ui.lbl_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_multi_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self._view.ui.tabWidget_2.currentIndex() == 1:
            self._view.ui.tabWidget_2.setCurrentIndex(0)
        self._view.ui.tabWidget_2.setTabEnabled(1, False)
        self._view.ui.tabWidget_2.setTabEnabled(0, True)
        self._view.ui.table_pred_analysis_multi_prot_to_predict.setEnabled(True)

    def show_prediction_configuration(self) -> None:
        """Opens the prediction configuration dialog window."""
        config = dialog_advanced_prediction_configurations.DialogAdvancedPredictionConfigurations(
            self.prediction_configuration,
        )
        config.exec_()
        self.prediction_configuration.amber_force_field = config.prediction_config.amber_force_field
        self.prediction_configuration.templates = config.prediction_config.templates

    def check_if_prediction_and_analysis_should_be_done(self):
        if self._view.ui.checkbox_add_analysis.isChecked():
            self._view.ui.btn_pred_analysis_multi_go_analysis_setup.setText("Go")
            self._view.ui.lbl_pred_analysis_multi_to_analysis_setup.setText("To Analysis Setup")
        else:
            self._view.ui.btn_pred_analysis_multi_go_analysis_setup.setText("Predict")
            self._view.ui.lbl_pred_analysis_multi_to_analysis_setup.setText("Protein Structure(s)")

    # <editor-fold desc="Prediction section">
    def multi_pred_analysis_validate_protein_name(self) -> None:
        """Validates the input of the protein name in real-time."""
        if safeguard.Safeguard.check_if_value_is_in_table_v_header(
                self._view.ui.txt_pred_analysis_multi_prot_name.text(),
                self._view.ui.table_pred_analysis_multi_prot_to_predict,
        ):
            self._view.ui.lbl_pred_analysis_multi_prot_name_status.setText("Protein name already used.")
            self._view.ui.btn_pred_analysis_multi_next.setEnabled(False)
            styles.color_button_not_ready(self._view.ui.btn_pred_analysis_multi_next)
        else:
            self._view.ui.btn_pred_analysis_multi_next.setEnabled(True)
            tools.validate_protein_name(
                self._view.ui.txt_pred_analysis_multi_prot_name,
                self._view.ui.lbl_pred_analysis_multi_prot_name_status,
                self._view.ui.btn_pred_analysis_multi_next,
            )

    def multi_pred_analysis_validate_protein_sequence(self) -> None:
        """Validates the input of the protein sequence in real-time."""
        tools.validate_protein_sequence(
            self._view.ui.txt_pred_analysis_multi_prot_seq,
            self._view.ui.lbl_pred_analysis_multi_prot_seq_status,
            self._view.ui.btn_pred_analysis_multi_prot_seq_add,
        )

    def multi_pred_analysis_check_if_list_is_empty(self) -> None:
        """Checks if the list of sequences of the protein is empty."""
        if self._view.ui.list_pred_analysis_multi_prot_seq_overview.count() == 0:
            styles.color_button_not_ready(self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2)
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2.setEnabled(False)
        else:
            styles.color_button_ready(self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2)
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2.setEnabled(True)

    def multi_pred_analysis_add_sequence_to_list(self) -> None:
        """Adds the entered sequence to the sequences of the protein."""
        self._view.ui.list_pred_analysis_multi_prot_seq_overview.addItem(
            QtWidgets.QListWidgetItem(self._view.ui.txt_pred_analysis_multi_prot_seq.toPlainText()),
        )
        self.multi_pred_analysis_check_if_list_is_empty()

    def multi_pred_analysis_remove_sequence_to_list(self) -> None:
        """Removes the entered sequence from the sequences of the protein."""
        self._view.ui.list_pred_analysis_multi_prot_seq_overview.takeItem(
            self._view.ui.list_pred_analysis_multi_prot_seq_overview.currentRow(),
        )
        self.multi_pred_analysis_check_if_list_is_empty()
        self._view.ui.btn_pred_analysis_multi_prot_seq_overview_remove.setEnabled(False)

    def multi_pred_analysis_check_if_table_is_empty(self) -> None:
        """Checks if the list of proteins to predict is empty."""
        if self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount() == 0:
            styles.color_button_not_ready(self._view.ui.btn_pred_analysis_multi_go_analysis_setup)
            gui_elements_to_show = [
                self._view.ui.lbl_pred_analysis_multi_prot_to_predict,
                self._view.ui.table_pred_analysis_multi_prot_to_predict,
                self._view.ui.btn_pred_analysis_multi_prot_to_predict_add,
            ]
            gui_elements_to_hide = [
                self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove,
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
                self._view.ui.lbl_pred_analysis_multi_advanced_config,
                self._view.ui.btn_pred_analysis_multi_advanced_config,
                self._view.ui.btn_pred_analysis_multi_go_analysis_setup,
                self._view.ui.lbl_pred_analysis_multi_to_analysis_setup,
                self._view.ui.checkbox_add_analysis,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.btn_pred_analysis_multi_go_analysis_setup.setEnabled(False)
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove.setEnabled(False)
        else:
            styles.color_button_ready(self._view.ui.btn_pred_analysis_multi_go_analysis_setup)
            self._view.ui.btn_pred_analysis_multi_start.setEnabled(True)
            gui_elements_to_show = [
                self._view.ui.lbl_pred_analysis_multi_prot_to_predict,
                self._view.ui.table_pred_analysis_multi_prot_to_predict,
                self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove,
                self._view.ui.btn_pred_analysis_multi_prot_to_predict_add,
                self._view.ui.lbl_pred_analysis_multi_advanced_config,
                self._view.ui.btn_pred_analysis_multi_advanced_config,
                self._view.ui.btn_pred_analysis_multi_go_analysis_setup,
                self._view.ui.lbl_pred_analysis_multi_to_analysis_setup,
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
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove.setEnabled(False)

    def multi_pred_analysis_add_protein_to_predict(self) -> None:
        """Adds the proteins to the list of proteins to predict."""
        for i in range(self._view.ui.list_pred_analysis_multi_prot_seq_overview.count()):
            self._view.ui.table_pred_analysis_multi_prot_to_predict.setRowCount(
                self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount() + 1,
            )
            self._view.ui.table_pred_analysis_multi_prot_to_predict.insertRow(
                self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount() + 1,
            )
            tmp_chain_seq = (
                constants.chain_dict.get(i),
                self._view.ui.list_pred_analysis_multi_prot_seq_overview.item(i).text(),
            )
            self._view.ui.table_pred_analysis_multi_prot_to_predict.setItem(
                self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount() - 1,
                0,
                QtWidgets.QTableWidgetItem(tmp_chain_seq[0]),
            )
            self._view.ui.table_pred_analysis_multi_prot_to_predict.setItem(
                self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount() - 1,
                1,
                QtWidgets.QTableWidgetItem(tmp_chain_seq[1]),
            )
            name_item = QtWidgets.QTableWidgetItem(self._view.ui.txt_pred_analysis_multi_prot_name.text())
            self._view.ui.table_pred_analysis_multi_prot_to_predict.setVerticalHeaderItem(
                self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount() - 1,
                name_item,
            )
        self._view.ui.table_pred_analysis_multi_prot_to_predict.resizeColumnsToContents()
        self.multi_pred_analysis_check_if_table_is_empty()
        self._view.ui.checkbox_add_analysis.show()
        self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove.setEnabled(False)

    def multi_pred_analysis_remove_protein_to_predict(self) -> None:
        """Removes the selected protein from the list of proteins to predict."""
        if self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount() == 1:
            self._view.ui.table_pred_analysis_multi_prot_to_predict.removeRow(0)
        else:
            self._view.ui.table_pred_analysis_multi_prot_to_predict.removeRow(
                self._view.ui.table_pred_analysis_multi_prot_to_predict.currentRow(),
            )
            prot_name = self._view.ui.table_pred_analysis_multi_prot_to_predict.verticalHeaderItem(
                self._view.ui.table_pred_analysis_multi_prot_to_predict.currentRow(),
            ).text()
            for i in range(self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount()):
                if self._view.ui.table_pred_analysis_multi_prot_to_predict.verticalHeaderItem(
                        i).text() == prot_name:
                    self._view.ui.table_pred_analysis_multi_prot_to_predict.setItem(
                        i,
                        0,
                        QtWidgets.QTableWidgetItem(constants.chain_dict.get(i)),
                    )
        self.multi_pred_analysis_check_if_table_is_empty()
        self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove.setEnabled(False)

    def multi_pred_analysis_add(self) -> None:
        """Shows the gui elements for the protein name."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_prot_to_predict,
            self._view.ui.table_pred_analysis_multi_prot_to_predict,
            self._view.ui.lbl_pred_analysis_multi_prot_name,
            self._view.ui.txt_pred_analysis_multi_prot_name,
            self._view.ui.lbl_pred_analysis_multi_prot_name_status,
            self._view.ui.btn_pred_analysis_multi_back,
            self._view.ui.btn_pred_analysis_multi_next,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove,
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_add,
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
            self._view.ui.lbl_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_multi_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.enable_text_box(
            self._view.ui.txt_pred_analysis_multi_prot_name,
            self._view.ui.lbl_pred_analysis_multi_prot_name,
        )
        gui_utils.disable_text_box(
            self._view.ui.txt_pred_analysis_multi_prot_seq,
            self._view.ui.lbl_pred_analysis_multi_prot_seq,
        )
        self._view.ui.btn_pred_analysis_multi_next.setEnabled(False)
        self._view.ui.txt_pred_analysis_multi_prot_name.clear()
        styles.color_button_not_ready(self._view.ui.btn_pred_analysis_multi_next)
        if self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount() > 0:
            try:
                self._view.ui.table_pred_analysis_multi_prot_to_predict.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection on Local Multimer Prediction in overview table.")

    def multi_pred_analysis_back(self) -> None:
        """Hides the gui elements for the protein name."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_prot_to_predict,
            self._view.ui.table_pred_analysis_multi_prot_to_predict,
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_add,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove,
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
            self._view.ui.lbl_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_multi_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self.multi_pred_analysis_check_if_table_is_empty()

    def multi_pred_analysis_next(self) -> None:
        """Shows the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_prot_to_predict,
            self._view.ui.table_pred_analysis_multi_prot_to_predict,
            self._view.ui.lbl_pred_analysis_multi_prot_name,
            self._view.ui.txt_pred_analysis_multi_prot_name,
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
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove,
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_add,
            self._view.ui.lbl_pred_analysis_multi_prot_name_status,
            self._view.ui.btn_pred_analysis_multi_back,
            self._view.ui.btn_pred_analysis_multi_next,
            self._view.ui.lbl_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_multi_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.enable_text_box(
            self._view.ui.txt_pred_analysis_multi_prot_seq,
            self._view.ui.lbl_pred_analysis_multi_prot_seq,
        )
        gui_utils.disable_text_box(
            self._view.ui.txt_pred_analysis_multi_prot_name,
            self._view.ui.lbl_pred_analysis_multi_prot_name,
        )
        self._view.ui.txt_pred_analysis_multi_prot_seq.clear()
        self._view.ui.list_pred_analysis_multi_prot_seq_overview.clear()
        self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2.setEnabled(False)
        self._view.ui.btn_pred_analysis_multi_prot_seq_overview_remove.setEnabled(False)
        styles.color_button_not_ready(self._view.ui.btn_pred_analysis_multi_prot_to_predict_add_2)

    def multi_pred_analysis_back_2(self) -> None:
        """Hides the gui elements for the protein sequence."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_prot_to_predict,
            self._view.ui.table_pred_analysis_multi_prot_to_predict,
            self._view.ui.lbl_pred_analysis_multi_prot_name_status,
            self._view.ui.btn_pred_analysis_multi_back,
            self._view.ui.btn_pred_analysis_multi_next,
            self._view.ui.lbl_pred_analysis_multi_prot_name,
            self._view.ui.txt_pred_analysis_multi_prot_name,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove,
            self._view.ui.btn_pred_analysis_multi_prot_to_predict_add,
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
            self._view.ui.lbl_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_advanced_config,
            self._view.ui.btn_pred_analysis_multi_go_analysis_setup,
            self._view.ui.lbl_pred_analysis_multi_to_analysis_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        gui_utils.enable_text_box(
            self._view.ui.txt_pred_analysis_multi_prot_name,
            self._view.ui.lbl_pred_analysis_multi_prot_name,
        )
        gui_utils.disable_text_box(
            self._view.ui.txt_pred_analysis_multi_prot_seq,
            self._view.ui.lbl_pred_analysis_multi_prot_seq,
        )

    def multi_pred_analysis_prot_seq_overview_item_changed(self) -> None:
        """Enables the remove button of the list of sequences of the protein."""
        self._view.ui.btn_pred_analysis_multi_prot_seq_overview_remove.setEnabled(True)

    def multi_pred_analysis_prot_to_predict_item_changed(self) -> None:
        """Enables the remove button of the list of proteins to predict."""
        self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove.setEnabled(True)

    # </editor-fold>

    def switch_multimer_pred_analysis_tab(self) -> None:
        """Switches the tabs from prediction to analysis and vice versa."""
        if self._view.ui.btn_pred_analysis_multi_go_analysis_setup.text() == "Predict":
            tmp_prediction_runs: list[
                prediction_protein_info.PredictionProteinInfo
            ] = prediction_util.get_prediction_name_and_seq_from_table(
                self._view.ui.table_pred_analysis_multi_prot_to_predict)
            self._view.close()
            self.job_input.emit(("job_input", tmp_prediction_runs, self.prediction_configuration, False))
        else:
            if self._view.ui.tabWidget_2.currentIndex() == 0:
                # goes from prediction to analysis
                self._view.ui.tabWidget_2.setCurrentIndex(1)
                gui_elements_to_show = [
                    self._view.ui.lbl_pred_analysis_multi_overview,
                    self._view.ui.list_pred_analysis_multi_overview,
                    self._view.ui.btn_pred_analysis_multi_add,
                    self._view.ui.btn_pred_analysis_multi_back_pred_setup,
                ]
                gui_elements_to_hide = [
                    self._view.ui.btn_pred_analysis_multi_remove,
                    self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
                    self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
                    self._view.ui.lbl_analysis_batch_vs_3,
                    self._view.ui.lbl_pred_analysis_multi_ref_chains,
                    self._view.ui.list_pred_analysis_multi_ref_chains,
                    self._view.ui.btn_pred_analysis_multi_back_4,
                    self._view.ui.btn_pred_analysis_multi_next_3,
                    self._view.ui.box_pred_analysis_multi_prot_struct_1,
                    self._view.ui.box_pred_analysis_multi_prot_struct_2,
                    self._view.ui.btn_pred_analysis_multi_back_3,
                    self._view.ui.btn_pred_analysis_multi_next_2,
                    self._view.ui.lbl_pred_analysis_multi_model_chains,
                    self._view.ui.list_pred_analysis_multi_model_chains,
                    self._view.ui.btn_pred_analysis_multi_back_5,
                    self._view.ui.btn_pred_analysis_multi_next_4,
                    self._view.ui.cb_pred_analysis_multi_images,
                    self._view.ui.btn_pred_analysis_multi_start,
                ]
                gui_utils.show_gui_elements(gui_elements_to_show)
                gui_utils.hide_gui_elements(gui_elements_to_hide)
                self._view.ui.tabWidget_2.setTabEnabled(1, True)
                self._view.ui.tabWidget_2.setTabEnabled(0, False)
                if self._view.ui.list_pred_analysis_multi_overview.count() > 0:
                    self._view.ui.btn_pred_analysis_multi_remove.show()
                    self._view.ui.btn_pred_analysis_multi_remove.setEnabled(False)
                    self._view.ui.btn_pred_analysis_multi_start.show()
                    self._view.ui.cb_pred_analysis_multi_images.hide()
                    styles.color_button_ready(self._view.ui.btn_pred_analysis_multi_start)
            else:
                # goes from analysis to prediction
                self._view.ui.tabWidget_2.setCurrentIndex(0)
                if self._view.ui.list_pred_analysis_multi_overview.count() > 0:
                    gui_elements_to_show = [
                        self._view.ui.lbl_pred_analysis_multi_prot_to_predict,
                        self._view.ui.table_pred_analysis_multi_prot_to_predict,
                        self._view.ui.btn_pred_analysis_multi_go_analysis_setup,
                        self._view.ui.lbl_pred_analysis_multi_to_analysis_setup,
                    ]
                    gui_elements_to_hide = [
                        self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove,
                        self._view.ui.btn_pred_analysis_multi_prot_to_predict_add,
                        self._view.ui.lbl_pred_analysis_multi_advanced_config,
                        self._view.ui.btn_pred_analysis_multi_advanced_config,
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
                else:
                    gui_elements_to_show = [
                        self._view.ui.lbl_pred_analysis_multi_prot_to_predict,
                        self._view.ui.table_pred_analysis_multi_prot_to_predict,
                        self._view.ui.btn_pred_analysis_multi_prot_to_predict_remove,
                        self._view.ui.btn_pred_analysis_multi_prot_to_predict_add,
                        self._view.ui.lbl_pred_analysis_multi_advanced_config,
                        self._view.ui.btn_pred_analysis_multi_advanced_config,
                        self._view.ui.btn_pred_analysis_multi_go_analysis_setup,
                        self._view.ui.lbl_pred_analysis_multi_to_analysis_setup,
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
                self._view.ui.tabWidget_2.setTabEnabled(0, True)
                self._view.ui.tabWidget_2.setTabEnabled(1, False)

    # <editor-fold desc="Analysis section">
    def multi_pred_analysis_structure_analysis_add(self) -> None:
        """Shows the gui elements to choose the two proteins."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_overview,
            self._view.ui.list_pred_analysis_multi_overview,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
            self._view.ui.box_pred_analysis_multi_prot_struct_1,
            self._view.ui.lbl_analysis_batch_vs_3,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
            self._view.ui.box_pred_analysis_multi_prot_struct_2,
            self._view.ui.btn_pred_analysis_multi_back_3,
            self._view.ui.btn_pred_analysis_multi_next_2,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_remove,
            self._view.ui.btn_pred_analysis_multi_add,
            self._view.ui.lbl_pred_analysis_multi_ref_chains,
            self._view.ui.list_pred_analysis_multi_ref_chains,
            self._view.ui.btn_pred_analysis_multi_back_4,
            self._view.ui.btn_pred_analysis_multi_next_3,
            self._view.ui.lbl_pred_analysis_multi_model_chains,
            self._view.ui.list_pred_analysis_multi_model_chains,
            self._view.ui.btn_pred_analysis_multi_back_5,
            self._view.ui.btn_pred_analysis_multi_next_4,
            self._view.ui.cb_pred_analysis_multi_images,
            self._view.ui.btn_pred_analysis_multi_start,
            self._view.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.lbl_pred_analysis_multi_prot_struct_1.clear()
        self._view.ui.lbl_pred_analysis_multi_prot_struct_2.clear()
        self._view.ui.lbl_pred_analysis_multi_prot_struct_1.setText("Protein structure 1")
        self._view.ui.lbl_pred_analysis_multi_prot_struct_2.setText("Protein structure 2")
        self.fill_multi_pred_analysis_protein_boxes()
        if self._view.ui.list_pred_analysis_multi_overview.count() > 0:
            try:
                self._view.ui.list_pred_analysis_multi_overview.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection in struction analysis overview.")

    def multi_pred_analysis_structure_analysis_back_3(self) -> None:
        """Hides the gui elements to choose the two proteins."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_overview,
            self._view.ui.list_pred_analysis_multi_overview,
            self._view.ui.btn_pred_analysis_multi_add,
            self._view.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_remove,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_3,
            self._view.ui.lbl_pred_analysis_multi_ref_chains,
            self._view.ui.list_pred_analysis_multi_ref_chains,
            self._view.ui.btn_pred_analysis_multi_back_4,
            self._view.ui.btn_pred_analysis_multi_next_3,
            self._view.ui.box_pred_analysis_multi_prot_struct_1,
            self._view.ui.box_pred_analysis_multi_prot_struct_2,
            self._view.ui.btn_pred_analysis_multi_back_3,
            self._view.ui.btn_pred_analysis_multi_next_2,
            self._view.ui.lbl_pred_analysis_multi_model_chains,
            self._view.ui.list_pred_analysis_multi_model_chains,
            self._view.ui.btn_pred_analysis_multi_back_5,
            self._view.ui.btn_pred_analysis_multi_next_4,
            self._view.ui.cb_pred_analysis_multi_images,
            self._view.ui.btn_pred_analysis_multi_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self._view.ui.list_pred_analysis_multi_overview.count() > 0:
            self._view.ui.btn_pred_analysis_multi_remove.show()
            self._view.ui.btn_pred_analysis_multi_remove.setEnabled(False)
            self._view.ui.btn_pred_analysis_multi_start.show()
            self._view.ui.cb_pred_analysis_multi_images.hide()
            styles.color_button_ready(self._view.ui.btn_pred_analysis_multi_start)

    def multi_pred_analysis_structure_analysis_next_3(self) -> None:
        """Shows the gui elements to select the chains in protein 2."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_overview,
            self._view.ui.list_pred_analysis_multi_overview,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_3,
            self._view.ui.lbl_pred_analysis_multi_ref_chains,
            self._view.ui.list_pred_analysis_multi_ref_chains,
            self._view.ui.lbl_pred_analysis_multi_model_chains,
            self._view.ui.list_pred_analysis_multi_model_chains,
            self._view.ui.btn_pred_analysis_multi_back_5,
            self._view.ui.btn_pred_analysis_multi_next_4,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_remove,
            self._view.ui.btn_pred_analysis_multi_add,
            self._view.ui.box_pred_analysis_multi_prot_struct_1,
            self._view.ui.box_pred_analysis_multi_prot_struct_2,
            self._view.ui.btn_pred_analysis_multi_back_3,
            self._view.ui.btn_pred_analysis_multi_next_2,
            self._view.ui.btn_pred_analysis_multi_back_4,
            self._view.ui.btn_pred_analysis_multi_next_3,
            self._view.ui.cb_pred_analysis_multi_images,
            self._view.ui.btn_pred_analysis_multi_start,
            self._view.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.list_pred_analysis_multi_model_chains.clear()
        self._view.ui.list_pred_analysis_multi_ref_chains.setEnabled(False)
        self._view.ui.btn_pred_analysis_multi_next_4.setEnabled(False)

        for i in range(self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount()):
            if (
                    self._view.ui.table_pred_analysis_multi_prot_to_predict.verticalHeaderItem(i).text()
                    == self._view.ui.box_pred_analysis_multi_prot_struct_2.currentText()
            ):
                self._view.ui.list_pred_analysis_multi_model_chains.addItem(
                    self._view.ui.table_pred_analysis_multi_prot_to_predict.item(i, 0).text(),
                )
        if self._view.ui.list_pred_analysis_multi_model_chains.count() == 0:
            tmp_protein = self._interface_manager.get_current_project().search_protein(
                self._view.ui.box_pred_analysis_multi_prot_struct_2.currentText(),
            )
            for tmp_chain in tmp_protein.chains:
                if tmp_chain.chain_type == "protein_chain":
                    self._view.ui.list_pred_analysis_multi_model_chains.addItem(tmp_chain.chain_letter)
        if len(self._view.ui.list_pred_analysis_multi_ref_chains.selectedItems()) == 1:
            self._view.ui.lbl_pred_analysis_multi_model_chains.setText(
                f"Select 1 chain in protein structure {self._view.ui.lbl_pred_analysis_multi_prot_struct_2.text()}.",
            )
        else:
            self._view.ui.lbl_pred_analysis_multi_model_chains.setText(
                f"Select {len(self._view.ui.list_pred_analysis_multi_ref_chains.selectedItems())} chains in "
                f"protein structure {self._view.ui.lbl_pred_analysis_multi_prot_struct_2.text()}.",
            )

    def multi_pred_analysis_structure_analysis_back_4(self) -> None:
        """Hides the gui elements to select the chains in protein 1."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_overview,
            self._view.ui.list_pred_analysis_multi_overview,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
            self._view.ui.box_pred_analysis_multi_prot_struct_1,
            self._view.ui.lbl_analysis_batch_vs_3,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
            self._view.ui.box_pred_analysis_multi_prot_struct_2,
            self._view.ui.btn_pred_analysis_multi_back_3,
            self._view.ui.btn_pred_analysis_multi_next_2,
            self._view.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_remove,
            self._view.ui.btn_pred_analysis_multi_add,
            self._view.ui.lbl_pred_analysis_multi_ref_chains,
            self._view.ui.list_pred_analysis_multi_ref_chains,
            self._view.ui.btn_pred_analysis_multi_back_4,
            self._view.ui.btn_pred_analysis_multi_next_3,
            self._view.ui.lbl_pred_analysis_multi_model_chains,
            self._view.ui.list_pred_analysis_multi_model_chains,
            self._view.ui.btn_pred_analysis_multi_back_5,
            self._view.ui.btn_pred_analysis_multi_next_4,
            self._view.ui.cb_pred_analysis_multi_images,
            self._view.ui.btn_pred_analysis_multi_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.lbl_pred_analysis_multi_prot_struct_1.setText("Protein structure 1")
        self._view.ui.lbl_pred_analysis_multi_prot_struct_2.setText("Protein structure 2")

    def multi_pred_analysis_structure_analysis_next_4(self) -> None:
        """Adds the protein pair to the list of protein pairs to analyze."""
        gui_elements_to_show = [
            self._view.ui.btn_pred_analysis_multi_remove,
            self._view.ui.btn_pred_analysis_multi_add,
            self._view.ui.lbl_pred_analysis_multi_overview,
            self._view.ui.list_pred_analysis_multi_overview,
            self._view.ui.btn_pred_analysis_multi_start,
            self._view.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self._view.ui.box_pred_analysis_multi_prot_struct_1,
            self._view.ui.box_pred_analysis_multi_prot_struct_2,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_3,
            self._view.ui.lbl_pred_analysis_multi_ref_chains,
            self._view.ui.list_pred_analysis_multi_ref_chains,
            self._view.ui.lbl_pred_analysis_multi_model_chains,
            self._view.ui.list_pred_analysis_multi_model_chains,
            self._view.ui.btn_pred_analysis_multi_back_3,
            self._view.ui.btn_pred_analysis_multi_next_2,
            self._view.ui.btn_pred_analysis_multi_back_4,
            self._view.ui.btn_pred_analysis_multi_next_3,
            self._view.ui.btn_pred_analysis_multi_back_5,
            self._view.ui.btn_pred_analysis_multi_next_4,
            self._view.ui.cb_pred_analysis_multi_images,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        prot_1_name = self._view.ui.lbl_pred_analysis_multi_prot_struct_1.text()
        prot_1_chains = []
        for chain in self._view.ui.list_pred_analysis_multi_ref_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self._view.ui.lbl_pred_analysis_multi_prot_struct_2.text()
        prot_2_chains = []
        for chain in self._view.ui.list_pred_analysis_multi_model_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        item = QtWidgets.QListWidgetItem(analysis_name)
        self._view.ui.list_pred_analysis_multi_overview.addItem(item)
        self._view.ui.btn_pred_analysis_multi_remove.setEnabled(False)
        styles.color_button_ready(self._view.ui.btn_pred_analysis_multi_start)

    def multi_pred_analysis_structure_analysis_back_5(self) -> None:
        """Hides the gui elements to select the chains in protein 2."""
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_overview,
            self._view.ui.list_pred_analysis_multi_overview,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_3,
            self._view.ui.lbl_pred_analysis_multi_ref_chains,
            self._view.ui.list_pred_analysis_multi_ref_chains,
            self._view.ui.btn_pred_analysis_multi_back_4,
            self._view.ui.btn_pred_analysis_multi_next_3,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_remove,
            self._view.ui.btn_pred_analysis_multi_add,
            self._view.ui.box_pred_analysis_multi_prot_struct_1,
            self._view.ui.box_pred_analysis_multi_prot_struct_2,
            self._view.ui.btn_pred_analysis_multi_back_3,
            self._view.ui.btn_pred_analysis_multi_next_2,
            self._view.ui.btn_pred_analysis_multi_back_5,
            self._view.ui.btn_pred_analysis_multi_next_4,
            self._view.ui.cb_pred_analysis_multi_images,
            self._view.ui.btn_pred_analysis_multi_start,
            self._view.ui.lbl_pred_analysis_multi_model_chains,
            self._view.ui.list_pred_analysis_multi_model_chains,
            self._view.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.list_pred_analysis_multi_ref_chains.setEnabled(True)

        # tmp_protein = self._interface_manager.get_current_project().search_protein(self._view.ui.box_pred_analysis_multi_prot_struct_2.currentText())
        # for tmp_chain in tmp_protein.chains:
        #     if tmp_chain.chain_type == "protein_chain":
        #         self._view.ui.list_pred_analysis_multi_ref_chains.addItem(tmp_chain.chain_letter)

    def multi_pred_analysis_structure_analysis_overview_clicked(self) -> None:
        """Enables the remove button."""
        self._view.ui.btn_pred_analysis_multi_remove.setEnabled(True)

    def fill_multi_pred_analysis_protein_boxes(self) -> None:
        """Fills the combo boxes with the protein names."""
        protein_names = []
        for i in range(self._view.ui.table_pred_analysis_multi_prot_to_predict.rowCount()):
            protein_names.append(
                self._view.ui.table_pred_analysis_multi_prot_to_predict.verticalHeaderItem(i).text())
        for tmp_protein in self._interface_manager.get_current_project().proteins:
            protein_names.append(tmp_protein.get_molecule_object())
        protein_names.insert(0, "")
        protein_names = list(set(protein_names))
        self._view.ui.box_pred_analysis_multi_prot_struct_1.clear()
        self._view.ui.box_pred_analysis_multi_prot_struct_2.clear()
        gui_utils.fill_combo_box(self._view.ui.box_pred_analysis_multi_prot_struct_1, protein_names)
        gui_utils.fill_combo_box(self._view.ui.box_pred_analysis_multi_prot_struct_2, protein_names)

    def remove_multi_pred_analysis_analysis_run(self) -> None:
        """Removes the selected protein pair from the list of protein pairs to analyze."""
        self._view.ui.list_pred_analysis_multi_overview.takeItem(
            self._view.ui.list_pred_analysis_multi_overview.currentRow(),
        )
        gui_elements_to_show = [
            self._view.ui.lbl_pred_analysis_multi_overview,
            self._view.ui.list_pred_analysis_multi_overview,
            self._view.ui.btn_pred_analysis_multi_add,
            self._view.ui.btn_pred_analysis_multi_back_pred_setup,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_pred_analysis_multi_remove,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_1,
            self._view.ui.lbl_pred_analysis_multi_prot_struct_2,
            self._view.ui.lbl_analysis_batch_vs_3,
            self._view.ui.lbl_pred_analysis_multi_ref_chains,
            self._view.ui.list_pred_analysis_multi_ref_chains,
            self._view.ui.btn_pred_analysis_multi_back_4,
            self._view.ui.btn_pred_analysis_multi_next_3,
            self._view.ui.box_pred_analysis_multi_prot_struct_1,
            self._view.ui.box_pred_analysis_multi_prot_struct_2,
            self._view.ui.btn_pred_analysis_multi_back_3,
            self._view.ui.btn_pred_analysis_multi_next_2,
            self._view.ui.lbl_pred_analysis_multi_model_chains,
            self._view.ui.list_pred_analysis_multi_model_chains,
            self._view.ui.btn_pred_analysis_multi_back_5,
            self._view.ui.btn_pred_analysis_multi_next_4,
            
            self._view.ui.cb_pred_analysis_multi_images,
            self._view.ui.btn_pred_analysis_multi_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self._view.ui.list_pred_analysis_multi_overview.count() > 0:
            self._view.ui.btn_pred_analysis_multi_remove.show()
            self._view.ui.btn_pred_analysis_multi_remove.setEnabled(False)
            self._view.ui.btn_pred_analysis_multi_start.show()
            self._view.ui.cb_pred_analysis_multi_images.hide()
            styles.color_button_ready(self._view.ui.btn_pred_analysis_multi_start)
        # if self._view.ui.list_pred_analysis_multi_overview.count() == 0:
        #
        #     self._view.ui.btn_pred_analysis_multi_back_pred_setup.show()
        #     self._view.ui.btn_pred_analysis_multi_remove.hide()

    def check_multi_pred_analysis_if_same_no_of_chains_selected(self) -> None:
        """Checks if the same number of chains were selected."""
        self._view.ui.btn_pred_analysis_multi_next_4.setEnabled(False)
        styles.color_button_not_ready(self._view.ui.btn_pred_analysis_multi_next_4)
        if self.no_of_selected_chains == len(self._view.ui.list_pred_analysis_multi_model_chains.selectedItems()):
            styles.color_button_ready(self._view.ui.btn_pred_analysis_multi_next_4)
            self._view.ui.btn_pred_analysis_multi_next_4.setEnabled(True)

        prot_1_name = self._view.ui.lbl_pred_analysis_multi_prot_struct_1.text()
        prot_1_chains = []
        for chain in self._view.ui.list_pred_analysis_multi_ref_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self._view.ui.lbl_pred_analysis_multi_prot_struct_2.text()
        prot_2_chains = []
        for chain in self._view.ui.list_pred_analysis_multi_model_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        for tmp_row in range(self._view.ui.list_pred_analysis_multi_overview.count()):
            if analysis_name == self._view.ui.list_pred_analysis_multi_overview.item(tmp_row).text():
                self._view.ui.btn_pred_analysis_multi_next_4.setEnabled(False)
                styles.color_button_not_ready(self._view.ui.btn_pred_analysis_multi_next_4)
                return

    def check_multi_pred_analysis_if_prot_structs_are_filled(self) -> None:
        """Checks if two proteins were selected."""
        prot_1 = self._view.ui.box_pred_analysis_multi_prot_struct_1.itemText(
            self._view.ui.box_pred_analysis_multi_prot_struct_1.currentIndex(),
        )
        prot_2 = self._view.ui.box_pred_analysis_multi_prot_struct_2.itemText(
            self._view.ui.box_pred_analysis_multi_prot_struct_2.currentIndex(),
        )
        if prot_1 != "" and prot_2 != "":
            self._view.ui.btn_pred_analysis_multi_next_2.setEnabled(True)
        else:
            self._view.ui.btn_pred_analysis_multi_next_2.setEnabled(False)

    def count_multi_pred_analysis_selected_chains_for_prot_struct_1(self) -> None:
        """Counts the number of chains in protein 1."""
        self.no_of_selected_chains = len(self._view.ui.list_pred_analysis_multi_ref_chains.selectedItems())
        if self.no_of_selected_chains > 0:
            self._view.ui.btn_pred_analysis_multi_next_3.setEnabled(True)
        else:
            self._view.ui.btn_pred_analysis_multi_next_3.setEnabled(False)

    # </editor-fold>
    # </editor-fold>
