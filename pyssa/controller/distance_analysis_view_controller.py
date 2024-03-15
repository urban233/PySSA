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
from pyssa.gui.ui.views import distance_analysis_view
from pyssa.internal.data_structures.data_classes import prediction_protein_info, prediction_configuration
from pyssa.internal.thread import tasks
from pyssa.internal.thread.async_pyssa import util_async
from pyssa.io_pyssa import safeguard
from pyssa.presenter import main_presenter_async
from pyssa.util import gui_utils, tools, constants, exit_codes, prediction_util


class DistanceAnalysisViewController(QtCore.QObject):

    job_input = pyqtSignal(tuple)

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager"):
        super().__init__()
        self._interface_manager = the_interface_manager
        self._view: "distance_analysis_view.DistanceAnalysisView" = the_interface_manager.get_distance_analysis_view()
        self.prediction_configuration = prediction_configuration.PredictionConfiguration(True, "pdb70")

        self._connect_all_ui_elements_to_slot_functions()
        self.display_distance_analysis()

    def open_help(self, a_page_name: str):
        """Opens the pyssa documentation window if it's not already open.

        Args:
            a_page_name (str): a name of a documentation page to display
        """
        self._interface_manager.status_bar_manager.show_temporary_message("Opening help center ...")
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
            self._interface_manager.status_bar_manager.show_temporary_message("Opening help center finished.")

    def _open_help_for_dialog(self):
        self.open_help("help/protein_structure_analysis/distance_analysis/")

    def _connect_all_ui_elements_to_slot_functions(self):
        self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)
        self._view.ui.btn_distance_analysis_add.clicked.connect(self.structure_analysis_add)
        self._view.ui.btn_distance_analysis_remove.clicked.connect(self.remove_analysis_run)
        self._view.ui.btn_distance_analysis_back.clicked.connect(self.structure_analysis_back)
        self._view.ui.btn_distance_analysis_next.clicked.connect(self.structure_analysis_next)
        self._view.ui.btn_distance_analysis_back_2.clicked.connect(self.structure_analysis_back_2)
        self._view.ui.btn_distance_analysis_next_2.clicked.connect(self.structure_analysis_next_2)
        self._view.ui.btn_distance_analysis_back_3.clicked.connect(self.structure_analysis_back_3)
        self._view.ui.btn_distance_analysis_next_3.clicked.connect(self.structure_analysis_next_3)
        self._view.ui.combobox_distance_analysis_prot_struct_1.currentIndexChanged.connect(
            self.check_if_prot_structs_are_filled_batch,
        )
        self._view.ui.combobox_distance_analysis_prot_struct_2.currentIndexChanged.connect(
            self.check_if_prot_structs_are_filled_batch,
        )
        self._view.ui.list_distance_analysis_prot_1_chains.itemSelectionChanged.connect(
            self.count_batch_selected_chains_for_prot_struct_1,
        )
        self._view.ui.list_distance_analysis_prot_2_chains.itemSelectionChanged.connect(
            self.check_if_same_no_of_chains_selected_batch,
        )
        self._view.ui.btn_distance_analysis_start.clicked.connect(self.start_process_batch)
        self._view.ui.list_distance_analysis_overview.clicked.connect(self.structure_analysis_overview_clicked)

    def display_distance_analysis(self) -> None:
        """Displays the job analysis work area."""
        self._view.ui.list_distance_analysis_prot_1_chains.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self._view.ui.list_distance_analysis_prot_2_chains.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        gui_elements_to_show = [
            self._view.ui.btn_distance_analysis_add,
            self._view.ui.lbl_distance_analysis_overview,
            self._view.ui.list_distance_analysis_overview,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_distance_analysis_remove,
            self._view.ui.combobox_distance_analysis_prot_struct_1,
            self._view.ui.combobox_distance_analysis_prot_struct_2,
            self._view.ui.lbl_distance_analysis_prot_struct_1,
            self._view.ui.lbl_distance_analysis_prot_struct_2,
            self._view.ui.lbl_distance_analysis_vs,
            self._view.ui.lbl_distance_analysis_prot_1_chains,
            self._view.ui.list_distance_analysis_prot_1_chains,
            self._view.ui.lbl_distance_analysis_prot_2_chains,
            self._view.ui.list_distance_analysis_prot_2_chains,
            self._view.ui.btn_distance_analysis_back,
            self._view.ui.btn_distance_analysis_next,
            self._view.ui.btn_distance_analysis_back_2,
            self._view.ui.btn_distance_analysis_next_2,
            self._view.ui.btn_distance_analysis_back_3,
            self._view.ui.btn_distance_analysis_next_3,
            self._view.ui.btn_distance_analysis_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)

        self._view.ui.list_distance_analysis_overview.clear()

    # <editor-fold desc="Logic">
    def start_process_batch(self):
        tmp_raw_analysis_run_names: list = []
        for row_no in range(self._view.ui.list_distance_analysis_overview.count()):
            tmp_raw_analysis_run_names.append(self._view.ui.list_distance_analysis_overview.item(row_no).text())
        self.job_input.emit(("job_input", tmp_raw_analysis_run_names, False))
        self._view.close()

    def structure_analysis_next(self) -> None:
        """Shows the gui elements to select the chains in protein 1."""
        self._view.wait_spinner.start()
        self._active_task = tasks.Task(
            target=main_presenter_async.check_chains_for_analysis,
            args=(
                self._view.ui.combobox_distance_analysis_prot_struct_1.currentText(),
                self._view.ui.combobox_distance_analysis_prot_struct_2.currentText(),
                self._interface_manager.get_current_project(),
            ),
            post_func=self.__await_structure_analysis_next,
        )
        self._active_task.start()

    def __await_structure_analysis_next(self, result: tuple) -> None:
        _, tmp_is_only_one_chain, tmp_analysis_run_name, tmp_protein_1, tmp_protein_2 = result

        if tmp_is_only_one_chain:
            gui_elements_to_show = [
                self._view.ui.btn_distance_analysis_remove,
                self._view.ui.btn_distance_analysis_add,
                self._view.ui.lbl_distance_analysis_overview,
                self._view.ui.list_distance_analysis_overview,
                
                
                self._view.ui.btn_distance_analysis_start,
            ]
            gui_elements_to_hide = [
                self._view.ui.combobox_distance_analysis_prot_struct_1,
                self._view.ui.combobox_distance_analysis_prot_struct_2,
                self._view.ui.lbl_distance_analysis_prot_struct_1,
                self._view.ui.lbl_distance_analysis_prot_struct_2,
                self._view.ui.lbl_distance_analysis_vs,
                self._view.ui.lbl_distance_analysis_prot_1_chains,
                self._view.ui.list_distance_analysis_prot_1_chains,
                self._view.ui.lbl_distance_analysis_prot_2_chains,
                self._view.ui.list_distance_analysis_prot_2_chains,
                self._view.ui.btn_distance_analysis_back,
                self._view.ui.btn_distance_analysis_next,
                self._view.ui.btn_distance_analysis_back_2,
                self._view.ui.btn_distance_analysis_next_2,
                self._view.ui.btn_distance_analysis_back_3,
                self._view.ui.btn_distance_analysis_next_3,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            item = QtWidgets.QListWidgetItem(tmp_analysis_run_name)
            self._view.ui.list_distance_analysis_overview.addItem(item)
            self._view.ui.btn_distance_analysis_remove.setEnabled(False)
        else:
            gui_elements_to_show = [
                self._view.ui.lbl_distance_analysis_overview,
                self._view.ui.list_distance_analysis_overview,
                self._view.ui.lbl_distance_analysis_prot_struct_1,
                self._view.ui.lbl_distance_analysis_prot_struct_2,
                self._view.ui.lbl_distance_analysis_vs,
                self._view.ui.lbl_distance_analysis_prot_1_chains,
                self._view.ui.list_distance_analysis_prot_1_chains,
                self._view.ui.btn_distance_analysis_back_2,
                self._view.ui.btn_distance_analysis_next_2,
            ]
            gui_elements_to_hide = [
                self._view.ui.btn_distance_analysis_remove,
                self._view.ui.btn_distance_analysis_add,
                self._view.ui.combobox_distance_analysis_prot_struct_1,
                self._view.ui.combobox_distance_analysis_prot_struct_2,
                self._view.ui.btn_distance_analysis_back,
                self._view.ui.btn_distance_analysis_next,
                self._view.ui.lbl_distance_analysis_prot_2_chains,
                self._view.ui.list_distance_analysis_prot_2_chains,
                self._view.ui.btn_distance_analysis_back_3,
                self._view.ui.btn_distance_analysis_next_3,
                
                
                self._view.ui.btn_distance_analysis_start,
            ]
            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.lbl_distance_analysis_prot_struct_1.setText(
                self._view.ui.combobox_distance_analysis_prot_struct_1.currentText(),
            )
            self._view.ui.lbl_distance_analysis_prot_struct_2.setText(
                self._view.ui.combobox_distance_analysis_prot_struct_2.currentText(),
            )
            self._view.ui.list_distance_analysis_prot_1_chains.clear()
            self._view.ui.btn_distance_analysis_next_2.setEnabled(False)
            self._view.ui.list_distance_analysis_prot_1_chains.setEnabled(True)

            for tmp_chain in tmp_protein_1.chains:
                if tmp_chain.chain_type == "protein_chain":
                    self._view.ui.list_distance_analysis_prot_1_chains.addItem(tmp_chain.chain_letter)

            if self._view.ui.list_distance_analysis_prot_1_chains.count() == 1:
                self._view.ui.lbl_distance_analysis_prot_1_chains.setText(
                    f"Select chain in protein structure {self._view.ui.lbl_distance_analysis_prot_struct_1.text()}.",
                )
            else:
                self._view.ui.lbl_distance_analysis_prot_1_chains.setText(
                    f"Select chains in protein structure {self._view.ui.lbl_distance_analysis_prot_struct_1.text()}.",
                )
        self._view.wait_spinner.stop()

    # </editor-fold>

    # <editor-fold desc="GUI">
    def structure_analysis_add(self) -> None:
        """Shows the gui elements to choose the two proteins."""
        gui_elements_to_show = [
            self._view.ui.lbl_distance_analysis_overview,
            self._view.ui.list_distance_analysis_overview,
            self._view.ui.lbl_distance_analysis_prot_struct_1,
            self._view.ui.combobox_distance_analysis_prot_struct_1,
            self._view.ui.lbl_distance_analysis_vs,
            self._view.ui.lbl_distance_analysis_prot_struct_2,
            self._view.ui.combobox_distance_analysis_prot_struct_2,
            self._view.ui.btn_distance_analysis_back,
            self._view.ui.btn_distance_analysis_next,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_distance_analysis_remove,
            self._view.ui.btn_distance_analysis_add,
            self._view.ui.lbl_distance_analysis_prot_1_chains,
            self._view.ui.list_distance_analysis_prot_1_chains,
            self._view.ui.btn_distance_analysis_back_2,
            self._view.ui.btn_distance_analysis_next_2,
            self._view.ui.lbl_distance_analysis_prot_2_chains,
            self._view.ui.list_distance_analysis_prot_2_chains,
            self._view.ui.btn_distance_analysis_back_3,
            self._view.ui.btn_distance_analysis_next_3,
            
            
            self._view.ui.btn_distance_analysis_start,
        ]

        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.lbl_distance_analysis_prot_struct_1.clear()
        self._view.ui.lbl_distance_analysis_prot_struct_2.clear()
        self._view.ui.lbl_distance_analysis_prot_struct_1.setText("Protein structure 1")
        self._view.ui.lbl_distance_analysis_prot_struct_2.setText("Protein structure 2")
        self.fill_protein_boxes_batch()
        if self._view.ui.list_distance_analysis_overview.count() > 0:
            try:
                self._view.ui.list_distance_analysis_overview.currentItem().setSelected(False)
            except AttributeError:
                constants.PYSSA_LOGGER.debug("No selection in struction analysis overview.")

    def structure_analysis_back(self) -> None:
        """Hides the gui elements to choose the two proteins."""
        gui_elements_to_show = [
            self._view.ui.lbl_distance_analysis_overview,
            self._view.ui.list_distance_analysis_overview,
            self._view.ui.btn_distance_analysis_add,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_distance_analysis_remove,
            self._view.ui.lbl_distance_analysis_prot_struct_1,
            self._view.ui.lbl_distance_analysis_prot_struct_2,
            self._view.ui.lbl_distance_analysis_vs,
            self._view.ui.lbl_distance_analysis_prot_1_chains,
            self._view.ui.list_distance_analysis_prot_1_chains,
            self._view.ui.btn_distance_analysis_back_2,
            self._view.ui.btn_distance_analysis_next_2,
            self._view.ui.combobox_distance_analysis_prot_struct_1,
            self._view.ui.combobox_distance_analysis_prot_struct_2,
            self._view.ui.btn_distance_analysis_back,
            self._view.ui.btn_distance_analysis_next,
            self._view.ui.lbl_distance_analysis_prot_2_chains,
            self._view.ui.list_distance_analysis_prot_2_chains,
            self._view.ui.btn_distance_analysis_back_3,
            self._view.ui.btn_distance_analysis_next_3,
            
            
            self._view.ui.btn_distance_analysis_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        if self._view.ui.list_distance_analysis_overview.count() > 0:
            self._view.ui.btn_distance_analysis_remove.show()
            self._view.ui.btn_distance_analysis_remove.setEnabled(False)
            self._view.ui.btn_distance_analysis_start.show()
            self._view.ui.lbl_distance_analysis_ray_trace_images.show()
            self._view.ui.cb_distance_analysis_ray_trace_images.show()

    def structure_analysis_next_2(self) -> None:
        """Shows the gui elements to select the chains in protein 2."""
        gui_elements_to_show = [
            self._view.ui.lbl_distance_analysis_overview,
            self._view.ui.list_distance_analysis_overview,
            self._view.ui.lbl_distance_analysis_prot_struct_1,
            self._view.ui.lbl_distance_analysis_prot_struct_2,
            self._view.ui.lbl_distance_analysis_vs,
            self._view.ui.lbl_distance_analysis_prot_1_chains,
            self._view.ui.list_distance_analysis_prot_1_chains,
            self._view.ui.lbl_distance_analysis_prot_2_chains,
            self._view.ui.list_distance_analysis_prot_2_chains,
            self._view.ui.btn_distance_analysis_back_3,
            self._view.ui.btn_distance_analysis_next_3,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_distance_analysis_remove,
            self._view.ui.btn_distance_analysis_add,
            self._view.ui.combobox_distance_analysis_prot_struct_1,
            self._view.ui.combobox_distance_analysis_prot_struct_2,
            self._view.ui.btn_distance_analysis_back,
            self._view.ui.btn_distance_analysis_next,
            self._view.ui.btn_distance_analysis_back_2,
            self._view.ui.btn_distance_analysis_next_2,
            
            
            self._view.ui.btn_distance_analysis_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.list_distance_analysis_prot_2_chains.clear()
        self._view.ui.list_distance_analysis_prot_1_chains.setEnabled(False)
        self._view.ui.btn_distance_analysis_next_3.setEnabled(False)

        tmp_protein = self._interface_manager.get_current_project().search_protein(
            self._view.ui.combobox_distance_analysis_prot_struct_2.currentText())
        for tmp_chain in tmp_protein.chains:
            if tmp_chain.chain_type == "protein_chain":
                self._view.ui.list_distance_analysis_prot_2_chains.addItem(tmp_chain.chain_letter)
        if len(self._view.ui.list_distance_analysis_prot_1_chains.selectedItems()) == 1:
            self._view.ui.lbl_distance_analysis_prot_2_chains.setText(
                f"Select 1 chain in protein structure {self._view.ui.lbl_distance_analysis_prot_struct_2.text()}.",
            )
        else:
            self._view.ui.lbl_distance_analysis_prot_2_chains.setText(
                f"Select {len(self._view.ui.list_distance_analysis_prot_1_chains.selectedItems())} chains in "
                f"protein structure {self._view.ui.lbl_distance_analysis_prot_struct_2.text()}.",
            )

    def structure_analysis_back_2(self) -> None:
        """Hides the gui elements to select the chains in protein 1."""
        gui_elements_to_show = [
            self._view.ui.lbl_distance_analysis_overview,
            self._view.ui.list_distance_analysis_overview,
            self._view.ui.lbl_distance_analysis_prot_struct_1,
            self._view.ui.combobox_distance_analysis_prot_struct_1,
            self._view.ui.lbl_distance_analysis_vs,
            self._view.ui.lbl_distance_analysis_prot_struct_2,
            self._view.ui.combobox_distance_analysis_prot_struct_2,
            self._view.ui.btn_distance_analysis_back,
            self._view.ui.btn_distance_analysis_next,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_distance_analysis_remove,
            self._view.ui.btn_distance_analysis_add,
            self._view.ui.lbl_distance_analysis_prot_1_chains,
            self._view.ui.list_distance_analysis_prot_1_chains,
            self._view.ui.btn_distance_analysis_back_2,
            self._view.ui.btn_distance_analysis_next_2,
            self._view.ui.lbl_distance_analysis_prot_2_chains,
            self._view.ui.list_distance_analysis_prot_2_chains,
            self._view.ui.btn_distance_analysis_back_3,
            self._view.ui.btn_distance_analysis_next_3,
            
            
            self._view.ui.btn_distance_analysis_start,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.lbl_distance_analysis_prot_struct_1.setText("Protein structure 1")
        self._view.ui.lbl_distance_analysis_prot_struct_2.setText("Protein structure 2")

    def structure_analysis_next_3(self) -> None:
        """Adds the protein pair to the list of protein pairs to analyze."""
        gui_elements_to_show = [
            self._view.ui.btn_distance_analysis_remove,
            self._view.ui.btn_distance_analysis_add,
            self._view.ui.lbl_distance_analysis_overview,
            self._view.ui.list_distance_analysis_overview,
            
            
            self._view.ui.btn_distance_analysis_start,
        ]
        gui_elements_to_hide = [
            self._view.ui.combobox_distance_analysis_prot_struct_1,
            self._view.ui.combobox_distance_analysis_prot_struct_2,
            self._view.ui.lbl_distance_analysis_prot_struct_1,
            self._view.ui.lbl_distance_analysis_prot_struct_2,
            self._view.ui.lbl_distance_analysis_vs,
            self._view.ui.lbl_distance_analysis_prot_1_chains,
            self._view.ui.list_distance_analysis_prot_1_chains,
            self._view.ui.lbl_distance_analysis_prot_2_chains,
            self._view.ui.list_distance_analysis_prot_2_chains,
            self._view.ui.btn_distance_analysis_back,
            self._view.ui.btn_distance_analysis_next,
            self._view.ui.btn_distance_analysis_back_2,
            self._view.ui.btn_distance_analysis_next_2,
            self._view.ui.btn_distance_analysis_back_3,
            self._view.ui.btn_distance_analysis_next_3,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        prot_1_name = self._view.ui.lbl_distance_analysis_prot_struct_1.text()
        prot_1_chains = []
        for chain in self._view.ui.list_distance_analysis_prot_1_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self._view.ui.lbl_distance_analysis_prot_struct_2.text()
        prot_2_chains = []
        for chain in self._view.ui.list_distance_analysis_prot_2_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        item = QtWidgets.QListWidgetItem(analysis_name)
        self._view.ui.list_distance_analysis_overview.addItem(item)
        self._view.ui.btn_distance_analysis_remove.setEnabled(False)

    def structure_analysis_back_3(self) -> None:
        """Hides the gui elements to select the chains in protein 2."""
        gui_elements_to_show = [
            self._view.ui.lbl_distance_analysis_overview,
            self._view.ui.list_distance_analysis_overview,
            self._view.ui.lbl_distance_analysis_prot_struct_1,
            self._view.ui.lbl_distance_analysis_prot_struct_2,
            self._view.ui.lbl_distance_analysis_vs,
            self._view.ui.lbl_distance_analysis_prot_1_chains,
            self._view.ui.list_distance_analysis_prot_1_chains,
            self._view.ui.btn_distance_analysis_back_2,
            self._view.ui.btn_distance_analysis_next_2,
        ]
        gui_elements_to_hide = [
            self._view.ui.btn_distance_analysis_remove,
            self._view.ui.btn_distance_analysis_add,
            self._view.ui.combobox_distance_analysis_prot_struct_1,
            self._view.ui.combobox_distance_analysis_prot_struct_2,
            self._view.ui.btn_distance_analysis_back,
            self._view.ui.btn_distance_analysis_next,
            self._view.ui.btn_distance_analysis_back_3,
            self._view.ui.btn_distance_analysis_next_3,
            
            
            self._view.ui.btn_distance_analysis_start,
            self._view.ui.lbl_distance_analysis_prot_2_chains,
            self._view.ui.list_distance_analysis_prot_2_chains,
        ]
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)
        self._view.ui.list_distance_analysis_prot_1_chains.setEnabled(True)

        # tmp_protein = self._interface_manager.get_current_project().search_protein(self._view.ui.combobox_distance_analysis_prot_struct_2.currentText())
        # for tmp_chain in tmp_protein.chains:
        #     if tmp_chain.chain_type == "protein_chain":
        #         self._view.ui.list_distance_analysis_prot_1_chains.addItem(tmp_chain.chain_letter)

    def structure_analysis_overview_clicked(self) -> None:
        """Enables the remove button."""
        self._view.ui.btn_distance_analysis_remove.setEnabled(True)

    def fill_protein_boxes_batch(self) -> None:
        """Fills the combo boxes with the protein names."""
        proteins = []
        for tmp_protein in self._interface_manager.get_current_project().proteins:
            proteins.append(tmp_protein.get_molecule_object())
        proteins.insert(0, "")
        self._view.ui.combobox_distance_analysis_prot_struct_1.clear()
        self._view.ui.combobox_distance_analysis_prot_struct_2.clear()
        gui_utils.fill_combo_box(self._view.ui.combobox_distance_analysis_prot_struct_1, proteins)
        gui_utils.fill_combo_box(self._view.ui.combobox_distance_analysis_prot_struct_2, proteins)

    def remove_analysis_run(self) -> None:
        """Removes the selected protein pair from the list of protein pairs to analyze."""
        self._view.ui.list_distance_analysis_overview.takeItem(self._view.ui.list_distance_analysis_overview.currentRow())
        if self._view.ui.list_distance_analysis_overview.count() == 0:
            gui_elements_to_show = [
                self._view.ui.lbl_distance_analysis_overview,
                self._view.ui.list_distance_analysis_overview,
                self._view.ui.btn_distance_analysis_add,
            ]

            gui_elements_to_hide = [
                self._view.ui.btn_distance_analysis_remove,
                self._view.ui.lbl_distance_analysis_prot_struct_1,
                self._view.ui.lbl_distance_analysis_prot_struct_2,
                self._view.ui.lbl_distance_analysis_vs,
                self._view.ui.lbl_distance_analysis_prot_1_chains,
                self._view.ui.list_distance_analysis_prot_1_chains,
                self._view.ui.btn_distance_analysis_back_2,
                self._view.ui.btn_distance_analysis_next_2,
                self._view.ui.combobox_distance_analysis_prot_struct_1,
                self._view.ui.combobox_distance_analysis_prot_struct_2,
                self._view.ui.btn_distance_analysis_back,
                self._view.ui.btn_distance_analysis_next,
                self._view.ui.lbl_distance_analysis_prot_2_chains,
                self._view.ui.list_distance_analysis_prot_2_chains,
                self._view.ui.btn_distance_analysis_back_3,
                self._view.ui.btn_distance_analysis_next_3,
                
                
                self._view.ui.btn_distance_analysis_start,
            ]

            gui_utils.show_gui_elements(gui_elements_to_show)
            gui_utils.hide_gui_elements(gui_elements_to_hide)
            self._view.ui.btn_distance_analysis_remove.hide()
        else:
            if self._view.ui.list_distance_analysis_overview.count() > 0:
                try:
                    self._view.ui.list_distance_analysis_overview.currentItem().setSelected(False)
                except AttributeError:
                    constants.PYSSA_LOGGER.debug("No selection in struction analysis overview.")
        self._view.ui.btn_distance_analysis_remove.setEnabled(False)

    def check_if_same_no_of_chains_selected_batch(self) -> None:
        """Checks if the same number of proteins were selected."""
        self._view.ui.btn_distance_analysis_next_3.setEnabled(False)

        if self.no_of_selected_chains == len(self._view.ui.list_distance_analysis_prot_2_chains.selectedItems()):
            self._view.ui.btn_distance_analysis_next_3.setEnabled(True)

        prot_1_name = self._view.ui.lbl_distance_analysis_prot_struct_1.text()
        prot_1_chains = []
        for chain in self._view.ui.list_distance_analysis_prot_1_chains.selectedItems():
            prot_1_chains.append(chain.text())
        prot_1_chains = ",".join([str(elem) for elem in prot_1_chains])
        prot_2_name = self._view.ui.lbl_distance_analysis_prot_struct_2.text()
        prot_2_chains = []
        for chain in self._view.ui.list_distance_analysis_prot_2_chains.selectedItems():
            prot_2_chains.append(chain.text())
        prot_2_chains = ",".join([str(elem) for elem in prot_2_chains])
        analysis_name = f"{prot_1_name};{prot_1_chains}_vs_{prot_2_name};{prot_2_chains}"
        for tmp_row in range(self._view.ui.list_distance_analysis_overview.count()):
            if analysis_name == self._view.ui.list_distance_analysis_overview.item(tmp_row).text():
                self._view.ui.btn_distance_analysis_next_3.setEnabled(False)
                return

    def check_if_prot_structs_are_filled_batch(self) -> None:
        """Checks if two proteins were selected."""
        prot_1 = self._view.ui.combobox_distance_analysis_prot_struct_1.itemText(
            self._view.ui.combobox_distance_analysis_prot_struct_1.currentIndex(),
        )
        prot_2 = self._view.ui.combobox_distance_analysis_prot_struct_2.itemText(
            self._view.ui.combobox_distance_analysis_prot_struct_2.currentIndex(),
        )
        if prot_1 != "" and prot_2 != "":
            self._view.ui.btn_distance_analysis_next.setEnabled(True)
        else:
            self._view.ui.btn_distance_analysis_next.setEnabled(False)

    def count_batch_selected_chains_for_prot_struct_1(self) -> None:
        """Counts the number of chains of protein 1."""
        self.no_of_selected_chains = len(self._view.ui.list_distance_analysis_prot_1_chains.selectedItems())
        if self.no_of_selected_chains > 0:
            self._view.ui.btn_distance_analysis_next_2.setEnabled(True)
        else:
            self._view.ui.btn_distance_analysis_next_2.setEnabled(False)

    # </editor-fold>