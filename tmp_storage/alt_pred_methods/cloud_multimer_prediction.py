#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2022
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
"""Module for storing code which can be used to implement the Google colab prediction method for a multimer"""

# -- Gui page management vars
self.cloud_pred_multimer_management: gui_page_management.GuiPageManagement

self.ui.list_cloud_pred_multi_seq_notebooks.addItem(constants.OFFICIAL_NOTEBOOK_NAME)
self.ui.list_cloud_pred_multi_seq_notebooks.currentItemChanged.connect(self.enable_cloud_predict_button)

# management
self._create_cloud_pred_multimer_management()

# setup defaults for pages
self._init_cloud_pred_multi_page()

# ----- Functions for GuiPageManagement obj creation
def _create_cloud_pred_multimer_management(self):
    # gui element management
    tmp_stages = [
        # protein name stage
        stage.Stage(
            {
                "label_protein_name": self.ui.lbl_cloud_pred_multi_protein_name,
                "text_field_protein_name": self.ui.txt_cloud_pred_multi_protein_name,
                "label_protein_name_status": self.ui.lbl_cloud_pred_multi_status_protein_name,
            },
            {
                "next_button": self.ui.btn_cloud_pred_multi_next,
            }
        ),
        # protein sequence stage
        stage.Stage(
            {
                "label_protein_sequence": self.ui.lbl_cloud_pred_multi_prot_seq_single,
                "text_field_protein_sequence": self.ui.txt_cloud_pred_multi_prot_seq,
                "label_protein_sequence_status": self.ui.lbl_cloud_pred_multi_status_prot_seq,
                "button_add_sequence": self.ui.btn_cloud_pred_multi_add_seq_single,
                "label_protein_sequence_overview": self.ui.lbl_cloud_pred_multi_prot_overview,
                "table_protein_sequence_overview": self.ui.table_cloud_pred_multi_prot_overview,
                "label_selected_protein_1": self.ui.lbl_cloud_pred_multi_selected_protein_1,
                "button_edit_selected_protein": self.ui.btn_cloud_pred_multi_selected_protein_edit,
                "label_selected_protein_2": self.ui.lbl_cloud_pred_multi_selected_protein_2,
                "button_remove_selected_protein": self.ui.btn_cloud_pred_multi_selected_protein_remove,
            },
            {
                "back_button": self.ui.btn_cloud_pred_multi_back,
                "next_button": self.ui.btn_cloud_pred_multi_next_2,
            }
        ),
        # prediction stage (with advanced configurations)
        stage.Stage(
            {
                "label_choose_notebook": self.ui.lbl_cloud_pred_multi_choose_notebook,
                "list_choose_notebook": self.ui.list_cloud_pred_multi_seq_notebooks,
                "label_advanced_config": self.ui.lbl_cloud_pred_multi_advanced_config,
                "button_advanced_config": self.ui.btn_cloud_pred_multi_advanced_config,
            },
            {
                "back_button": self.ui.btn_cloud_pred_multi_back_2,
                "predict_button": self.ui.btn_cloud_pred_multi_predict,
            }
        ),
    ]
    self.cloud_pred_multimer_management = gui_page_management.GuiPageManagement(tmp_stages)

# connections
self.ui.btn_pred_cloud_multimer_page.clicked.connect(self.display_cloud_pred_multi)

# multimer cloud prediction page
self.ui.btn_cloud_pred_multi_next.clicked.connect(self.show_cloud_pred_multi_stage_1)
self.ui.btn_cloud_pred_multi_back.clicked.connect(self.show_cloud_pred_multi_stage_0)
self.ui.btn_cloud_pred_multi_next_2.clicked.connect(self.show_cloud_pred_multi_stage_2)
self.ui.btn_cloud_pred_multi_back_2.clicked.connect(self.show_cloud_pred_multi_stage_1)
self.ui.txt_cloud_pred_multi_prot_seq.textChanged.connect(self.validate_cloud_pred_multi)
self.ui.txt_cloud_pred_multi_protein_name.textChanged.connect(self.validate_cloud_pred_multi)
self.ui.btn_cloud_pred_multi_add_seq_single.clicked.connect(self.add_sequence_cloud_pred_multi)
self.ui.btn_cloud_pred_multi_selected_protein_edit.clicked.connect(self.edit_protein_cloud_pred_multi)
self.ui.btn_cloud_pred_multi_selected_protein_remove.clicked.connect(self.remove_protein_cloud_pred_multi)
self.ui.btn_cloud_pred_multi_predict.clicked.connect(self.cloud_multi_predict)


def _init_cloud_pred_multi_page(self):
    self.cloud_pred_multimer_management.show_stage_x(0)
    self.cloud_pred_multimer_management.disable_all_next_buttons()
    self.cloud_pred_multimer_management.set_empty_string_in_label()

def display_cloud_pred_multi(self):
    tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 17, "Cloud Multimer Prediction")
    item = self.ui.list_cloud_pred_multi_seq_notebooks.findItems("AlphaFold",
                                                                 Qt.QtCore.Qt.MatchContains |
                                                                 Qt.QtCore.Qt.MatchExactly
                                                                 )
    self.ui.list_cloud_pred_multi_seq_notebooks.setCurrentItem(item[0])

# ----- Functions for Multimer Cloud Prediction
def show_cloud_pred_multi_stage_0(self):
    self.cloud_pred_multimer_management.show_stage_x(0)

def show_cloud_pred_multi_stage_1(self):
    self.cloud_pred_multimer_management.show_stage_x(1)
    first_table_item = self.ui.table_cloud_pred_multi_prot_overview.item(0, 0)
    if first_table_item is not None:
        for row in range(self.ui.table_cloud_pred_multi_prot_overview.rowCount()):
            new_item = QTableWidgetItem(self.ui.txt_cloud_pred_multi_protein_name.text())
            self.ui.table_cloud_pred_multi_prot_overview.setItem(row, 0, new_item)

def show_cloud_pred_multi_stage_2(self):
    self.cloud_pred_multimer_management.show_stage_x(2)

def validate_cloud_pred_multi(self):
    self.cloud_pred_multimer_management.create_validation()

def add_sequence_cloud_pred_multi(self):
    no_of_rows = self.ui.table_cloud_pred_multi_prot_overview.rowCount()
    if no_of_rows < 26:
        self.ui.table_cloud_pred_multi_prot_overview.setRowCount(no_of_rows + 1)
        name = QTableWidgetItem(self.ui.txt_cloud_pred_multi_protein_name.text())
        chain = QTableWidgetItem(constants.chain_dict[no_of_rows])
        sequence = QTableWidgetItem(self.ui.txt_cloud_pred_multi_prot_seq.text())
        self.ui.table_cloud_pred_multi_prot_overview.setItem(no_of_rows, 0, name)
        self.ui.table_cloud_pred_multi_prot_overview.setItem(no_of_rows, 1, chain)
        self.ui.table_cloud_pred_multi_prot_overview.setItem(no_of_rows, 2, sequence)
    else:
        self.ui.btn_cloud_pred_multi_add_seq_single.setEnabled(False)

def edit_protein_cloud_pred_multi(self):
    dialog = dialog_edit_prot_table.DialogEditProtTable()
    dialog.protein_name = self.ui.table_cloud_pred_multi_prot_overview.item(
        self.ui.table_cloud_pred_multi_prot_overview.currentRow(), 0
    )
    dialog.chain = self.ui.table_cloud_pred_multi_prot_overview.item(
        self.ui.table_cloud_pred_multi_prot_overview.currentRow(), 1
    )
    dialog.sequence = self.ui.table_cloud_pred_multi_prot_overview.item(
        self.ui.table_cloud_pred_multi_prot_overview.currentRow(), 2
    )
    dialog.set_all_parameters(self.ui.table_cloud_pred_multi_prot_overview,
                              self.ui.table_cloud_pred_multi_prot_overview.currentRow())
    dialog.exec_()

def remove_protein_cloud_pred_multi(self):
    self.ui.table_cloud_pred_multi_prot_overview.removeRow(self.ui.table_cloud_pred_multi_prot_overview.currentRow())

def enable_cloud_predict_button(self):
    self.ui.btn_cloud_pred_multi_predict.setEnabled(True)
    styles.color_button_ready(self.ui.btn_prediction_only_start)

def cloud_multi_predict(self):
    """This function is used to predict with any google colab notebook.

    """
    self.ui.action_toggle_notebook_visibility.setVisible(True)
    self.web_gui = web_interface.WebInterface()
    # prepare sequences
    seqs = []
    for row in range(self.ui.table_cloud_pred_multi_prot_overview.rowCount()):
        tmp_seq = self.ui.table_cloud_pred_multi_prot_overview.item(row, 2).text()
        seqs.append(tmp_seq)
    complete_sequence = ':'.join(seqs)

    self.web_gui.set_protein_sequence(complete_sequence)
    self.web_gui.set_job_name(self.ui.lbl_current_project_name.text())
    self.web_gui.show_interface()
    if self.web_gui.exit_code != 0:
        print("An error occurred!")
        return
    self.ui.action_toggle_notebook_visibility.setVisible(False)
    # colabfold: AlphaFold2_mmseqs2 notebook specific process
    archive = f"{constants.NOTEBOOK_RESULTS_ZIP_NAME}.zip"
    shutil.unpack_archive(f"{self.scratch_path}/{archive}",
                          f"{self.scratch_path}/{constants.NOTEBOOK_RESULTS_ZIP_NAME}", "zip")
    prediction_results: list[str] = os.listdir(
        pathlib.Path(f"{constants.SCRATCH_DIR}/{constants.NOTEBOOK_RESULTS_ZIP_NAME}"))
    for filename in prediction_results:
        check = filename.find("_relaxed_rank_1")
        if check != -1:
            src = pathlib.Path(f"{constants.SCRATCH_DIR}/{constants.NOTEBOOK_RESULTS_ZIP_NAME}/{filename}")
            dest = pathlib.Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{filename}")
            shutil.copy(src, dest)
            os.rename(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{filename}",
                      f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{self.ui.txt_prediction_only_protein_name.text()}.pdb")
            break
    shutil.rmtree(f"{self.scratch_path}/{constants.NOTEBOOK_RESULTS_ZIP_NAME}")
    os.remove(f"{self.scratch_path}/{archive}")
    try:
        cmd.load(
            f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{self.ui.txt_prediction_only_protein_name.text()}.pdb")
    except pymol.CmdException:
        print("Loading the model failed.")
        return
    self._project_watcher.show_valid_options(self.ui)
    self.display_view_page()

@staticmethod
def activate_second_python():
        subprocess.run(["python.exe", f"{os.path.expanduser('~')}\\anaconda3\\envs\\pyssa\\lib\\site-packages\\pmg_tk\\startup\\tmpPySSA\\pyssa\\colab_notebook.py"])