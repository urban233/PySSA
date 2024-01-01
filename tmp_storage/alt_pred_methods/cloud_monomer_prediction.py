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
"""Module for storing code which can be used to implement the Google colab prediction method for a monomer"""

# -- Gui page management vars
self.cloud_pred_monomer_management: gui_page_management.GuiPageManagement

self.ui.list_new_seq_notebooks.addItem(constants.OFFICIAL_NOTEBOOK_NAME)
self.ui.list_new_seq_notebooks.currentItemChanged.connect(self.enable_predict_button)

# management
self._create_cloud_pred_monomer_management()
# setup defaults for pages
self._init_new_sequence_page()

# ----- Functions for GuiPageManagement obj creation
def _create_cloud_pred_monomer_management(self):
    # gui element management
    tmp_stages = [
        # protein name stage
        stage.Stage(
            {
                "label_protein_name": self.ui.lbl_prediction_only_protein_name,
                "text_field_protein_name": self.ui.txt_prediction_only_protein_name,
                "label_protein_name_status": self.ui.lbl_prediction_only_status_protein_name,
            },
            {
                "next_button": self.ui.btn_prediction_only_next,
            },
        ),
        # protein sequence stage
        stage.Stage(
            {
                "label_protein_sequence": self.ui.lbl_cloud_pred_mono_prot_seq,
                "text_field_protein_sequence": self.ui.txt_cloud_pred_mono_prot_seq,
                "label_protein_sequence_status": self.ui.lbl_cloud_pred_mono_status_prot_seq,
            },
            {
                "back_button": self.ui.btn_cloud_pred_mono_back,
                "next_button": self.ui.btn_cloud_pred_mono_next_2,
            },
        ),
        # prediction stage (with advanced configurations)
        stage.Stage(
            {
                "label_choose_notebook": self.ui.lbl_prediction_only_choose_notebook,
                "list_choose_notebook": self.ui.list_new_seq_notebooks,
                "label_advanced_config": self.ui.lbl_cloud_pred_mono_advanced_config,
                "button_advanced_config": self.ui.btn_cloud_pred_mono_advanced_config,
            },
            {
                "back_button": self.ui.btn_prediction_only_back,
                "predict_button": self.ui.btn_prediction_only_start,
            },
        ),
    ]
    self.cloud_pred_monomer_management = gui_page_management.GuiPageManagement(tmp_stages)


# connections
self.ui.btn_pred_cloud_monomer_page.clicked.connect(self.display_prediction_only_page)
# monomer cloud prediction page
self.ui.txt_prediction_only_protein_name.textChanged.connect(self.validate_cloud_mono)
self.ui.btn_prediction_only_next.clicked.connect(self.show_cloud_pred_mono_stage_1)
self.ui.btn_cloud_pred_mono_back.clicked.connect(self.show_cloud_pred_mono_stage_0)
self.ui.txt_cloud_pred_mono_prot_seq.textChanged.connect(self.validate_cloud_mono)
self.ui.btn_cloud_pred_mono_next_2.clicked.connect(self.show_cloud_pred_mono_stage_2)
self.ui.btn_prediction_only_back.clicked.connect(self.show_cloud_pred_mono_stage_1)
self.ui.btn_cloud_pred_mono_advanced_config.clicked.connect(self.show_prediction_configuration)
self.ui.btn_prediction_only_start.clicked.connect(self.predict_only)


def _init_new_sequence_page(self):
    """This function clears all text fields and hides everything which is needed"""
    self.cloud_pred_monomer_management.show_stage_x(0)
    self.cloud_pred_monomer_management.disable_all_next_buttons()
    self.cloud_pred_monomer_management.set_empty_string_in_label()

    # # clears everything
    # self.ui.txt_prediction_only_protein_name.clear()
    # self.ui.txt_cloud_pred_mono_prot_seq.clear()
    # # sets up defaults: Prediction
    # self.ui.btn_prediction_only_start.setEnabled(False)
    # self.ui.btn_prediction_only_next.setEnabled(False)
    # self.ui.btn_cloud_pred_mono_next_2.setEnabled(False)
    # self.ui.lbl_cloud_pred_mono_status_prot_seq.setText("")
    #
    # gui_elements = [
    #     self.ui.lbl_cloud_pred_mono_prot_seq,
    #     self.ui.txt_cloud_pred_mono_prot_seq,
    #     self.ui.lbl_cloud_pred_mono_status_prot_seq,
    #     self.ui.btn_cloud_pred_mono_back,
    #     self.ui.btn_cloud_pred_mono_next_2,
    #     self.ui.lbl_prediction_only_choose_notebook,
    #     self.ui.list_new_seq_notebooks,
    #     self.ui.lbl_cloud_pred_mono_advanced_config,
    #     self.ui.btn_cloud_pred_mono_advanced_config,
    #     self.ui.btn_prediction_only_back,
    #     self.ui.btn_prediction_only_start,
    # ]
    # gui_utils.hide_gui_elements(gui_elements)


def display_prediction_only_page(self):
    """This function displays the prediction only work area"""
    tools.switch_page(self.ui.stackedWidget, self.ui.lbl_page_title, 1, "New Sequence")
    item = self.ui.list_new_seq_notebooks.findItems("AlphaFold", Qt.QtCore.Qt.MatchContains | Qt.QtCore.Qt.MatchExactly)
    self.ui.list_new_seq_notebooks.setCurrentItem(item[0])


# ----- Functions for Monomer Cloud Prediction
def show_cloud_pred_mono_stage_0(self):
    self.cloud_pred_monomer_management.show_stage_x(0)
    self.ui.btn_cloud_pred_mono_advanced_config.hide()


def show_cloud_pred_mono_stage_1(self):
    self.cloud_pred_monomer_management.show_stage_x(1)


def show_cloud_pred_mono_stage_2(self):
    self.cloud_pred_monomer_management.show_stage_x(2)


def validate_cloud_mono(self):
    self.cloud_pred_monomer_management.create_validation()


@staticmethod
def show_prediction_configuration():
    config = dialog_advanced_prediction_configurations.DialogAdvancedPredictionConfigurations()
    config.exec_()


# def validate_protein_name(self):
#     """This function validates the input of the project name in real-time
#
#     """
#     tools.validate_protein_name(self.ui.txt_prediction_only_protein_name,
#                                 self.ui.lbl_prediction_only_status_protein_name,
#                                 self.ui.btn_prediction_only_next)
#
# def validate_protein_sequence(self):
#     """This function validates the input of the protein sequence in real-time
#
#     """
#     tools.validate_protein_sequence(self.ui.txt_cloud_pred_mono_prot_seq,
#                                     self.ui.lbl_cloud_pred_mono_status_prot_seq,
#                                     self.ui.btn_cloud_pred_mono_next_2)
#
# def show_cloud_prediction_mono_protein_sequence(self):
#     gui_elements_hide = [
#         self.ui.btn_prediction_only_next,
#     ]
#     gui_elements_show = [
#         self.ui.lbl_cloud_pred_mono_prot_seq,
#         self.ui.lbl_cloud_pred_mono_status_prot_seq,
#         self.ui.txt_cloud_pred_mono_prot_seq,
#         self.ui.btn_cloud_pred_mono_back,
#         self.ui.btn_cloud_pred_mono_next_2,
#     ]
#     gui_utils.manage_gui_visibility(gui_elements_show, gui_elements_hide)
#     gui_utils.disable_text_box(self.ui.txt_prediction_only_protein_name,
#                                self.ui.lbl_prediction_only_protein_name)
#
# def hide_cloud_prediction_mono_protein_sequence(self):
#     gui_elements_hide = [
#         self.ui.lbl_cloud_pred_mono_prot_seq,
#         self.ui.lbl_cloud_pred_mono_status_prot_seq,
#         self.ui.txt_cloud_pred_mono_prot_seq,
#         self.ui.btn_cloud_pred_mono_back,
#         self.ui.btn_cloud_pred_mono_next_2,
#     ]
#     gui_elements_show = [
#         self.ui.btn_prediction_only_next,
#     ]
#     gui_utils.manage_gui_visibility(gui_elements_show, gui_elements_hide)
#     gui_utils.enable_text_box(self.ui.txt_prediction_only_protein_name,
#                               self.ui.lbl_prediction_only_protein_name)
#
# def show_prediction_only_choose_notebook(self):
#     gui_elements_show = [
#         self.ui.lbl_prediction_only_choose_notebook,
#         self.ui.list_new_seq_notebooks,
#         self.ui.btn_prediction_only_back,
#         self.ui.btn_prediction_only_start,
#         self.ui.list_new_seq_notebooks,
#         self.ui.lbl_cloud_pred_mono_advanced_config,
#         self.ui.btn_cloud_pred_mono_advanced_config,
#     ]
#     gui_elements_hide = [
#         self.ui.btn_cloud_pred_mono_back,
#         self.ui.btn_cloud_pred_mono_next_2,
#     ]
#     gui_utils.manage_gui_visibility(gui_elements_show, gui_elements_hide)
#     gui_utils.disable_text_box(self.ui.txt_cloud_pred_mono_prot_seq,
#                                self.ui.lbl_cloud_pred_mono_prot_seq)
#     self.ui.btn_prediction_only_next.setEnabled(True)
#     styles.color_button_ready(self.ui.btn_prediction_only_next)
#
# def hide_prediction_only_choose_notebook(self):
#     gui_elements_hide = [
#         self.ui.lbl_prediction_only_choose_notebook,
#         self.ui.list_new_seq_notebooks,
#         self.ui.btn_prediction_only_back,
#         self.ui.btn_prediction_only_start,
#         self.ui.lbl_cloud_pred_mono_advanced_config,
#         self.ui.btn_cloud_pred_mono_advanced_config,
#     ]
#     gui_elements_show = [
#         self.ui.btn_cloud_pred_mono_back,
#         self.ui.btn_cloud_pred_mono_next_2,
#     ]
#     gui_utils.manage_gui_visibility(gui_elements_show, gui_elements_hide)
#     gui_utils.enable_text_box(self.ui.txt_cloud_pred_mono_prot_seq,
#                               self.ui.lbl_cloud_pred_mono_prot_seq)
#     self.ui.btn_prediction_only_next.setEnabled(True)
#     styles.color_button_ready(self.ui.btn_prediction_only_next)
#
# def check_prediction_only_if_txt_notebook_url_is_filled(self):
#     """This function checks if a reference pdb file is selected.
#
#     """
#     self.__check_start_possibility_prediction_only()
# def validate_project_name(self):
#     """This function validates the input of the project name in real-time
#
#     """
#     if self.ui.list_widget_projects.currentItem() is not None:
#         self.ui.list_widget_projects.currentItem().setSelected(False)
#     # set color for lineEdit
#     self.ui.txt_prediction_project_name.setStyleSheet("background-color: #FC5457")
#     if len(self.ui.txt_prediction_project_name.text()) == 0:
#         self.ui.btn_prediction_next_1.setEnabled(False)
#         styles.color_button_not_ready(self.ui.btn_prediction_next_1)
#         return
#     elif len(self.ui.txt_prediction_project_name.text()) > 20:
#         self.ui.btn_prediction_next_1.setEnabled(False)
#         styles.color_button_not_ready(self.ui.btn_prediction_next_1)
#         self.ui.lbl_prediction_status_project_name.setText("Project name is too long (max. 20 characters).")
#         return
#     else:
#         regex = Qt.QtCore.QRegularExpression()
#         # TO-DO: has no dash in regex!
#         regex.setPattern("(([a-z])|([A-Z])|([0-9])|(-)|(_)){0,20}")
#         validator = QtGui.QRegularExpressionValidator(regex)
#         for i in range(len(self.ui.txt_prediction_project_name.text())):
#             result = validator.validate(self.ui.txt_prediction_project_name.text(), i)
#             if result[0] > 0:
#                 self.ui.txt_prediction_project_name.setStyleSheet("background-color: #33C065")
#                 self.ui.btn_prediction_next_1.setEnabled(True)
#                 styles.color_button_ready(self.ui.btn_prediction_next_1)
#             else:
#                 self.ui.txt_prediction_project_name.setStyleSheet("background-color: #FC5457")
#                 self.ui.btn_prediction_next_1.setEnabled(False)
#                 styles.color_button_not_ready(self.ui.btn_prediction_next_1)
#                 self.ui.lbl_prediction_status_project_name.setText("Invalid character.")
#                 return
#         item = self.ui.list_widget_projects.findItems(self.ui.txt_prediction_project_name.text(),
#                                                       Qt.QtCore.Qt.MatchContains |
#                                                       Qt.QtCore.Qt.MatchExactly
#                                                       )
#         if len(item) != 0:
#             self.ui.list_widget_projects.setCurrentItem(item[0])
#             self.ui.txt_prediction_project_name.setStyleSheet("background-color: #FC5457")
#             self.ui.btn_prediction_next_1.setEnabled(False)
#             styles.color_button_not_ready(self.ui.btn_prediction_next_1)
#             self.ui.lbl_prediction_status_project_name.setText("Project name already exists.")
#         # else:
#         #     self.ui.list_widget_projects.currentItem().setSelected(False)
#         #     self.ui.txt_prediction_project_name.setStyleSheet("background-color: green")
#         #     self.ui.btn_prediction_next_1.setEnabled(True)
#         #     styles.color_button_ready(self.ui.btn_prediction_next_1)
#         print("Check successful.")


def predict_only(self):
    """This function is used to predict with any google colab notebook."""
    # create prediction_params file for colab notebook
    seq = self.ui.txt_cloud_pred_mono_prot_seq.text()
    job_name = self.ui.lbl_current_project_name.text()
    prediction_params = {
        "seq": seq,
        "job_name": job_name,
        "amber": "false",
        "templates": ("pdb70", "1"),
    }
    params_file = open(f"{constants.SCRATCH_DIR}\\prediction_params.json", "w", encoding="utf-8")
    json.dump(prediction_params, params_file, indent=4)
    params_file.close()
    self.activate_second_python()
    archive = f"{constants.NOTEBOOK_RESULTS_ZIP_NAME}.zip"
    if not os.path.exists(f"{self.scratch_path}/{archive}"):
        print("Prediction was aborted.")
        os.remove(f"{constants.SCRATCH_DIR}\\prediction_params.json")
        return
    # colabfold: AlphaFold2_mmseqs2 notebook specific process
    shutil.unpack_archive(
        f"{self.scratch_path}/{archive}", f"{self.scratch_path}/{constants.NOTEBOOK_RESULTS_ZIP_NAME}", "zip"
    )
    prediction_results: list[str] = os.listdir(
        pathlib.Path(f"{constants.SCRATCH_DIR}/{constants.NOTEBOOK_RESULTS_ZIP_NAME}")
    )
    for filename in prediction_results:
        check = filename.find("_relaxed_rank_1")
        if check != -1:
            src = pathlib.Path(f"{constants.SCRATCH_DIR}/{constants.NOTEBOOK_RESULTS_ZIP_NAME}/{filename}")
            dest = pathlib.Path(f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{filename}")
            shutil.copy(src, dest)
            os.rename(
                f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{filename}",
                f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{self.ui.txt_prediction_only_protein_name.text()}.pdb",
            )
            break
    shutil.rmtree(f"{self.scratch_path}/{constants.NOTEBOOK_RESULTS_ZIP_NAME}")
    os.remove(f"{self.scratch_path}/{archive}")
    try:
        cmd.load(
            f"{self.workspace_path}/{self.ui.lbl_current_project_name.text()}/pdb/{self.ui.txt_prediction_only_protein_name.text()}.pdb"
        )
    except pymol.CmdException:
        print("Loading the model failed.")
        return
    self._project_watcher.show_valid_options(self.ui)
    self.display_view_page()


@staticmethod
def activate_second_python():
    subprocess.run(
        [
            "python.exe",
            f"{os.path.expanduser('~')}\\anaconda3\\envs\\pyssa\\lib\\site-packages\\pmg_tk\\startup\\tmpPySSA\\pyssa\\colab_notebook.py",
        ]
    )
