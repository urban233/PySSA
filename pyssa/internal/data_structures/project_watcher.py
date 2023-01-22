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
"""Module for the project watcher class which is hard coded based on pyssa"""
import os
from internal.data_structures import project
from util import gui_utils


class ProjectWatcher:
    """A watcher class that observes the project status

    Attributes:
        current_project:
            project object of the current active project
        no_of_pdb_files:
            number of how many pdb files are stored in the project
        no_of_results:
            number of how many results are in the project

    """
    def __init__(self,
                 current_project: project.Project,
                 no_of_pdb_files=0,
                 no_of_results=0,
                 ):
        self.current_project: project.Project = current_project
        self.no_of_pdb_files: int = no_of_pdb_files
        self.no_of_results: int = no_of_results

    def show_valid_options(self, ui):
        """This function shows all valid options based on the number of pdb files and results in the project

        Args:
            ui:
                ui variable of the Qt main window

        Notes:
            This function counts the pdb files and results automatically!
        """
        self.count_pdb_files()
        self.count_results()

        if self.no_of_pdb_files is None:
            gui_elements_to_show = [
                ui.btn_new_page,
                ui.btn_open_page,
                ui.btn_delete_page,
            ]
            gui_elements_to_hide = [
                ui.btn_close_project,
                ui.btn_save_project,
                ui.btn_view_page,
                ui.btn_use_page,
                ui.btn_edit_page,
                ui.btn_hotspots_page,
                ui.lbl_pred_cloud,
                ui.btn_pred_cloud_monomer_page,
                ui.btn_pred_cloud_monomer_vs_pdb_page,
                ui.btn_pred_cloud_multimer_page,
                ui.btn_pred_cloud_multimer_vs_pdb_page,
                ui.lbl_pred_local,
                ui.btn_pred_local_monomer_page,
                ui.btn_pred_local_monomer_vs_pdb_page,
                ui.btn_pred_local_multimer_page,
                ui.btn_pred_local_multimer_vs_pdb_page,
                ui.lbl_analysis,
                ui.btn_single_analysis_page,
                ui.btn_batch_analysis_page,
                ui.btn_results_page,
                ui.lbl_handle_pymol_session,
                ui.btn_image_page,
            ]
            gui_utils.manage_gui_visibility(gui_elements_to_show, gui_elements_to_hide)
        elif self.no_of_pdb_files == 0:
            gui_elements_to_show = [
                ui.btn_close_project,
                ui.btn_use_page,
                ui.btn_edit_page,
                ui.btn_save_project,
                ui.btn_view_page,
                ui.lbl_pred_local,
                ui.btn_pred_local_monomer_page,
                ui.btn_pred_local_multimer_page,
            ]
            gui_elements_to_hide = [
                ui.btn_new_page,
                ui.btn_open_page,
                ui.btn_delete_page,
                ui.lbl_pred_cloud,
                ui.btn_pred_cloud_monomer_page,
                ui.btn_pred_cloud_multimer_page,
                ui.btn_pred_cloud_monomer_vs_pdb_page,
                ui.btn_pred_cloud_multimer_vs_pdb_page,
                ui.btn_pred_local_monomer_vs_pdb_page,
                ui.btn_pred_local_multimer_vs_pdb_page,
                ui.lbl_analysis,
                ui.btn_single_analysis_page,
                ui.btn_batch_analysis_page,
                ui.btn_results_page,
                ui.lbl_handle_pymol_session,
                ui.btn_image_page,
                ui.btn_hotspots_page,
            ]
            gui_utils.manage_gui_visibility(gui_elements_to_show, gui_elements_to_hide)
        elif self.no_of_pdb_files == 1:
            gui_elements_to_show = [
                ui.btn_close_project,
                ui.btn_save_project,
                ui.btn_use_page,
                ui.btn_edit_page,
                ui.btn_view_page,
                ui.lbl_pred_local,
                ui.btn_pred_local_monomer_vs_pdb_page,
                ui.btn_pred_local_multimer_vs_pdb_page,
                ui.lbl_handle_pymol_session,
                ui.btn_image_page,
                ui.btn_hotspots_page,
            ]
            gui_elements_to_hide = [
                ui.btn_new_page,
                ui.btn_open_page,
                ui.btn_delete_page,
                ui.btn_pred_cloud_monomer_page,
                ui.btn_pred_cloud_multimer_page,
                ui.btn_pred_local_monomer_page,
                ui.btn_pred_local_multimer_page,
                ui.lbl_analysis,
                ui.btn_single_analysis_page,
                ui.btn_batch_analysis_page,
                ui.btn_results_page,
                ui.lbl_pred_cloud,
                ui.btn_pred_cloud_monomer_vs_pdb_page,
                ui.btn_pred_cloud_multimer_vs_pdb_page,
            ]
            gui_utils.manage_gui_visibility(gui_elements_to_show, gui_elements_to_hide)
        # elif self.no_of_pdb_files == 2:
        #     gui_elements_to_show = [
        #         ui.btn_close_project,
        #         ui.btn_save_project,
        #         ui.btn_use_page,
        #         ui.btn_edit_page,
        #         ui.btn_view_page,
        #         ui.lbl_analysis,
        #         ui.btn_single_analysis_page,
        #         ui.lbl_handle_pymol_session,
        #         ui.btn_image_page,
        #     ]
        #     gui_elements_to_hide = [
        #         ui.btn_new_page,
        #         ui.btn_open_page,
        #         ui.btn_delete_page,
        #         ui.lbl_pred_cloud,
        #         ui.btn_pred_cloud_monomer_page,
        #         ui.btn_pred_cloud_monomer_vs_pdb_page,
        #         ui.btn_pred_cloud_multimer_page,
        #         ui.btn_pred_cloud_multimer_vs_pdb_page,
        #         ui.lbl_pred_local,
        #         ui.btn_pred_local_monomer_vs_pdb_page,
        #         ui.btn_pred_local_multimer_vs_pdb_page,
        #         ui.btn_pred_local_monomer_page,
        #         ui.btn_pred_local_multimer_page,
        #         ui.btn_batch_analysis_page,
        #         ui.btn_hotspots_page,
        #     ]
        #     if self.no_of_results > 0:
        #         gui_elements_to_show.append(ui.btn_results_page)
        #     else:
        #         gui_elements_to_hide.append(ui.btn_results_page)
        #     gui_utils.manage_gui_visibility(gui_elements_to_show, gui_elements_to_hide)
        elif self.no_of_pdb_files >= 2:
            gui_elements_to_show = [
                ui.btn_close_project,
                ui.btn_save_project,
                ui.btn_use_page,
                ui.btn_edit_page,
                ui.btn_view_page,
                ui.lbl_analysis,
                ui.btn_batch_analysis_page,
                ui.lbl_handle_pymol_session,
                ui.btn_image_page,
                ui.btn_hotspots_page,
            ]
            gui_elements_to_hide = [
                ui.btn_new_page,
                ui.btn_open_page,
                ui.btn_delete_page,
                ui.btn_single_analysis_page,
                ui.lbl_pred_cloud,
                ui.btn_pred_cloud_monomer_page,
                ui.btn_pred_cloud_monomer_vs_pdb_page,
                ui.btn_pred_cloud_multimer_page,
                ui.btn_pred_cloud_multimer_vs_pdb_page,
                ui.lbl_pred_local,
                ui.btn_pred_local_monomer_vs_pdb_page,
                ui.btn_pred_local_multimer_vs_pdb_page,
                ui.btn_pred_local_monomer_page,
                ui.btn_pred_local_multimer_page,
            ]
            if self.no_of_results > 0:
                gui_elements_to_show.append(ui.btn_results_page)
            else:
                gui_elements_to_hide.append(ui.btn_results_page)
            gui_utils.manage_gui_visibility(gui_elements_to_show, gui_elements_to_hide)

    def count_pdb_files(self):
        """This function counts the number of pdb files in the project pdb directory

        Notes:
            You do NOT need to use this function before using "show_valid_options"!!
        """
        if self.current_project.get_pdb_path() == "pdb":
            self.no_of_pdb_files = None
        elif not os.path.exists(self.current_project.get_pdb_path()):
            raise ValueError(f"The pdb path {self.current_project.get_pdb_path()} does not exist!")
        else:
            self.no_of_pdb_files = len(os.listdir(self.current_project.get_pdb_path()))

    def count_results(self):
        """This function counts the number of results in the project

        Notes:
            You do NOT need to use this function before using "show_valid_options"!!
        """
        if self.current_project.get_results_path() == "results":
            print("Empty project.")
        elif not os.path.exists(self.current_project.get_results_path()):
            raise ValueError(f"The results path {self.current_project.get_results_path()} does not exist!")
        else:
            self.no_of_results = len(os.listdir(self.current_project.get_results_path()))
