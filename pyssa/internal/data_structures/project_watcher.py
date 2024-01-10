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
"""Module for the project watcher class which is hard coded based on pyssa."""
import pathlib

from pyssa.internal.data_structures import project
from pyssa.util import gui_utils
from pyssa.util import tools


class ProjectWatcher:
    """A watcher class that observes the project status.

    Attributes:
        current_project:
            project object of the current active project
        no_of_pdb_files:
            number of how many pdb files are stored in the project
        no_of_results:
            number of how many results are in the project

    """

    def __init__(
        self,
        current_project: "project.Project",
        no_of_pdb_files: int = 0,
        no_of_results: int = 0,
        make_images: bool = False,
    ) -> None:
        """Constructor.

        Args:
            current_project: the current project.
            no_of_pdb_files: the number of pdb files in the current project.
            no_of_results: the number of results in the current project.
            make_images: a flag whether images should be generated or not.
        """
        self.current_project: project.Project = current_project
        self.no_of_pdb_files: int = no_of_pdb_files
        self.no_of_results: int = no_of_results
        self.make_images: bool = make_images
        self.on_home_page: bool = True

    def show_valid_options(self, ui) -> None:  # noqa: ANN001 #TODO: needs to be checked
        """This function shows all valid options based on the number of pdb files and results in the project.

        Args:
            ui: an ui variable of the Qt main window

        Notes:
            This function counts the pdb files and results automatically!
        """
        self.count_proteins_in_project()
        self.count_results()
        self.check_images()

        if self.no_of_pdb_files is None:
            gui_elements_to_show = [
                ui.btn_new_page,
                ui.btn_open_page,
                ui.btn_delete_page,
                ui.btn_import_project,
            ]
            gui_elements_to_hide = [
                ui.btn_close_project,
                ui.btn_view_page,
                ui.btn_use_page,
                ui.btn_export_project,
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
                ui.btn_prediction_abort,
                ui.lbl_analysis,
                ui.btn_single_analysis_page,
                ui.btn_batch_analysis_page,
                ui.btn_image_analysis_page,
                ui.lbl_pred_analysis,
                ui.btn_pred_analysis_monomer_page,
                ui.btn_pred_analysis_multimer_page,
                ui.btn_results_page,
                ui.btn_analysis_abort,
                ui.lbl_handle_pymol_session,
                ui.btn_manage_session,
                ui.btn_image_page,
            ]
            gui_utils.manage_gui_visibility(gui_elements_to_show, gui_elements_to_hide)
        elif self.no_of_pdb_files == 0:
            gui_elements_to_show = [
                ui.btn_close_project,
                ui.btn_use_page,
                ui.btn_edit_page,
                ui.btn_view_page,
                ui.lbl_pred_local,
                ui.btn_pred_local_monomer_page,
                ui.btn_pred_local_multimer_page,
            ]
            gui_elements_to_hide = [
                ui.btn_new_page,
                ui.btn_open_page,
                ui.btn_delete_page,
                ui.btn_import_project,
                ui.btn_export_project,
                ui.lbl_pred_cloud,
                ui.btn_pred_cloud_monomer_page,
                ui.btn_pred_cloud_multimer_page,
                ui.btn_pred_cloud_monomer_vs_pdb_page,
                ui.btn_pred_cloud_multimer_vs_pdb_page,
                ui.btn_pred_local_monomer_vs_pdb_page,
                ui.btn_pred_local_multimer_vs_pdb_page,
                ui.btn_prediction_abort,
                ui.lbl_analysis,
                ui.btn_single_analysis_page,
                ui.btn_batch_analysis_page,
                ui.btn_image_analysis_page,
                ui.btn_results_page,
                ui.btn_analysis_abort,
                ui.lbl_pred_analysis,
                ui.btn_pred_analysis_monomer_page,
                ui.btn_pred_analysis_multimer_page,
                ui.lbl_handle_pymol_session,
                ui.btn_manage_session,
                ui.btn_image_page,
                ui.btn_hotspots_page,
            ]
            gui_utils.manage_gui_visibility(gui_elements_to_show, gui_elements_to_hide)
        elif self.no_of_pdb_files == 1:
            gui_elements_to_show = [
                ui.btn_close_project,
                ui.btn_use_page,
                ui.btn_export_project,
                ui.btn_edit_page,
                ui.btn_view_page,
                ui.lbl_pred_analysis,
                ui.btn_pred_analysis_monomer_page,
                ui.btn_pred_analysis_multimer_page,
                ui.lbl_handle_pymol_session,
                ui.btn_manage_session,
                ui.btn_image_page,
                ui.btn_hotspots_page,
            ]
            gui_elements_to_hide = [
                ui.btn_new_page,
                ui.btn_open_page,
                ui.btn_delete_page,
                ui.btn_import_project,
                ui.btn_pred_cloud_monomer_page,
                ui.btn_pred_cloud_multimer_page,
                ui.btn_pred_local_monomer_page,
                ui.btn_pred_local_multimer_page,
                ui.btn_prediction_abort,
                ui.lbl_analysis,
                ui.btn_single_analysis_page,
                ui.btn_batch_analysis_page,
                ui.btn_image_analysis_page,
                ui.btn_results_page,
                ui.btn_analysis_abort,
                ui.lbl_pred_cloud,
                ui.btn_pred_cloud_monomer_vs_pdb_page,
                ui.btn_pred_cloud_multimer_vs_pdb_page,
                ui.lbl_pred_local,
                ui.btn_pred_local_monomer_vs_pdb_page,
                ui.btn_pred_local_multimer_vs_pdb_page,
            ]
            gui_utils.manage_gui_visibility(gui_elements_to_show, gui_elements_to_hide)
        elif self.no_of_pdb_files >= 2:
            gui_elements_to_show = [
                ui.btn_close_project,
                ui.btn_use_page,
                ui.btn_export_project,
                ui.btn_edit_page,
                ui.btn_view_page,
                ui.lbl_analysis,
                ui.btn_batch_analysis_page,
                ui.lbl_handle_pymol_session,
                ui.btn_manage_session,
                ui.btn_image_page,
                ui.btn_hotspots_page,
            ]
            gui_elements_to_hide = [
                ui.btn_new_page,
                ui.btn_open_page,
                ui.btn_delete_page,
                ui.btn_import_project,
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
                ui.lbl_pred_analysis,
                ui.btn_pred_analysis_monomer_page,
                ui.btn_pred_analysis_multimer_page,
                ui.btn_prediction_abort,
                ui.btn_analysis_abort,
            ]
            if self.no_of_results > 0:
                gui_elements_to_show.append(ui.btn_results_page)
            else:
                gui_elements_to_hide.append(ui.btn_results_page)
            if self.make_images is False:
                gui_elements_to_hide.append(ui.btn_image_analysis_page)
            else:
                gui_elements_to_show.append(ui.btn_image_analysis_page)
            gui_utils.manage_gui_visibility(gui_elements_to_show, gui_elements_to_hide)

    def check_workspace_for_projects(self, workspace: pathlib.Path, ui) -> None:  # noqa: ANN001
        """Checks the workspace for any projects.

        Args:
            workspace: the path to the workspace.
            ui: an instance of the PySSA ui class.
        """
        if len(tools.scan_workspace_for_valid_projects(workspace, ui.list_delete_projects)) > 0:
            ui.btn_delete_page.show()
            ui.btn_open_page.show()
        else:
            ui.btn_delete_page.hide()
            ui.btn_open_page.hide()

    def count_proteins_in_project(self) -> None:
        """Counts the number of pdb files in the project pdb directory.

        Notes:
            You do NOT need to use this function before using "show_valid_options"!!
        """
        if not self.on_home_page:
            self.no_of_pdb_files = len(self.current_project.proteins)
        else:
            self.no_of_pdb_files = None

    def count_results(self) -> None:
        """Counts the number of results in the project.

        Notes:
            You do NOT need to use this function before using "show_valid_options"!!
        """
        if not self.on_home_page:
            self.no_of_results = len(self.current_project.protein_pairs)
        else:
            self.no_of_results = None

    def check_images(self) -> None:
        """Checks if images are within the project."""
        for tmp_protein_pair in self.current_project.protein_pairs:
            if len(tmp_protein_pair.distance_analysis.analysis_results.structure_aln_image) == 0:
                self.make_images = True
                return
        self.make_images = False

    def count_images(self) -> None:
        """This function counts the number of images in the project.

        Notes:
            You do NOT need to use this function before using "show_valid_options"!!
        """
        if not self.on_home_page:
            no_of_protein_pair_with_images = 0
            for tmp_protein_pair in self.current_project.protein_pairs:
                if len(tmp_protein_pair.distance_analysis.analysis_results.interesting_regions_images) > 0:
                    no_of_protein_pair_with_images += 1
            self.no_of_images = no_of_protein_pair_with_images
        else:
            self.no_of_images = None
