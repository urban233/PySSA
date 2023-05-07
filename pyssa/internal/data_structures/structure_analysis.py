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
"""Module for structure analysis class"""
import logging
import os.path
import pathlib
import shutil
from typing import TYPE_CHECKING
from pyssa.io_pyssa import path_util
from pyssa.logging_pyssa import log_handlers
from pyssa.internal.data_structures import protein_pair
from pyssa.util import constants

if TYPE_CHECKING:
    from pyssa.internal.data_structures import protein_pair
    from pyssa.internal.data_structures import project

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class Analysis:
    """This class contains information about the type of analysis"""

    analysis_list: list[protein_pair.ProteinPair] = []
    app_project: 'project.Project'

    def __init__(self, app_project):
        self.app_project = app_project

    def run_analysis(self):
        logger.debug(self.analysis_list)
        for tmp_protein_pair in self.analysis_list:
            tmp_protein_pair.distance_analysis.do_analysis_in_pymol(self.app_project)
            # TODO: Comment back in!!!
            if not os.path.exists(constants.SCRATCH_DIR_IMAGES):
                os.mkdir(constants.SCRATCH_DIR_IMAGES)
            if not os.path.exists(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR):
                os.mkdir(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR)
            if not os.path.exists(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR):
                os.mkdir(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR)
            tmp_protein_pair.distance_analysis.take_image_of_protein_pair(filename=f"structure_aln_{tmp_protein_pair.name}",
                                                                          representation="cartoon")
            tmp_protein_pair.distance_analysis.analysis_results.set_structure_aln_image(
                path_util.FilePath(
                    pathlib.Path(
                        f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR}/structure_aln_{tmp_protein_pair.name}.png")
                )
            )
            tmp_protein_pair.distance_analysis.take_image_of_interesting_regions(tmp_protein_pair.distance_analysis.cutoff,
                                                                                 f"interesting_reg_{tmp_protein_pair.name}")
            interesting_region_filepaths = []
            for tmp_filename in os.listdir(constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR):
                interesting_region_filepaths.append(
                    path_util.FilePath(
                        pathlib.Path(f"{constants.SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR}/{tmp_filename}")
                    )
                )
            tmp_protein_pair.distance_analysis.analysis_results.set_interesting_region_images(interesting_region_filepaths)
            self.app_project.add_protein_pair(tmp_protein_pair)
            shutil.rmtree(constants.SCRATCH_DIR_IMAGES)
        self.analysis_list.clear()

# class StructureAnalysis:
#     """This class is used for the structure analysis
#
#     """
#     reference_protein = []
#     model_proteins = []
#     ref_chains = []
#     model_chains = []
#     response_create_images = False
#     export_dir: str = ""
#     _cycles: int = 1
#     _cutoff: float = 1.0
#     _figure_size: tuple[float, float] = (11.0, 6.0)
#     _distance_label: str = "$\\alpha$-C"
#     _alignment_file_name: str = "alignment_file_model_s"
#     _session_file_name: str = "session_file_model_s"  # little Tesla easter egg ;)
#
#     def __init__(self, reference_protein: list[protein.Protein], model_proteins: list[protein.Protein],
#                  ref_chains, model_chains, export_dir, cycles=1, cutoff=1.0) -> None:
#         self.reference_protein = reference_protein
#         self.model_proteins = model_proteins
#         self.ref_chains = ref_chains
#         self.model_chains = model_chains
#         self.export_dir = export_dir
#         self._cycles = cycles
#         self._cutoff = cutoff
#
#     def get_export_dir(self) -> str:
#         """This function returns the value for export_dir
#
#         Returns:
#             export_dir
#         """
#         return self.export_dir
#
#     def set_export_dir(self, value) -> None:
#         """This function sets the value for export_dir
#
#         Args:
#             value:
#                 value which should get set
#         """
#         self.export_dir = value
#
#     def get_cycles(self) -> int:
#         """This function returns the value for get_cycles
#
#         Returns:
#             get_cycles
#         """
#         return self._cycles
#
#     def set_cycles(self, value) -> None:
#         """This function sets the value for get_cycles
#
#         Args:
#             value:
#                 value which should get set
#         """
#         self._cycles = value
#
#     def get_cutoff(self) -> float:
#         """This function returns the value for cutoff
#
#         Returns:
#             cutoff
#         """
#         return self._cutoff
#
#     def set_cutoff(self, value) -> None:
#         """This function sets the value for cutoff
#
#         Args:
#             value:
#                 value which should get set
#         """
#         self._cutoff = value
#
#     def get_figure_size(self) -> tuple[float, float]:
#         """This function returns the value for figure_size
#
#         Returns:
#             figure_size
#         """
#         return self._figure_size
#
#     def set_figure_size(self, value) -> None:
#         """This function sets the value for figure_size
#
#         Args:
#             value:
#                 value which should get set
#         """
#         self._figure_size = value
#
#     def get_distance_label(self) -> str:
#         """This function returns the value for distance_label
#
#         Returns:
#             distance_label
#         """
#         return self._distance_label
#
#     def set_distance_label(self, value) -> None:
#         """This function sets the value for distance_label
#
#         Args:
#             value:
#                 value which should get set
#         """
#         self._distance_label = value
#
#     def get_alignment_file_name(self) -> str:
#         """This function returns the value for alignment_file_name
#
#         Returns:
#             alignment_file_name
#         """
#         return self._alignment_file_name
#
#     def set_alignment_file_name(self, value) -> None:
#         """This function sets the value for alignment_file_name
#
#         Args:
#             value:
#                 value which should get set
#         """
#         self._alignment_file_name = value
#
#     def get_session_file_name(self) -> str:
#         """This function returns the value for session_file_name
#
#         Returns:
#             session_file_name
#         """
#         return self._session_file_name
#
#     def set_session_file_name(self, value) -> None:
#         """This function sets the value for session_file_name
#
#         Args:
#             value:
#                 value which should get set
#         """
#         self._session_file_name = value
#
#     @staticmethod
#     def prepare_chains_for_selection(txt_box_chain_ref):
#         """This function prepares the raw chain input into a format which can be used by the
#         function create_selection_for_proteins()
#
#         Returns:
#             a string which is in the format for the function create_selection_for_proteins()
#         """
#         return txt_box_chain_ref.text().split(",")
#
#     @staticmethod
#     def create_selection_for_proteins(chains: list, proteins: list[protein.Protein]) -> None:
#         """This function creates the selections for the proteins and sets these directly into
#         the Protein objects
#
#         Args:
#             chains:
#                 list which contains the raw chain information from the input box
#             proteins:
#                 list which contains the Protein objects
#         """
#         pymol_io.load_protein(proteins[0].filepath, proteins[0].filename, proteins[0].molecule_object)
#         d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
#                  'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
#                  'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
#                  'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
#         i = 0
#         for chain in chains:
#             atoms_of_chain = cmd.get_model(f"chain {chain}")
#             first_aa = atoms_of_chain.atom[0].resn
#             if first_aa not in d3to1:
#                 chains.pop(i)
#             i += 1
#
#         # sets selection for any chain combination of reference
#         selection = ""
#         seperator = ", "
#         tmp_list = []
#         # creates list containing all possible selections for the Protein
#         for chain in chains:
#             tmp_selection = f"/{proteins[0].molecule_object}//{chain}//CA"
#             tmp_list.append(tmp_selection)
#             selection = seperator.join(tmp_list)
#         # adds the selections into the Protein object
#         for tmp_protein in proteins:
#             if selection == "":
#                 selection = f"/{tmp_protein.molecule_object}////CA"
#             tmp_protein.set_selection(selection)
#
#     def create_protein_pairs(self) -> list['protein_pair.ProteinPair']:
#         """This function creates Protein pairs based on the proteins in the object
#
#         Returns:
#             protein_pairs: list of Protein pairs
#         """
#         tmp_protein_pairs = []
#         for model in self.model_proteins:
#             tmp_protein_pair = protein_pair.ProteinPair(self.reference_protein[0], model, self.export_dir)
#             tmp_protein_pairs.append(tmp_protein_pair)
#         return tmp_protein_pairs
#
#     def do_analysis_in_pymol(self, protein_pairs: list['protein_pair.ProteinPair'],
#                              status_bar_obj: Qt.QtWidgets.QStatusBar) -> None:
#         """This function does the actual structure analysis in pymol based on the Protein pairs
#
#         Args:
#             protein_pairs:
#                 list of Protein pairs which were created with the function create_protein_pairs
#             status_bar_obj:
#                 the statusbar object from the main window
#         """
#         for tmp_protein_pair in protein_pairs:
#             # reinitialize pymol session
#             msg = "Reinitialize pymol session ..."
#
#             # load both proteins in the pymol session
#             msg = "Finished reinitializing pymol session. | Load proteins ..."
#             tmp_protein_pair.load_protein_pair()
#
#             # color Protein with default colors; ref: green, model: blue
#             msg = "Finished loading proteins. | Color proteins ..."
#             tmp_protein_pair.color_protein_pair()
#
#             # do the structure alignment
#             msg = "Finished coloring proteins. | Align proteins ..."
#             try:
#                 align_results = tmp_protein_pair.align_protein_pair(self._cycles, self._cutoff, self._alignment_file_name)
#             except pymol.CmdException:
#                 if len(self.ref_chains) != 0:
#                     tools.quick_log_and_display("critical", "There is a problem with the chain selection!",
#                                                 status_bar_obj, "Critical Error: Please check your chain information!")
#                     gui_utils.critical_message("Please check your chain information!",
#                                                "There is a problem with the chain selection!"
#                                                "You should check if you have entered the correct chains.")
#                 return
#             rmsd_dict = {
#                 "rmsd": str(round(align_results[0], 2)),
#                 "aligned_residues": str(align_results[1]),
#             }
#             rmsd_file = open(pathlib.Path(f"{self.export_dir}/rmsd.json"), "w", encoding="utf-8")
#             json.dump(rmsd_dict, rmsd_file, indent=4)
#             rmsd_file.close()
#
#             # do the distance calculation
#             msg = "Finished aligning proteins. | Calculate distances ..."
#             distance_results = tmp_protein_pair.calculate_distance_between_ca_atoms(self._alignment_file_name)
#             tmp_protein_pair.export_distance_between_ca_atoms(distance_results)
#
#             # create an instance of the Graphics class
#             # msg = "Finished calculating distances. | Create distance plot ..."
#             # tools.quick_log_and_display(GLOBAL_VAR_LOG_TYPE, msg, status_bar_obj, msg)
#             # graphics_instance: graphics.Graphics = graphics.Graphics(tmp_protein_pair, distance_results, self._figure_size)
#
#             # TODO: is this needed due to the pyqtgraph plots?
#             # # create distance plot
#             # fig = graphics_instance.create_distance_plot(self._distance_label, self._cutoff)
#             # # save distance plot
#             # if not os.path.exists(f"{tmp_protein_pair.results_dir}/plots"):
#             #     os.mkdir(f"{tmp_protein_pair.results_dir}/plots")
#             # if not os.path.exists(f"{tmp_protein_pair.results_dir}/plots/distance_plot"):
#             #     os.mkdir(f"{tmp_protein_pair.results_dir}/plots/distance_plot")
#             # plt.savefig(f"{tmp_protein_pair.results_dir}/plots/distance_plot/"
#             #             f"distance_plot_{tmp_protein_pair.model_obj.molecule_object}.svg")
#             # plt.close(fig)
#
#             # # create distance histogram
#             # msg = "Finished creating distance plot. | Create distance histogram ..."
#             # tools.quick_log_and_display(GLOBAL_VAR_LOG_TYPE, msg, status_bar_obj, msg)
#             # graphics_instance.create_distance_histogram()
#             # # save distance histogram
#             # if not os.path.exists(f"{tmp_protein_pair.results_dir}/plots"):
#             #     os.mkdir(f"{tmp_protein_pair.results_dir}/plots")
#             # if not os.path.exists(
#             #         f"{tmp_protein_pair.results_dir}/plots/distance_histogram"):
#             #     os.mkdir(f"{tmp_protein_pair.results_dir}/plots/distance_histogram")
#             # plt.savefig(f"{tmp_protein_pair.results_dir}/plots/distance_histogram"
#             #             f"/distance_histogram_{tmp_protein_pair.model_obj.molecule_object}.svg")
#
#             # if self.response_create_images is True:
#             #     # take image of whole structure alignment
#             #     msg = "Finished creating distance histogram. | Take image of structure alignment ..."
#             #     tools.quick_log_and_display(GLOBAL_VAR_LOG_TYPE, msg, status_bar_obj, msg)
#             #     graphics_instance.take_image_of_protein_pair(self._alignment_file_name, "cartoon", "structure_alignment")
#             #     # take image of interesting regions
#             #     msg = f"Finished taking image of structure alignment. | Take images "\
#             #           f"of interesting regions (within {self._cutoff} angstrom)"
#             #     tools.quick_log_and_display(GLOBAL_VAR_LOG_TYPE, msg, status_bar_obj, msg)
#             #     graphics_instance.take_image_of_interesting_regions(3.0, "interesting_region", opaque_background=1)
#             # else:
#             #     # take image of whole structure alignment
#             #     msg = "Finished creating distance histogram. | Create scene of structure alignment ..."
#             #     tools.quick_log_and_display(GLOBAL_VAR_LOG_TYPE, msg, status_bar_obj, msg)
#             #     graphics_instance.take_image_of_protein_pair(self._alignment_file_name, "cartoon",
#             #                                                  "structure_alignment", take_images=False)
#             #     # take image of interesting regions
#             #     msg = f"Finished creating scene of structure alignment. | Create scenes " \
#             #           f"of interesting regions (within {self._cutoff} angstrom)"
#             #     tools.quick_log_and_display(GLOBAL_VAR_LOG_TYPE, msg, status_bar_obj, msg)
#             #     graphics_instance.take_image_of_interesting_regions(3.0, "interesting_region", opaque_background=1, take_images=False)
#
#             # NOT THIS
#             # # color residues by distance
#             # graphics_instance.color_by_distance(ALIGNMENT_FILE_NAME)
#             # print(f"Finished coloring of prediction with color_by_distance.")
#             # graphics_instance.take_image_of_protein_pair(ALIGNMENT_FILE_NAME, "cartoon",
#             #                                              "coloredByRMSD")
#             # graphics_instance.create_gif()
#             # print(f"Finished creation of gif which shows the whole predicted "
#             #       f"structure colored by distance.")
#
#             # save pymol session
#             msg = "Finished taking images of interesting regions. | Save PyMOL session"
#             tmp_protein_pair.save_session_of_protein_pair(self._session_file_name)
#
#             msg = "Finished saving PyMOL session. | Protein structure analysis has successfully finished."
