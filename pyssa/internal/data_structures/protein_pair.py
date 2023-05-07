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
"""Module for the protein pair class"""
import os
import pathlib
import logging
from pyssa.logging_pyssa import log_handlers
from pyssa.io_pyssa import path_util
from pyssa.internal.portal import pymol_io
from pyssa.internal.portal import protein_pair_operations
from pyssa.io_pyssa import filesystem_io
from pyssa.util import protein_util
from pyssa.util import pyssa_keys
from xml.etree import ElementTree
from pyssa.io_pyssa.xml_pyssa import element_names
from pyssa.io_pyssa.xml_pyssa import attribute_names
from typing import TYPE_CHECKING
from pyssa.util import constants
from pyssa.io_pyssa import binary_data

if TYPE_CHECKING:
    from pyssa.internal.data_structures import protein
    from pyssa.internal.analysis_types import distance_analysis

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class ProteinPair:
    """This class consists of two Protein objects. It is used to have a better workflow for the analysis."""

    # <editor-fold desc="Class attributes">
    """
    the first protein of the protein pair
    """
    protein_1: 'protein.Protein'
    """
    the second protein of the protein pair
    """
    protein_2: 'protein.Protein'
    """
    a directory where all results related to the protein will be stored
    """
    distance_analysis: 'distance_analysis.DistanceAnalysis' = None
    """
    the full filepath where the session file is stored
    """
    pymol_session_filepath: path_util.FilePath
    """
    a base64 string of the pymol session
    """
    pymol_session: str

    def __init__(self, protein_1: 'protein.Protein', protein_2: 'protein.Protein') -> None:
        """Constructor.

        Args:
            protein_1 (core.protein.Protein):
                reference Protein object
            protein_2 (core.Protein):
                model Protein object

        Raises:
            NotADirectoryError: If directory not found.
        """
        self.protein_1: protein.Protein = protein_1
        self.protein_2: protein.Protein = protein_2
        self.name = f"{self.protein_1.get_molecule_object()}_with_{self.protein_2.get_molecule_object()}"
        self.pymol_session = pymol_io.convert_pymol_session_to_base64_string(self.name)
        # protein_pair_dirname = f"{protein_pairs_dirname}/{self.name}"
        #
        # self.protein_pair_subdirs = {
        #     pyssa_keys.PROTEIN_PAIR_SUBDIR: pathlib.Path(f"{protein_pair_dirname}"),
        #     pyssa_keys.PROTEIN_PAIR_SESSION_SUBDIR: pathlib.Path(f"{protein_pair_dirname}/session"),
        #     pyssa_keys.PROTEIN_PAIR_RESULTS_SUBDIR: pathlib.Path(f"{protein_pair_dirname}/results"),
        #     pyssa_keys.PROTEIN_PAIR_OBJECTS_SUBDIR: pathlib.Path(f"{protein_pair_dirname}/.objects"),
        # }
        # for key in self.protein_pair_subdirs:
        #     if not os.path.exists(self.protein_pair_subdirs.get(key)):
        #         os.mkdir(self.protein_pair_subdirs.get(key))
        # self.export_dirname = self.protein_pair_subdirs.get(pyssa_keys.PROTEIN_PAIR_RESULTS_SUBDIR)
        # self.pymol_session_filepath = path_util.FilePath(f"{self.protein_pair_subdirs.get(pyssa_keys.PROTEIN_PAIR_SESSION_SUBDIR)}/{self.name}_session.pse")

    def load_protein_pair_in_pymol(self):
        """This function loads to proteins into a NEW pymol session.

        """
        self.protein_1.load_protein_in_pymol()
        self.protein_2.load_protein_in_pymol()

    def load_pymol_session(self):
        """This function loads the existing pymol session of the pair.

        """
        session_filepath = pathlib.Path(f"{constants.CACHE_PYMOL_SESSION_DIR}/{self.name}_session.pse")
        if not os.path.exists(constants.CACHE_PYMOL_SESSION_DIR):
            os.mkdir(constants.CACHE_PYMOL_SESSION_DIR)
        binary_data.write_binary_file_from_base64_string(filepath=session_filepath, base64_data=self.pymol_session)
        pymol_io.load_pymol_session(session_filepath)

    def color_protein_pair(self) -> None:
        """This function colors both the reference and the model Protein.

        Note:
            Only the official colors from PyMOL are supported. These can
            be looked up under the `color values`_ page.

        Args:
            color_ref (str, optional):
                defines color for the reference Protein
            color_model (str, optional):
                defines color for the model Protein

        Raises:
            pymol.CmdException:
                Exception is raised if one or both proteins
                does not exist as pymol objects.

        """
        protein_pair_operations.color_protein_pair(self.protein_1.get_molecule_object(),
                                                   self.protein_2.get_molecule_object())

    def save_session_of_protein_pair(self) -> None:
        """This function saves the pymol session of the Protein pair.

        Note:
            The pse file will be saved under the relative path
            (if export_data_dir = "data/results"):
            ``data/results/sessions``

            The file name (filename) MUST NOT have the file extension .pse!

        """
        protein_pair_operations.save_session_of_protein_pair(constants.CACHE_PYMOL_SESSION_DIR)

    def set_distance_analysis(self, value):
        self.distance_analysis = value

    # def align_protein_pair(self, cycle_number: int, cutoff_value: float,
    #                        alignment_filename: str = None) -> tuple[float, int]:
    #     """This function aligns the model to the reference Protein, with the
    #     `align`_ command from PyMOL.
    #
    #     Note:
    #         Before this function can be used, the load_protein_pair
    #         needs to be executed. Moreover, two selections are needed
    #         which have to be set through the set_selection function,
    #         for each Protein.
    #
    #         If an alignment file should be saved, it will be stored
    #         under the relative path (if export_data_dir =
    #         "data/results"): ``data/results/alignment_files``
    #
    #         The filename of the alignment file (alignment_filename)
    #         MUST NOT have the file extension .aln!
    #
    #     Args:
    #         cycle_number (int):
    #             defines the number of refinement cycles
    #         cutoff_value (float):
    #             defines the value when residues are not aligned
    #         alignment_filename (str, optional):
    #             string which functions as filename for the align object
    #
    #     Returns:
    #         tuple of the results
    #
    #     Raises:
    #         ValueError:
    #             Exception is raised if
    #             the align_reference or align_model string is empty or
    #             the cycle_number or cutoff_value are equal or below zero.
    #
    #     Example:
    #         The ``alignment_filename`` variable should be as follows::
    #
    #             Yes:
    #                 # This is a correct file name, because it has NO
    #                 # file extension
    #                 f"aln_results_Bmp2_0_CA"
    #
    #             No:
    #                 # This file name is wrong, because it has a
    #                 # trailing .csv which is not compatible with
    #                 # the function
    #                 f"aln_results_Bmp2_0_CA.aln"
    #
    #     .. _align:
    #         https://pymolwiki.org/index.php/Align
    #     """
    #     # argument test
    #     if not safeguard.Safeguard.check_if_pymol_selection_is_valid(self.ref_obj.selection):
    #         raise ValueError(f"The reference has an invalid selection: {self.ref_obj.selection}")
    #     if not safeguard.Safeguard.check_if_pymol_selection_is_valid(self.model_obj.selection):
    #         raise ValueError(f"The model has an invalid selection: {self.model_obj.selection}")
    #     if not safeguard.Safeguard.check_if_number_is_positive(cycle_number):
    #         raise ValueError("Number of cycles must be greater or equal than zero.")
    #     if not safeguard.Safeguard.check_if_number_is_positive(cutoff_value):
    #         raise ValueError("The cutoff needs to be greater than zero.")
    #
    #     # This block runs if an alignObject should be created.
    #     if alignment_filename is not None:
    #         results = cmd.align(target=self.ref_obj.selection,
    #                             mobile=self.model_obj.selection,
    #                             object=alignment_filename,
    #                             cycles=cycle_number,
    #                             cutoff=cutoff_value,
    #                             quiet=0)
    #
    #         if not os.path.exists(pathlib.Path(f"{self.results_dir}/alignment_files")):
    #             os.mkdir(pathlib.Path(f"{self.results_dir}/alignment_files"))
    #
    #         # save the align object from pymol as alignment file
    #         cmd.save(pathlib.Path(f"{self.results_dir}/alignment_files/{alignment_filename}.aln"))
    #         return results[0], results[1]
    #
    #     # This block runs if no alignObject should be created.
    #     else:
    #         results = cmd.align(target=self.ref_obj.selection,
    #                             mobile=self.model_obj.selection,
    #                             cycles=cycle_number,
    #                             cutoff=cutoff_value,
    #                             quiet=1)
    #
    #         return results[0], results[1]
    #
    # def calculate_distance_between_ca_atoms(self, alignment_filename: str,
    #                                         cutoff: float = 20.0) -> dict[str, np.ndarray]:
    #     # def calculate_distance_between_ca_atoms(self, alignment_filename: str,
    #     #                                         cutoff: float = 20.0) -> Dict[str, np.ndarray]:
    #     """This function calculates the distances between the aligned alpha-C
    #     atoms. It uses the alignment file generated by the alignProteinPair
    #     function as base.
    #
    #     Note:
    #         The filename of the alignment file (alignment_filename) MUST NOT
    #         have the file extension .aln!
    #
    #     Args:
    #         alignment_filename (str):
    #             filename of the alignment file, generated from the
    #             alignProteinPair function
    #         cutoff (float, optional):
    #             all distances which are above the cutoff will NOT be in the
    #             results dataframe.
    #             The default value is 20 angstrom.
    #
    #     Returns:
    #         result_hashtable (dict[str, np.ndarray]):
    #             a hash table which has the following keys:
    #             'index', 'ref_chain', 'ref_pos', 'ref_resi', 'model_chain',
    #             'model_pos', 'model_resi', 'distance'
    #
    #     Example:
    #         The ``alignment_filename`` variable should be as follows::
    #
    #             Yes:
    #                 # This is a correct file name, because it has NO
    #                 # file extension
    #                 f"aln_results_Bmp2_0_CA"
    #
    #             No:
    #                 # This file name is wrong, because it has a
    #                 # trailing .csv which is not compatible with
    #                 # the function
    #                 f"aln_results_Bmp2_0_CA.aln"
    #     """
    #     # argument test
    #     if not safeguard.Safeguard.check_if_file_is_readable(pathlib.Path(f"{self.results_dir}/alignment_files/{alignment_filename}.aln")):
    #         print(f"File not found, in {self.results_dir}.")
    #
    #     cmd.create(f"{self.ref_obj.molecule_object}_CA", f"/{self.ref_obj.molecule_object}////CA")
    #     ref_ca_obj = cmd.get_model(f"{self.ref_obj.molecule_object}_CA")
    #     cmd.create(f"{self.model_obj.molecule_object}_CA", f"/{self.model_obj.molecule_object}////CA")
    #     model_ca_obj = cmd.get_model(f"{self.model_obj.molecule_object}_CA")
    #     # read in alignment file from alignProteinPair function
    #     align = AlignIO.read(pathlib.Path(f"{self.results_dir}/alignment_files/{alignment_filename}.aln"), "clustal")
    #     index_list: [int] = []
    #     ref_chain_list: [str] = []
    #     ref_pos_list: [int] = []
    #     ref_resi_list: [str] = []
    #     model_chain_list: [str] = []
    #     model_pos_list: [int] = []
    #     model_resi_list: [str] = []
    #     distance_list: [float] = []
    #
    #     i = 0  # i for alignment file
    #     j = 0  # j for reference
    #     k = 0  # k for model
    #     index = 0
    #     while i < int(align.get_alignment_length()):
    #         # gets executed if the reference contains a "-" in the alignment
    #         if align[0, i] == "-":
    #             i += 1
    #             k += 1
    #             j = j
    #         # gets executed if the model contains a "-" in the alignment
    #         elif align[1, i] == "-":
    #             i += 1
    #             j += 1
    #             k = k
    #         # gets executed if no "-" is found in the alignment
    #         else:
    #             # create var for chain, position and residue name for
    #             # the reference
    #             chain_ref: str = ref_ca_obj.atom[j].chain
    #             pos_ref: str = ref_ca_obj.atom[j].resi
    #             resi_ref: str = ref_ca_obj.atom[j].resn
    #
    #             # create var for chain, position and residue name for
    #             # the model
    #             chain_model: str = model_ca_obj.atom[k].chain
    #             pos_model: str = model_ca_obj.atom[k].resi
    #             resi_model: str = model_ca_obj.atom[k].resn
    #
    #             # calculate the distance between the alpha-C atoms
    #             atom1 = f"/{self.ref_obj.molecule_object}_CA//{chain_ref}/{pos_ref}/"
    #             atom2 = f"/{self.model_obj.molecule_object}_CA//{chain_model}/{pos_model}/"
    #             distance = round(cmd.get_distance(atom1, atom2), 2)
    #
    #             # TODO: should there be a cutoff?
    #             # gets executed if the distance is greater than the
    #             # pre-defined cutoff
    #             # if distance > cutoff:
    #             #     i += 1
    #             #     j += 1
    #             #     k += 1
    #             # else:
    #             #     # append calculated data to each separate list
    #             #     index_list.append(index)
    #             #     ref_chain_list.append(chain_ref)
    #             #     ref_pos_list.append(pos_ref)
    #             #     ref_resi_list.append(resi_ref)
    #             #     model_chain_list.append(chain_model)
    #             #     model_pos_list.append(pos_model)
    #             #     model_resi_list.append(resi_model)
    #             #     distance_list.append(distance)
    #             #     # increment all indices
    #             #     i += 1
    #             #     j += 1
    #             #     k += 1
    #             #     index += 1
    #
    #             # append calculated data to each separate list
    #             index_list.append(index)
    #             ref_chain_list.append(chain_ref)
    #             ref_pos_list.append(pos_ref)
    #             ref_resi_list.append(resi_ref)
    #             model_chain_list.append(chain_model)
    #             model_pos_list.append(pos_model)
    #             model_resi_list.append(resi_model)
    #             distance_list.append(distance)
    #             # increment all indices
    #             i += 1
    #             j += 1
    #             k += 1
    #             index += 1
    #
    #     index_array: np.ndarray = np.array(index_list)
    #     ref_chain_array: np.ndarray = np.array(ref_chain_list)
    #     ref_pos_array: np.ndarray = np.array(ref_pos_list)
    #     ref_resi_array: np.ndarray = np.array(ref_resi_list)
    #     model_chain_array: np.ndarray = np.array(model_chain_list)
    #     model_pos_array: np.ndarray = np.array(model_pos_list)
    #     model_resi_array: np.ndarray = np.array(model_resi_list)
    #     distance_array: np.ndarray = np.array(distance_list)
    #
    #     result_hashtable: dict[str, np.ndarray] = {'index': index_array,
    #                                                'ref_chain':
    #                                                    ref_chain_array,
    #                                                'ref_pos':
    #                                                    ref_pos_array,
    #                                                'ref_resi':
    #                                                    ref_resi_array,
    #                                                'model_chain':
    #                                                    model_chain_array,
    #                                                'model_pos':
    #                                                    model_pos_array,
    #                                                'model_resi':
    #                                                    model_resi_array,
    #                                                'distance':
    #                                                    distance_array}
    #
    #     # return the hast table with all results
    #     return result_hashtable
    #
    #     # general exception if something goes wrong
    #     # except:
    #     #     print(f"Error. The index i is {i}, j is {j} and k is {k}.")
    #     #     print(f"The length of the alignment is "
    #     #           f"{int(align.get_alignment_length())}.")
    #     #     print(f"{align[0, i]},{align[1, i]}")
    #     #     print(result_hashtable)
    #
    # def export_distance_between_ca_atoms(
    #         self, distance_results: dict[str, np.ndarray]) -> None:
    #     """This function exports the results from the function
    #     calculate_distance_between_ca_atoms as CSV file
    #
    #     Args:
    #         distance_results (dict[str, np.ndarray]):
    #             hash table which contains the full information of the distance
    #             calculations
    #     """
    #     if not safeguard.Safeguard.check_if_dict_is_empty(distance_results):
    #         raise ValueError("The dictionary of the distances is empty!")
    #     distance_data = {
    #         "index": distance_results.get("index"),
    #         "ref_chain": distance_results.get("ref_chain"),
    #         "ref_pos": distance_results.get("ref_pos"),
    #         "ref_resi": distance_results.get("ref_resi"),
    #         "model_chain": distance_results.get("model_chain"),
    #         "model_pos": distance_results.get("model_pos"),
    #         "model_resi": distance_results.get("model_resi"),
    #         "distance": distance_results.get("distance")
    #     }
    #     tmp_frame = pd.DataFrame(distance_data)
    #     if not os.path.exists(f"{self.results_dir}/distance_csv"):
    #         os.mkdir(f"{self.results_dir}/distance_csv")
    #     # exports dataframe to csv file
    #     tmp_frame.to_csv(f"{self.results_dir}/distance_csv/distances.csv")
    #
    # def get_rmsd_from_cealign(self) -> float:
    #     """This function aligns the reference to the model Protein and returns
    #     the RMSD of the structure alignment.
    #
    #     Returns:
    #         float:
    #             results generated from the cealign function from PyMol
    #     """
    #     # aligning the two proteins with the PyMol command cealign
    #     results = cmd.cealign(
    #         target=self.ref_obj.molecule_object,
    #         mobile=self.model_obj.molecule_object,
    #         quiet=1)
    #     value_rmsd: float = results.get('RMSD')
    #     return value_rmsd
    #
    # @staticmethod
    # def print_results_of_align_protein_pair(result: tuple[float, int],
    #                                         model_name: str) -> None:
    #     """This function prints the results of the function alignProteinPair in
    #     a specific format.
    #
    #     Args:
    #         result (pd.DataFrame):
    #             dataframe which holds the RMSD value and
    #             the number of aligned amino acids
    #         model_name (str):
    #             the original name of the aligned model
    #
    #     Raises:
    #         ValueError:
    #             Exception is raised if the dataframe is NULL, the size of the
    #             dataframe is zero or the model_name variable is empty.
    #     """
    #     # argument test
    #     if result is None:
    #         raise ValueError("The DataFrame is None.")
    #     if len(result) == 0:
    #         raise ValueError("The DataFrame has a size of zero.")
    #     if model_name == "":
    #         raise ValueError("The model name is empty.")
    #
    #     print("----------------------------------")
    #     print(model_name)
    #     print(result)
    #     print("----------------------------------\n")
    #
    # # def export_multiple_align_protein_pair_results(self,
    # #                                                align_results: LinkedList,
    # #                                                model_names: tuple,
    # #                                                filename: str) -> None:
    # #     """This function exports the rmsd and aligned amino acids from the
    # #     function align_protein_pair as CSV file.
    # #
    # #     Args:
    # #         align_results (LinkedList):
    # #             Linked list of tuples from align_protein_pair
    # #         model_names (tuple):
    # #             tuple which contains all model names from the structure
    # #             alignments
    # #         filename (str):
    # #             name of the exported csv file
    # #     """
    # #     export_data_array = np.ndarray((len(align_results), 2), dtype=float)
    # #
    # #     i = 0
    # #     for tmp_tuple in align_results:
    # #         export_data_array[i, 0] = tmp_tuple[0]
    # #         export_data_array[i, 1] = tmp_tuple[1]
    # #         i += 1
    # #
    # #     index = model_names
    # #
    # #     # check if path exists where the data will be exported,
    # #     # if not the directory will be created
    # #     if not os.path.exists(f"{self.results_dir}/rmsds"):
    # #         os.mkdir(f"{self.results_dir}/rmsds")
    # #
    # #     df_export = pd.DataFrame(export_data_array,
    # #                              columns=["RMSD", "Aligned AA"], index=index)
    # #     df_export.to_csv(f"{self.results_dir}/rmsds/{filename}.csv")
    #
    # # def export_distance_as_csv_file(self, results: dict[str, np.ndarray],
    # #                                 filename: str) -> None:
    # #     """This function exports the results dataframe from
    # #     calculate_distance_between_ca_atoms() as csv format.
    # #
    # #     Note:
    # #         The csv file will be saved under the relative path
    # #         (if export_data_dir = "data/results"):
    # #         ``data/results/distances``
    # #
    # #         The file name (filename) MUST NOT have the file extension .csv!
    # #
    # #     Args:
    # #         results (dict[str, np.ndarray]):
    # #             hash table which contains all results from
    # #             calculate_distance_between_ca_atoms()
    # #         filename (str):
    # #             name of the csv file
    # #
    # #     Example:
    # #         The ``filename`` variable should be as follows::
    # #
    # #             Yes:
    # #                 # This is a correct file name, because it has NO
    # #                 # file extension
    # #                 f"results_Bmp2_0_CA_Distance"
    # #
    # #             No:
    # #                 # This file name is wrong, because it has a
    # #                 # trailing .csv which is not compatible with
    # #                 # the function
    # #                 f"results_Bmp2_0_CA_Distance.csv"
    # #     """
    # #     # check if path exists where the data will be exported,
    # #     # if not the directory will be created
    # #     if not os.path.exists(f"{self.results_dir}/distances"):
    # #         os.mkdir(f"{self.results_dir}/distances")
    # #
    # #     # convert hash table in pandas DataFrame
    # #     df_results: pd.DataFrame = pd.DataFrame.from_dict(results)
    # #     # export data to csv file
    # #     df_results.to_csv(f"{self.results_dir}/distances/{filename}.csv")
    #
    # def take_image_of_protein_pair(self,
    #                                color_ref: str,
    #                                color_model: str,
    #                                representation: str,
    #                                view_point,
    #                                filename: str,
    #                                cycles: int,
    #                                cutoff: float,
    #                                selection: str = "",
    #                                ray_shadows: bool = False,
    #                                opaque_background: int = 0) -> None:
    #     """This function takes an image of the whole Protein/Protein pair.
    #
    #     Note:
    #         The png file will be saved under the relative path
    #         (if export_data_dir = "data/results"):
    #         ``data/results/images``
    #
    #     Args:
    #         color_ref (str):
    #             defines color for the reference Protein
    #         color_model (str):
    #             defines color for the model Protein
    #         representation (str):
    #             defines the type of molecular representation
    #             like cartoon or ribbon
    #         view_point:
    #             position of camera
    #         filename (str):
    #             name of the png image file
    #         cycles (int):
    #             defines the number of refinement cycles
    #         cutoff (float):
    #             defines the value when residues are not aligned
    #         selection (str, optional):
    #             the atoms which MUST NOT displayed in the image
    #         ray_shadows (bool, optional):
    #             false if no shadows, true if shadows should be displayed
    #         opaque_background (int, optional):
    #             0 for a transparent background and 1 for a white background
    #
    #     Raises:
    #         ValueError: If opaque_background is not 0 or 1.
    #     """
    #     # argument test
    #     # if opaque_background != 0 or opaque_background != 1:
    #     #     raise Exception(
    #     #         "ValueError: The value for opaque_background MUST be 0 or 1!")
    #
    #     # determine the option for ray_shadows
    #     opt_ray_shadows: str
    #     if not ray_shadows:
    #         opt_ray_shadows = "off"
    #     else:
    #         opt_ray_shadows = "on"
    #
    #     # begin of image creation
    #     cmd.reinitialize()
    #     self.load_protein_pair()
    #     self.color_protein_pair(color_ref, color_model)
    #     cmd.show(representation)
    #
    #     if selection != "":
    #         cmd.hide(representation, selection)
    #
    #     self.align_protein_pair(cycle_number=cycles,
    #                             cutoff_value=cutoff)
    #
    #     cmd.orient()
    #     cmd.center()
    #     cmd.set_view(view_point)
    #     cmd.set("ray_trace_mode", 3)
    #     cmd.bg_color("white")
    #     cmd.set("antialias", 2)
    #     cmd.set("ray_shadows", opt_ray_shadows)
    #     cmd.set('ray_opaque_background', opaque_background)
    #     cmd.ray(2400, 2400, renderer=0)
    #
    #     if not os.path.exists(f"{self.results_dir}/images"):
    #         os.mkdir(f"{self.results_dir}/images")
    #
    #     # save image as 300 dpi png image
    #     cmd.png(f'{self.results_dir}/images/{filename}.png', dpi=300)

    def serialize_protein_pair(self, xml_protein_pairs_element: ElementTree.Element):
        """This function serialize the protein pair object

        """
        tmp_protein_pair = ElementTree.SubElement(xml_protein_pairs_element, element_names.PROTEIN_PAIR)
        tmp_protein_pair.set(attribute_names.PROTEIN_PAIR_NAME, str(self.name))
        tmp_protein_pair.set(attribute_names.PROTEIN_PAIR_PROT_1_MOLECULE_OBJECT,
                             str(self.protein_1.get_molecule_object()))
        tmp_protein_pair.set(attribute_names.PROTEIN_PAIR_PROT_1_ID, str(self.protein_1.get_id()))
        tmp_protein_pair.set(attribute_names.PROTEIN_PAIR_PROT_2_MOLECULE_OBJECT,
                             str(self.protein_2.get_molecule_object()))
        tmp_protein_pair.set(attribute_names.PROTEIN_PAIR_PROT_2_ID, str(self.protein_2.get_id()))
        tmp_session_data = ElementTree.SubElement(tmp_protein_pair, element_names.PROTEIN_PAIR_SESSION)
        tmp_session_data.set(attribute_names.PROTEIN_PAIR_SESSION, self.pymol_session)
        if self.distance_analysis is not None:
            self.distance_analysis.serialize_distance_analysis(tmp_protein_pair)
        # tmp_protein.set(attribute_names.ID, str(self._id))
        # tmp_protein.set(attribute_names.PROTEIN_MOLECULE_OBJECT, self._pymol_molecule_object)
        # tmp_protein.set(attribute_names.PROTEIN_SELECTION, self.pymol_selection.selection_string)
        # tmp_protein.append(bio_data.convert_pdb_data_list_to_xml_string(self._pdb_data))
        # tmp_session_data = ElementTree.SubElement(tmp_protein, element_names.PROTEIN_SESSION)
        # tmp_session_data.set(attribute_names.PROTEIN_SESSION, self.pymol_session)
        #
        # protein_pair_serializer = filesystem_io.ObjectSerializer(self, self.protein_pair_subdirs.get(pyssa_keys.PROTEIN_PAIR_OBJECTS_SUBDIR), self.name)
        # protein_pair_dict = {
        #     "prot_1_molecule_object": self.protein_1.get_molecule_object(),
        #     "prot_1_import_data_dir": str(self.protein_1.pdb_filepath.get_filepath()),
        #     "prot_1_export_data_dir": str(self.protein_1.export_dirname),
        #     "prot_1_filename": str(self.protein_1.pdb_filepath.get_filename()),
        #     "prot_1_selection": self.protein_1.pymol_selection.selection_string,
        #     "prot_1_chains": protein_util.get_chains_as_list_of_tuples(self.protein_1.chains),
        #     "prot_2_molecule_object": self.protein_2.get_molecule_object(),
        #     "prot_2_import_data_dir": str(self.protein_2.pdb_filepath.get_filepath()),
        #     "prot_2_export_data_dir": str(self.protein_2.export_dirname),
        #     "prot_2_filename": str(self.protein_2.pdb_filepath.get_filename()),
        #     "prot_2_selection": self.protein_2.pymol_selection.selection_string,
        #     "prot_2_chains": protein_util.get_chains_as_list_of_tuples(self.protein_2.chains),
        #     "export_dirname": str(self.SCRATCH_DIR),
        #     "pymol_session_file": str(self.pymol_session_filepath.get_filepath()),
        #     "name": self.name,
        #     "protein_pairs_dirname": str(self.protein_pair_subdirs.get(pyssa_keys.PROTEIN_PAIR_SUBDIR)),
        # }
        # protein_pair_serializer.set_custom_object_dict(protein_pair_dict)
        # protein_pair_serializer.serialize_object()

    @staticmethod
    def deserialize_protein_pair(protein_obj_json_file: path_util.FilePath):
        """This function constructs the protein pair object from
        the json file

        Returns:
            two complete protein objects and a protein pair object deserialized from a json file
        """
        return filesystem_io.ObjectDeserializer(protein_obj_json_file.get_dirname(), protein_obj_json_file.get_filename()).deserialize_protein_pair()

    def create_plain_text_memory_mirror(self):
        mirror = [
            self.protein_1.get_molecule_object(),
            str(self.protein_1.pdb_filepath.get_filepath()),
            str(self.protein_1.export_dirname),
            str(self.protein_1.pdb_filepath.get_filename()),
            self.protein_1.pymol_selection.selection_string,
            protein_util.get_chains_as_list_of_tuples(self.protein_1.chains),
            self.protein_2.get_molecule_object(),
            str(self.protein_2.pdb_filepath.get_filepath()),
            str(self.protein_2.export_dirname),
            str(self.protein_2.pdb_filepath.get_filename()),
            self.protein_2.pymol_selection.selection_string,
            protein_util.get_chains_as_list_of_tuples(self.protein_2.chains),
            str(self.SCRATCH_DIR),
            str(self.pymol_session_filepath.get_filepath()),
            self.name,
            str(self.protein_pair_subdirs.get(pyssa_keys.PROTEIN_PAIR_SUBDIR)),
        ]
        return mirror