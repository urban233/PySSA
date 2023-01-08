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
import json
import os
import pathlib

import pymol
import pandas as pd
import numpy as np
from pymol import cmd
from Bio import AlignIO
from pyssa.pymol_protein_tools import protein


class ProteinPair:
    """This class consists of two Protein objects. It is used to have a better workflow for
    the analysis.

    """

    def __init__(self, reference_obj: protein.Protein, model_obj: protein.Protein,
                 results_dir: pathlib.Path) -> None:
        """Constructor.

        Args:
            reference_obj (core.protein.Protein):
                reference Protein object
            model_obj (core.Protein):
                model Protein object
            results_dir (str):
                directory where all results will be stored.
                All subdirectories like ``images``, ``alignment_files``
                and ``distances`` will be created automatically.
                The results_dir will then function as parent directory.

        Raises:
            NotADirectoryError: If directory not found.
        """
        self.ref_obj: protein.Protein = reference_obj
        self.model_obj: protein.Protein = model_obj
        self.results_dir: pathlib.Path = results_dir
        self.name = "generic"
        self.cutoff = 0

        # argument test
        if not os.path.exists(f"{results_dir}"):
            raise NotADirectoryError

    def load_protein_pair(self) -> None:
        """This function loads to proteins with the `load`_ command from
        PyMOL. The first one is the reference which comes from the constructor.
        The second one is the model. This text will be automatically synced.

        Raises:
            pymol.CmdException:
                Exception is raised if load command fails.

        .. _load:
            https://pymolwiki.org/index.php/Load
        """
        # check if the protein is a duplicate
        if self.ref_obj.molecule_object.find("_1"):
            tmp_ref_molecule_object = self.ref_obj.molecule_object.replace("_1", "")
        else:
            tmp_ref_molecule_object = self.ref_obj.molecule_object
        if self.model_obj.molecule_object.find("_2"):
            tmp_model_molecule_object = self.model_obj.molecule_object.replace("_2", "")
        else:
            tmp_model_molecule_object = self.model_obj.molecule_object
        # loading the reference in the active PyMol session
        cmd.load(f"{self.ref_obj.filepath}/"
                 f"{tmp_ref_molecule_object}.pdb", object=self.ref_obj.molecule_object)
        # loading the model in the active PyMol session
        cmd.load(f"{self.model_obj.filepath}/"
                 f"{tmp_model_molecule_object}.pdb", object=self.model_obj.molecule_object)
        print(cmd.get_object_list())

    def color_protein_pair(self, color_ref="green", color_model="blue") -> None:
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

        .. _color values:
            https://pymolwiki.org/index.php/Color_Values
        """
        # argument test
        # checks if either the reference or the model is an actual object
        # in the memory
        # if cmd.get_model(self.ref_obj.molecule_object) is None \
        #         or cmd.get_model(self.model_obj.molecule_object) is None:
        #     raise pymol.CmdException(f"Either the reference or the model is "
        #                              f"not in the pymol session as object.")
        # # checks if both the reference and the model are actual objects
        # # in the memory
        # if cmd.get_model(self.ref_obj.molecule_object) is None \
        #         and cmd.get_model(self.model_obj.molecule_object) is None:
        #     raise pymol.CmdException(f"Both, the reference and the model are "
        #                              f"not in the pymol session as objects.")
        # actual color cmd command
        cmd.color(color_ref, self.ref_obj.molecule_object)
        cmd.color(color_model, self.model_obj.molecule_object)

    def align_protein_pair(self, cycle_number: int, cutoff_value: float,
                           alignment_filename: str = None) -> tuple[float, int]:
        """This function aligns the model to the reference Protein, with the
        `align`_ command from PyMOL.

        Note:
            Before this function can be used, the load_protein_pair
            needs to be executed. Moreover, two selections are needed
            which have to be set through the set_selection function,
            for each Protein.

            If an alignment file should be saved, it will be stored
            under the relative path (if export_data_dir =
            "data/results"): ``data/results/alignment_files``

            The filename of the alignment file (alignment_filename)
            MUST NOT have the file extension .aln!

        Args:
            cycle_number (int):
                defines the number of refinement cycles
            cutoff_value (float):
                defines the value when residues are not aligned
            alignment_filename (str, optional):
                string which functions as filename for the align object

        Returns:
            tuple of the results

        Raises:
            ValueError:
                Exception is raised if
                the align_reference or align_model string is empty or
                the cycle_number or cutoff_value are equal or below zero.

        Example:
            The ``alignment_filename`` variable should be as follows::

                Yes:
                    # This is a correct file name, because it has NO
                    # file extension
                    f"aln_results_Bmp2_0_CA"

                No:
                    # This file name is wrong, because it has a
                    # trailing .csv which is not compatible with
                    # the function
                    f"aln_results_Bmp2_0_CA.aln"

        .. _align:
            https://pymolwiki.org/index.php/Align
        """
        # argument test
        if self.ref_obj.selection is None or self.model_obj.selection is None:
            raise ValueError("Either the reference or the model "
                             "has an empty selection.")
        if self.ref_obj.selection is None and self.model_obj.selection is None:
            raise ValueError("Both, the reference and the model "
                             "have an empty selection.")
        if cycle_number < 0:
            raise ValueError("Number of cycles must be greater or equal than zero.")
        if cutoff_value <= 0:
            raise ValueError("The cutoff needs to be greater than zero.")

        # This block runs if an alignObject should be created.
        if alignment_filename is not None:
            results = cmd.align(target=self.ref_obj.selection,
                                mobile=self.model_obj.selection,
                                object=alignment_filename,
                                cycles=cycle_number,
                                cutoff=cutoff_value,
                                quiet=0)

            if not os.path.exists(f"{self.results_dir}/alignment_files"):
                os.mkdir(f"{self.results_dir}/alignment_files")

            # save the align object from pymol as alignment file
            cmd.save(f"{self.results_dir}/alignment_files/"
                     f"{alignment_filename}.aln")

            return results[0], results[1]

        # This block runs if no alignObject should be created.
        else:
            results = cmd.align(target=self.ref_obj.selection,
                                mobile=self.model_obj.selection,
                                cycles=cycle_number,
                                cutoff=cutoff_value,
                                quiet=1)

            return results[0], results[1]
            # tup_results = (results[0], results[1])
            # return tup_results

    def calculate_distance_between_ca_atoms(self, alignment_filename: str,
                                            cutoff: float = 20.0) -> dict[str, np.ndarray]:
        # def calculate_distance_between_ca_atoms(self, alignment_filename: str,
        #                                         cutoff: float = 20.0) -> Dict[str, np.ndarray]:
        """This function calculates the distances between the aligned alpha-C
        atoms. It uses the alignment file generated by the alignProteinPair
        function as base.

        Note:
            The filename of the alignment file (alignment_filename) MUST NOT
            have the file extension .aln!

        Args:
            alignment_filename (str):
                filename of the alignment file, generated from the
                alignProteinPair function
            cutoff (float, optional):
                all distances which are above the cutoff will NOT be in the
                results dataframe.
                The default value is 20 angstrom.

        Returns:
            result_hashtable (dict[str, np.ndarray]):
                a hash table which has the following keys:
                'index', 'ref_chain', 'ref_pos', 'ref_resi', 'model_chain',
                'model_pos', 'model_resi', 'distance'

        Example:
            The ``alignment_filename`` variable should be as follows::

                Yes:
                    # This is a correct file name, because it has NO
                    # file extension
                    f"aln_results_Bmp2_0_CA"

                No:
                    # This file name is wrong, because it has a
                    # trailing .csv which is not compatible with
                    # the function
                    f"aln_results_Bmp2_0_CA.aln"
        """
        # argument test
        try:
            file = open(f"{self.results_dir}/alignment_files/"
                        f"{alignment_filename}.aln", "r")
            file.close()
        except FileNotFoundError:
            print(f"File not found, in {self.results_dir}.")

        cmd.create(f"{self.ref_obj.molecule_object}_CA",
                   f"/{self.ref_obj.molecule_object}////CA")
        ref_ca_obj = cmd.get_model(f"{self.ref_obj.molecule_object}_CA")

        cmd.create(f"{self.model_obj.molecule_object}_CA",
                   f"/{self.model_obj.molecule_object}////CA")
        model_ca_obj = cmd.get_model(f"{self.model_obj.molecule_object}_CA")

        # read in alignment file from alignProteinPair function
        align = AlignIO.read(f"{self.results_dir}/alignment_files/"
                             f"{alignment_filename}.aln", "clustal")

        index_list: [int] = []
        ref_chain_list: [str] = []
        ref_pos_list: [int] = []
        ref_resi_list: [str] = []
        model_chain_list: [str] = []
        model_pos_list: [int] = []
        model_resi_list: [str] = []
        distance_list: [float] = []

        i = 0  # i for alignment file
        j = 0  # j for reference
        k = 0  # k for model
        index = 0
        #try:
        while i < int(align.get_alignment_length()):
            # gets executed if the reference contains a "-" in the alignment
            if align[0, i] == "-":
                i += 1
                k += 1
                j = j
            # gets executed if the model contains a "-" in the alignment
            elif align[1, i] == "-":
                i += 1
                j += 1
                k = k
            # gets executed if no "-" is found in the alignment
            else:
                # create var for chain, position and residue name for
                # the reference
                chain_ref: str = ref_ca_obj.atom[j].chain
                pos_ref: str = ref_ca_obj.atom[j].resi
                resi_ref: str = ref_ca_obj.atom[j].resn

                # create var for chain, position and residue name for
                # the model
                chain_model: str = model_ca_obj.atom[k].chain
                pos_model: str = model_ca_obj.atom[k].resi
                resi_model: str = model_ca_obj.atom[k].resn

                # calculate the distance between the alpha-C atoms
                atom1 = f"/{self.ref_obj.molecule_object}_CA//{chain_ref}/{pos_ref}/"
                atom2 = f"/{self.model_obj.molecule_object}_CA//{chain_model}/{pos_model}/"
                distance = round(cmd.get_distance(atom1, atom2), 2)

                # TODO: should there be a cutoff?
                # gets executed if the distance is greater than the
                # pre-defined cutoff
                # if distance > cutoff:
                #     i += 1
                #     j += 1
                #     k += 1
                # else:
                #     # append calculated data to each separate list
                #     index_list.append(index)
                #     ref_chain_list.append(chain_ref)
                #     ref_pos_list.append(pos_ref)
                #     ref_resi_list.append(resi_ref)
                #     model_chain_list.append(chain_model)
                #     model_pos_list.append(pos_model)
                #     model_resi_list.append(resi_model)
                #     distance_list.append(distance)
                #     # increment all indices
                #     i += 1
                #     j += 1
                #     k += 1
                #     index += 1

                # append calculated data to each separate list
                index_list.append(index)
                ref_chain_list.append(chain_ref)
                ref_pos_list.append(pos_ref)
                ref_resi_list.append(resi_ref)
                model_chain_list.append(chain_model)
                model_pos_list.append(pos_model)
                model_resi_list.append(resi_model)
                distance_list.append(distance)
                # increment all indices
                i += 1
                j += 1
                k += 1
                index += 1

        index_array: np.ndarray = np.array(index_list)
        ref_chain_array: np.ndarray = np.array(ref_chain_list)
        ref_pos_array: np.ndarray = np.array(ref_pos_list)
        ref_resi_array: np.ndarray = np.array(ref_resi_list)
        model_chain_array: np.ndarray = np.array(model_chain_list)
        model_pos_array: np.ndarray = np.array(model_pos_list)
        model_resi_array: np.ndarray = np.array(model_resi_list)
        distance_array: np.ndarray = np.array(distance_list)

        result_hashtable: dict[str, np.ndarray] = {'index': index_array,
                                                   'ref_chain':
                                                       ref_chain_array,
                                                   'ref_pos':
                                                       ref_pos_array,
                                                   'ref_resi':
                                                       ref_resi_array,
                                                   'model_chain':
                                                       model_chain_array,
                                                   'model_pos':
                                                       model_pos_array,
                                                   'model_resi':
                                                       model_resi_array,
                                                   'distance':
                                                       distance_array}

        # return the hast table with all results
        return result_hashtable

        # general exception if something goes wrong
        # except:
        #     print(f"Error. The index i is {i}, j is {j} and k is {k}.")
        #     print(f"The length of the alignment is "
        #           f"{int(align.get_alignment_length())}.")
        #     print(f"{align[0, i]},{align[1, i]}")
        #     print(result_hashtable)

    def export_distance_between_ca_atoms(
            self, distance_results: dict[str, np.ndarray]) -> None:
        """This function exports the results from the function
        calculate_distance_between_ca_atoms as CSV file

        Args:
            distance_results (dict[str, np.ndarray]):
                hash table which contains the full information of the distance
                calculations
        """
        distance_data = {
            "index": distance_results.get("index"),
            "ref_chain": distance_results.get("ref_chain"),
            "ref_pos": distance_results.get("ref_pos"),
            "ref_resi": distance_results.get("ref_resi"),
            "model_chain": distance_results.get("model_chain"),
            "model_pos": distance_results.get("model_pos"),
            "model_resi": distance_results.get("model_resi"),
            "distance": distance_results.get("distance")
        }
        tmp_frame = pd.DataFrame(distance_data)

        # check if path exists where the data will be exported,
        # if not the directory will be created
        if not os.path.exists(f"{self.results_dir}/distance_csv"):
            os.mkdir(f"{self.results_dir}/distance_csv")
        # exports dataframe to csv file
        tmp_frame.to_csv(f"{self.results_dir}/distance_csv/distances.csv")

    def save_session_of_protein_pair(self, filename: str) -> None:
        """This function saves the pymol session of the Protein pair.

        Note:
            The pse file will be saved under the relative path
            (if export_data_dir = "data/results"):
            ``data/results/sessions``

            The file name (filename) MUST NOT have the file extension .pse!

        Args:
            filename (str):
                name of the session file

        Example:
            The ``filename`` variable should be as follows::

                Yes:
                    # This is a correct file name, because it has NO
                    # file extension
                    f"Bmp2_0_CA_session"

                No:
                    # This file name is wrong, because it has a
                    # trailing .pse which is not compatible with
                    # the function
                    f"Bmp2_0_CA_session.pse"
        """
        if not os.path.exists(f"{self.results_dir}/sessions"):
            os.mkdir(f"{self.results_dir}/sessions")
        cmd.save(f"{self.results_dir}/sessions/{filename}.pse")

    def get_rmsd_from_cealign(self) -> float:
        """This function aligns the reference to the model Protein and returns
        the RMSD of the structure alignment.

        Returns:
            float:
                results generated from the cealign function from PyMol
        """
        # aligning the two proteins with the PyMol command cealign
        results = cmd.cealign(
            target=self.ref_obj.molecule_object,
            mobile=self.model_obj.molecule_object,
            quiet=1)
        value_rmsd: float = results.get('RMSD')
        return value_rmsd

    @staticmethod
    def print_results_of_align_protein_pair(result: tuple[float, int],
                                            model_name: str) -> None:
        """This function prints the results of the function alignProteinPair in
        a specific format.

        Args:
            result (pd.DataFrame):
                dataframe which holds the RMSD value and
                the number of aligned amino acids
            model_name (str):
                the original name of the aligned model

        Raises:
            ValueError:
                Exception is raised if the dataframe is NULL, the size of the
                dataframe is zero or the model_name variable is empty.
        """
        # argument test
        if result is None:
            raise ValueError("The DataFrame is None.")
        if len(result) == 0:
            raise ValueError("The DataFrame has a size of zero.")
        if model_name == "":
            raise ValueError("The model name is empty.")

        print("----------------------------------")
        print(model_name)
        print(result)
        print("----------------------------------\n")

    # def export_multiple_align_protein_pair_results(self,
    #                                                align_results: LinkedList,
    #                                                model_names: tuple,
    #                                                filename: str) -> None:
    #     """This function exports the rmsd and aligned amino acids from the
    #     function align_protein_pair as CSV file.
    #
    #     Args:
    #         align_results (LinkedList):
    #             Linked list of tuples from align_protein_pair
    #         model_names (tuple):
    #             tuple which contains all model names from the structure
    #             alignments
    #         filename (str):
    #             name of the exported csv file
    #     """
    #     export_data_array = np.ndarray((len(align_results), 2), dtype=float)
    #
    #     i = 0
    #     for tmp_tuple in align_results:
    #         export_data_array[i, 0] = tmp_tuple[0]
    #         export_data_array[i, 1] = tmp_tuple[1]
    #         i += 1
    #
    #     index = model_names
    #
    #     # check if path exists where the data will be exported,
    #     # if not the directory will be created
    #     if not os.path.exists(f"{self.results_dir}/rmsds"):
    #         os.mkdir(f"{self.results_dir}/rmsds")
    #
    #     df_export = pd.DataFrame(export_data_array,
    #                              columns=["RMSD", "Aligned AA"], index=index)
    #     df_export.to_csv(f"{self.results_dir}/rmsds/{filename}.csv")

    # def export_distance_as_csv_file(self, results: dict[str, np.ndarray],
    #                                 filename: str) -> None:
    #     """This function exports the results dataframe from
    #     calculate_distance_between_ca_atoms() as csv format.
    #
    #     Note:
    #         The csv file will be saved under the relative path
    #         (if export_data_dir = "data/results"):
    #         ``data/results/distances``
    #
    #         The file name (filename) MUST NOT have the file extension .csv!
    #
    #     Args:
    #         results (dict[str, np.ndarray]):
    #             hash table which contains all results from
    #             calculate_distance_between_ca_atoms()
    #         filename (str):
    #             name of the csv file
    #
    #     Example:
    #         The ``filename`` variable should be as follows::
    #
    #             Yes:
    #                 # This is a correct file name, because it has NO
    #                 # file extension
    #                 f"results_Bmp2_0_CA_Distance"
    #
    #             No:
    #                 # This file name is wrong, because it has a
    #                 # trailing .csv which is not compatible with
    #                 # the function
    #                 f"results_Bmp2_0_CA_Distance.csv"
    #     """
    #     # check if path exists where the data will be exported,
    #     # if not the directory will be created
    #     if not os.path.exists(f"{self.results_dir}/distances"):
    #         os.mkdir(f"{self.results_dir}/distances")
    #
    #     # convert hash table in pandas DataFrame
    #     df_results: pd.DataFrame = pd.DataFrame.from_dict(results)
    #     # export data to csv file
    #     df_results.to_csv(f"{self.results_dir}/distances/{filename}.csv")

    def take_image_of_protein_pair(self,
                                   color_ref: str,
                                   color_model: str,
                                   representation: str,
                                   view_point,
                                   filename: str,
                                   cycles: int,
                                   cutoff: float,
                                   selection: str = "",
                                   ray_shadows: bool = False,
                                   opaque_background: int = 0) -> None:
        """This function takes an image of the whole Protein/Protein pair.

        Note:
            The png file will be saved under the relative path
            (if export_data_dir = "data/results"):
            ``data/results/images``

        Args:
            color_ref (str):
                defines color for the reference Protein
            color_model (str):
                defines color for the model Protein
            representation (str):
                defines the type of molecular representation
                like cartoon or ribbon
            view_point:
                position of camera
            filename (str):
                name of the png image file
            cycles (int):
                defines the number of refinement cycles
            cutoff (float):
                defines the value when residues are not aligned
            selection (str, optional):
                the atoms which MUST NOT displayed in the image
            ray_shadows (bool, optional):
                false if no shadows, true if shadows should be displayed
            opaque_background (int, optional):
                0 for a transparent background and 1 for a white background

        Raises:
            ValueError: If opaque_background is not 0 or 1.
        """
        # argument test
        # if opaque_background != 0 or opaque_background != 1:
        #     raise Exception(
        #         "ValueError: The value for opaque_background MUST be 0 or 1!")

        # determine the option for ray_shadows
        opt_ray_shadows: str
        if not ray_shadows:
            opt_ray_shadows = "off"
        else:
            opt_ray_shadows = "on"

        # begin of image creation
        cmd.reinitialize()
        self.load_protein_pair()
        self.color_protein_pair(color_ref, color_model)
        cmd.show(representation)

        if selection != "":
            cmd.hide(representation, selection)

        self.align_protein_pair(cycle_number=cycles,
                                cutoff_value=cutoff)

        cmd.orient()
        cmd.center()
        cmd.set_view(view_point)
        cmd.set("ray_trace_mode", 3)
        cmd.bg_color("white")
        cmd.set("antialias", 2)
        cmd.set("ray_shadows", opt_ray_shadows)
        cmd.set('ray_opaque_background', opaque_background)
        cmd.ray(2400, 2400, renderer=0)

        # check if path exists where the data will be exported,
        # if not the directory will be created
        if not os.path.exists(f"{self.results_dir}/images"):
            os.mkdir(f"{self.results_dir}/images")

        # save image as 300 dpi png image
        cmd.png(f'{self.results_dir}/images/{filename}.png', dpi=300)

    def serialize_protein_pair(self, filepath) -> None:
        """This function serialize the protein pair object

        """
        if not os.path.exists(filepath):
            print(f"The filepath: {filepath} does not exists!")
            return

        protein_structures_dict = {
            "prot_1_molecule_object": self.ref_obj.molecule_object,
            "prot_1_import_data_dir": str(self.ref_obj.filepath),
            "prot_1_export_data_dir": str(self.ref_obj.export_data_dir),
            "prot_1_selection": self.ref_obj.selection,
            "prot_1_sequence": self.ref_obj.sequence,
            "prot_1_chains": self.ref_obj.chains,
            "prot_2_molecule_object": self.model_obj.molecule_object,
            "prot_2_import_data_dir": str(self.model_obj.filepath),
            "prot_2_export_data_dir": str(self.model_obj.export_data_dir),
            "prot_2_selection": self.model_obj.selection,
            "prot_2_sequence": self.model_obj.sequence,
            "prot_2_chains": self.model_obj.chains,
            "results_dir": str(self.results_dir),
            "cutoff": self.cutoff,
            "name": self.name,
        }
        # protein_structures_dict = {
        #     "ref_obj": self.ref_obj,
        #     "prot_1_molecule_object": protein_structure_1.molecule_object,
        #     "prot_1_import_data_dir": protein_structure_1.filepath,
        #     "prot_1_export_data_dir": protein_structure_1.export_data_dir,
        #     "prot_1_selection": protein_structure_1.selection,
        #     "prot_1_sequence": protein_structure_1.sequence,
        #     "prot_1_chains": protein_structure_1.chains,
        #     "model_obj": self.model_obj,
        #     "prot_2_molecule_object": protein_structure_2.molecule_object,
        #     "prot_2_import_data_dir": protein_structure_2.filepath,
        #     "prot_2_export_data_dir": protein_structure_2.export_data_dir,
        #     "prot_2_selection": protein_structure_2.selection,
        #     "prot_2_sequence": protein_structure_2.sequence,
        #     "prot_2_chains": protein_structure_2.chains,
        #     "results_dir": self.results_dir,
        #     "cutoff": self.cutoff,
        #     "name": self.name,
        # }

        protein_file = open(pathlib.Path(f"{filepath}/{self.name}.json"), "w", encoding="utf-8")
        json.dump(protein_structures_dict, protein_file, indent=4)

    @staticmethod
    def deserialize_protein_pair(protein_obj_json_file):
        """This function constructs the protein pair object from
        the json file

        Returns:
            two complete protein objects and a protein pair object deserialized from a json file
        """
        try:
            protein_obj_file = open(protein_obj_json_file, "r", encoding="utf-8")
        except FileNotFoundError:
            print(f"There is no valid protein pair json file under: {protein_obj_json_file}")
            return
        protein_dict = json.load(protein_obj_file)

        tmp_protein_pair = ProteinPair(protein_dict.get("prot_1_molecule_object"),
                                       protein_dict.get("prot_2_molecule_object"),
                                       pathlib.Path(protein_dict.get("results_dir")),
                                       )
        tmp_protein_pair.cutoff = protein_dict.get("cutoff")
        tmp_protein_pair.name = protein_dict.get("name")

        if protein_dict.get("prot_1_export_data_dir") == "None":
            export_data_dir = None
        else:
            export_data_dir = protein_dict.get("prot_1_export_data_dir")
        tmp_protein_1 = protein.Protein(protein_dict.get("prot_1_molecule_object"),
                                        protein_dict.get("prot_1_import_data_dir"),
                                        export_data_dir=export_data_dir,
                                        )
        tmp_protein_1.set_sequence(protein_dict.get("prot_1_sequence"))
        tmp_protein_1.set_selection(protein_dict.get("prot_1_selection"))
        tmp_protein_1.set_chains(protein_dict.get("prot_1_chains"))

        if protein_dict.get("prot_2_export_data_dir") == "None":
            export_data_dir = None
        else:
            export_data_dir = protein_dict.get("prot_2_export_data_dir")
        tmp_protein_2 = protein.Protein(protein_dict.get("prot_2_molecule_object"),
                                        protein_dict.get("prot_2_import_data_dir"),
                                        export_data_dir=export_data_dir,
                                        )
        tmp_protein_2.set_sequence(protein_dict.get("prot_2_sequence"))
        tmp_protein_2.set_selection(protein_dict.get("prot_2_selection"))
        tmp_protein_2.set_chains(protein_dict.get("prot_2_chains"))
        return tmp_protein_1, tmp_protein_2, tmp_protein_pair
