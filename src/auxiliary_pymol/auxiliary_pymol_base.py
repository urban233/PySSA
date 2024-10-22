#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
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
"""Module for the auxiliary pymol class."""
import os
import pathlib
import pymol2

import numpy as np

import local_constants
import utils

__docformat__ = "google"

import logger

module_logger = logger.setup_logger(__file__, log_path=f"{local_constants.LOG_PATH}\\auxiliary_pymol")


class AuxiliaryPyMOL:
    """A container for static methods that use an independent pymol instance, for each static method."""

    @staticmethod
    def do_distance_analysis(
        the_protein_pair_name: str,
        a_protein_1_pdb_cache_filepath: str,
        a_protein_2_pdb_cache_filepath: str,
        a_protein_1_pymol_selection_string: str,
        a_protein_2_pymol_selection_string: str,
        a_cutoff: float,
        the_cycles: int,
    ) -> tuple[tuple, str]:
        """Does the actual distance analysis in pymol.

        Args:
            the_protein_pair_name (str): The name of the protein pair.
            a_protein_1_pdb_cache_filepath (str): The filepath of the PDB cache file for protein 1.
            a_protein_2_pdb_cache_filepath (str): The filepath of the PDB cache file for protein 2.
            a_protein_1_pymol_selection_string (str): The PyMOL selection string for protein 1.
            a_protein_2_pymol_selection_string (str): The PyMOL selection string for protein 2.
            a_cutoff (float): The distance cutoff value.
            the_cycles (int): The number of cycles to perform for alignment.

        Returns:
            A tuple containing:
                A tuple of distance analysis results including the index, protein chain, position, residue label,
                  and distance values between the C-alpha atoms of the proteins.
                A string representing the alignment results including the RMSD value and aligned residue count.
        """
        # <editor-fold desc="Checks">
        if the_protein_pair_name is None or the_protein_pair_name == "":
            raise ValueError("the_protein_pair_name is either None or empty.")
        if a_protein_1_pdb_cache_filepath is None or a_protein_1_pdb_cache_filepath == "":
            raise ValueError("a_protein_1_pdb_cache_filepath is either None or empty.")
        if not os.path.exists(a_protein_1_pdb_cache_filepath):
            raise FileNotFoundError("a_protein_1_pdb_cache_filepath could not be found.")
        if a_protein_2_pdb_cache_filepath is None or a_protein_2_pdb_cache_filepath == "":
            raise ValueError("a_protein_2_pdb_cache_filepath is either None or empty.")
        if not os.path.exists(a_protein_2_pdb_cache_filepath):
            raise FileNotFoundError("a_protein_2_pdb_cache_filepath could not be found.")
        if a_protein_1_pymol_selection_string is None:
            raise ValueError("a_protein_1_pymol_selection_string is None.")
        if a_cutoff is None or a_cutoff < 0:
            raise ValueError("a_cutoff is either None or less than 0.")
        if the_cycles is None or the_cycles < 0:
            raise ValueError("the_cycles is either None or less than 0.")

        # </editor-fold>

        try:
            tmp_protein_1_name = pathlib.Path(a_protein_1_pdb_cache_filepath).name.replace(".pdb", "")
            tmp_protein_2_name = pathlib.Path(a_protein_2_pdb_cache_filepath).name.replace(".pdb", "")

            with pymol2.PyMOL() as auxiliary_pymol:
                # <editor-fold desc="Prepare session">

                # Load both proteins in a new session
                auxiliary_pymol.cmd.load(filename=str(a_protein_1_pdb_cache_filepath), object=tmp_protein_1_name)
                auxiliary_pymol.cmd.load(filename=str(a_protein_2_pdb_cache_filepath), object=tmp_protein_2_name)

                auxiliary_pymol.cmd.color(local_constants.DEFAULT_COLOR_PROTEIN_1, tmp_protein_1_name)
                auxiliary_pymol.cmd.color(local_constants.DEFAULT_COLOR_PROTEIN_2, tmp_protein_2_name)

                auxiliary_pymol.cmd.zoom("all")
                auxiliary_pymol.cmd.scene("base", action="store")
                # auxiliary_pymol.cmd.reset()
                # auxiliary_pymol.cmd.scene("base", action="update")

                auxiliary_pymol.cmd.bg_color(local_constants.PYMOL_DEFAULT_BACKGROUND_COLOR)
                auxiliary_pymol.cmd.set("valence", 0)
                auxiliary_pymol.cmd.set("scene_buttons", 0)
                auxiliary_pymol.cmd.set("ray_trace_mode", local_constants.PYMOL_DEFAULT_RAY_TRACE_MODE)
                auxiliary_pymol.cmd.set("antialias", local_constants.PYMOL_DEFAULT_ANTIALIAS)
                auxiliary_pymol.cmd.set("ambient", local_constants.PYMOL_DEFAULT_AMBIENT)
                auxiliary_pymol.cmd.set("cartoon_fancy_helices", local_constants.PYMOL_DEFAULT_FANCY_HELICES)
                auxiliary_pymol.cmd.set(
                    "cartoon_discrete_colors", local_constants.PYMOL_DEFAULT_CARTOON_DISCRETE_COLORS,
                )
                auxiliary_pymol.cmd.set("cartoon_sampling", local_constants.PYMOL_DEFAULT_CARTOON_SAMPLING)
                auxiliary_pymol.cmd.set("spec_power", local_constants.PYMOL_DEFAULT_SPEC_POWER)
                auxiliary_pymol.cmd.set("spec_reflect", local_constants.PYMOL_DEFAULT_SPEC_REFLECT)
                auxiliary_pymol.cmd.set(
                    "ray_transparency_contrast", local_constants.PYMOL_DEFAULT_RAY_TRANSPARENCY_CONTRAST,
                )
                auxiliary_pymol.cmd.set(
                    "ray_transparency_oblique", local_constants.PYMOL_DEFAULT_RAY_TRANSPARENCY_OBLIQUE,
                )  # noqa: E501
                auxiliary_pymol.cmd.set(
                    "ray_transparency_oblique_power", local_constants.PYMOL_DEFAULT_RAY_OBLIQUE_POWER,
                )
                auxiliary_pymol.cmd.set("ray_trace_color", local_constants.PYMOL_DEFAULT_RAY_COLOR)
                auxiliary_pymol.cmd.unset("depth_cue")
                # </editor-fold>

                # <editor-fold desc="Align proteins">
                tmp_align_results = auxiliary_pymol.cmd.align(
                    target=a_protein_1_pymol_selection_string,
                    mobile=a_protein_2_pymol_selection_string,
                    cutoff=a_cutoff,
                    cycles=the_cycles,
                    object="aln",
                    quiet=0,
                )

                fasta_prot_1 = auxiliary_pymol.cmd.get_fastastr(a_protein_1_pymol_selection_string)
                seq_len_protein_1 = len(fasta_prot_1[fasta_prot_1.find("\n"):])
                rmsd_dict = {
                    "rmsd": str(round(tmp_align_results[0], 2)),
                    "aligned_residues": f"{str(tmp_align_results[1])} / {seq_len_protein_1}",
                }
                # </editor-fold>

                # <editor-fold desc="Calculate distances of C-alpha atoms">
                index_list = []
                ref_chain_list: [str] = []
                ref_pos_list: [int] = []
                ref_resi_list: [str] = []
                model_chain_list: [str] = []
                model_pos_list: [int] = []
                model_resi_list: [str] = []
                distance_list = []

                idx2resi = []
                auxiliary_pymol.cmd.iterate(
                    "aln", "idx2resi.append((model, chain, resi, resn))", space={"idx2resi": idx2resi},
                )
                prot_1_indices = []
                prot_2_indices = []
                for tmp_prot_atom in idx2resi:
                    module_logger.debug(idx2resi)
                    if tmp_prot_atom[0] == tmp_protein_1_name:
                        try:
                            tmp_residue_number = int(tmp_prot_atom[2])
                        except ValueError:
                            tmp_residue_number = int(tmp_prot_atom[2][:-1])
                        prot_1_indices.append((tmp_prot_atom[1], tmp_residue_number, tmp_prot_atom[3]))
                    if tmp_prot_atom[0] == tmp_protein_2_name:
                        try:
                            tmp_residue_number = int(tmp_prot_atom[2])
                        except ValueError:
                            tmp_residue_number = int(tmp_prot_atom[2][:-1])
                        prot_2_indices.append((tmp_prot_atom[1], tmp_residue_number, tmp_prot_atom[3]))
                # calculate the distance between the alpha-C atoms
                for resi_no in range(len(prot_1_indices)):
                    atom1 = f"/{tmp_protein_1_name}//{prot_1_indices[resi_no][0]}/{prot_1_indices[resi_no][1]}/CA"
                    atom2 = f"/{tmp_protein_2_name}//{prot_2_indices[resi_no][0]}/{prot_2_indices[resi_no][1]}/CA"
                    # Raises Exception if "Error: Selection 1: More than one atom found"
                    distance = round(auxiliary_pymol.cmd.get_distance(atom1, atom2, state=1), 2)
                    module_logger.debug(f"Distance of atom_1 {atom1} and atom_2 {atom2} is {distance}")

                    ref_chain_list.append(prot_1_indices[resi_no][0])
                    ref_pos_list.append(int(prot_1_indices[resi_no][1]))
                    ref_resi_list.append(prot_1_indices[resi_no][2])
                    model_chain_list.append(prot_2_indices[resi_no][0])
                    model_pos_list.append(int(prot_2_indices[resi_no][1]))
                    model_resi_list.append(prot_2_indices[resi_no][2])
                    distance_list.append(distance)
                    index_list.append(resi_no)

                index_array: np.ndarray = np.array(index_list)
                ref_chain_array: np.ndarray = np.array(ref_chain_list)
                ref_pos_array: np.ndarray = np.array(ref_pos_list)
                ref_resi_array: np.ndarray = np.array(ref_resi_list)
                model_chain_array: np.ndarray = np.array(model_chain_list)
                model_pos_array: np.ndarray = np.array(model_pos_list)
                model_resi_array: np.ndarray = np.array(model_resi_list)
                distance_array: np.ndarray = np.array(distance_list)

                result_hashtable: dict[str, np.ndarray] = {
                    local_constants.ARRAY_DISTANCE_INDEX: index_array,
                    local_constants.ARRAY_DISTANCE_PROT_1_CHAIN: ref_chain_array,
                    local_constants.ARRAY_DISTANCE_PROT_1_POSITION: ref_pos_array,
                    local_constants.ARRAY_DISTANCE_PROT_1_RESI: ref_resi_array,
                    local_constants.ARRAY_DISTANCE_PROT_2_CHAIN: model_chain_array,
                    local_constants.ARRAY_DISTANCE_PROT_2_POSITION: model_pos_array,
                    local_constants.ARRAY_DISTANCE_PROT_2_RESI: model_resi_array,
                    local_constants.ARRAY_DISTANCE_DISTANCES: distance_array,
                }
                distances = result_hashtable

                # </editor-fold>

                # <editor-fold desc="Create alignment scene">
                auxiliary_pymol.cmd.show("cartoon")
                # fixme: this could become useful if certain objects should be hidden from the scene
                # if selection != "":
                #     auxiliary_pymol.cmd.hide("cartoon", selection)
                auxiliary_pymol.cmd.hide("cgo", "all")
                auxiliary_pymol.cmd.reset()
                auxiliary_pymol.cmd.zoom("all")
                auxiliary_pymol.cmd.scene(
                    key=f"{tmp_protein_1_name}-{tmp_protein_2_name}",
                    action="store",
                )
                # </editor-fold>

                # <editor-fold desc="Create scenes for interesting regions">
                j: int = 0
                i: int = 0
                for distance_value in distances.get("distance"):
                    if float(distance_value) > float(a_cutoff):
                        # here algorithm for image
                        ref_pos_array: np.ndarray = distances.get("ref_pos")
                        ref_pos: int = ref_pos_array[i]
                        ref_chain_array: np.ndarray = distances.get("ref_chain")
                        ref_chain: str = ref_chain_array[i]
                        model_pos_array: np.ndarray = distances.get("model_pos")
                        model_pos: int = model_pos_array[i]
                        model_chain_array: np.ndarray = distances.get("model_chain")
                        model_chain: str = model_chain_array[i]

                        measurement_obj: str = f"measure{j}"
                        # create two atoms for the get_distance command
                        atom1: str = f"/{tmp_protein_1_name}//{ref_chain}/{ref_pos}/CA"
                        atom2: str = f"/{tmp_protein_2_name}//{model_chain}/{model_pos}/CA"
                        # zoom to reference amino acid
                        auxiliary_pymol.cmd.zoom(f"/{tmp_protein_1_name}//{ref_chain}/{ref_pos}", 10)
                        # create distance object with get_distance command
                        auxiliary_pymol.cmd.distance(measurement_obj, atom1, atom2)
                        auxiliary_pymol.cmd.label(atom1, "'%s-%s' % (resn, resi)")
                        auxiliary_pymol.cmd.label(atom2, "'%s-%s' % (resn, resi)")
                        auxiliary_pymol.cmd.set("label_position", (0, 0, 10))
                        # set scene
                        auxiliary_pymol.cmd.scene(key=f"{ref_pos}-{model_pos}", action="store")
                        # hide created labels
                        auxiliary_pymol.cmd.hide("labels", atom1)
                        auxiliary_pymol.cmd.hide("labels", atom2)
                        auxiliary_pymol.cmd.hide("labels", measurement_obj)
                        auxiliary_pymol.cmd.hide("dashes", measurement_obj)
                        i += 1
                        j += 1
                    else:
                        i += 1
                # </editor-fold>

                # <editor-fold desc="Save session">
                session_filepath = pathlib.Path(f"{local_constants.SCRATCH_DIR}/{the_protein_pair_name}_session.pse")
                auxiliary_pymol.cmd.save(session_filepath)
                base64_string = utils.create_base64_string_from_file(str(session_filepath))
                os.remove(session_filepath)
                # </editor-fold>

                # Convert numpy arrays to lists
                result_lists_hashtable: dict[str, list] = {
                    local_constants.ARRAY_DISTANCE_INDEX: index_list,
                    local_constants.ARRAY_DISTANCE_PROT_1_CHAIN: ref_chain_list,
                    local_constants.ARRAY_DISTANCE_PROT_1_POSITION: ref_pos_list,
                    local_constants.ARRAY_DISTANCE_PROT_1_RESI: ref_resi_list,
                    local_constants.ARRAY_DISTANCE_PROT_2_CHAIN: model_chain_list,
                    local_constants.ARRAY_DISTANCE_PROT_2_POSITION: model_pos_list,
                    local_constants.ARRAY_DISTANCE_PROT_2_RESI: model_resi_list,
                    local_constants.ARRAY_DISTANCE_DISTANCES: distance_list,
                }
                distance_analysis_results_object_values = (
                    result_lists_hashtable,
                    base64_string,
                    float(rmsd_dict["rmsd"]),
                    rmsd_dict["aligned_residues"],
                )
        except Exception as e:
            if str(e).find("Error: Selection 1: More than one atom found") != -1:
                raise AttributeError("Error: Selection 1: More than one atom found")
            else:
                tmp_msg: str = f"The analysis in PyMOL failed with the error: {e}"
                print(tmp_msg)
        else:
            return distance_analysis_results_object_values, base64_string

    @staticmethod
    def create_ray_traced_image(
        an_image_filepath: str,
        the_cached_session: str,
        an_image_ray_trace_mode: str,
        an_image_ray_texture: str,
        an_image_renderer: int,
    ) -> None:
        """Creates a ray traced image using PyMOL.

        Args:
            an_image_filepath: The file path to save the image.
            the_cached_session: The path to the cached session.
            an_image_ray_trace_mode: The ray trace mode for the image.
            an_image_ray_texture: The ray texture for the image.
            an_image_renderer: The renderer for the image.

        Raises:
            Exception: If any unexpected exception occurs during the process.
        """
        # <editor-fold desc="Checks">
        if an_image_filepath is None or an_image_filepath == "":
            raise ValueError("an_image_filepath is either None or empty.")
        if the_cached_session is None or the_cached_session == "":
            raise ValueError("the_cached_session is either None or empty.")
        if an_image_ray_trace_mode is None or an_image_ray_trace_mode == "":
            raise ValueError("an_image_ray_trace_mode is either None or empty.")
        if an_image_ray_texture is None or an_image_ray_texture == "":
            raise ValueError("an_image_ray_texture is either None or empty.")
        if an_image_renderer is None:
            raise ValueError("an_image_renderer is None.")

        # </editor-fold>

        try:
            with pymol2.PyMOL() as auxiliary_pymol:
                auxiliary_pymol.cmd.load(the_cached_session)
                auxiliary_pymol.cmd.set("ray_trace_mode", an_image_ray_trace_mode)
                auxiliary_pymol.cmd.set("ray_texture", an_image_ray_texture)
                auxiliary_pymol.cmd.set("max_threads", str(os.cpu_count() - 3))
                auxiliary_pymol.cmd.ray(2400, 2400, renderer=int(an_image_renderer))
                auxiliary_pymol.cmd.png(an_image_filepath, dpi=300)
        except Exception as e:
            print(f"Unexpected exception. {e}")

    @staticmethod
    def create_pymol_session_for_protein(a_pdb_filepath: str) -> str:
        """Creates a new pymol session for a given protein.

        Args:
            a_pdb_filepath (str): The file path of the protein PDB file.

        Returns:
            A base64 encoded string representation of the PyMOL session for the protein
            or an empty string if the argument is illegal or an error was caught.

        Note:
            This method creates a PyMOL session for the given protein PDB file. It loads the PDB file into PyMOL,
            sets various PyMOL settings to customize the visualization, saves the session, converts it to a
            base64 encoded string, and finally returns the string representation of the PyMOL session.
        """
        # <editor-fold desc="Checks">
        if a_pdb_filepath is None or a_pdb_filepath == "":
            return ""
        if not os.path.exists(a_pdb_filepath):
            return ""

        # </editor-fold>

        try:
            tmp_protein_name: str = pathlib.Path(a_pdb_filepath).name.replace(".pdb", "")

            with pymol2.PyMOL() as auxiliary_pymol:
                auxiliary_pymol.cmd.load(filename=str(a_pdb_filepath), object=tmp_protein_name)
                auxiliary_pymol.cmd.bg_color(local_constants.PYMOL_DEFAULT_BACKGROUND_COLOR)
                auxiliary_pymol.cmd.set("valence", 0)
                auxiliary_pymol.cmd.set("scene_buttons", 0)
                auxiliary_pymol.cmd.set("ray_trace_mode", local_constants.PYMOL_DEFAULT_RAY_TRACE_MODE)
                auxiliary_pymol.cmd.set("antialias", local_constants.PYMOL_DEFAULT_ANTIALIAS)
                auxiliary_pymol.cmd.set("ambient", local_constants.PYMOL_DEFAULT_AMBIENT)
                auxiliary_pymol.cmd.set("cartoon_fancy_helices", local_constants.PYMOL_DEFAULT_FANCY_HELICES)
                auxiliary_pymol.cmd.set(
                    "cartoon_discrete_colors", local_constants.PYMOL_DEFAULT_CARTOON_DISCRETE_COLORS,
                )
                auxiliary_pymol.cmd.set("cartoon_sampling", local_constants.PYMOL_DEFAULT_CARTOON_SAMPLING)
                auxiliary_pymol.cmd.set("spec_power", local_constants.PYMOL_DEFAULT_SPEC_POWER)
                auxiliary_pymol.cmd.set("spec_reflect", local_constants.PYMOL_DEFAULT_SPEC_REFLECT)
                auxiliary_pymol.cmd.set(
                    "ray_transparency_contrast", local_constants.PYMOL_DEFAULT_RAY_TRANSPARENCY_CONTRAST,
                )
                auxiliary_pymol.cmd.set(
                    "ray_transparency_oblique", local_constants.PYMOL_DEFAULT_RAY_TRANSPARENCY_OBLIQUE,
                )  # noqa: E501
                auxiliary_pymol.cmd.set(
                    "ray_transparency_oblique_power", local_constants.PYMOL_DEFAULT_RAY_OBLIQUE_POWER,
                )
                auxiliary_pymol.cmd.set("ray_trace_color", local_constants.PYMOL_DEFAULT_RAY_COLOR)
                auxiliary_pymol.cmd.unset("depth_cue")
                auxiliary_pymol.cmd.color("green", tmp_protein_name)
                auxiliary_pymol.cmd.reset()
                auxiliary_pymol.cmd.scene("base", action="store")
                session_filepath: pathlib.Path = pathlib.Path(
                    f"{local_constants.SCRATCH_DIR}/{tmp_protein_name}_session.pse",
                )
                auxiliary_pymol.cmd.save(session_filepath)
                base64_string: str = utils.create_base64_string_from_file(str(session_filepath))
                os.remove(session_filepath)
        except Exception as e:
            return ""
        else:
            return base64_string

    @staticmethod
    def get_protein_chains(a_pdb_filepath: str) -> list[tuple]:
        """Divides the chains from a protein, into protein and non-protein chains.

        Args:
            a_pdb_filepath: A file path to a .pdb file.

        Returns:
            a list of all chains (as tuples) from the protein, divided into protein and non-protein chains
            or an empty list if the argument is illegal or an error was caught.
        """
        # <editor-fold desc="Checks">
        if a_pdb_filepath is None or a_pdb_filepath == "":
            return []
        if not os.path.exists(a_pdb_filepath):
            return []

        # </editor-fold>

        try:
            with pymol2.PyMOL() as auxiliary_pymol:
                auxiliary_pymol.cmd.load(filename=str(a_pdb_filepath))
                tmp_protein_name: str = pathlib.Path(a_pdb_filepath).name.replace(".pdb", "")
                tmp_chains: list[str] = auxiliary_pymol.cmd.get_chains()
                i: int = 0
                chains_of_protein: list[tuple] = []
                for tmp_chain in tmp_chains:
                    sequence_of_chain = auxiliary_pymol.cmd.get_model(f"chain {tmp_chain}")
                    if sequence_of_chain.atom[0].resn in local_constants.AMINO_ACID_CODE:
                        fasta_sequence_of_chain = auxiliary_pymol.cmd.get_fastastr(f"chain {tmp_chain}")
                        fasta_sequence_of_chain_without_header = fasta_sequence_of_chain[
                            fasta_sequence_of_chain.find("\n") :
                        ]
                        complete_sequence_of_chain = (
                            tmp_protein_name,
                            fasta_sequence_of_chain_without_header.replace("\n", ""),
                        )
                        chains_of_protein.append(
                            (tmp_chain, complete_sequence_of_chain, local_constants.CHAIN_TYPE_PROTEIN),
                        )
                    else:
                        fasta_sequence_of_chain = auxiliary_pymol.cmd.get_fastastr(f"chain {tmp_chain}")
                        fasta_sequence_of_chain_without_header = fasta_sequence_of_chain[
                            fasta_sequence_of_chain.find("\n") :
                        ]
                        complete_sequence_of_chain = (
                            tmp_protein_name,
                            fasta_sequence_of_chain_without_header.replace("\n", ""),
                        )
                        chains_of_protein.append(
                            (tmp_chain, complete_sequence_of_chain, local_constants.CHAIN_TYPE_NON_PROTEIN),
                        )
                    i += 1
        except Exception as e:
            return []
        else:
            return chains_of_protein

    @staticmethod
    def consolidate_molecule_object_to_first_state(a_pdb_filepath: str) -> str:
        """Consolidate molecule object states.

        This method takes a PDB file path as input and consolidates the molecule object states to the first state.
        It uses PyMOL to perform the consolidation.

        Args:
            a_pdb_filepath (str): The file path of the PDB file.

        Returns:
            The filepath of the consolidated PDB file if successful, or an empty string if an error occurred.
        """
        # <editor-fold desc="Checks">
        if a_pdb_filepath is None or a_pdb_filepath == "":
            return ""
        if not os.path.exists(a_pdb_filepath):
            return ""

        # </editor-fold>

        with pymol2.PyMOL() as auxiliary_pymol:
            try:
                auxiliary_pymol.cmd.load(filename=str(a_pdb_filepath))
                tmp_protein_name = pathlib.Path(a_pdb_filepath).name.replace(".pdb", "")
                tmp_protein_name = tmp_protein_name.replace(" ", "_")
            except Exception as e:
                print(f"Protein states could not be consolidated! Ran into error: {e}")
                return ""

            tmp_states = auxiliary_pymol.cmd.count_states(tmp_protein_name)
            print(f"Number of states: {tmp_states}")
            if tmp_states > 1:
                try:
                    # Create a new object with only the first state
                    auxiliary_pymol.cmd.create("new_object", tmp_protein_name, 1, 1)
                    # Delete the original molecule to keep only the new object
                    auxiliary_pymol.cmd.delete(tmp_protein_name)
                    # Rename the new object to the original name if needed
                    auxiliary_pymol.cmd.set_name("new_object", tmp_protein_name)
                    if auxiliary_pymol.cmd.count_states(tmp_protein_name) > 1:
                        print("Molecule object has still more than one state.")
                    if auxiliary_pymol.cmd.select(f"/{tmp_protein_name}//A/1/CA") > 1:
                        print(
                            f"Molecule object has still more than one state. {auxiliary_pymol.cmd.select(f'/{tmp_protein_name}//A/1/CA')}",
                        )

                    tmp_pdb_cache_filepath = pathlib.Path(
                        f"{local_constants.CACHE_PROTEIN_DIR}/{tmp_protein_name}.pdb",
                    )
                    auxiliary_pymol.cmd.save(tmp_pdb_cache_filepath)
                except Exception as e:
                    print(f"Protein states could not be consolidated! Ran into error: {e}")
                    return ""
                else:
                    print("Protein states successfully consolidated!")
                    return str(tmp_pdb_cache_filepath)
            else:
                print("No states to consolidate.")
                return a_pdb_filepath

    @staticmethod
    def get_all_scenes_of_session(a_pymol_session: str) -> list:
        """Gets all scenes from a given pymol session.

        Args:
            a_pymol_session: A base64 string representing the PyMOL session data.

        Returns:
            A list of all scenes in the PyMOL session or an empty list if the argument is illegal or an error was caught.
        """
        # <editor-fold desc="Checks">
        if a_pymol_session == "":
            return []

        # </editor-fold>

        try:
            tmp_session_filepath: pathlib.Path = pathlib.Path(f"{local_constants.SCRATCH_DIR}/temp_session.pse")
            utils.write_binary_file_from_base64_string(tmp_session_filepath, a_pymol_session)
            with pymol2.PyMOL() as auxiliary_pymol:
                auxiliary_pymol.cmd.load(str(tmp_session_filepath))
                tmp_all_scenes: list = auxiliary_pymol.cmd.get_scene_list()
            os.remove(tmp_session_filepath)
        except Exception as e:
            return []
        else:
            return tmp_all_scenes

    @staticmethod
    def clean_protein_update(a_pymol_session: str, a_protein_name: str) -> tuple[str, str]:
        """Cleans a protein from all sugar and solvent molecules.

        Args:
            a_pymol_session: A base64 pymol session string.
            a_protein_name: A name of a protein.

        Returns:
            A pymol session base64 string with the cleaned protein or a tuple with two empty strings if
            an argument is illegal or an error was caught.
        """
        # <editor-fold desc="Checks">
        if a_pymol_session == "":
            return "", ""
        if a_protein_name == "":
            return "", ""

        # </editor-fold>

        try:
            tmp_session_filepath: pathlib.Path = pathlib.Path(f"{local_constants.SCRATCH_DIR}/temp_session.pse")
            utils.write_binary_file_from_base64_string(tmp_session_filepath, a_pymol_session)
            with pymol2.PyMOL() as auxiliary_pymol:
                auxiliary_pymol.cmd.load(str(tmp_session_filepath))
                tmp_all_object_names = auxiliary_pymol.cmd.get_names()
                if len(tmp_all_object_names) == 0:
                    # return "ERROR: No proteins found", f"Value of tmp_all_object_names: {tmp_all_object_names}"
                    return "", ""
                if a_protein_name not in tmp_all_object_names:
                    # return f"ERROR: Protein with the name {a_protein_name} not found!", f"Value of tmp_all_object_names: {tmp_all_object_names}"
                    return "", ""
                auxiliary_pymol.cmd.remove("solvent")
                auxiliary_pymol.cmd.remove("organic")
                tmp_export_session_filepath: pathlib.Path = pathlib.Path(
                    f"{local_constants.SCRATCH_DIR}/export_temp_session.pse",
                )
                tmp_export_pdb_filepath: pathlib.Path = pathlib.Path(
                    f"{local_constants.SCRATCH_DIR}/export_temp_clean.pdb",
                )
                auxiliary_pymol.cmd.save(str(tmp_export_session_filepath))
                auxiliary_pymol.cmd.save(str(tmp_export_pdb_filepath))
                base64_string: str = utils.create_base64_string_from_file(str(tmp_export_session_filepath))
                # os.remove(tmp_export_session_filepath)
            os.remove(tmp_session_filepath)
        except Exception as e:
            return "", ""
        else:
            return base64_string, str(tmp_export_pdb_filepath)

    @staticmethod
    def clean_protein_new() -> None:
        """clean_protein_new method is used to clean a protein by removing solvent molecules and organic molecules.
        
        However, this method is not yet fully implemented and currently raises a NotImplementedError.

        Raises:
            NotImplementedError: If cleaning a protein that generates a new protein structure is attempted, this exception is raised.
        """
        # fixme: Old code below that might be useful
        # if new_protein is False:
        #     try:
        #         self.load_protein_pymol_session()
        #         remove_solvent_molecules_in_protein()
        #         remove_organic_molecules_in_protein()
        #     except pymol.CmdException:
        #         return  # noqa: RET502 #TODO: needs more thoughts
        #     tmp_was_successful, tmp_pdb_filepath = pymol_io.save_protein_to_pdb_file(constants.CACHE_PROTEIN_DIR,
        #                                                                              str(self._id))
        #     if tmp_was_successful:
        #         self._pdb_data = bio_data.parse_pdb_file(tmp_pdb_filepath)
        #         logger.debug(self._pdb_data)
        #     else:
        #         logger.error("The protein could not be cleaned, because the new pdb file could not be found!")
        #         raise RuntimeError("The protein could not be cleaned, because the new pdb file could not be found!")
        raise NotImplementedError("Cleaning a protein that generates a new protein structure is not yet implemented!")
