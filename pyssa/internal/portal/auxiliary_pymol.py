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
"""Module for the auxiliary pymol class."""
import logging
import os
import pathlib

import pymol2
import pymol
import numpy as np
from pyssa.internal.data_structures import project, structure_prediction, job, protein, results, protein_pair
from pyssa.internal.portal import pymol_io, graphic_operations
from pyssa.io_pyssa import bio_data, path_util, binary_data
from pyssa.logging_pyssa import log_handlers
from pyssa.util import enums, exception, exit_codes, constants, constant_messages, pyssa_keys

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class AuxiliaryPyMOL:

    def __init__(self):
        pass

    @staticmethod
    def create_pymol_session_for_protein(a_protein: "protein.Protein") -> str:
        # <editor-fold desc="Checks">
        if not os.path.exists(constants.CACHE_PROTEIN_DIR):
            os.mkdir(constants.CACHE_PROTEIN_DIR)
        if len(a_protein.get_pdb_data()) == 0:
            raise ValueError("No pdb data in current object!")
        # </editor-fold>

        pdb_filepath = pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{a_protein.get_molecule_object()}.pdb")
        try:
            bio_data.build_pdb_file(a_protein.get_pdb_data(),
                                    pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{a_protein.get_molecule_object()}.pdb"))
        except exception.IllegalArgumentError:
            logger.error(f"The argument pdb data is not usable: {a_protein.get_pdb_data}.")
            raise exception.UnableToCreatePdbFileError("")
        except exception.DirectoryNotFoundError:
            logger.error(f"The argument pdb_filepath is illegal: {pdb_filepath}!")
            raise exception.UnableToCreatePdbFileError("")
        except PermissionError:
            logger.error(f"The argument pdb_filepath is illegal: {pdb_filepath}!")
            raise exception.UnableToCreatePdbFileError("")
        except exception.UnableToOpenFileError:
            logger.error("pdb file could not be opened for writing.")
            raise exception.UnableToOpenFileError("")

        with pymol2.PyMOL() as auxiliary_pymol:
            auxiliary_pymol.cmd.load(
                filename=str(pathlib.Path(f"{constants.CACHE_PROTEIN_DIR}/{a_protein.get_molecule_object()}.pdb")),
                object=a_protein.get_molecule_object()
            )
            auxiliary_pymol.cmd.bg_color(constants.PYMOL_DEFAULT_BACKGROUND_COLOR)
            auxiliary_pymol.cmd.set("valence", 0)
            auxiliary_pymol.cmd.set("scene_buttons", 0)
            auxiliary_pymol.cmd.set("ray_trace_mode", constants.PYMOL_DEFAULT_RAY_TRACE_MODE)
            auxiliary_pymol.cmd.set("antialias", constants.PYMOL_DEFAULT_ANTIALIAS)
            auxiliary_pymol.cmd.set("ambient", constants.PYMOL_DEFAULT_AMBIENT)
            auxiliary_pymol.cmd.set("cartoon_fancy_helices", constants.PYMOL_DEFAULT_FANCY_HELICES)
            auxiliary_pymol.cmd.set("cartoon_discrete_colors", constants.PYMOL_DEFAULT_CARTOON_DISCRETE_COLORS)
            auxiliary_pymol.cmd.set("cartoon_sampling", constants.PYMOL_DEFAULT_CARTOON_SAMPLING)
            auxiliary_pymol.cmd.set("spec_power", constants.PYMOL_DEFAULT_SPEC_POWER)
            auxiliary_pymol.cmd.set("spec_reflect", constants.PYMOL_DEFAULT_SPEC_REFLECT)
            auxiliary_pymol.cmd.set("ray_transparency_contrast", constants.PYMOL_DEFAULT_RAY_TRANSPARENCY_CONTRAST)
            auxiliary_pymol.cmd.set("ray_transparency_oblique", constants.PYMOL_DEFAULT_RAY_TRANSPARENCY_OBLIQUE)  # noqa: E501
            auxiliary_pymol.cmd.set("ray_transparency_oblique_power", constants.PYMOL_DEFAULT_RAY_OBLIQUE_POWER)
            auxiliary_pymol.cmd.set("ray_trace_color", constants.PYMOL_DEFAULT_RAY_COLOR)
            auxiliary_pymol.cmd.unset("depth_cue")
            auxiliary_pymol.cmd.color("green", a_protein.get_molecule_object())
            auxiliary_pymol.cmd.scene("base", action="store")
            session_filepath = pathlib.Path(f"{constants.SCRATCH_DIR}/{a_protein.get_molecule_object()}_session.pse")
            auxiliary_pymol.cmd.save(session_filepath)
            base64_string = binary_data.create_base64_string_from_file(path_util.FilePath(session_filepath))
            os.remove(session_filepath)
        return base64_string

    @staticmethod
    def create_ray_traced_image(
            an_image_filepath,
            the_cached_session,
            an_image_ray_trace_mode,
            an_image_ray_texture,
            an_image_renderer
    ):
        with pymol2.PyMOL() as auxiliary_pymol:
            auxiliary_pymol.cmd.load(the_cached_session)
            auxiliary_pymol.cmd.set("ray_trace_mode", an_image_ray_trace_mode)
            auxiliary_pymol.cmd.set("ray_texture", an_image_ray_texture)
            try:
                auxiliary_pymol.cmd.set("max_threads", str(os.cpu_count() - 3))
                auxiliary_pymol.cmd.ray(2400, 2400, renderer=int(an_image_renderer))
                auxiliary_pymol.cmd.png(an_image_filepath, dpi=300)
            except pymol.CmdException:
                logger.warning("Unexpected exception.")

    @staticmethod
    def do_distance_analysis(
            tmp_protein_pair: "protein_pair.ProteinPair",
    ):
        if tmp_protein_pair is None:
            logger.error(constant_messages.ARGUMENT_IS_ILLEGAL)
            raise exception.IllegalArgumentError(constant_messages.ARGUMENT_IS_ILLEGAL)

        tmp_protein_1_name = tmp_protein_pair.protein_1.get_molecule_object()
        tmp_protein_1_pdb_cache_filepath = pathlib.Path(
            f"{constants.CACHE_PROTEIN_DIR}/{tmp_protein_1_name}.pdb"
        )
        tmp_protein_2_name = tmp_protein_pair.protein_2.get_molecule_object()
        tmp_protein_2_pdb_cache_filepath = pathlib.Path(
            f"{constants.CACHE_PROTEIN_DIR}/{tmp_protein_2_name}.pdb"
        )
        try:
            bio_data.build_pdb_file(
                tmp_protein_pair.protein_1.get_pdb_data(),
                tmp_protein_1_pdb_cache_filepath
            )
            bio_data.build_pdb_file(
                tmp_protein_pair.protein_2.get_pdb_data(),
                tmp_protein_2_pdb_cache_filepath
            )
        except Exception as e:
            logger.error(f"PDB file could not be built. Error: {e}")

        try:
            with pymol2.PyMOL() as auxiliary_pymol:
                # <editor-fold desc="Prepare session">
                # Load both proteins in a new session
                auxiliary_pymol.cmd.load(
                    filename=str(tmp_protein_1_pdb_cache_filepath),
                    object=tmp_protein_1_name
                )
                auxiliary_pymol.cmd.load(
                    filename=str(tmp_protein_2_pdb_cache_filepath),
                    object=tmp_protein_2_name
                )
                logger.info(f"Loaded protein pair: {tmp_protein_pair.name} in pymol session.")

                auxiliary_pymol.cmd.color(constants.DEFAULT_COLOR_PROTEIN_1, tmp_protein_1_name)
                auxiliary_pymol.cmd.color(constants.DEFAULT_COLOR_PROTEIN_2, tmp_protein_2_name)
                logger.info(f"Colored protein pair: {tmp_protein_pair.name} in pymol session.")

                auxiliary_pymol.cmd.scene("base", action="store")

                auxiliary_pymol.cmd.bg_color(constants.PYMOL_DEFAULT_BACKGROUND_COLOR)
                auxiliary_pymol.cmd.set("valence", 0)
                auxiliary_pymol.cmd.set("scene_buttons", 0)
                auxiliary_pymol.cmd.set("ray_trace_mode", constants.PYMOL_DEFAULT_RAY_TRACE_MODE)
                auxiliary_pymol.cmd.set("antialias", constants.PYMOL_DEFAULT_ANTIALIAS)
                auxiliary_pymol.cmd.set("ambient", constants.PYMOL_DEFAULT_AMBIENT)
                auxiliary_pymol.cmd.set("cartoon_fancy_helices", constants.PYMOL_DEFAULT_FANCY_HELICES)
                auxiliary_pymol.cmd.set("cartoon_discrete_colors", constants.PYMOL_DEFAULT_CARTOON_DISCRETE_COLORS)
                auxiliary_pymol.cmd.set("cartoon_sampling", constants.PYMOL_DEFAULT_CARTOON_SAMPLING)
                auxiliary_pymol.cmd.set("spec_power", constants.PYMOL_DEFAULT_SPEC_POWER)
                auxiliary_pymol.cmd.set("spec_reflect", constants.PYMOL_DEFAULT_SPEC_REFLECT)
                auxiliary_pymol.cmd.set("ray_transparency_contrast", constants.PYMOL_DEFAULT_RAY_TRANSPARENCY_CONTRAST)
                auxiliary_pymol.cmd.set("ray_transparency_oblique",
                                        constants.PYMOL_DEFAULT_RAY_TRANSPARENCY_OBLIQUE)  # noqa: E501
                auxiliary_pymol.cmd.set("ray_transparency_oblique_power", constants.PYMOL_DEFAULT_RAY_OBLIQUE_POWER)
                auxiliary_pymol.cmd.set("ray_trace_color", constants.PYMOL_DEFAULT_RAY_COLOR)
                auxiliary_pymol.cmd.unset("depth_cue")
                # </editor-fold>

                # <editor-fold desc="Align proteins">
                tmp_align_results = auxiliary_pymol.cmd.align(
                    target=tmp_protein_pair.protein_1.pymol_selection.selection_string,
                    mobile=tmp_protein_pair.protein_2.pymol_selection.selection_string,
                    cutoff=tmp_protein_pair.distance_analysis.cutoff,
                    cycles=tmp_protein_pair.distance_analysis.cycles,
                    object="aln",
                    quiet=0,
                )
                logger.info(f"Aligned protein pair: {tmp_protein_pair.name} in pymol session.")

                fasta_prot_1 = auxiliary_pymol.cmd.get_fastastr(tmp_protein_pair.protein_1.pymol_selection.selection_string)
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
                auxiliary_pymol.cmd.iterate("aln", "idx2resi.append((model, chain, resi, resn))", space={"idx2resi": idx2resi})
                prot_1_indices = []
                prot_2_indices = []
                for tmp_prot_atom in idx2resi:
                    if tmp_prot_atom[0] == tmp_protein_pair.protein_1.get_molecule_object():
                        prot_1_indices.append((tmp_prot_atom[1], tmp_prot_atom[2], tmp_prot_atom[3]))
                    if tmp_prot_atom[0] == tmp_protein_pair.protein_2.get_molecule_object():
                        prot_2_indices.append((tmp_prot_atom[1], tmp_prot_atom[2], tmp_prot_atom[3]))
                # calculate the distance between the alpha-C atoms
                for resi_no in range(len(prot_1_indices)):
                    atom1 = f"/{tmp_protein_pair.protein_1.get_molecule_object()}//{prot_1_indices[resi_no][0]}/{prot_1_indices[resi_no][1]}/CA"
                    atom2 = f"/{tmp_protein_pair.protein_2.get_molecule_object()}//{prot_2_indices[resi_no][0]}/{prot_2_indices[resi_no][1]}/CA"
                    distance = round(auxiliary_pymol.cmd.get_distance(atom1, atom2, state=-1), 2)

                    ref_chain_list.append(prot_1_indices[resi_no][0])
                    ref_pos_list.append(prot_1_indices[resi_no][1])
                    ref_resi_list.append(prot_1_indices[resi_no][2])
                    model_chain_list.append(prot_2_indices[resi_no][0])
                    model_pos_list.append(prot_2_indices[resi_no][1])
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
                    pyssa_keys.ARRAY_DISTANCE_INDEX: index_array,
                    pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN: ref_chain_array,
                    pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION: ref_pos_array,
                    pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI: ref_resi_array,
                    pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN: model_chain_array,
                    pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION: model_pos_array,
                    pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI: model_resi_array,
                    pyssa_keys.ARRAY_DISTANCE_DISTANCES: distance_array,
                }
                logger.info(f"Calculated distances of protein pair: {tmp_protein_pair.name}.")
                distances = result_hashtable
                analysis_results = results.DistanceAnalysisResults(
                    distances,
                    pymol_io.convert_pymol_session_to_base64_string(tmp_protein_pair.name),
                    float(rmsd_dict["rmsd"]),
                    rmsd_dict["aligned_residues"],
                )
                # </editor-fold>

                # <editor-fold desc="Create alignment scene">
                auxiliary_pymol.cmd.show("cartoon")
                # fixme: this could become useful if certain objects should be hidden from the scene
                # if selection != "":
                #     auxiliary_pymol.cmd.hide("cartoon", selection)
                auxiliary_pymol.cmd.hide("cgo", "all")
                auxiliary_pymol.cmd.orient()
                auxiliary_pymol.cmd.center()
                auxiliary_pymol.cmd.scene(
                    key=f"{tmp_protein_pair.protein_1.get_molecule_object()}-{tmp_protein_pair.protein_2.get_molecule_object()}",
                    action="store",
                )
                # </editor-fold>

                # <editor-fold desc="Create scenes for interesting regions">
                j: int = 0
                i: int = 0
                for distance_value in distances.get("distance"):
                    if float(distance_value) > float(tmp_protein_pair.distance_analysis.cutoff):
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
                        # set image parameters
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
                session_filepath = pathlib.Path(
                    f"{constants.SCRATCH_DIR}/{tmp_protein_pair.name}_session.pse")
                auxiliary_pymol.cmd.save(session_filepath)
                base64_string = binary_data.create_base64_string_from_file(path_util.FilePath(session_filepath))
                os.remove(session_filepath)
                # </editor-fold>
        except Exception as e:
            tmp_msg: str = f"The analysis in PyMOL failed with the error: {e}"
            logger.error(tmp_msg)
            raise exception.UnableToDoAnalysisError(tmp_msg)
        logger.info(
            f"Packed results into the DistanceAnalysisResults object for the protein pair: {tmp_protein_pair.name}.",
        )
        return analysis_results, base64_string


