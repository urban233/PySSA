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
"""This module contains helper function for the protein pair class."""
import logging
import pathlib
import numpy as np
from pyssa.io_pyssa import safeguard
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants
from pyssa.util import types
from pyssa.internal.data_structures.data_classes import atom_for_distance_calc
from Bio import AlignIO
from pymol import cmd

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def calculate_distance_between_ca_atoms(protein_pair: types.PROTEIN_PAIR, alignment_filename: str,
                                        cutoff: float = 20.0) -> dict[str, np.ndarray]:
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
    if not safeguard.Safeguard.check_if_file_is_readable(
            pathlib.Path(f"{protein_pair.results_dir}/alignment_files/{alignment_filename}.aln")):
        print(f"File not found, in {protein_pair.results_dir}.")

    cmd.create(f"{protein_pair.ref_obj.molecule_object}_CA", f"/{protein_pair.ref_obj.molecule_object}////CA")
    ref_ca_obj = cmd.get_model(f"{protein_pair.ref_obj.molecule_object}_CA")
    cmd.create(f"{protein_pair.model_obj.molecule_object}_CA", f"/{protein_pair.model_obj.molecule_object}////CA")
    model_ca_obj = cmd.get_model(f"{protein_pair.model_obj.molecule_object}_CA")
    # read in alignment file from alignProteinPair function
    align = AlignIO.read(pathlib.Path(f"{protein_pair.results_dir}/alignment_files/{alignment_filename}.aln"), "clustal")
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
            atom1 = f"/{protein_pair.ref_obj.molecule_object}_CA//{chain_ref}/{pos_ref}/"
            atom2 = f"/{protein_pair.model_obj.molecule_object}_CA//{chain_model}/{pos_model}/"
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
