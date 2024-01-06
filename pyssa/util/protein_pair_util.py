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
import numpy as np
from pymol import cmd

from pyssa.logging_pyssa import log_handlers
from pyssa.util import pyssa_keys

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def calculate_distance_between_ca_atoms(ref_prot_name: str, model_prot_name: str) -> dict[str, np.ndarray]:
    """Calculates the distance between two alpha carbon-atoms.

    Args:
        ref_prot_name: The name of the reference protein.
        model_prot_name: The name of the model protein.
    """
    index_list = []
    ref_chain_list: [str] = []
    ref_pos_list: [int] = []
    ref_resi_list: [str] = []
    model_chain_list: [str] = []
    model_pos_list: [int] = []
    model_resi_list: [str] = []
    distance_list = []

    idx2resi = []
    cmd.iterate("aln", "idx2resi.append((model, chain, resi, resn))", space={"idx2resi": idx2resi})
    prot_1_indices = []
    prot_2_indices = []
    for tmp_prot_atom in idx2resi:
        if tmp_prot_atom[0] == ref_prot_name:
            prot_1_indices.append((tmp_prot_atom[1], tmp_prot_atom[2], tmp_prot_atom[3]))
        if tmp_prot_atom[0] == model_prot_name:
            prot_2_indices.append((tmp_prot_atom[1], tmp_prot_atom[2], tmp_prot_atom[3]))
    # calculate the distance between the alpha-C atoms
    for resi_no in range(len(prot_1_indices)):
        atom1 = f"/{ref_prot_name}//{prot_1_indices[resi_no][0]}/{prot_1_indices[resi_no][1]}/CA"
        atom2 = f"/{model_prot_name}//{prot_2_indices[resi_no][0]}/{prot_2_indices[resi_no][1]}/CA"
        distance = round(cmd.get_distance(atom1, atom2, state=-1), 2)

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

    # return the hast table with all results
    return result_hashtable


def get_chain_and_position(results_hashtable: dict[str, np.ndarray], i: int) -> tuple:
    """This function gets the chain and the residue postion from the results hash table.

    Args:
        results_hashtable (dict(str, np.ndarray)):
            hash table which contains all information from the
            distance calculation
        i (int):
            interator for the results hash table index

    Returns:
        tuple
    """
    ref_chain_array = results_hashtable.get(pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN)
    ref_chain = ref_chain_array[i]
    ref_pos_array = results_hashtable.get(pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION)
    ref_pos = ref_pos_array[i]
    ref_resi_array = results_hashtable.get(pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI)
    ref_resi = ref_resi_array[i]
    model_chain_array = results_hashtable.get(pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN)
    model_chain = model_chain_array[i]
    model_pos_array = results_hashtable.get(pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION)
    model_pos = model_pos_array[i]
    model_resi_array = results_hashtable.get(pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI)
    model_resi = model_resi_array[i]
    return (ref_chain, ref_pos, ref_resi, model_chain, model_pos, model_resi)
