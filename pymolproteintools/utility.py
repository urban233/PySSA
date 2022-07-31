# ----------------------------------------------------------------------- #
# Python Package: pymolproteintools
# Version 1.0 for Python 3.9
# ----------------------------------------------------------------------- #
# Authors: Martin Urban & Hannah Kullik
# Westfaelische Hochschule
# Recklinghausen, Germany
#
# Copyright 2022 Martin Urban & Hannah Kullik
#
# Citation:
# Martin Urban & Hannah Kullik, pymolproteintools (PPT),
# Version 1.0, Recklinghausen, Germany, 2022.
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License (LGPL) as
# published by the Free Software Foundation, version 3 of the License.
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License (LGPL) for more details.
# You should have received a copy of the GNU Lesser General Public
# License along with this program. If not, see
# <http://www.gnu.org/licenses/>.

# ----------------------------------------------------------------------- #
import pandas as pd
from pymol import cmd
import numpy as np


def get_fasta_sequence(protein_object: str) -> str:
    """Creates the FASTA sequence of the proteinObject in PyMOL.

    Args:
        protein_object (str):
            name of the protein object in PyMOL
    Returns:
        str:
            FASTA sequence of the proteinObject
    """

    fasta_seq: str = cmd.get_fastastr(protein_object)
    # indexOfEscapeChar: int = fasta_seq.find("\n")
    # sequence = fasta_seq[indexOfEscapeChar + 1:len(fasta_seq)]
    sequence = fasta_seq.replace("\n", "")
    return sequence


def get_chain_and_position(results_hashtable, i):
    """This function gets the chain and the residue postion from the
    results hash table

    Args:
        results_hashtable (dict(str, np.ndarray)):
            hash table which contains all information from the
            distance calculation
        i (int):
            interator for the results hash table index

    Returns:
        tuple
    """
    ref_chain_array = results_hashtable.get("ref_chain")
    ref_chain = ref_chain_array[i]
    ref_pos_array = results_hashtable.get("ref_pos")
    ref_pos = ref_pos_array[i]
    ref_resi_array = results_hashtable.get("ref_resi")
    ref_resi = ref_resi_array[i]
    model_chain_array = results_hashtable.get("model_chain")
    model_chain = model_chain_array[i]
    model_pos_array = results_hashtable.get("model_pos")
    model_pos = model_pos_array[i]
    model_resi_array = results_hashtable.get("model_resi")
    model_resi = model_resi_array[i]
    tmp_tuple = (ref_chain, ref_pos, ref_resi, model_chain, model_pos,
                 model_resi)
    return tmp_tuple


# def convert_linked_list_in_1d_array(linked_list, datatype=None) -> np.ndarray:
#     """This function converts a linked list from the structlinks.DataStructures
#     package into a numpy array.
#
#     Args:
#         linked_list:
#             a linked list which should be converted
#         datatype:
#             defines the data type of the array
#     Returns:
#         tmp_array:
#             one dimensional array with the data from the linked list
#     """
#     if datatype != object:
#         tmp_array: np.ndarray = np.arange(0, len(linked_list), dtype=datatype)
#     else:
#         tmp_array: np.ndarray = np.array(range(len(linked_list)), dtype='|S3')
#
#     i: int = 0
#     for j in linked_list:
#         tmp_array.flat[i] = j
#         i += 1
#     return tmp_array
#
#
# def generate_float_numbers(start_number: float, end_number: float,
#                            step: float) -> np.ndarray:
#     """This function generates rising float numbers between the start_number and
#     end_number, in a pre-defined step.
#
#     Args:
#         start_number (float):
#              first number of the float list
#         end_number (float):
#              last number of the float list
#         step (float):
#              step in which all numbers between start_number and end_number
#              are generates
#
#     Returns:
#         np.ndarray:
#             an array of rising float numbers
#
#     Raises:
#         ValueError: If end_number is greater or equal than start_number.
#     Todo:
#         * create unittest
#     """
#     if start_number >= end_number:
#         raise ValueError
#
#     numbers_ll: list = [start_number]
#     number = start_number + step
#     while number <= end_number:
#         numbers_ll.append(number)
#         number += step
#
#     numbers_array: np.ndarray = convert_linked_list_in_1d_array(numbers_ll,
#                                                                 float)
#     return numbers_array