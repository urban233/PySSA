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
"""Module contains helpful functions for the other modules"""

from pymol import cmd


def get_fasta_sequence(protein_object: str) -> str:
    """Creates the FASTA sequence of the proteinObject in PyMOL.

    Args:
        protein_object (str):
            name of the Protein object in PyMOL
    Returns:
        str:
            FASTA sequence of the proteinObject
    """

    fasta_seq: str = cmd.get_fastastr(protein_object)
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
    tmp_tuple = (ref_chain, ref_pos, ref_resi, model_chain, model_pos, model_resi)
    return tmp_tuple
