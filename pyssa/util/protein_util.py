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
"""This module contains helper function for the protein class."""
import logging
from pyssa.io_pyssa import safeguard
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants
from pyssa.util import types

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def check_if_protein_is_from_file_or_id(molecule_object) -> tuple:
    """This function checks if a protein is from a .pdb file or from a PDB ID.

    Args:
        molecule_object:
            the name of the protein which is also used within pymol
    Returns:
        a tuple with the molecule object and the filename
    """
    # <editor-fold desc="Checks">
    if not safeguard.Safeguard.check_if_value_is_not_none(molecule_object) or molecule_object == "":
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")

    # </editor-fold>

    if molecule_object.find(".pdb") != -1:
        # pdb file is given
        tmp_molecule_object = molecule_object.replace(".pdb", "")
        tmp_filename = molecule_object
    else:
        # PDB ID is given
        tmp_filename = f"{molecule_object}.pdb"
        tmp_molecule_object = molecule_object

    return tmp_molecule_object, tmp_filename


def filter_chains_for_protein_chains(chains: list[types.CHAIN]) -> list[types.CHAIN]:
    """This function filters the chains for protein chains only

    Args:
        chains:
            a list of chains which occur in the protein
    Returns:
        a list of protein chains only
    """
    # <editor-fold desc="Checks">
    if not safeguard.Safeguard.check_if_value_is_not_none(chains) or not safeguard.Safeguard.check_if_list_is_empty(chains):
        logger.error("An argument is illegal.")
        raise ValueError("An argument is illegal.")
    # </editor-fold>

    protein_chains = []
    for tmp_chain in chains:
        if tmp_chain.chain_type == constants.CHAIN_TYPE_PROTEIN:
            protein_chains.append(tmp_chain)
    return protein_chains
