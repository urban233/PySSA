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
"""This module contains the selection class."""
import logging
from pyssa.io_pyssa import safeguard
from pyssa.logging_pyssa import log_handlers
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pyssa.internal.data_structures import chain

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class Selection:
    """This class contains information about a selection."""

    # <editor-fold desc="Class attribute">
    """
    a pymol conform selection string
    """
    selection_string: str
    """
    the name of the protein which is also used within pymol
    """
    molecule_object: str
    """
    a list of all current chains from the selection
    """
    selection_chain_letters: list[str]

    # </editor-fold>

    def __init__(self, molecule_object: str) -> None:
        """Constructor.

        Args:
            molecule_object: the name of the protein which is also used within pymol
        """
        # <editor-fold desc="Checks">
        if not safeguard.Safeguard.check_if_value_is_not_none(molecule_object) or molecule_object == "":
            logger.error("An argument is illegal.")
            raise ValueError("An argument is illegal.")

        # </editor-fold>

        self.molecule_object = molecule_object

    def set_selections_from_chains_ca(self, chains: list["chain.Chain"]) -> None:
        """This function sets a selection based on the chains of the protein. The selection selects only the alpha-C's.

        Args:
            chains:
                a list of chains
        """
        seperator = ", "
        tmp_list = []
        self.selection_chain_letters = []
        for tmp_chain in chains:
            tmp_selection = f"/{self.molecule_object}//{tmp_chain.chain_letter}//CA"
            self.selection_chain_letters.append(tmp_chain.chain_letter)
            tmp_list.append(tmp_selection)
            self.selection_string = seperator.join(tmp_list)

    def set_selections_without_chains_ca(self) -> None:
        """Sets a selection without any chains of the protein.

        Notes:
            The selection selects only the alpha-C's.
        """
        self.selection_string = f"/{self.molecule_object}////CA"

    def set_single_selection(self, segi: str, chain: str, resi: str, atom_name: str) -> None:
        """Creates a single pymol selection with only one chain and one resi.

        Args:
            segi:
                a segment identifier
            chain:
                a chain identifier
            resi:
                a residue name or position
            atom_name:
                a type of atom
        """
        self.selection_string = f"/{self.molecule_object}/{segi}/{chain}/{resi}/{atom_name}"

    def set_custom_selection(self, sele_string: str) -> None:
        """Sets a custom selection as selection string.

        Args:
            sele_string: a custom pymol selection string.
        """
        self.selection_string = sele_string
