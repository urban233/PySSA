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
"""Module for all asynchronous functions that are related to the protein pair object."""
from pyssa.controller import pymol_session_manager
from pyssa.internal.portal import graphic_operations


def color_protein_pair_by_rmsd_value(a_protein_pair: "protein_pair.ProteinPair",
                                     the_pymol_session_manager: "pymol_session_manager.PymolSessionManager") -> tuple:
    """Colors a given protein pair by their rmsd value.

    Args:
        a_project: a project object containing the protein pair.
        a_results_name: the name of the protein pair.

    Returns:
        a tuple with ("result", an_existing_protein_pair_object)
    """
    the_pymol_session_manager.color_protein_pair_by_rmsd(a_protein_pair)
    return ("result", a_protein_pair)
