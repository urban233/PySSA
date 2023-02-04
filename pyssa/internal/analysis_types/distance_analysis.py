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
"""This module contains the distance analysis class."""
import os
import logging
import pathlib
import json
from pymol import cmd
from typing import TYPE_CHECKING
from pyssa.logging_pyssa import log_handlers
from pyssa.util import protein_pair_util
from pyssa.internal.portal import protein_pair_operations

if TYPE_CHECKING:
    from pyssa.internal.data_structures import protein_pair
    from pyssa.internal.data_structures import settings

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class DistanceAnalysis:
    """This class contains information about the distance analysis"""

    # <editor-fold desc="Class attributes">
    """
    a pair of proteins which get analyzed 
    """
    protein_pair_for_analysis: protein_pair.ProteinPair
    """
    the number of refinement cycles for the align command from PyMOL
    """
    cycles: int
    """
    the cutoff value for the align command from PyMOL
    """
    cutoff: float
    """
    the filename of the alignement file which gets created during the align command
    """
    alignment_file_name: str
    """
    the session filename for the PyMOL session of the analysis
    """
    session_file_name: str
    """
    the settings from pyssa
    """
    app_settings: settings.Settings

    # </editor-fold>

    def __init__(self, protein_pair_for_analysis: 'protein_pair.ProteinPair', app_settings: 'settings.Settings'):
        """Constructor

        Args:
            protein_pair:
                a pair of proteins which get analyzed
            app_settings:
                the settings from pyssa
        """
        self.protein_pair_for_analysis: protein_pair.ProteinPair = protein_pair_for_analysis
        self.app_settings: settings.Settings = app_settings
        self.cutoff: float = app_settings.cutoff
        self.cycles: int = app_settings.cycles

    def create_align_selections(self) -> None:
        """This function creates the selection which are needed for the align command.

        """
        if len(self.protein_pair_for_analysis.ref_obj.chains) > 0:
            self.protein_pair_for_analysis.ref_obj.selection.set_selections_from_chains_ca(self.protein_pair_for_analysis.ref_obj.chains)
        else:
            self.protein_pair_for_analysis.model_obj.selection.set_selections_without_chains_ca()
        if len(self.protein_pair_for_analysis.model_obj.chains) > 0:
            self.protein_pair_for_analysis.model_obj.selection.set_selections_from_chains_ca(self.protein_pair_for_analysis.model_obj.chains)
        else:
            self.protein_pair_for_analysis.model_obj.selection.set_selections_without_chains_ca()

    def align_protein_pair_for_analysis(self) -> tuple:
        """This function aligns the protein pair with the PyMOL align command.

        Returns:
            a tuple with the rmsd and aligned amino acids
        """
        results = protein_pair_operations.align_protein_pair(self.protein_pair_for_analysis.ref_obj.selection.selection_string,
                                                             self.protein_pair_for_analysis.model_obj.selection.selection_string,
                                                             self.alignment_file_name)
        # save the align object from pymol as alignment file
        if not os.path.exists(pathlib.Path(f"{self.protein_pair_for_analysis.results_dir}/alignment_files")):
            os.mkdir(pathlib.Path(f"{self.protein_pair_for_analysis.results_dir}/alignment_files"))
        cmd.save(pathlib.Path(f"{self.protein_pair_for_analysis.results_dir}/alignment_files/{self.alignment_file_name}.aln"))
        return results[0], results[1]

    def do_analysis_in_pymol(self):
        self.protein_pair_for_analysis.load_protein_pair()
        self.protein_pair_for_analysis.color_protein_pair()
        align_results = self.align_protein_pair_for_analysis()
        rmsd_dict = {
            "rmsd": str(round(align_results[0], 2)),
            "aligned_residues": str(align_results[1]),
        }
        rmsd_file = open(pathlib.Path(f"{self.protein_pair_for_analysis.results_dir}/rmsd.json"), "w", encoding="utf-8")
        json.dump(rmsd_dict, rmsd_file, indent=4)
        rmsd_file.close()
        protein_pair_util.calculate_distance_between_ca_atoms(self.protein_pair_for_analysis, self.alignment_file_name)
        self.protein_pair_for_analysis.save_session_of_protein_pair(self.session_file_name)

