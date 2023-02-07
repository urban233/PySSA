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

from io_pyssa import filesystem_io
from pyssa.logging_pyssa import log_handlers
from pyssa.util import protein_pair_util
from pyssa.internal.portal import protein_pair_operations
from pyssa.io_pyssa import path_util
from pyssa.util import pyssa_keys

if TYPE_CHECKING:
    from pyssa.internal.data_structures import protein_pair
    from pyssa.internal.data_structures import settings
    from pyssa.internal.data_structures import selection

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class DistanceAnalysis:
    """This class contains information about the distance analysis"""

    # <editor-fold desc="Class attributes">
    """
    a pair of proteins which get analyzed 
    """
    protein_pair_for_analysis: 'protein_pair.ProteinPair'
    """
    the settings from pyssa
    """
    app_settings: 'settings.Settings'
    """
    the cutoff value for the align command from PyMOL
    """
    cutoff: float
    """
    the number of refinement cycles for the align command from PyMOL
    """
    cycles: int
    """
    the filename of the alignement file which gets created during the align command
    """
    alignment_file_name: str
    """
    a directory where all results related to the protein will be stored
    """
    export_dirname: pathlib.Path
    """
    the session filepath for the PyMOL session of the analysis
    """
    pymol_session_filepath: path_util.FilePath

    # </editor-fold>

    def __init__(self, protein_pair_for_analysis: 'protein_pair.ProteinPair', app_settings: 'settings.Settings',
                 distance_analysis_dirname: pathlib.Path):
        """Constructor

        Args:
            protein_pair_for_analysis:
                a pair of proteins which get analyzed
            app_settings:
                the settings from pyssa
            export_dirname:
                a directory where all results will be stored
        """
        self.protein_pair_for_analysis: protein_pair.ProteinPair = protein_pair_for_analysis
        self.app_settings: settings.Settings = app_settings
        self.cutoff: float = app_settings.cutoff
        self.cycles: int = app_settings.cycles

        self.distance_analysis_subdirs = {
            pyssa_keys.DISTANCE_ANALYSIS_SUBDIR: pathlib.Path(f"{distance_analysis_dirname}"),
            pyssa_keys.DISTANCE_ANALYSIS_SESSION_SUBDIR: pathlib.Path(f"{distance_analysis_dirname}/session"),
            pyssa_keys.DISTANCE_ANALYSIS_RESULTS_SUBDIR: pathlib.Path(f"{distance_analysis_dirname}/results"),
            pyssa_keys.DISTANCE_ANALYSIS_OBJECTS_SUBDIR: pathlib.Path(f"{distance_analysis_dirname}/.objects"),
        }
        for key in self.distance_analysis_subdirs:
            if not os.path.exists(self.distance_analysis_subdirs.get(key)):
                os.mkdir(self.distance_analysis_subdirs.get(key))

        self.export_dirname = self.distance_analysis_subdirs.get(pyssa_keys.DISTANCE_ANALYSIS_RESULTS_SUBDIR)
        self.pymol_session_filepath = path_util.FilePath(f"{self.distance_analysis_subdirs.get(pyssa_keys.DISTANCE_ANALYSIS_SESSION_SUBDIR)}/{self.protein_pair_for_analysis.name}_analysis_session.pse")

    def save_distance_analysis_session(self) -> None:
        """This function saves the pymol session of the Protein pair distance analysis.

        """
        protein_pair_operations.save_session_of_protein_pair(self.pymol_session_filepath)

    def create_align_selections(self,
                                protein_1_selection: 'selection.Selection',
                                protein_2_selection: 'selection.Selection') -> None:
        """This function creates the selection which are needed for the align command.

        """
        # <editor-fold desc="Checks">
        if protein_1_selection.molecule_object != self.protein_pair_for_analysis.protein_1.get_molecule_object():
            logger.error("Selection is illegal.")
            raise ValueError("Selection is illegal.")
        if protein_2_selection.molecule_object != self.protein_pair_for_analysis.protein_2.get_molecule_object():
            logger.error("Selection is illegal.")
            raise ValueError("Selection is illegal.")

        # </editor-fold>

        self.protein_pair_for_analysis.protein_1.pymol_selection = protein_1_selection
        self.protein_pair_for_analysis.protein_2.pymol_selection = protein_2_selection

    def align_protein_pair_for_analysis(self) -> tuple:
        """This function aligns the protein pair with the PyMOL align command.

        Returns:
            a tuple with the rmsd and aligned amino acids
        """
        results = protein_pair_operations.align_protein_pair(self.protein_pair_for_analysis.protein_1.pymol_selection.selection_string,
                                                             self.protein_pair_for_analysis.protein_2.pymol_selection.selection_string,
                                                             self.alignment_file_name)
        # save the align object from pymol as alignment file
        if not os.path.exists(pathlib.Path(f"{self.protein_pair_for_analysis.export_dirname}/alignment_files")):
            os.mkdir(pathlib.Path(f"{self.protein_pair_for_analysis.export_dirname}/alignment_files"))
        cmd.save(pathlib.Path(f"{self.protein_pair_for_analysis.export_dirname}/alignment_files/{self.alignment_file_name}.aln"))
        return results[0], results[1]

    def do_analysis_in_pymol(self) -> None:
        """This function does the distance analysis of the protein pair.

        """
        self.protein_pair_for_analysis.load_protein_pair_in_pymol() # This creates a new pymol session
        self.protein_pair_for_analysis.color_protein_pair()
        align_results = self.align_protein_pair_for_analysis()
        rmsd_dict = {
            "rmsd": str(round(align_results[0], 2)),
            "aligned_residues": str(align_results[1]),
        }
        rmsd_file = open(pathlib.Path(f"{self.protein_pair_for_analysis.export_dirname}/rmsd.json"), "w", encoding="utf-8")
        json.dump(rmsd_dict, rmsd_file, indent=4)
        rmsd_file.close()
        protein_pair_util.calculate_distance_between_ca_atoms(self.protein_pair_for_analysis, self.alignment_file_name)
        self.save_distance_analysis_session()

    def serialize_distance_analysis(self):
        """This function serialize the protein pair object

                """
        distance_analysis_serializer = filesystem_io.ObjectSerializer(self, self.distance_analysis_subdirs.get(
            pyssa_keys.DISTANCE_ANALYSIS_OBJECTS_SUBDIR), f"{self.protein_pair_for_analysis.name}_analysis")
        distance_analysis_dict = {
            "protein_pair_for_analysis_name": str(self.protein_pair_for_analysis.name),
            "distance_analysis_dirname": str(self.distance_analysis_subdirs.get(pyssa_keys.DISTANCE_ANALYSIS_SUBDIR)),
            "cutoff": str(self.cutoff),
            "cycles": str(self.cycles),
            "export_dirname": str(self.export_dirname),
            "pymol_session_filepath": str(self.pymol_session_filepath.get_filepath()),
        }
        distance_analysis_serializer.set_custom_object_dict(distance_analysis_dict)
        distance_analysis_serializer.serialize_object()
