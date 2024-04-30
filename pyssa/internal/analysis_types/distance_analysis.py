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
from typing import TYPE_CHECKING
from xml.etree import ElementTree
import numpy as np
from pyssa.io_pyssa.xml_pyssa import element_names
from pyssa.io_pyssa.xml_pyssa import attribute_names
from pyssa.logging_pyssa import log_handlers
from pyssa.util import protein_pair_util
from pyssa.io_pyssa import path_util
from pyssa.util import constants
from pyssa.util import exception
from pyssa.internal.data_structures import results

if TYPE_CHECKING:
    from pyssa.internal.data_structures import protein_pair
    from pyssa.internal.data_structures import settings
    from pyssa.internal.data_structures import selection

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class DistanceAnalysis:
    """This class contains information about the distance analysis."""

    # <editor-fold desc="Class attributes">
    """
    a pair of proteins which get analyzed
    """
    _protein_pair_for_analysis: "protein_pair.ProteinPair"
    """
    the settings from pyssa
    """
    app_settings: "settings.Settings"
    """
    the cutoff value for the align command from PyMOL
    """
    cutoff: float
    """
    the number of refinement cycles for the align command from PyMOL
    """
    cycles: int
    """
    the size of the figure
    """
    figure_size: tuple[float, float]
    """
    the filename of the alignment file which gets created during the align command
    """
    alignment_file_name: str = "aln"
    """
    a directory where all results related to the protein will be stored
    """
    export_dirname: pathlib.Path
    """
    the session filepath for the PyMOL session of the analysis
    """
    pymol_session_filepath: path_util.FilePath
    """
    """
    distance_analysis_data: dict[str, np.ndarray] = {}
    """
    """
    analysis_results: "results.DistanceAnalysisResults" = None
    """
    """
    rmsd_dict: dict = {
        "rmsd": 0.0,
        "aligned_aa": 0,
    }

    # </editor-fold>

    def __init__(
        self,
        protein_pair_for_analysis: "protein_pair.ProteinPair",
        a_cutoff: float, cycles: int,
    ) -> None:
        """Constructor.

        Args:
            protein_pair_for_analysis: A pair of proteins which get analyzed.
            app_settings: The settings from PySSA.
        """
        self._protein_pair_for_analysis: protein_pair.ProteinPair = protein_pair_for_analysis
        self.name = f"dist_analysis_{self._protein_pair_for_analysis.name}"
        self.cutoff: float = a_cutoff
        self.cycles: int = cycles
        self.figure_size = (11.0, 6.0)
        self.alignment_file_name = "aln"

    def create_align_selections(
        self,
        protein_1_selection: "selection.Selection",
        protein_2_selection: "selection.Selection",
    ) -> None:
        """This function creates the selection which are needed for the align command."""
        logger.debug(
            f"1st argument of <create_align_selections>: "
            f"{protein_1_selection.selection_string} {protein_1_selection}",
        )
        logger.debug(
            f"2nd argument of <create_align_selections>: "
            f"{protein_2_selection.selection_string} {protein_2_selection}",
        )
        # <editor-fold desc="Checks">
        if (
            protein_1_selection.molecule_object != self._protein_pair_for_analysis.protein_1.get_molecule_object()
        ):  # pylint: disable=line-too-long
            logger.error("Selection is illegal.")
            raise ValueError("Selection is illegal.")
        if (
            protein_2_selection.molecule_object != self._protein_pair_for_analysis.protein_2.get_molecule_object()
        ):  # pylint: disable=line-too-long
            logger.error("Selection is illegal.")
            raise ValueError("Selection is illegal.")

        # </editor-fold>

        self._protein_pair_for_analysis.protein_1.pymol_selection = protein_1_selection
        self._protein_pair_for_analysis.protein_2.pymol_selection = protein_2_selection
        logger.debug(
            f"Prot 1 sele in <create_align_selections>: "
            f"{self._protein_pair_for_analysis.protein_1.pymol_selection.selection_string} "
            f"{self._protein_pair_for_analysis.protein_1.pymol_selection}",
        )
        logger.debug(
            f"Prot 2 sele in <create_align_selections>: "
            f"{self._protein_pair_for_analysis.protein_2.pymol_selection.selection_string} "
            f"{self._protein_pair_for_analysis.protein_2.pymol_selection}",
        )