#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/zielesny/PySSA>
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
"""Module for the protein pair class."""
import copy
import logging
from typing import TYPE_CHECKING
from pyssa.logging_pyssa import log_handlers
from pyssa.io_pyssa import path_util
from pyssa.util import exception

if TYPE_CHECKING:
    from pyssa.internal.data_structures import protein
    from pyssa.internal.data_structures import structure_analysis
    from pyssa.internal.data_structures import sequence

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class ProteinPair:
    """This class consists of two Protein objects. It is used to have a better workflow for the analysis."""

    # <editor-fold desc="Class attributes">
    """
    the unique identifier of the protein
    """
    _id: int

    """
    the first protein of the protein pair
    """
    protein_1: "protein.Protein"

    """
    the second protein of the protein pair
    """
    protein_2: "protein.Protein"

    """
    a directory where all results related to the protein will be stored
    """
    distance_analysis: "structure_analysis.DistanceAnalysis" = None

    """
    the full filepath where the session file is stored
    """
    pymol_session_filepath: path_util.FilePath

    """
    a base64 string of the pymol session
    """
    pymol_session: str

    """
    the project id from the database
    """
    db_project_id: int

    def __init__(self, protein_1: "protein.Protein", protein_2: "protein.Protein") -> None:
        """Constructor.

        Args:
            protein_1 (core.protein.Protein):
                reference Protein object
            protein_2 (core.Protein):
                model Protein object

        Raises:
            IllegalArgumentError: If an argument is None.
        """
        if protein_1 is None:
            logger.error(f"The argument 'protein_1' is illegal: {protein_1}!")
            raise exception.IllegalArgumentError("An argument is illegal.")
        if protein_2 is None:
            logger.error(f"The argument 'protein_2' is illegal: {protein_2}!")
            raise exception.IllegalArgumentError("An argument is illegal.")

        self.protein_1: "protein.Protein" = copy.deepcopy(protein_1)
        self.protein_2: "protein.Protein" = copy.deepcopy(protein_2)
        self.name = f"{self.protein_1.get_molecule_object()}_with_{self.protein_2.get_molecule_object()}"
        self.pymol_session = ""

    def set_distance_analysis(self, a_value: "distance_analysis.DistanceAnalysis") -> None:
        """Sets the distance analysis object.

        Args:
            a_value (distance_analysis.DistanceAnalysis): a distance analysis object

        Raises:
            IllegalArgumentError: If the argument is None.
        """
        if a_value is None:
            logger.error(f"The argument 'a_value' is illegal: {a_value}!")
            raise exception.IllegalArgumentError("")
        self.distance_analysis = a_value

    def get_id(self):
        return self._id

    def set_id(self, value):
        self._id = value

    def get_protein_name_without_chains(self) -> str:
        return f"{self.protein_1.get_molecule_object()}-{self.protein_2.get_molecule_object()}"
