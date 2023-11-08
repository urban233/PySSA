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
"""Module for the protein pair class."""
import os
import pathlib
import logging
from pyssa.logging_pyssa import log_handlers
from pyssa.io_pyssa import path_util
from pyssa.internal.portal import pymol_io
from pyssa.internal.portal import protein_pair_operations
from pyssa.io_pyssa import filesystem_io
from pyssa.util import protein_util
from pyssa.util import pyssa_keys
from xml.etree import ElementTree
from pyssa.io_pyssa.xml_pyssa import element_names
from pyssa.io_pyssa.xml_pyssa import attribute_names
from typing import TYPE_CHECKING
from pyssa.util import constants
from pyssa.io_pyssa import binary_data

if TYPE_CHECKING:
    from pyssa.internal.data_structures import protein
    from pyssa.internal.analysis_types import distance_analysis

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class ProteinPair:
    """This class consists of two Protein objects. It is used to have a better workflow for the analysis."""

    # <editor-fold desc="Class attributes">
    """
    the first protein of the protein pair
    """
    protein_1: 'protein.Protein'
    """
    the second protein of the protein pair
    """
    protein_2: 'protein.Protein'
    """
    a directory where all results related to the protein will be stored
    """
    distance_analysis: 'distance_analysis.DistanceAnalysis' = None
    """
    the full filepath where the session file is stored
    """
    pymol_session_filepath: path_util.FilePath
    """
    a base64 string of the pymol session
    """
    pymol_session: str

    def __init__(self, protein_1: 'protein.Protein', protein_2: 'protein.Protein') -> None:
        """Constructor.

        Args:
            protein_1 (core.protein.Protein):
                reference Protein object
            protein_2 (core.Protein):
                model Protein object

        Raises:
            NotADirectoryError: If directory not found.
        """
        self.protein_1: protein.Protein = protein_1
        self.protein_2: protein.Protein = protein_2
        self.name = f"{self.protein_1.get_molecule_object()}_with_{self.protein_2.get_molecule_object()}"
        self.pymol_session = pymol_io.convert_pymol_session_to_base64_string(self.name)
        # protein_pair_dirname = f"{protein_pairs_dirname}/{self.name}"
        #
        # self.protein_pair_subdirs = {
        #     pyssa_keys.PROTEIN_PAIR_SUBDIR: pathlib.Path(f"{protein_pair_dirname}"),
        #     pyssa_keys.PROTEIN_PAIR_SESSION_SUBDIR: pathlib.Path(f"{protein_pair_dirname}/session"),
        #     pyssa_keys.PROTEIN_PAIR_RESULTS_SUBDIR: pathlib.Path(f"{protein_pair_dirname}/results"),
        #     pyssa_keys.PROTEIN_PAIR_OBJECTS_SUBDIR: pathlib.Path(f"{protein_pair_dirname}/.objects"),
        # }
        # for key in self.protein_pair_subdirs:
        #     if not os.path.exists(self.protein_pair_subdirs.get(key)):
        #         os.mkdir(self.protein_pair_subdirs.get(key))
        # self.export_dirname = self.protein_pair_subdirs.get(pyssa_keys.PROTEIN_PAIR_RESULTS_SUBDIR)
        # self.pymol_session_filepath = path_util.FilePath(f"{self.protein_pair_subdirs.get(pyssa_keys.PROTEIN_PAIR_SESSION_SUBDIR)}/{self.name}_session.pse")

    def load_protein_pair_in_pymol(self):
        """This function loads to proteins into a NEW pymol session."""
        self.protein_1.load_protein_in_pymol()
        self.protein_2.load_protein_in_pymol()

    def load_pymol_session(self):
        """This function loads the existing pymol session of the pair."""
        session_filepath = pathlib.Path(f"{constants.CACHE_PYMOL_SESSION_DIR}/{self.name}_session.pse")
        if not os.path.exists(constants.CACHE_PYMOL_SESSION_DIR):
            os.mkdir(constants.CACHE_PYMOL_SESSION_DIR)
        binary_data.write_binary_file_from_base64_string(filepath=session_filepath, base64_data=self.pymol_session)
        pymol_io.load_pymol_session(session_filepath)

    def color_protein_pair(self) -> None:
        """This function colors both the reference and the model Protein.

        Note:
            Only the official colors from PyMOL are supported. These can
            be looked up under the `color values`_ page.

        Args:
            color_ref (str, optional):
                defines color for the reference Protein
            color_model (str, optional):
                defines color for the model Protein

        Raises:
            pymol.CmdException:
                Exception is raised if one or both proteins
                does not exist as pymol objects.

        """
        protein_pair_operations.color_protein_pair(self.protein_1.get_molecule_object(),
                                                   self.protein_2.get_molecule_object())

    def save_session_of_protein_pair(self) -> None:
        """This function saves the pymol session of the Protein pair.

        Note:
            The pse file will be saved under the relative path
            (if export_data_dir = "data/results"):
            ``data/results/sessions``

            The file name (filename) MUST NOT have the file extension .pse!

        """
        self.pymol_session = protein_pair_operations.save_session_of_protein_pair(self.name)

    def set_distance_analysis(self, value):
        self.distance_analysis = value

    def serialize_protein_pair(self, xml_protein_pairs_element: ElementTree.Element):
        """This function serialize the protein pair object."""
        tmp_protein_pair = ElementTree.SubElement(xml_protein_pairs_element, element_names.PROTEIN_PAIR)
        tmp_protein_pair.set(attribute_names.PROTEIN_PAIR_NAME, str(self.name))
        tmp_protein_pair.set(attribute_names.PROTEIN_PAIR_PROT_1_MOLECULE_OBJECT,
                             str(self.protein_1.get_molecule_object()))
        tmp_protein_pair.set(attribute_names.PROTEIN_PAIR_PROT_1_ID, str(self.protein_1.get_id()))
        tmp_protein_pair.set(attribute_names.PROTEIN_PAIR_PROT_2_MOLECULE_OBJECT,
                             str(self.protein_2.get_molecule_object()))
        tmp_protein_pair.set(attribute_names.PROTEIN_PAIR_PROT_2_ID, str(self.protein_2.get_id()))
        tmp_session_data = ElementTree.SubElement(tmp_protein_pair, element_names.PROTEIN_PAIR_SESSION)
        tmp_session_data.set(attribute_names.PROTEIN_PAIR_SESSION, self.pymol_session)
        if self.distance_analysis is not None:
            self.distance_analysis.serialize_distance_analysis(tmp_protein_pair)

    @staticmethod
    def deserialize_protein_pair(protein_obj_json_file: path_util.FilePath):
        """This function constructs the protein pair object from
        the json file.

        Returns:
            two complete protein objects and a protein pair object deserialized from a json file
        """
        return filesystem_io.ObjectDeserializer(protein_obj_json_file.get_dirname(), protein_obj_json_file.get_filename()).deserialize_protein_pair()

    def create_plain_text_memory_mirror(self):
        mirror = [
            self.protein_1.get_molecule_object(),
            str(self.protein_1.pdb_filepath.get_filepath()),
            str(self.protein_1.export_dirname),
            str(self.protein_1.pdb_filepath.get_filename()),
            self.protein_1.pymol_selection.selection_string,
            protein_util.get_chains_as_list_of_tuples(self.protein_1.chains),
            self.protein_2.get_molecule_object(),
            str(self.protein_2.pdb_filepath.get_filepath()),
            str(self.protein_2.export_dirname),
            str(self.protein_2.pdb_filepath.get_filename()),
            self.protein_2.pymol_selection.selection_string,
            protein_util.get_chains_as_list_of_tuples(self.protein_2.chains),
            str(self.SCRATCH_DIR),
            str(self.pymol_session_filepath.get_filepath()),
            self.name,
            str(self.protein_pair_subdirs.get(pyssa_keys.PROTEIN_PAIR_SUBDIR)),
        ]
        return mirror