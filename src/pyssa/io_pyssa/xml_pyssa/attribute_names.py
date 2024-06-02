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
"""XML attribute names which are used in PySSA."""
ID = "id"

PROJECT_NAME = "name"
PROJECT_WORKSPACE_PATH = "workspace_path"
PROJECT_CREATION_DATE = "creation_date"
PROJECT_OS = "os"

PROTEIN_MOLECULE_OBJECT = "pymol_molecule_object"
PROTEIN_SELECTION = "pymol_selection"
PROTEIN_PDB_FILEPATH = "pdb_filepath"
PROTEIN_FASTA_FILEPATH = "fasta_filepath"
PROTEIN_SESSION_FILEPATH = "pymol_session_filepath"
PROTEIN_SESSION = "session"

PROTEIN_PAIR_NAME = "name"
PROTEIN_PAIR_PROT_1_MOLECULE_OBJECT = "protein_1_molecule_object"
PROTEIN_PAIR_PROT_1_ID = "protein_1_id"
PROTEIN_PAIR_PROT_2_MOLECULE_OBJECT = "protein_2_molecule_object"
PROTEIN_PAIR_PROT_2_ID = "protein_2_id"
PROTEIN_PAIR_SESSION = "session"

DISTANCE_ANALYSIS_NAME = "name"
DISTANCE_ANALYSIS_CUTOFF = "cutoff"
DISTANCE_ANALYSIS_CYCLES = "cycles"
DISTANCE_ANALYSIS_RMSD = "rmsd"
DISTANCE_ANALYSIS_ALIGNED_AA = "aligned_amino_acids"


DISTANCE_ANALYSIS_SESSION = "session"

DISTANCE_ANALYSIS_STRUCTURE_ALN_IMAGE_BASENAME = "basename"
DISTANCE_ANALYSIS_ALN_IMAGES_INTERESTING_REGIONS_BASENAME = "basename"

SEQUENCE_NAME = "name"
SEQUENCE_SEQ = "seq"
SEQUENCE_ID = "id"
