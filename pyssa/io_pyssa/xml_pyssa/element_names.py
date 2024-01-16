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
"""XML element names which are used in PySSA."""
PROJECT = "project"
PROJECT_INFO = "project_info"
PROTEINS = "proteins"
PROTEIN_PAIRS = "protein_pairs"
PROTEIN = "protein"
PROTEIN_SESSION = "session_data"

PROTEIN_PAIR = "protein_pair"
PROTEIN_PAIR_SESSION = "session_data"

SEQUENCE = "sequence"
SEQUENCES = "sequences"

DISTANCE_ANALYSIS = "distance_analysis"
DISTANCE_ANALYSIS_RESULTS = "results"
DISTANCE_ANALYSIS_DISTANCE_RESULTS = "results_distance"

DISTANCE_ANALYSIS_INDEX_LIST = "index"
DISTANCE_ANALYSIS_PROT_1_CHAIN_LIST = "protein_1_chains"
DISTANCE_ANALYSIS_PROT_1_POSITION_LIST = "protein_1_positions"
DISTANCE_ANALYSIS_PROT_1_RESIDUE_LIST = "protein_1_residues"
DISTANCE_ANALYSIS_PROT_2_CHAIN_LIST = "protein_2_chains"
DISTANCE_ANALYSIS_PROT_2_POSITION_LIST = "protein_2_positions"
DISTANCE_ANALYSIS_PROT_2_RESIDUE_LIST = "protein_2_residues"
DISTANCE_ANALYSIS_DISTANCES_LIST = "distances"

DISTANCE_ANALYSIS_IMAGES = "auto_images"
DISTANCE_ANALYSIS_STRUCTURE_ALN_IMAGE = "structure_aln_image"
DISTANCE_ANALYSIS_ALN_IMAGES_INTERESTING_REGIONS = "interesting_reg_image"

DISTANCE_ANALYSIS_SESSION = "session_data"
