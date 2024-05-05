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
"""A module to hold the global constants from the PySSA in the auxiliary scope."""
import pathlib
import os

SETTINGS_DIR = str(pathlib.Path(f"{os.path.expanduser('~')}/.pyssa/"))

SCRATCH_DIR = pathlib.Path(f"{SETTINGS_DIR}/scratch")
SCRATCH_DIR_ANALYSIS = pathlib.Path(f"{SCRATCH_DIR}/analysis")
SCRATCH_DIR_IMAGES = pathlib.Path(f"{SCRATCH_DIR_ANALYSIS}/images")
SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR = pathlib.Path(f"{SCRATCH_DIR_IMAGES}/structure_alignment")
SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR = pathlib.Path(
    f"{SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR}/interesting_regions",
)
CACHE_DIR = pathlib.Path(f"{SETTINGS_DIR}/.cache")
CACHE_PROTEIN_DIR = pathlib.Path(f"{CACHE_DIR}/pdb_files")
CACHE_PYMOL_SESSION_DIR = pathlib.Path(f"{CACHE_DIR}/sessions")
CACHE_CSV_DIR = pathlib.Path(f"{CACHE_DIR}/csv")
CACHE_IMAGES = pathlib.Path(f"{CACHE_DIR}/images")
CACHE_STRUCTURE_ALN_IMAGES_DIR = pathlib.Path(f"{CACHE_DIR}/images/structure_alignment")
CACHE_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR = pathlib.Path(f"{CACHE_STRUCTURE_ALN_IMAGES_DIR}/interesting_regions")

PYMOL_DEFAULT_BACKGROUND_COLOR = "black"
PYMOL_DEFAULT_RAY_TRACE_MODE = 1
PYMOL_DEFAULT_ANTIALIAS = 4
PYMOL_DEFAULT_AMBIENT = 0.5
PYMOL_DEFAULT_FANCY_HELICES = 1
PYMOL_DEFAULT_CARTOON_DISCRETE_COLORS = "on"
PYMOL_DEFAULT_CARTOON_SAMPLING = 50
PYMOL_DEFAULT_SPEC_POWER = 1500
PYMOL_DEFAULT_SPEC_REFLECT = 2.0
PYMOL_DEFAULT_RAY_TRANSPARENCY_CONTRAST = 0.20
PYMOL_DEFAULT_RAY_TRANSPARENCY_OBLIQUE = 1.0
PYMOL_DEFAULT_RAY_OBLIQUE_POWER = 20
PYMOL_DEFAULT_RAY_COLOR = "black"

DEFAULT_COLOR_PROTEIN_1 = "green"
DEFAULT_COLOR_PROTEIN_2 = "blue"

ARRAY_DISTANCE_INDEX = "index"
ARRAY_DISTANCE_PROT_1_CHAIN = "ref_chain"
ARRAY_DISTANCE_PROT_1_POSITION = "ref_pos"
ARRAY_DISTANCE_PROT_1_RESI = "ref_resi"
ARRAY_DISTANCE_PROT_2_CHAIN = "model_chain"
ARRAY_DISTANCE_PROT_2_POSITION = "model_pos"
ARRAY_DISTANCE_PROT_2_RESI = "model_resi"
ARRAY_DISTANCE_DISTANCES = "distance"

CHAIN_TYPE_PROTEIN = "protein_chain"
CHAIN_TYPE_NON_PROTEIN = "non_protein_chain"

AMINO_ACID_CODE = {
    "CYS": "C",
    "ASP": "D",
    "SER": "S",
    "GLN": "Q",
    "LYS": "K",
    "ILE": "I",
    "PRO": "P",
    "THR": "T",
    "PHE": "F",
    "ASN": "N",
    "GLY": "G",
    "HIS": "H",
    "LEU": "L",
    "ARG": "R",
    "TRP": "W",
    "ALA": "A",
    "VAL": "V",
    "GLU": "E",
    "TYR": "Y",
    "MET": "M",
}