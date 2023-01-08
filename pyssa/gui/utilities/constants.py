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
"""Module for storing project-wide constants"""

import os
import pathlib
from pathlib import Path


VERSION_NUMBER = "v0.9.2"
# important PATHs
# settings path: /home/$USER/.pyssa/settings.xml
SETTINGS_DIR = Path(f"{os.path.expanduser('~')}/.pyssa/")
SETTINGS_FILE = "settings.xml"
SETTINGS_FILENAME = "settings.json"
SETTINGS_FULL_FILEPATH = pathlib.Path(f"{SETTINGS_DIR}/{SETTINGS_FILENAME}")

DEFAULT_WORKSPACE_PATH = pathlib.Path(f"{os.path.expanduser('~')}/.pyssa/default_workspace")

SCRATCH_DIR = Path(f"{SETTINGS_DIR}/scratch")
PREDICTION_FASTA_DIR = Path(f"{SCRATCH_DIR}/local_predictions/fasta")
PREDICTION_PDB_DIR = Path(f"{SCRATCH_DIR}/local_predictions/pdb")

# Constants for the structure prediction
OFFICIAL_NOTEBOOK_NAME = "AlphaFold Colab"
OFFICIAL_NOTEBOOK_URL = "https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb#scrollTo=rowN0bVYLe9n"
NOTEBOOK_URL = "https://colab.research.google.com/drive/1bJXKZ9Fva7Rk0E4z5nS2wPdwwdnEevxb#scrollTo=CcOzpV-SHPrS"
NOTEBOOK_RESULTS_ZIP_NAME = "prediction"

# linux, macOS, Windows path
path_list = [
    f"{os.path.expanduser('~')}/anaconda3/envs/pymol_plugin/lib/python3.9/site-packages/pmg_tk/startup/tmpPySSA/",
    f"{os.path.expanduser('~')}/opt/anaconda3/envs/pyssa/lib/python3.9/site-packages/pmg_tk/startup/tmpPySSA/",
    f"{os.path.expanduser('~')}\\anaconda3\\envs\\pyssa\\lib\\site-packages\\pmg_tk\\startup\\pyssa\\",
]

chain_dict = {
    0: "A",
    1: "B",
    2: "C",
    3: "D",
    4: "E",
    5: "F",
    6: "G",
    7: "H",
    8: "I",
    9: "J",
    10: "K",
    11: "L",
    12: "M",
    13: "N",
    14: "O",
    15: "P",
    16: "Q",
    17: "R",
    18: "S",
    19: "T",
    20: "U",
    21: "V",
    22: "W",
    23: "X",
    24: "Y",
    25: "Z",
}