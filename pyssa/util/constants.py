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
import logging
import os
import pathlib
from pathlib import Path

PLUGIN_NAME = "tmpPySSA"
PLUGIN_PATH = pathlib.Path(f"{os.path.expanduser('~')}\\AppData\\Roaming\\pymol\\startup\\{PLUGIN_NAME}")
PLUGIN_ROOT_PATH = os.path.realpath(os.path.join(os.path.dirname(__file__), "..", ".."))
VERSION_NUMBER = "v0.9.3"
# important PATHs
# settings path: /home/$USER/.pyssa/settings.xml
SETTINGS_DIR = Path(f"{os.path.expanduser('~')}/.pyssa/")
SETTINGS_FILE_NAME = "settings"
SETTINGS_FILENAME = "settings.json"
SETTINGS_FULL_FILEPATH = pathlib.Path(f"{SETTINGS_DIR}/{SETTINGS_FILENAME}")

DEFAULT_WORKSPACE_PATH = pathlib.Path(f"{os.path.expanduser('~')}/.pyssa/default_workspace")

SCRATCH_DIR = Path(f"{SETTINGS_DIR}/scratch")
CACHE_DIR = Path(f"{SETTINGS_DIR}/.cache")
CACHE_PROTEIN_DIR = Path(f"{CACHE_DIR}/pdb_files")
PREDICTION_FASTA_DIR = Path(f"{SCRATCH_DIR}/local_predictions/fasta")
PREDICTION_PDB_DIR = Path(f"{SCRATCH_DIR}/local_predictions/pdb")

POWERSHELL_EXE = "C:\\Windows\\System32\\WindowsPowerShell\\v1.0\\powershell.exe"
# Constants for structure prediction
CONVERT_DOS_TO_UNIX = pathlib.Path(f"{PLUGIN_ROOT_PATH}/scripts/powershell/convert_dos_to_unix.ps1")
# TODO: original paths, please uncomment before deployment!!!
# COLABFOLD_PREDICT_SCRIPT = f"/mnt/c/Users/{os.getlogin()}/AppData/Roaming/pymol/startup/{PLUGIN_NAME}/scripts/unix/colabfold_predict.sh"
# COLABFOLD_PREDICT_NO_TEMPLATES_SCRIPT = f"/mnt/c/Users/{os.getlogin()}/AppData/Roaming/pymol/startup/{PLUGIN_NAME}/scripts/unix/colabfold_predict_no_templates.sh"
# INSTALLATION_COLABFOLD_SCRIPT = f"/mnt/c/Users/{os.getlogin()}/AppData/Roaming/pymol/startup/{PLUGIN_NAME}/scripts/unix/installation_colabfold.sh"
COLABFOLD_PREDICT_SCRIPT = f"/mnt/c/Users/{os.getlogin()}/github_repos/{PLUGIN_NAME}/scripts/unix/colabfold_predict.sh"
COLABFOLD_PREDICT_NO_TEMPLATES_SCRIPT = f"/mnt/c/Users/{os.getlogin()}/github_repos/{PLUGIN_NAME}/scripts/unix/colabfold_predict_no_templates.sh"
INSTALLATION_COLABFOLD_SCRIPT = f"/mnt/c/Users/{os.getlogin()}/github_repos/{PLUGIN_NAME}/scripts/unix/installation_colabfold.sh"
NOTEBOOK_RESULTS_ZIP_NAME = "prediction"

# OFFICIAL_NOTEBOOK_NAME = "AlphaFold Colab"
# OFFICIAL_NOTEBOOK_URL = "https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb#scrollTo=rowN0bVYLe9n"
# NOTEBOOK_URL = "https://colab.research.google.com/drive/1bJXKZ9Fva7Rk0E4z5nS2wPdwwdnEevxb#scrollTo=CcOzpV-SHPrS"
REMOVE_WSL_POWERSHELL = pathlib.Path(f"{PLUGIN_ROOT_PATH}/scripts/powershell/remove_wsl_env.ps1")

# Constants for config file
# TODO: uncomment constant below before deployment
# WSL_CONF_PATH = f"/mnt/c/Users/{os.getlogin()}/AppData/Roaming/pymol/startup/{PLUGIN_NAME}/config/wsl/wsl.conf"
WSL_CONF_PATH = f"/mnt/c/Users/{os.getlogin()}/github_repos/{PLUGIN_NAME}/config/wsl/wsl.conf"

WSL_DISTRO_NAME = "UbuntuColabfold"
WSL_STORAGE_PATH = pathlib.Path(f"C:/Users/{os.getlogin()}/.pyssa/wsl/{WSL_DISTRO_NAME}")
WSL_DISTRO_IMPORT_PATH = pathlib.Path(f"C:/Users/{os.getlogin()}/.pyssa/{WSL_DISTRO_NAME}.tar")

# thread main tasks
PREDICTION_TASK = "Structure Prediction"
ANALYSIS_TASK = "Structure Analysis"

# all loggers used in the pyssa project
PYSSA_LOGGER = logging.getLogger("PySSA-Logger")
PREDICTION_WORKER_LOGGER = logging.getLogger("PredictionWorker")
ANALYSIS_WORKER_LOGGER = logging.getLogger("AnalysisWorker")

AMINO_ACID_CODE = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                   'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                   'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                   'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

CHAIN_TYPE_PROTEIN = "protein_chain"
CHAIN_TYPE_NON_PROTEIN = "non_protein_chain"

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

# Protein subdir dict keys


