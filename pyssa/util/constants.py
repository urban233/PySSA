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
"""Module for storing project-wide constants."""
import logging
import os
import datetime
import pathlib
from pathlib import Path
from pyssa.util import globals

PLUGIN_NAME = "PySSA"
PLUGIN_PATH = globals.g_plugin_path
PLUGIN_ROOT_PATH = globals.g_plugin_root_path
PLUGIN_LOGO_FILEPATH = str(pathlib.Path(f"{PLUGIN_ROOT_PATH}/assets/images/pyssa_logo.png"))
VERSION_NUMBER = "v0.9.70"
PLUGIN_PATH_WSL_NOTATION = "/mnt/c/ProgramData/pyssa/mambaforge_pyssa/pyssa-mamba-env/Lib/site-packages/pymol/pymol_path/data/startup/PySSA"
# important PATHs
# settings path: /home/$USER/.pyssa/settings.xml
SETTINGS_DIR = str(pathlib.Path(f"{os.path.expanduser('~')}/.pyssa/"))
SETTINGS_DIR_UNIX_NOTATION = SETTINGS_DIR.replace("\\", "/")
SETTINGS_FILE_NAME = "settings"
SETTINGS_FILENAME = "settings.json"
SETTINGS_FULL_FILEPATH = pathlib.Path(f"{SETTINGS_DIR}/{SETTINGS_FILENAME}")

DEFAULT_WORKSPACE_PATH = pathlib.Path(f"{os.path.expanduser('~')}/.pyssa/default_workspace")

CONTAINER_NAME = "localcolabfold-container"
IMAGE_NAME = "localhost/localcolabfold-ubuntu2204:1.5.1.2"
SCRATCH_DIR = Path(f"{SETTINGS_DIR}/scratch")
SCRATCH_DIR_ANALYSIS = Path(f"{SCRATCH_DIR}/analysis")
SCRATCH_DIR_IMAGES = Path(f"{SCRATCH_DIR_ANALYSIS}/images")
SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR = Path(f"{SCRATCH_DIR_IMAGES}/structure_alignment")
SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR = Path(f"{SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR}/interesting_regions")
CACHE_DIR = Path(f"{SETTINGS_DIR}/.cache")
CACHE_PROTEIN_DIR = Path(f"{CACHE_DIR}/pdb_files")
CACHE_PYMOL_SESSION_DIR = Path(f"{CACHE_DIR}/sessions")
CACHE_CSV_DIR = Path(f"{CACHE_DIR}/csv")
CACHE_IMAGES = Path(f"{CACHE_DIR}/images")
CACHE_STRUCTURE_ALN_IMAGES_DIR = Path(f"{CACHE_DIR}/images/structure_alignment")
CACHE_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR = Path(f"{CACHE_STRUCTURE_ALN_IMAGES_DIR}/interesting_regions")
PREDICTION_FASTA_DIR = Path(f"{SCRATCH_DIR}/local_predictions/fasta")
WSL_SCRATCH_DIR = r"\\wsl.localhost\almaColabfold9\home\rhel_user\scratch"
WSL_PREDICTION_FASTA_DIR = r"\\wsl.localhost\almaColabfold9\home\rhel_user\scratch\local_predictions\fasta"
PREDICTION_PDB_DIR = Path(f"{SCRATCH_DIR}/local_predictions/pdb")
WSL_PREDICTION_PDB_DIR = r"\\wsl.localhost\almaColabfold9\home\rhel_user\scratch\local_predictions\pdb"

ESMFOLD_DIR = Path(f"{SCRATCH_DIR}/esmfold")
ESMFOLD_PDB_DIR = Path(f"{SCRATCH_DIR}/esmfold/pdb")

UNIX_SCRIPTS_SCIEBO_URL = "https://w-hs.sciebo.de/s/X3L7pnr4wfqy6gu/download"
UNIX_SCRIPTS_USER_DIR = pathlib.Path(f"{SETTINGS_DIR}/scripts/unix")

DEMO_PROJECT_SCIEBO_URL = "https://w-hs.sciebo.de/s/ZHJa6XB9SKWtqGi/download"

POWERSHELL_EXE = "C:\\Windows\\System32\\WindowsPowerShell\\v1.0\\powershell.exe"
# Constants for structure prediction
CONVERT_DOS_TO_UNIX = pathlib.Path(f"{PLUGIN_PATH}/scripts/batch/convert.bat")
CONVERT_DOS_TO_UNIX_EXE = pathlib.Path(f"{PLUGIN_PATH}/externals/dos2unix.exe")
UNIX_SCRIPTS_WIN = pathlib.Path(f"{PLUGIN_PATH}/scripts/unix/")
COLABFOLD_PREDICT_SCRIPT_WIN = pathlib.Path(f"{UNIX_SCRIPTS_WIN}/colabfold_predict.sh")

INSTALL_WSL_PS1 = pathlib.Path(f"{PLUGIN_ROOT_PATH}/scripts/powershell/install_wsl.ps1")
PREDICTION_PS1 = pathlib.Path(f"{PLUGIN_ROOT_PATH}/scripts/powershell/run_prediction.ps1")
INSTALL_WSL = pathlib.Path(f"{PLUGIN_ROOT_PATH}/scripts/batch/install_wsl.bat")
INSTALL_LOCAL_COLABFOLD_DISTRO = pathlib.Path(f"{PLUGIN_ROOT_PATH}/scripts/batch/import_distro.bat")
UNINSTALL_LOCAL_COLABFOLD_DISTRO = pathlib.Path(f"{PLUGIN_ROOT_PATH}/scripts/batch/uninstall_distro.bat")
# TODO: original paths, please uncomment before deployment!!!
COLABFOLD_PREDICT_SCRIPT_OLD = f"/mnt/c/Users/{os.getlogin()}/AppData/Roaming/pymol/startup/{PLUGIN_NAME}/scripts/unix/colabfold_predict.sh"
COLABFOLD_PREDICT_NO_TEMPLATES_SCRIPT_OLD = f"/mnt/c/Users/{os.getlogin()}/AppData/Roaming/pymol/startup/{PLUGIN_NAME}/scripts/unix/colabfold_predict_no_templates.sh"
INSTALLATION_COLABFOLD_SCRIPT_OLD = f"/mnt/c/Users/{os.getlogin()}/AppData/Roaming/pymol/startup/{PLUGIN_NAME}/scripts/unix/installation_colabfold.sh"

COLABFOLD_PREDICT_SCRIPT = f"/mnt/c/Users/{os.getlogin()}/.pyssa/scripts/unix/colabfold_predict.sh"
COLABFOLD_PREDICT_NO_TEMPLATES_SCRIPT = f"/mnt/c/Users/{os.getlogin()}/.pyssa/scripts/unix/colabfold_predict_no_templates.sh"
INSTALLATION_COLABFOLD_SCRIPT = f"/mnt/c/Users/{os.getlogin()}/.pyssa/scripts/unix/installation_colabfold.sh"
COLABFOLD_PREDICT_NO_AMBER_SCRIPT = f"/mnt/c/Users/{os.getlogin()}/.pyssa/scripts/unix/colabfold_predict_no_amber.sh"
COLABFOLD_PREDICT_NO_AMBER_AND_TEMPLATES_SCRIPT = f"/mnt/c/Users/{os.getlogin()}/.pyssa/scripts/unix/colabfold_predict_no_amber_and_templates.sh"
COLABFOLD_LOG_FILE_PATH = pathlib.Path(f"{PREDICTION_PDB_DIR}/log.txt")
#COLABFOLD_PREDICT_SCRIPT = f"/mnt/c/Users/{os.getlogin()}/github_repos/{PLUGIN_NAME}/scripts/unix/colabfold_predict.sh"
#COLABFOLD_PREDICT_NO_TEMPLATES_SCRIPT = f"/mnt/c/Users/{os.getlogin()}/github_repos/{PLUGIN_NAME}/scripts/unix/colabfold_predict_no_templates.sh"
#COLABFOLD_PREDICT_NO_AMBER_SCRIPT_OLD = f"/mnt/c/Users/{os.getlogin()}/github_repos/{PLUGIN_NAME}/scripts/unix/colabfold_predict_no_amber.sh"
#COLABFOLD_PREDICT_NO_AMBER_AND_TEMPLATES_SCRIPT_OLD = f"/mnt/c/Users/{os.getlogin()}/github_repos/{PLUGIN_NAME}/scripts/unix/colabfold_predict_no_amber_and_templates.sh"
#INSTALLATION_COLABFOLD_SCRIPT = f"/mnt/c/Users/{os.getlogin()}/github_repos/{PLUGIN_NAME}/scripts/unix/installation_colabfold.sh"
NOTEBOOK_RESULTS_ZIP_NAME = "prediction"

current_time = datetime.datetime.now()
LOG_FILENAME = f"{current_time.year}-{current_time.month:02d}-{current_time.day:02d}_{current_time.hour:02d}-{current_time.minute:02d}.log"
LOG_FILEPATH = pathlib.Path(f"{SETTINGS_DIR}/logs/{LOG_FILENAME}")
LOG_PATH = pathlib.Path(f"{SETTINGS_DIR}/logs")

TUTORIAL_PATH = "C:\\ProgramData\\pyssa\\tutorials"
DOCS_PATH = "C:\\ProgramData\\pyssa\\user_guide.pdf"

# OFFICIAL_NOTEBOOK_NAME = "AlphaFold Colab"
# OFFICIAL_NOTEBOOK_URL = "https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb#scrollTo=rowN0bVYLe9n"
# NOTEBOOK_URL = "https://colab.research.google.com/drive/1bJXKZ9Fva7Rk0E4z5nS2wPdwwdnEevxb#scrollTo=CcOzpV-SHPrS"
REMOVE_WSL_POWERSHELL = pathlib.Path(f"{PLUGIN_ROOT_PATH}/scripts/powershell/remove_wsl_env.ps1")
ADD_WSL_POWERSHELL = pathlib.Path(f"{PLUGIN_ROOT_PATH}/scripts/powershell/add_wsl_env.ps1")
# Constants for config file
# TODO: uncomment constant below before deployment
# WSL_CONF_PATH = f"/mnt/c/Users/{os.getlogin()}/AppData/Roaming/pymol/startup/{PLUGIN_NAME}/config/wsl/wsl.conf"
WSL_CONF_PATH = f"/mnt/c/ProgramData/pyssa/plugin/Miniconda3/envs/pyssa_colab/Lib/site-packages/pymol/pymol_path/data/startup/{PLUGIN_NAME}/config/wsl/wsl.conf"
WSL_DISTRO_NAME = "UbuntuColabfold"
WSL_STORAGE_PATH = pathlib.Path(f"C:\\ProgramData\\pyssa\\wsl\\UbuntuColabfold")
#WSL_DISTRO_IMPORT_PATH = pathlib.Path(f"C:/Users/{os.getlogin()}/.pyssa/{WSL_DISTRO_NAME}.tar")
WSL_DISTRO_IMPORT_PATH = pathlib.Path(f"C:/Users/{os.getlogin()}/.pyssa/")
DISTRO_DOWNLOAD_URL = "https://mega.nz/file/tz9wlLIQ#1qRxBdslCnOuUmLk2ytYHhSkItBsbuet3PTkZuvo-to"
WSL_DISK_PATH = pathlib.Path(f"{WSL_STORAGE_PATH}/ext4.vhdx")
# thread main tasks
PREDICTION_TASK = "Structure Prediction"
ANALYSIS_TASK = "Structure Analysis"

PREDICTION_TYPE_PRED = 1
PREDICTION_TYPE_PRED_MONO_ANALYSIS = 2
PREDICTION_TYPE_PRED_MULTI_ANALYSIS = 3

# all loggers used in the pyssa project
PYSSA_LOGGER = logging.getLogger("PySSA-Logger")
PREDICTION_WORKER_LOGGER = logging.getLogger("PredictionWorker")
ANALYSIS_WORKER_LOGGER = logging.getLogger("AnalysisWorker")

# docs paths
# TODO: get correct paths
DOCS_PDF = ""
DOCS_HTML = pathlib.Path(f"{PLUGIN_PATH}/docs/html/index.html")

# Paths of help html files
HELP_HOME_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/home.html")
HELP_ANALYSIS_IMAGES_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/analysis_images.html")
HELP_GLOBAL_SETTINGS_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/global_settings.html")
HELP_HOTSPOTS_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/hotspots.html")
HELP_IMAGE_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/image.html")
HELP_LOCAL_MONOMER_PREDICTION_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/local_monomer_prediction.html")
HELP_LOCAL_MULTIMER_PREDICTION_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/local_multimer_prediction.html")
HELP_MANAGE_PYMOL_SESSION_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/manage_pymol_session.html")
HELP_MONOMER_PREDICTION_ANALYSIS_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/monomer_prediction_analysis.html")
HELP_MULTIMER_PREDICTION_ANALYSIS_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/multimer_prediction_analysis.html")
HELP_VIEW_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/view_proteins_of_current_project.html")
HELP_CREATE_NEW_PROJECT_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/create_new_project.html")
HELP_DELETE_PROJECT_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/delete_project.html")
HELP_EDIT_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/edit_proteins_of_current_project.html")
HELP_OPEN_EXISTING_PROJECT_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/open_existing_project.html")
HELP_STRUCTURE_ANALYSIS_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/structure_analysis.html")
HELP_USE_EXISTING_PROJECT_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/use_existing_project.html")
HELP_RESULTS_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/results.html")
CHANGELOG_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/changelog/changelog.html")

# page header names
PAGE_HOME = "Home"
PAGE_ANALYSIS_IMAGES = "Analysis Images"
PAGE_HOTSPOTS = "Hotspots"
PAGE_IMAGE = "Image"
PAGE_LOCAL_MONOMER_PREDICTION = "Local Monomer Prediction"
PAGE_LOCAL_MULTIMER_PREDICTION = "Local Multimer Prediction"
PAGE_MANAGE_PYMOL_SESSION = "Manage PyMOL session"
PAGE_MONOMER_PREDICTION_ANALYSIS = "Monomer Prediction + Analysis"
PAGE_MULTIMER_PREDICTION_ANALYSIS = "Multimer Prediction + Analysis"
PAGE_VIEW = "View proteins of current project"
PAGE_CREATE_NEW_PROJECT = "Create new project"
PAGE_DELETE_PROJECT = "Delete existing project"
PAGE_EDIT = "Edit proteins of current project"
PAGE_OPEN_EXISTING_PROJECT = "Open existing project"
PAGE_STRUCTURE_ANALYSIS = "Structure Analysis"
PAGE_USE_EXISTING_PROJECT = "Use existing project"
PAGE_RESULTS = "Results"

# help paths
PAGE_HELP_PATHS_DICT = {
    PAGE_HOME: HELP_HOME_HTML_PATH,
    PAGE_ANALYSIS_IMAGES: HELP_ANALYSIS_IMAGES_HTML_PATH,
    PAGE_HOTSPOTS: HELP_HOTSPOTS_HTML_PATH,
    PAGE_IMAGE: HELP_IMAGE_HTML_PATH,
    PAGE_LOCAL_MONOMER_PREDICTION: HELP_LOCAL_MONOMER_PREDICTION_HTML_PATH,
    PAGE_LOCAL_MULTIMER_PREDICTION: HELP_LOCAL_MULTIMER_PREDICTION_HTML_PATH,
    PAGE_MANAGE_PYMOL_SESSION: HELP_MANAGE_PYMOL_SESSION_HTML_PATH,
    PAGE_MONOMER_PREDICTION_ANALYSIS: HELP_MONOMER_PREDICTION_ANALYSIS_HTML_PATH,
    PAGE_MULTIMER_PREDICTION_ANALYSIS: HELP_MULTIMER_PREDICTION_ANALYSIS_HTML_PATH,
    PAGE_VIEW: HELP_VIEW_HTML_PATH,
    PAGE_CREATE_NEW_PROJECT: HELP_CREATE_NEW_PROJECT_HTML_PATH,
    PAGE_DELETE_PROJECT: HELP_DELETE_PROJECT_HTML_PATH,
    PAGE_EDIT: HELP_EDIT_HTML_PATH,
    PAGE_OPEN_EXISTING_PROJECT: HELP_OPEN_EXISTING_PROJECT_HTML_PATH,
    PAGE_STRUCTURE_ANALYSIS: HELP_STRUCTURE_ANALYSIS_HTML_PATH,
    PAGE_USE_EXISTING_PROJECT: HELP_USE_EXISTING_PROJECT_HTML_PATH,
    PAGE_RESULTS: HELP_RESULTS_HTML_PATH,
}

# pymol parameters
PYMOL_COLORS = [
            "",
            "red",
            "green",
            "limegreen",
            "blue",
            "skyblue",
            "yellow",
            "limon",
            "magenta",
            "hotpink",
            "violet",
            "cyan",
            "greencyan",
            "orange",
            "lightorange",
            "white",
        ]

AMINO_ACID_CODE = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                   'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                   'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                   'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

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


CVM_NORMAL = "normal"
CVM_NORMAL_PROT_1_COLOR = "green"
CVM_NORMAL_PROT_2_COLOR = "blue"
CVM_DEUTERANOPIA = "Red-green (green weak, deuteranopia)"
CVM_DEUTERANOPIA_PROT_1_COLOR = "orange"
CVM_DEUTERANOPIA_PROT_2_COLOR = "blue"

CVM_PROTANOPIA = "Red-green (red weak, protanopia)"
CVM_PROTANOPIA_PROT_1_COLOR = "orange"
CVM_PROTANOPIA_PROT_2_COLOR = "blue"

CVM_TRITANOPIA = "Blue-yellow (tritanopia)"
CVM_TRITANOPIA_PROT_1_COLOR = "red"
CVM_TRITANOPIA_PROT_2_COLOR = "blue"


CHAIN_TYPE_PROTEIN = "protein_chain"
CHAIN_TYPE_NON_PROTEIN = "non_protein_chain"

ONLY_ANALYSIS = "only-analysis"
PREDICTION_ANALYSIS = "prediction-analysis"

IMAGES_ALL = "all images made"
IMAGES_STRUCT_ALN_ONLY = "only structure alignment image made"
IMAGES_NONE = "none images made"

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
