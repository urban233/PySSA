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
from pyssa.util import globals, enums

PLUGIN_NAME = "PySSA"
PLUGIN_PATH = globals.g_plugin_path
PLUGIN_ROOT_PATH = globals.g_plugin_root_path
PLUGIN_EXTRA_TOOLS_PATH = str(pathlib.Path(f"{PLUGIN_ROOT_PATH}/extra_tools/"))
PLUGIN_DOCS_PATH = str(pathlib.Path(f"{PLUGIN_ROOT_PATH}/docs/pyssa-documentation"))
PLUGIN_LOGO_FILEPATH = str(pathlib.Path(f"{PLUGIN_ROOT_PATH}/assets/images/pyssa_logo.png"))
PLUGIN_LOGO_WITH_FONT_FILEPATH = str(pathlib.Path(f"{PLUGIN_ROOT_PATH}/assets/images/logo_type_2.tiff"))
VERSION_NUMBER = "v0.10.15"
PLUGIN_PATH_WSL_NOTATION = (
    "/mnt/c/ProgramData/pyssa/mambaforge_pyssa/pyssa-mamba-env/Lib/site-packages/pymol/pymol_path/data/startup/PySSA"
)
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
SCRATCH_DIR_STRUCTURE_ALN_IMAGES_INTERESTING_REGIONS_DIR = Path(
    f"{SCRATCH_DIR_STRUCTURE_ALN_IMAGES_DIR}/interesting_regions",
)
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
COLABFOLD_PREDICT_SCRIPT_OLD = (
    f"/mnt/c/Users/{os.getlogin()}/AppData/Roaming/pymol/startup/{PLUGIN_NAME}/scripts/unix/colabfold_predict.sh"
)
COLABFOLD_PREDICT_NO_TEMPLATES_SCRIPT_OLD = f"/mnt/c/Users/{os.getlogin()}/AppData/Roaming/pymol/startup/{PLUGIN_NAME}/scripts/unix/colabfold_predict_no_templates.sh"  # noqa: E501
INSTALLATION_COLABFOLD_SCRIPT_OLD = (
    f"/mnt/c/Users/{os.getlogin()}/AppData/Roaming/pymol/startup/{PLUGIN_NAME}/scripts/unix/installation_colabfold.sh"
)

COLABFOLD_PREDICT_SCRIPT = f"/mnt/c/Users/{os.getlogin()}/.pyssa/scripts/unix/colabfold_predict.sh"
COLABFOLD_PREDICT_NO_TEMPLATES_SCRIPT = (
    f"/mnt/c/Users/{os.getlogin()}/.pyssa/scripts/unix/colabfold_predict_no_templates.sh"
)
INSTALLATION_COLABFOLD_SCRIPT = f"/mnt/c/Users/{os.getlogin()}/.pyssa/scripts/unix/installation_colabfold.sh"
COLABFOLD_PREDICT_NO_AMBER_SCRIPT = f"/mnt/c/Users/{os.getlogin()}/.pyssa/scripts/unix/colabfold_predict_no_amber.sh"
COLABFOLD_PREDICT_NO_AMBER_AND_TEMPLATES_SCRIPT = (
    f"/mnt/c/Users/{os.getlogin()}/.pyssa/scripts/unix/colabfold_predict_no_amber_and_templates.sh"
)
COLABFOLD_LOG_FILE_PATH = pathlib.Path(f"{PREDICTION_PDB_DIR}/log.txt")
NOTEBOOK_RESULTS_ZIP_NAME = "prediction"

current_time = datetime.datetime.now()
LOG_FILENAME = f"{current_time.year}-{current_time.month:02d}-{current_time.day:02d}_{current_time.hour:02d}-{current_time.minute:02d}.log"  # noqa: E501
LOG_FILEPATH = pathlib.Path(f"{SETTINGS_DIR}/logs/{LOG_FILENAME}")
LOG_PATH = pathlib.Path(f"{SETTINGS_DIR}/logs")

TUTORIAL_PATH = "C:\\ProgramData\\pyssa\\tutorials"
DOCS_PATH = "C:\\ProgramData\\pyssa\\user_guide.pdf"

REMOVE_WSL_POWERSHELL = pathlib.Path(f"{PLUGIN_ROOT_PATH}/scripts/powershell/remove_wsl_env.ps1")
ADD_WSL_POWERSHELL = pathlib.Path(f"{PLUGIN_ROOT_PATH}/scripts/powershell/add_wsl_env.ps1")
# Constants for config file
# TODO: uncomment constant below before deployment
# WSL_CONF_PATH = f"/mnt/c/Users/{os.getlogin()}/AppData/Roaming/pymol/startup/{PLUGIN_NAME}/config/wsl/wsl.conf"
WSL_CONF_PATH = f"/mnt/c/ProgramData/pyssa/plugin/Miniconda3/envs/pyssa_colab/Lib/site-packages/pymol/pymol_path/data/startup/{PLUGIN_NAME}/config/wsl/wsl.conf"  # noqa: E501
WSL_DISTRO_NAME = "UbuntuColabfold"
WSL_STORAGE_PATH = pathlib.Path("C:\\ProgramData\\pyssa\\wsl\\UbuntuColabfold")
# WSL_DISTRO_IMPORT_PATH = pathlib.Path(f"C:/Users/{os.getlogin()}/.pyssa/{WSL_DISTRO_NAME}.tar")
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

WINDOW_TITLE_OF_HELP_CENTER = "PySSA - Documentation Center"
WINDOW_TITLE_OF_PYSSA = "PySSA"
WINDOW_TITLE_OF_PYMOL_PART = "PyMOL"
HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH = f"{PLUGIN_PATH}\\winbatch\\bring_docs_to_front.exe"
ARRANGE_WINDOWS_EXE_FILEPATH = f"{PLUGIN_PATH}\\winbatch\\arrange_windows.exe"

# Paths of help html files
HELP_HOME_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/home.html")
HELP_ANALYSIS_IMAGES_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/analysis_images.html")
HELP_GLOBAL_SETTINGS_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/global_settings.html")
HELP_HOTSPOTS_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/hotspots.html")
HELP_IMAGE_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/image.html")
HELP_LOCAL_MONOMER_PREDICTION_HTML_PATH = pathlib.Path(
    f"{PLUGIN_PATH}/docs/internal_help/html/local_monomer_prediction.html",
)
HELP_LOCAL_MULTIMER_PREDICTION_HTML_PATH = pathlib.Path(
    f"{PLUGIN_PATH}/docs/internal_help/html/local_multimer_prediction.html",
)
HELP_MANAGE_PYMOL_SESSION_HTML_PATH = pathlib.Path(f"{PLUGIN_PATH}/docs/internal_help/html/manage_pymol_session.html")
HELP_MONOMER_PREDICTION_ANALYSIS_HTML_PATH = pathlib.Path(
    f"{PLUGIN_PATH}/docs/internal_help/html/monomer_prediction_analysis.html",
)
HELP_MULTIMER_PREDICTION_ANALYSIS_HTML_PATH = pathlib.Path(
    f"{PLUGIN_PATH}/docs/internal_help/html/multimer_prediction_analysis.html",
)
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

DEFAULT_HISTOGRAM_PROPERTIES = {
    enums.HistogramPropertiesEnum.X_AXIS_UNITS: 10,
    enums.HistogramPropertiesEnum.DISTANCE_INTERVAL: 1.0,
}

# possible setting parameters
STATUS_MESSAGE_TIMEOUT = 5000  # value in msec

# pymol parameters
PYMOL_COLORS = [
    "red", "tv_red", "salmon", "raspberry",
    "green", "tv_green", "palegreen", "forest",
    "blue", "tv_blue", "lightblue", "skyblue",
    "yellow", "tv_yellow", "paleyellow", "sand",
    "magenta", "purple", "pink", "hotpink",
    "cyan", "aquamarine", "palecyan", "teal",
    "orange", "tv_orange", "lightorange", "olive",
    "white", "grey70", "grey30", "black", "By RMSD"
]

PYMOL_COLORS_WITH_INDICES = {
    4: "red",
    32: "tv_red",
    9: "salmon",
    5268: "raspberry",
    3: "green",
    33: "tv_green",
    5259: "palegreen",
    22: "forest",
    2: "blue",
    34: "tv_blue",
    5263: "lightblue",
    5277: "skyblue",
    6: "yellow",
    35: "tv_yellow",
    5256: "paleyellow",
    5269: "sand",
    8: "magenta",
    19: "purple",
    48: "pink",
    12: "hotpink",
    5: "cyan",
    5257: "aquamarine",
    5265: "palecyan",
    20: "teal",
    13: "orange",
    37: "tv_orange",
    5264: "lightorange",
    18: "olive",
    0: "white",
    124: "grey70",
    84: "grey30",
    1: "black",
    "By RMSD": "By RMSD"
}

PYMOL_COLORS_WITH_INDICES_2 = {
    4: "red",
    3: "green",
    15: "limegreen",
    2: "blue",
    5277: "skyblue",
    6: "yellow",
    5276: "limon",
    8: "magenta",
    12: "hotpink",
    53: "violet",
    5: "cyan",
    5275: "greencyan",
    13: "orange",
    5264: "lightorange",
    0: "white",
}

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

chain_dict_reverse = {
    "A": 0,
    "B": 1,
    "C": 2,
    "D": 3,
    "E": 4,
    "F": 5,
    "G": 6,
    "H": 7,
    "I": 8,
    "J": 9,
    "K": 10,
    "L": 11,
    "M": 12,
    "N": 13,
    "O": 14,
    "P": 15,
    "Q": 16,
    "R": 17,
    "S": 18,
    "T": 19,
    "U": 20,
    "V": 21,
    "W": 22,
    "X": 23,
    "Y": 24,
    "Z": 25,
}

PYMOL_REPR_STATES_WITH_INDICES = {0: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 0}, 512: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 1}, 256: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 0}, 768: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 1}, 128: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 0}, 640: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 1}, 384: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 0}, 896: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 1}, 64: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 0}, 576: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 1}, 320: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 0}, 832: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 1}, 192: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 0}, 704: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 1}, 448: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 0}, 960: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 1}, 32: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 0}, 544: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 1}, 288: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 0}, 800: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 1}, 160: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 0}, 672: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 1}, 416: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 0}, 928: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 1}, 96: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 0}, 608: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 1}, 352: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 0}, 864: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 1}, 224: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 0}, 736: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 1}, 480: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 0}, 992: {'sticks': 0, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 1}, 4: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 0}, 516: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 1}, 260: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 0}, 772: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 1}, 132: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 0}, 644: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 1}, 388: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 0}, 900: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 1}, 68: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 0}, 580: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 1}, 324: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 0}, 836: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 1}, 196: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 0}, 708: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 1}, 452: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 0}, 964: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 1}, 36: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 0}, 548: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 1}, 292: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 0}, 804: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 1}, 164: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 0}, 676: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 1}, 420: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 0}, 932: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 1}, 100: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 0}, 612: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 1}, 356: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 0}, 868: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 1}, 228: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 0}, 740: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 1}, 484: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 0}, 996: {'sticks': 0, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 1}, 2: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 0}, 514: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 1}, 258: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 0}, 770: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 1}, 130: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 0}, 642: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 1}, 386: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 0}, 898: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 1}, 66: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 0}, 578: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 1}, 322: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 0}, 834: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 1}, 194: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 0}, 706: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 1}, 450: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 0}, 962: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 1}, 34: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 0}, 546: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 1}, 290: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 0}, 802: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 1}, 162: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 0}, 674: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 1}, 418: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 0}, 930: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 1}, 98: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 0}, 610: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 1}, 354: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 0}, 866: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 1}, 226: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 0}, 738: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 1}, 482: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 0}, 994: {'sticks': 0, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 1}, 6: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 0}, 518: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 1}, 262: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 0}, 774: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 1}, 134: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 0}, 646: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 1}, 390: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 0}, 902: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 1}, 70: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 0}, 582: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 1}, 326: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 0}, 838: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 1}, 198: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 0}, 710: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 1}, 454: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 0}, 966: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 1}, 38: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 0}, 550: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 1}, 294: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 0}, 806: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 1}, 166: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 0}, 678: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 1}, 422: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 0}, 934: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 1}, 102: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 0}, 614: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 1}, 358: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 0}, 870: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 1}, 230: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 0}, 742: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 1}, 486: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 0}, 998: {'sticks': 0, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 1}, 1: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 0}, 513: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 1}, 257: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 0}, 769: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 1}, 129: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 0}, 641: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 1}, 385: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 0}, 897: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 1}, 65: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 0}, 577: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 1}, 321: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 0}, 833: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 1}, 193: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 0}, 705: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 1}, 449: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 0}, 961: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 1}, 33: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 0}, 545: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 1}, 289: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 0}, 801: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 1}, 161: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 0}, 673: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 1}, 417: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 0}, 929: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 1}, 97: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 0}, 609: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 1}, 353: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 0}, 865: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 1}, 225: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 0}, 737: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 1}, 481: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 0}, 993: {'sticks': 1, 'spheres': 0, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 1}, 5: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 0}, 517: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 1}, 261: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 0}, 773: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 1}, 133: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 0}, 645: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 1}, 389: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 0}, 901: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 1}, 69: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 0}, 581: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 1}, 325: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 0}, 837: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 1}, 197: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 0}, 709: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 1}, 453: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 0}, 965: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 1}, 37: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 0}, 549: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 1}, 293: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 0}, 805: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 1}, 165: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 0}, 677: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 1}, 421: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 0}, 933: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 1}, 101: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 0}, 613: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 1}, 357: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 0}, 869: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 1}, 229: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 0}, 741: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 1}, 485: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 0}, 997: {'sticks': 1, 'spheres': 0, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 1}, 3: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 0}, 515: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 1}, 259: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 0}, 771: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 1}, 131: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 0}, 643: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 1}, 387: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 0}, 899: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 1}, 67: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 0}, 579: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 1}, 323: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 0}, 835: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 1}, 195: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 0}, 707: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 1}, 451: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 0}, 963: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 1}, 35: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 0}, 547: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 1}, 291: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 0}, 803: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 1}, 163: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 0}, 675: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 1}, 419: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 0}, 931: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 1}, 99: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 0}, 611: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 1}, 355: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 0}, 867: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 1}, 227: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 0}, 739: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 1}, 483: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 0}, 995: {'sticks': 1, 'spheres': 1, 'surface': 0, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 1}, 7: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 0}, 519: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 1}, 263: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 0}, 775: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 1}, 135: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 0}, 647: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 1}, 391: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 0}, 903: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 1}, 71: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 0}, 583: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 1}, 327: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 0}, 839: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 1}, 199: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 0}, 711: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 1}, 455: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 0}, 967: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 0, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 1}, 39: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 0}, 551: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 0, 'dots': 1}, 295: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 0}, 807: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 0, 'mesh': 1, 'dots': 1}, 167: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 0}, 679: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 0, 'dots': 1}, 423: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 0}, 935: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 0, 'lines': 1, 'mesh': 1, 'dots': 1}, 103: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 0}, 615: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 0, 'dots': 1}, 359: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 0}, 871: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 0, 'mesh': 1, 'dots': 1}, 231: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 0}, 743: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 0, 'dots': 1}, 487: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 0}, 999: {'sticks': 1, 'spheres': 1, 'surface': 1, 'cartoon': 1, 'ribbon': 1, 'lines': 1, 'mesh': 1, 'dots': 1}}

PYMOL_REPS_WITH_INDICES = {
        1: 'sticks',
        2: 'spheres',
        3: 'ball and stick',  # Nothing is known on how to hide this repr
        4: 'surface',
        32: 'cartoon',
        64: 'ribbon',
        128: 'lines',
        256: 'mesh',
        512: 'dots'
}