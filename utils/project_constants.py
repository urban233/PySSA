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

import os
from pathlib import Path

# important PATHs
# settings path: /home/$USER/.pyssa/settings.xml
SETTINGS_DIR = Path(f"{os.path.expanduser('~')}/.pyssa/")
SETTINGS_FILE = "settings.xml"
# SETTINGS = SETTINGS_DIR / SETTINGS_FILE
#
# # XML tags & attributes
# WORKSPACE_PATH_TAG = "workspacePath"
# PDB_STORAGE_PATH_TAG = "pdbPath"
# ZIP_STORAGE_PATH_TAG = "zipPath"
# CYCLES_VALUE_TAG = "cyclesValue"
# CUTOFF_VALUE_TAG = "cutoffValue"
# ATTRIBUTE = "value"

# Constants for the structure prediction
#FULL_FILENAME_PREDICTION_ZIP = Path(f"{os.path.expanduser('~')}/Downloads/prediction.zip")
OFFICIAL_NOTEBOOK_URL = "https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb#scrollTo=rowN0bVYLe9n"
NOTEBOOK_URL = "https://colab.research.google.com/drive/1bJXKZ9Fva7Rk0E4z5nS2wPdwwdnEevxb#scrollTo=CcOzpV-SHPrS"
# Constants for the structure analysis
