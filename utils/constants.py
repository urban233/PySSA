import os

# important PATHs
# settings path for linux local testing
# SETTINGS_DIR = "/home/matt/Documents/settings.xml"
# settings path for mac testing
# SETTINGS_DIR = "/Users/matt/Documents/settings.xml"
# settings path for windows
SETTINGS_DIR = f"{os.path.expanduser('~')}/.pyssa/settings.xml"
# settings path for pymol testing
# SETTINGS_DIR = f"{os.path.expanduser('~')}/anaconda3/envs/pymol_plugin/lib/python3.9/site-packages/pmg_tk/startup/pymol_plugin/settings/settings.xml"

# XML nodes & attributes
WORKSPACE_PATH_NODE = "workspacePath"
PDB_STORAGE_PATH_NODE = "pdbPath"
ZIP_STORAGE_PATH_NODE = "zipPath"
CYCLES_VALUE_NODE = "cyclesValue"
CUTOFF_VALUE_NODE = "cutoffValue"
ATTRIBUTE = "value"

# Constants for the structure prediction
FULL_FILENAME_PREDICTION_ZIP = f"{os.path.expanduser('~')}/Downloads/prediction.zip"
OFFICIAL_NOTEBOOK_URL = "https://colab.research.google.com/github/deepmind/alphafold/blob/main/notebooks/AlphaFold.ipynb#scrollTo=rowN0bVYLe9n"
NOTEBOOK_URL = "https://colab.research.google.com/drive/1bJXKZ9Fva7Rk0E4z5nS2wPdwwdnEevxb#scrollTo=CcOzpV-SHPrS"
# Constants for the structure analysis
