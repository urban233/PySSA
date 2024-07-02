#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
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
"""Module that holds enumerations for the auxiliary PyMOL which are the same as "enums"."""
import enum


class JobType(enum.Enum):
    """An enum for all job types."""
    PREDICTION = "structure prediction"
    DISTANCE_ANALYSIS = "distance analysis"
    PREDICTION_AND_DISTANCE_ANALYSIS = "prediction and distance analysis"
    RAY_TRACING = "ray-tracing"
    GENERAL_PURPOSE = "general-purpose"


class JobProgress(enum.Enum):
    """An enum for all job progress states."""
    WAITING = "waiting"
    RUNNING = "running"
    FINISHED = "finished"
    FAILED = "failed"


class JobDescriptionKeys(enum.Enum):
    """An enum for all job description keys.

    Notes:
        IMPORTANT: always use the .value of the enum! Otherwise, the connection to the auxiliary pymol would fail!
    """
    JOB_TYPE = "job_type"
    # Prediction keys
    PDB_FILEPATH = "pdb_filepath"
    # Distance analysis keys
    PROTEIN_PAIR_NAME = "the_protein_pair_name"
    PROTEIN_1_PDB_CACHE_FILEPATH = "a_protein_1_pdb_cache_filepath"
    PROTEIN_2_PDB_CACHE_FILEPATH = "a_protein_2_pdb_cache_filepath"
    PROTEIN_1_PYMOL_SELECTION_STRING = "a_protein_1_pymol_selection_string"
    PROTEIN_2_PYMOL_SELECTION_STRING = "a_protein_2_pymol_selection_string"
    CUTOFF = "a_cutoff"
    CYCLES = "the_cycles"
    # Ray-tracing keys
    IMAGE_DESTINATION_FILEPATH = "dest"
    CACHED_SESSION_FILEPATH = "cached"
    RAY_TRACE_MODE = "mode"
    RAY_TEXTURE = "texture"
    RAY_TRACING_RENDERER = "renderer"
    # General purpose
    JOB_SHORT_DESCRIPTION = "job_short_description"
    PYMOL_SESSION = "pymol_session"
    PROTEIN_NAME = "protein_name"


class JobShortDescription(enum.Enum):
    """An enum for all job short descriptions.

    Notes:
        IMPORTANT: always use the .value of the enum! Otherwise, the connection to the auxiliary pymol would fail!
    """
    RUN_STRUCTURE_PREDICTION = "Run ColabFold structure prediction."
    RUN_DISTANCE_ANALYSIS = "Run distance analysis."
    CREATE_RAY_TRACED_IMAGE = "Create new ray traced image."
    CREATE_NEW_PROTEIN_PYMOL_SESSION = "Create new protein pymol session."
    GET_ALL_CHAINS_OF_GIVEN_PROTEIN = "Get all chains of a given protein."
    CONSOLIDATE_MOLECULE_OBJECT_TO_FIRST_STATE = "Consolidate molecule object to first state."
    GET_ALL_SCENES_OF_SESSION = "Get all scenes of a given pymol session."
    CLEAN_PROTEIN_UPDATE_STRUCTURE = "Clean the existing protein structure."
