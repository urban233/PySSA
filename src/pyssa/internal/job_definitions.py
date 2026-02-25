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
"""Pure worker functions for job execution.

This module contains pure functions extracted from legacy job classes.
Each function accepts a set of inputs and produces outputs without external
side effects (no database updates, no UI signals, no logging).

These functions are designed to be used with the JobScheduler and are
picklable for execution in child processes via ProcessRuntime.
"""
import pathlib
import time
from typing import Any, TYPE_CHECKING
import subprocess

from src.pyssa.internal.data_structures import structure_prediction, structure_analysis, protein
from src.pyssa.internal.pymol import pml_worker, pml_enums
from src.pyssa.io_pyssa import bio_data
from src.pyssa.util import analysis_util, constants

if TYPE_CHECKING:
  from src.pyssa.internal.data_structures.data_classes import prediction_protein_info
  from src.pyssa.internal.data_structures.data_classes import prediction_configuration
  from src.pyssa.internal.data_structures import project
  from src.pyssa.internal.data_structures import protein_pair


def run_prediction_job(
    a_list_of_prediction_protein_infos: list["prediction_protein_info.PredictionProteinInfo"],
    a_prediction_configuration: "prediction_configuration.PredictionConfiguration",
) -> dict[str, Any]:
  """Run structure prediction job.

  Args:
      a_list_of_prediction_protein_infos: Protein information for prediction.
      a_prediction_configuration: Configuration for prediction.

  Returns:
      Dictionary containing:
          - "success": bool indicating if job succeeded
          - "best_prediction_models": list of best prediction model paths
          - "error_message": str with error details if failed
  """
  structure_prediction_obj = structure_prediction.StructurePrediction(
      a_list_of_prediction_protein_infos,
      a_prediction_configuration,
  )

  structure_prediction_obj.create_tmp_directories()
  structure_prediction_obj.create_fasta_files_for_prediction()
  structure_prediction_obj.run_prediction()
  best_prediction_models = structure_prediction_obj.move_best_prediction_models()

  tmp_proteins = []
  for tmp_prediction_model in best_prediction_models:
    tmp_protein = protein.Protein(tmp_prediction_model[0].name)
    tmp_pdb_filepath = pathlib.Path(
      f"{pathlib.Path(constants.PREDICTION_PDB_DIR)}/{tmp_prediction_model[0].name}.pdb"
    )
    tmp_protein.add_protein_structure_data_from_local_pdb_file(tmp_pdb_filepath)
    bio_data.build_pdb_file(tmp_protein.get_pdb_data(), str(tmp_pdb_filepath))
    tmp_reply_data = pml_worker.PmlWorker.one_shot_do(
      pml_enums.PmlCommand.CREATE_NEW_SESSION,
      args=(str(tmp_pdb_filepath),)
    )
    tmp_protein.pymol_session = tmp_reply_data
    tmp_proteins.append(tmp_protein)

  subprocess.run(["wsl", "--shutdown"], creationflags=subprocess.CREATE_NO_WINDOW)

  return {
      "success": True,
      "predicted_proteins": tmp_proteins,
  }


def run_distance_analysis_job(
    frozen_project: Any,
    list_with_analysis_names: list,
    cutoff: float,
    cycles: int,
) -> dict[str, bool | list["protein_pair.ProteinPair"]]:
  """Run distance analysis job.

  Args:
      frozen_project: Frozen snapshot of project state.
      list_with_analysis_names: List of analysis names to perform.
      cutoff: Distance cutoff value.
      cycles: Number of analysis cycles.

  Returns:
      Dictionary containing:
          - "success": bool indicating if job succeeded
          - "protein_pairs": list of analyzed protein pairs
          - "error_message": str with error details if failed
  """
  analysis_runs = structure_analysis.Analysis(frozen_project)
  analysis_runs.analysis_list = (
      analysis_util.transform_gui_input_to_practical_data(
          list_with_analysis_names,
          frozen_project,
          cutoff,
          cycles,
      )
  )

  analysis_runs.run_analysis("distance", False)

  protein_pairs = []
  for tmp_protein_pair in analysis_runs.analysis_list:
    tmp_protein_pair.db_project_id = frozen_project.get_id()
    protein_pairs.append(tmp_protein_pair)

  return {
      "success": True,
      "protein_pairs": protein_pairs,
  }


def run_prediction_and_distance_analysis_job(
    a_list_of_prediction_protein_infos: list["prediction_protein_info.PredictionProteinInfo"],
    a_prediction_configuration: "prediction_configuration.PredictionConfiguration",
    frozen_project: "project.Project",
    list_with_analysis_names: list,
    cutoff: float,
    cycles: int,
) -> dict[str, Any]:
  """Run combined prediction and distance analysis job.

  Args:
      a_list_of_prediction_protein_infos: Protein information for prediction.
      a_prediction_configuration: Configuration for prediction.
      frozen_project: Frozen snapshot of project state.
      list_with_analysis_names: List of analysis names to perform.
      cutoff: Distance cutoff value.
      cycles: Number of analysis cycles.

  Returns:
      Dictionary containing:
          - "success": bool indicating if job succeeded
          - "best_prediction_models": list of best prediction model paths
          - "protein_pairs": list of analyzed protein pairs
          - "error_message": str with error details if failed
  """
  prediction_result = run_prediction_job(
      a_list_of_prediction_protein_infos,
      a_prediction_configuration,
  )

  if not prediction_result["success"]:
    return prediction_result

  # Add newly predicted proteins to the frozen project so they can be analyzed
  for tmp_protein in prediction_result["predicted_proteins"]:
    frozen_project.add_existing_protein(tmp_protein)

  analysis_result = run_distance_analysis_job(
      frozen_project,
      list_with_analysis_names,
      cutoff,
      cycles,
  )

  if not analysis_result["success"]:
    return analysis_result

  return {
      "success": True,
      "predicted_proteins": prediction_result["predicted_proteins"],
      "protein_pairs": analysis_result["protein_pairs"],
  }


def run_ray_tracing_job(
    dest_image_filepath: str,
    cached_session_filepath: str,
    image_ray_trace_mode: int,
    image_ray_texture: int,
    image_renderer: str,
) -> dict[str, Any]:
  """Run ray tracing job.

  Note:
      This functionality is currently deprecated and does nothing.

  Args:
      dest_image_filepath: Destination path for rendered image.
      cached_session_filepath: Path to cached session file.
      image_ray_trace_mode: Ray trace mode setting.
      image_ray_texture: Ray texture setting.
      image_renderer: Renderer name.

  Returns:
      Dictionary containing:
          - "success": bool (currently always False)
          - "error_message": str indicating deprecation
  """
  with pml_worker.PmlWorker.session(cached_session_filepath) as worker:
    worker.do(
      pml_enums.PmlCommand.SET,
      args=("ray_trace_mode", image_ray_trace_mode)
    )
    worker.do(
      pml_enums.PmlCommand.SET,
      args=("ray_texture", image_ray_texture)
    )
    worker.do(
      pml_enums.PmlCommand.RAY,
      args=(1920, 1080, image_renderer),
      sync=True
    )
    worker.do(pml_enums.PmlCommand.PNG, args=(dest_image_filepath, 1920, 1080, 300, 1), sync=True)
  return {
      "success": False,
      "error_message": "RayTracingJob functionality is no longer supported.",
  }


def run_simple_image_job(
        dest_image_filepath: str,
        cached_session_filepath: str,
) -> dict[str, Any]:
  """Run a simple image job.

  Note:
      This functionality is currently deprecated and does nothing.

  Args:
      dest_image_filepath: Destination path for rendered image.
      cached_session_filepath: Path to cached session file.

  Returns:
      Dictionary containing:
          - "success": bool (currently always False)
          - "error_message": str indicating deprecation
  """
  with pml_worker.PmlWorker.session(cached_session_filepath) as worker:
    worker.do(pml_enums.PmlCommand.DRAW, args=(1920, 1080), sync=True)
    worker.do(pml_enums.PmlCommand.PNG, args=(dest_image_filepath, 1920, 1080), sync=True)
  return {
    "success": False,
    "error_message": "RayTracingJob functionality is no longer supported.",
  }
