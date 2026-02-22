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
from typing import Any
import subprocess

from src.pyssa.internal.data_structures import structure_prediction, structure_analysis
from src.pyssa.util import analysis_util


def run_prediction_job(
    prediction_protein_infos: Any,
    prediction_configuration: Any,
) -> dict[str, Any]:
  """Run structure prediction job.

  Args:
      prediction_protein_infos: Protein information for prediction.
      prediction_configuration: Configuration for prediction.

  Returns:
      Dictionary containing:
          - "success": bool indicating if job succeeded
          - "best_prediction_models": list of best prediction model paths
          - "error_message": str with error details if failed
  """
  structure_prediction_obj = structure_prediction.StructurePrediction(
      prediction_protein_infos,
      prediction_configuration,
  )

  structure_prediction_obj.create_tmp_directories()
  structure_prediction_obj.create_fasta_files_for_prediction()
  structure_prediction_obj.run_prediction()
  best_prediction_models = structure_prediction_obj.move_best_prediction_models()

  subprocess.run(["wsl", "--shutdown"], creationflags=subprocess.CREATE_NO_WINDOW)

  return {
      "success": True,
      "best_prediction_models": best_prediction_models,
  }


def run_distance_analysis_job(
    frozen_project: Any,
    list_with_analysis_names: list,
    cutoff: float,
    cycles: int,
) -> dict[str, Any]:
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

  analysis_runs.run_analysis("distance", False, None, None)

  protein_pairs = []
  for tmp_protein_pair in analysis_runs.analysis_list:
    tmp_protein_pair.db_project_id = frozen_project.get_id()
    protein_pairs.append(tmp_protein_pair)

  return {
      "success": True,
      "protein_pairs": protein_pairs,
  }


def run_prediction_and_distance_analysis_job(
    prediction_protein_infos: Any,
    prediction_configuration: Any,
    frozen_project: Any,
    list_with_analysis_names: list,
    cutoff: float,
    cycles: int,
) -> dict[str, Any]:
  """Run combined prediction and distance analysis job.

  Args:
      prediction_protein_infos: Protein information for prediction.
      prediction_configuration: Configuration for prediction.
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
      prediction_protein_infos,
      prediction_configuration,
  )

  if not prediction_result["success"]:
    return prediction_result

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
      "best_prediction_models": prediction_result["best_prediction_models"],
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
  return {
      "success": False,
      "error_message": "RayTracingJob functionality is no longer supported.",
  }
