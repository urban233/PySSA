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
"""Module for all asynchronous functions used in the main presenter."""
import logging
import os
import pathlib
import shutil
import subprocess
from typing import TYPE_CHECKING

from pymol import cmd

from pyssa.internal.data_structures import project, protein, structure_analysis, structure_prediction
from pyssa.internal.data_structures.data_classes import (
    prediction_protein_info,
    prediction_configuration,
    current_session,
)
from pyssa.internal.portal import pymol_io, graphic_operations
from pyssa.io_pyssa import path_util, filesystem_io, bio_data
from pyssa.logging_pyssa import log_handlers
from pyssa.util import analysis_util, exception, exit_codes, constants

if TYPE_CHECKING:
    from pyssa.internal.data_structures import settings
    from pyssa.internal.data_structures import protein_pair

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def open_project(
    the_workspace_path: pathlib.Path,
    the_project_name: "project.Project",
    the_application_settings: "settings.Settings",
) -> tuple:
    """Opens a project and creates a new project object.

    Args:
        the_workspace_path: the current workspace path.
        the_project_name: the name of the project to open.
        the_application_settings: the settings of the application.

    Returns:
        a tuple with ("result", a_new_project_obj)
    """
    # TODO: needs checks
    # TODO: needs tests!
    tmp_project_path = pathlib.Path(f"{str(the_workspace_path)}/{the_project_name}")
    return ("result", project.Project.deserialize_project(tmp_project_path, the_application_settings))


def create_new_project(
    the_project_name: str,
    the_workspace_path: pathlib.Path,
    an_add_protein_flag: bool = False,
    protein_source_information: str = "",
) -> tuple:
    """Creates a new project object.

    Args:
        the_project_name: the project name to use for the new project.
        the_workspace_path: the current workspace path.
        an_add_protein_flag: a boolean value if a protein should be added or not.
        protein_source_information: either a filepath or a pdb id

    Returns:
        a tuple with ("results", a_new_project_obj)
    """
    # TODO: checks are needed
    tmp_project: "project.Project" = project.Project(the_project_name, the_workspace_path)
    if an_add_protein_flag:
        if len(protein_source_information) == 4:
            tmp_ref_protein = pymol_io.get_protein_from_pdb(protein_source_information.upper())
        else:
            # local pdb file as input
            pdb_filepath = path_util.FilePath(pathlib.Path(protein_source_information))
            graphic_operations.setup_default_session_graphic_settings()
            tmp_ref_protein = protein.Protein(
                molecule_object=pdb_filepath.get_filename(),
                pdb_filepath=pdb_filepath,
            )
        # fixme: should the disulfid-bonds be displayed globally
        # cmd.select(name="disulfides", selection="byres (resn CYS and name SG) within 2 of (resn CYS and name SG)")
        # cmd.color(color="atomic", selection="disulfides and not elem C")
        # cmd.set("valence", 0)  # this needs to be better implemented
        # cmd.show("sticks", "disulfides")
        tmp_project.add_existing_protein(tmp_ref_protein)
    tmp_project.serialize_project(
        pathlib.Path(f"{the_workspace_path}/{tmp_project.get_project_name()}.xml"),
    )
    return ("result", tmp_project)


def save_project(a_project: "project.Project", placeholder: int) -> tuple:
    """Saves the project through serialization.

    Args:
        a_project: a project that should be saved.
        placeholder: a placeholder argument so that python accepts the two arguments as tuple.

    Returns:
        a tuple with ("result", the_existing_project_object)
    """
    # TODO: needs checks
    # TODO: needs tests!
    a_project.serialize_project(a_project.get_project_xml_path())
    return ("result", a_project)


def clean_protein_new(
    a_protein_name: str,
    a_project: "project.Project",
) -> tuple:
    """Cleans a protein by creating a duplicate and removing all solvent and sugar molecules.

    Args:
        a_protein_name: the name of the protein to clean.
        a_project: the current project.

    Returns:
        a tuple with ("result", a_updates_project_object)
    """
    # TODO: needs checks
    tmp_protein = a_project.search_protein(
        a_protein_name,
    )
    clean_tmp_protein = tmp_protein.clean_protein(new_protein=True)
    constants.PYSSA_LOGGER.info("The protein %s has been cleaned.", clean_tmp_protein.get_molecule_object())
    a_project.add_existing_protein(clean_tmp_protein)
    a_project.serialize_project(a_project.get_project_xml_path())
    return ("result", a_project)


def clean_protein_update(a_protein_name: str, a_project: "project.Project") -> tuple:
    """Cleans a protein by removing all solvent and sugar molecules in the current molecule object.

    Args:
        a_protein_name: the name of the protein to clean.
        a_project: the current project.

    Returns:
        a tuple with ("result", a_updates_project_object)
    """
    # TODO: needs checks
    # TODO: needs tests!
    tmp_protein = a_project.search_protein(
        a_protein_name,
    )
    tmp_protein.clean_protein()
    constants.PYSSA_LOGGER.info("The protein %s has been cleaned.", tmp_protein.get_molecule_object())
    a_project.serialize_project(a_project.get_project_xml_path())
    return ("result", a_project)


def delete_protein(a_protein_name: str, a_project: "project.Project") -> tuple:
    """Deletes a certain protein from a project.

    Args:
        a_protein_name: the name of the protein to remove from the project.
        a_project: the current project.

    Returns:
        a tuple with ("result", an_existing_project_object)
    """
    # TODO: needs checks
    # TODO: needs tests!
    try:
        a_project.delete_specific_protein(a_protein_name)
    except ValueError:
        constants.PYSSA_LOGGER.error(
            "The protein %s could not be deleted, because it is not in the project.",
            a_protein_name,
        )
    return ("result", a_project)


def check_for_cleaning(a_protein_name: str, a_project: "project.Project") -> tuple:
    """Deletes a certain protein from a project.

    Args:
        a_protein_name: the name of the protein to remove from the project.
        a_project: the current project.

    Returns:
        a tuple with ("result", is_cleanable, is_in_protein_pair)
    """
    # TODO: needs checks
    # TODO: needs tests!
    is_cleanable: bool = False
    is_in_protein_pair: bool = False
    tmp_protein = a_project.search_protein(a_protein_name.replace(".pdb", ""))
    cmd.reinitialize()
    tmp_protein.load_protein_in_pymol()
    if cmd.select("organic") > 0 or cmd.select("solvent") > 0:
        is_cleanable = True
    if a_project.check_if_protein_is_in_any_protein_pair(a_protein_name) is True:
        is_in_protein_pair = True
    return ("result", is_cleanable, is_in_protein_pair)


def add_existing_protein_to_project(the_protein_information: tuple, a_project: "project.Project") -> tuple:
    """Adds a protein based on the filepath/ PDB id to a project.

    Args:
        the_protein_information (tuple[str, int]): a pair of either filepath or id and the length of the first one.
        a_project: a project where the protein should be added.

    Note:
        If a PDB id is used than the first element of the tuple contains the id and the second the length = 4.
        If a filepath is used than the first element of the tuple contains the id and the second the length of
        the entire filepath.
    """
    # TODO: checks are needed
    if the_protein_information[1] == 4:
        # PDB ID is used
        tmp_protein = pymol_io.get_protein_from_pdb(the_protein_information[0])
    else:
        # Filepath is used
        pdb_filepath: "path_util.FilePath" = path_util.FilePath(pathlib.Path(the_protein_information[0]))
        graphic_operations.setup_default_session_graphic_settings()
        tmp_protein_name: str = pdb_filepath.get_filename().replace(" ", "_")
        tmp_protein = protein.Protein(
            molecule_object=tmp_protein_name,
            pdb_filepath=pdb_filepath,
        )
    a_project.add_existing_protein(tmp_protein)
    return ("result", a_project)


def save_selected_protein_structure_as_pdb_file(
    a_protein_name: str,
    a_project: "project.Project",
    a_filepath: str,
) -> tuple:
    """Saves a given protein structure to a pdb file."""
    tmp_protein = a_project.search_protein(a_protein_name)
    try:
        bio_data.convert_xml_string_to_pdb_file(
            bio_data.convert_pdb_data_list_to_xml_string(tmp_protein.get_pdb_data()),
            pathlib.Path(a_filepath),
        )
    except Exception as e:
        logger.error(f"Saving protein to pdb file ended with error: {e}")
        return (exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0], exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[1])
    else:
        return (exit_codes.EXIT_CODE_ZERO[0], exit_codes.EXIT_CODE_ZERO[1])


def rename_selected_protein_structure(
    a_protein_name: str,
    the_new_protein_name: str,
    a_project: "project.Project",
) -> tuple:
    """Deletes a certain protein from a project.

    Args:
        a_protein_name: the name of the protein to rename.
        the_new_protein_name: the new name for the given protein.
        a_project: the current project.

    Returns:
        a tuple with ("result", an_existing_protein_object)
    """
    tmp_protein = a_project.search_protein(
        a_protein_name,
    )
    tmp_protein.set_molecule_object(the_new_protein_name)
    return ("result", tmp_protein)


def predict_protein_with_colabfold(
    the_prediction_protein_infos: list["prediction_protein_info.PredictionProteinInfo"],
    the_prediction_configuration: "prediction_configuration.PredictionConfiguration",
    a_project: "project.Project",
) -> tuple:
    """Runs structure prediction for a monomeric protein.

    Args:
        the_prediction_protein_infos: a list with protein names and sequences to predict.
        the_prediction_configuration: a prediction configuration for the monomeric protein prediction.
        a_project: a project object that contains all the proteins.

    Returns:
        a tuple with the exit code and exit code description
    """
    # TODO: needs checks
    # TODO: needs tests!
    structure_prediction_obj = structure_prediction.StructurePrediction(
        the_prediction_protein_infos,
        the_prediction_configuration,
        a_project,
    )
    structure_prediction_obj.create_tmp_directories()
    logger.info("Tmp directories were created.")

    # Create fasta files for prediction
    try:
        structure_prediction_obj.create_fasta_files_for_prediction()
    except exception.FastaFilesNotCreatedError:
        logger.error("Fasta files were not created.")
        return (exit_codes.ERROR_WRITING_FASTA_FILES[0], exit_codes.ERROR_WRITING_FASTA_FILES[1])
    except exception.FastaFilesNotFoundError:
        logger.error("Fasta files were not found.")
        return (exit_codes.ERROR_FASTA_FILES_NOT_FOUND[0], exit_codes.ERROR_FASTA_FILES_NOT_FOUND[1])
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return (exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0], exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[1])
    else:
        logger.info("Fasta files were successfully created.")

    # Run structure prediction
    try:
        structure_prediction_obj.run_prediction()
    except exception.PredictionEndedWithError:
        logger.error("Prediction ended with error.")
        return (exit_codes.ERROR_PREDICTION_FAILED[0], exit_codes.ERROR_PREDICTION_FAILED[1])
    else:
        logger.info("Prediction process finished.")

    try:
        structure_prediction_obj.move_best_prediction_models()
        logger.info("Saved predicted pdb file into XML file.")
    except exception.UnableToFindColabfoldModelError:
        logger.error("Could not move rank 1 model, because it does not exists.")
        return (
            exit_codes.ERROR_COLABFOLD_MODEL_NOT_FOUND[0],
            exit_codes.ERROR_COLABFOLD_MODEL_NOT_FOUND[1],
        )
    except FileNotFoundError:
        logger.error("Could not move rank 1 model, because it does not exists.")
        return (
            exit_codes.ERROR_COLABFOLD_MODEL_NOT_FOUND[0],
            exit_codes.ERROR_COLABFOLD_MODEL_NOT_FOUND[1],
        )
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        logger.error("Could not move rank 1 model, because it does not exists.")
        return (exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0], exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[1])
    else:
        subprocess.run(["wsl", "--shutdown"])
        logger.info("WSL gets shutdown.")
        return (exit_codes.EXIT_CODE_ZERO[0], exit_codes.EXIT_CODE_ZERO[1])


def run_distance_analysis(
    a_list_with_analysis_names: list,
    a_project: "project.Project",
    the_settings: "settings.Settings",
    an_make_images_flag: bool,
) -> tuple:
    """Runs the distance analysis for all protein pairs in the given job.

    Args:
        a_list_with_analysis_names: a list of the raw QListWidget analysis names.
        a_project: the current project.
        the_settings: the application settings.
        an_make_images_flag: True if images should be generated and False if not.

    Returns:
        a tuple containing exit codes.
    """
    # TODO: checks are needed
    logger.info("Running distance analysis in QThread using the Task class.")
    try:
        analysis_runs = structure_analysis.Analysis(a_project)
        analysis_runs.analysis_list = analysis_util.transform_gui_input_to_practical_data(
            a_list_with_analysis_names,
            a_project,
            the_settings,
        )
        analysis_runs.run_analysis("distance", an_make_images_flag)
    except exception.UnableToSetupAnalysisError:
        logger.error("Setting up the analysis runs failed therefore the distance analysis failed.")
        return (
            exit_codes.ERROR_DISTANCE_ANALYSIS_FAILED[0],
            exit_codes.ERROR_DISTANCE_ANALYSIS_FAILED[1],
        )
    except Exception as e:
        logger.error(f"Unknown error: {e}")
        return (exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[0], exit_codes.EXIT_CODE_ONE_UNKNOWN_ERROR[1])
    else:
        return (exit_codes.EXIT_CODE_ZERO[0], exit_codes.EXIT_CODE_ZERO[1])


def load_results(a_project: "project.Project", a_results_name: str) -> tuple:
    """Loads the results of a given protein pair name.

    Args:
        a_project: a project object containing the protein pair.
        a_results_name: the name of the protein pair.
    """
    shutil.rmtree(constants.CACHE_DIR)
    os.mkdir(constants.CACHE_DIR)
    os.mkdir(constants.CACHE_IMAGES)
    tmp_protein_pair = a_project.search_protein_pair(a_results_name)
    filesystem_io.XmlDeserializer(pathlib.Path(a_project.get_project_xml_path())).deserialize_analysis_images(
        tmp_protein_pair.name,
        tmp_protein_pair.distance_analysis.analysis_results,
    )
    if (
        len(tmp_protein_pair.distance_analysis.analysis_results.structure_aln_image) != 0
        and len(tmp_protein_pair.distance_analysis.analysis_results.interesting_regions_images) != 0
    ):
        # if both image types were made during analysis
        tmp_protein_pair.distance_analysis.analysis_results.create_image_png_files_from_base64()
        image_type = constants.IMAGES_ALL
    elif (
        len(tmp_protein_pair.distance_analysis.analysis_results.structure_aln_image) != 0
        and len(tmp_protein_pair.distance_analysis.analysis_results.interesting_regions_images) == 0
    ):
        # only struct align image were made
        tmp_protein_pair.distance_analysis.analysis_results.create_image_png_files_from_base64()
        image_type = constants.IMAGES_STRUCT_ALN_ONLY
    else:
        # no images were made
        image_type = constants.IMAGES_NONE

    cmd.reinitialize()
    tmp_protein_pair.load_pymol_session()
    tmp_current_session = current_session.CurrentSession(
        "protein_pair",
        tmp_protein_pair.name,
        tmp_protein_pair.pymol_session,
    )
    cmd.scene(
        f"{tmp_protein_pair.protein_1.get_molecule_object()}-{tmp_protein_pair.protein_2.get_molecule_object()}",
        action="recall",
    )
    return (
        "result",
        image_type,
        tmp_current_session,
        tmp_protein_pair.distance_analysis.analysis_results.rmsd,
        tmp_protein_pair.distance_analysis.analysis_results.aligned_aa,
    )


def color_protein_pair_by_rmsd_value(a_project: "project.Project", a_results_name: str) -> tuple:
    """Colors a given protein pair by their rmsd value.

    Args:
        a_project: a project object containing the protein pair.
        a_results_name: the name of the protein pair.

    Returns:
        a tuple with ("result", an_existing_protein_pair_object)
    """
    tmp_protein_pair: "protein_pair.ProteinPair" = a_project.search_protein_pair(a_results_name)
    graphic_operations.color_protein_pair_by_rmsd(tmp_protein_pair)
    return ("result", tmp_protein_pair)
