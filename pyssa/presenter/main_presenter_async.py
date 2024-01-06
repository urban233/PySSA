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
import pathlib
from typing import TYPE_CHECKING

from pyssa.internal.data_structures import project, protein, structure_analysis
from pyssa.internal.portal import pymol_io, graphic_operations
from pyssa.io_pyssa import path_util
from pyssa.logging_pyssa import log_handlers
from pyssa.util import analysis_util, exception, exit_codes

if TYPE_CHECKING:
    from pyssa.internal.data_structures import settings

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
