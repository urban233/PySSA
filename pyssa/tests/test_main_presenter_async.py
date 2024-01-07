import os.path
import pathlib

from pyssa.internal.data_structures import project, settings, protein
from pyssa.tests.util import basic_objects
from pyssa.tests.util import io_test
from pyssa.presenter import main_presenter_async
from pyssa.util import constants, exit_codes


class TestMainPresenterAsync:
    """A collection of tests for the main_presenter_async module."""

    def test_async_create_new_project_without_protein(self) -> None:
        """Tests if a project can be created without adding a protein."""
        main_presenter_async.create_new_project(
            "pytest-1",
            constants.DEFAULT_WORKSPACE_PATH,
        )
        tmp_project_filepath: pathlib.Path = pathlib.Path(f"{constants.DEFAULT_WORKSPACE_PATH}/pytest-1.xml")
        if os.path.exists(tmp_project_filepath):
            os.remove(tmp_project_filepath)
            assert True
        else:
            assert False

    def test_async_create_new_project_with_pdb_id(self) -> None:
        """Tests if a project can be created if a pdb id is provided."""
        tmp_result = main_presenter_async.create_new_project(
            "pytest-1",
            constants.DEFAULT_WORKSPACE_PATH,
            True,
            "1pt4",
        )
        tmp_project_filepath: pathlib.Path = pathlib.Path(f"{constants.DEFAULT_WORKSPACE_PATH}/pytest-1.xml")
        if os.path.exists(tmp_project_filepath) and len(tmp_result[1].proteins) == 1:
            os.remove(tmp_project_filepath)
            assert True
        else:
            assert False

    def test_async_create_new_project_with_local_pdb_file(self) -> None:
        """Tests if a project can be created if a local pdb file is provided."""
        tmp_pdb_filepath: str = io_test.download_1pt4_pdb_file()
        tmp_result = main_presenter_async.create_new_project(
            "pytest-1",
            constants.DEFAULT_WORKSPACE_PATH,
            True,
            tmp_pdb_filepath,
        )
        tmp_project_filepath: pathlib.Path = pathlib.Path(f"{constants.DEFAULT_WORKSPACE_PATH}/pytest-1.xml")
        if os.path.exists(tmp_pdb_filepath):
            os.remove(tmp_pdb_filepath)
        if os.path.exists(tmp_project_filepath) and len(tmp_result[1].proteins) == 1:
            os.remove(tmp_project_filepath)
            assert True
        else:
            assert False

    def test_run_distance_analysis_with_cyclotide(self) -> None:
        """Tests if the distance analysis runs successfully with a cyclotide."""
        tmp_project: "project.Project" = project.Project("pytest-1", constants.DEFAULT_WORKSPACE_PATH)
        # add the two proteins to the project
        tmp_prot_1: "protein.Protein" = basic_objects.create_protein_obj_from_pdb_id("1nb1")
        tmp_prot_2: "protein.Protein" = basic_objects.create_protein_obj_from_pdb_file(
            "https://w-hs.sciebo.de/s/oF7XbNNo3fJrCpf/download", "kB1_pred.pdb",
        )
        tmp_project.add_existing_protein(tmp_prot_1)
        tmp_project.add_existing_protein(tmp_prot_2)
        tmp_settings: "settings.Settings" = settings.Settings(constants.SETTINGS_DIR, constants.SETTINGS_FILENAME)
        tmp_analysis_names: list = ["kB1_pred;A_vs_1nb1;A"]
        result: tuple = main_presenter_async.run_distance_analysis(
            tmp_analysis_names,
            tmp_project,
            tmp_settings,
            False,
        )
        if os.path.exists(tmp_project.get_project_xml_path()):
            os.remove(tmp_project.get_project_xml_path())
        if result[0] == exit_codes.EXIT_CODE_ZERO[0]:
            assert True
        else:
            print(result)
            assert False

    def test_run_distance_analysis_with_bmp2_protein(self) -> None:
        """Tests if the distance analysis runs successfully with the bmp2 proteins."""
        tmp_project: "project.Project" = project.Project("pytest-1", constants.DEFAULT_WORKSPACE_PATH)
        # add the two proteins to the project
        tmp_prot_1: "protein.Protein" = basic_objects.create_protein_obj_from_pdb_id("3BMP")
        tmp_prot_2: "protein.Protein" = basic_objects.create_protein_obj_from_pdb_id("6OMN")
        tmp_project.add_existing_protein(tmp_prot_1)
        tmp_project.add_existing_protein(tmp_prot_2)
        tmp_settings: "settings.Settings" = settings.Settings(constants.SETTINGS_DIR, constants.SETTINGS_FILENAME)
        tmp_analysis_names: list = ["6OMN;E_vs_3BMP;A"]
        result: tuple = main_presenter_async.run_distance_analysis(
            tmp_analysis_names,
            tmp_project,
            tmp_settings,
            False,
        )
        if os.path.exists(tmp_project.get_project_xml_path()):
            os.remove(tmp_project.get_project_xml_path())
        if result[0] == exit_codes.EXIT_CODE_ZERO[0]:
            assert True
        else:
            assert False

    def test_clean_protein_new(self) -> None:
        """Tests the cleaning process if a new protein should be created."""
        tmp_project: "project.Project" = project.Project("pytest-1", constants.DEFAULT_WORKSPACE_PATH)
        # add the two proteins to the project
        tmp_prot_1: "protein.Protein" = basic_objects.create_protein_obj_from_pdb_id("6OMN")
        tmp_project.add_existing_protein(tmp_prot_1)
        result = main_presenter_async.clean_protein_new("6OMN", tmp_project)
        if os.path.exists(tmp_project.get_project_xml_path()):
            os.remove(tmp_project.get_project_xml_path())
        if len(result[1].proteins) == 2:
            assert True
        else:
            assert False

