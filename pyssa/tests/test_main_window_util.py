import os
import requests

from pyssa.util import main_window_util
from pyssa.internal.data_structures import project
from pyssa.util import constants


class TestMainWindowUtil:
    """A collection of tests from the main window util module."""

    @staticmethod
    def test_add_protein_with_pdb_id_to_project() -> None:
        """Tests if a pdb id and filepath can be added to an existing project."""
        return_tuple: tuple[str, int] = ("3bmp", 4)
        tmp_project: "project.Project" = project.Project("pytest-1", constants.DEFAULT_WORKSPACE_PATH)
        tmp_project = main_window_util.add_protein_to_project(return_tuple, tmp_project)
        assert len(tmp_project.proteins) == 1

    @staticmethod
    def test_add_protein_from_filesystem_to_project() -> None:
        """Tests if a pdb id and filepath can be added to an existing project."""
        local_filepath = os.path.join(os.path.join(os.path.expanduser("~"), "Downloads"),
                                      'pytest protein kalata family with a lot of spaces.pdb')
        url = 'https://files.rcsb.org/download/1PT4.pdb'
        response = requests.get(url)
        if response.status_code == 200:
            with open(local_filepath, 'wb') as file:
                file.write(response.content)
        else:
            print(f"Failed to download file. Status code: {response.status_code}")
            assert False
        return_tuple: tuple[str, int] = (local_filepath, len(local_filepath))
        tmp_project: "project.Project" = project.Project("pytest-1", constants.DEFAULT_WORKSPACE_PATH)
        tmp_project = main_window_util.add_protein_to_project(return_tuple, tmp_project)
        # Clean up
        os.remove(local_filepath)
        # assertion
        assert len(tmp_project.proteins) == 1
