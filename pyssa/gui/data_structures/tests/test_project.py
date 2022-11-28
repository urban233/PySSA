import unittest
from pyssa.gui.data_structures import project
from pyssa.gui.utilities import global_variables


class MyTestCase(unittest.TestCase):
    # basic variables for all tests
    workspace_path = global_variables.global_var_settings_obj.get_workspace_path()

    def test_create_project(self) -> None:
        """Tests the creation of a project object

        """
        test_project_obj = project.Project("test-project", self.workspace_path)
        if test_project_obj is not None:
            self.assertTrue(True)
        else:
            self.assertTrue(False)


if __name__ == '__main__':
    unittest.main()
