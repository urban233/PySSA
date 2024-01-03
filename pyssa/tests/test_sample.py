import os.path
import pathlib

from pyssa.util import constants
from pyssa.main import MainWindow



def func(x):
    return x + 1


def test_answer():
    tmp_main_window = MainWindow()
    tmp_main_window.ui.txt_new_project_name.setText("pytest-test")
    tmp_main_window.create_new_project()
    assert os.path.exists(pathlib.Path(f"{constants.DEFAULT_WORKSPACE_PATH}/pytest-test.xml"))
