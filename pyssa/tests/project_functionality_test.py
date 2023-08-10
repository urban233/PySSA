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
"""Module which tests the functionality of the project system."""
import os.path
import pathlib
import sys
import time
import unittest
from PyQt5.QtWidgets import QApplication
from pyssa.main import MainWindow


class TestMainWindow(unittest.TestCase):
    """Test class."""
    def setUp(self) -> None:
        """Sets up the gui for the tests."""
        self.app = QApplication(sys.argv)
        self.window = MainWindow()
        self.window.show()

    def tearDown(self) -> None:
        """Cleanes all rests of the tests."""
        self.window.hide()
        del self.window
        del self.app

    def test_create_new_project(self) -> None:
        """Tests if a new project can be created."""
        # parameters for the test
        project_name = "Test-project"
        project_path = pathlib.Path(f"{self.window.workspace_path}/{project_name}.xml")
        # gui process
        self.window.ui.btn_new_page.animateClick()
        self.window.ui.txt_new_project_name.setText(project_name)
        self.window.ui.btn_new_create_project.click()
        time.sleep(5)
        # test
        check = os.path.exists(project_path)
        if check:
            os.remove(project_path)
        self.assertTrue(check)
    
    def test_create_open_project(self) -> None:
        """Tests if an existing project can be opened."""
        # parameters for the test
        project_name = "Test-project"
        project_path = pathlib.Path(f"{self.window.workspace_path}/{project_name}.xml")
        # gui process
        self.window.ui.btn_new_page.animateClick()
        self.window.ui.txt_new_project_name.setText(project_name)
        self.window.ui.btn_new_create_project.click()
        time.sleep(5)
        self.window.ui.btn_close_project.click()
        self.window.ui.btn_open_page.click()
        self.window.ui.txt_open_search.setText(project_name)
        self.window.ui.btn_open_open_project.click()
        # test
        if self.window.ui.lbl_current_project_name.text() == project_name:
            check = True
        else:
            check = False 
        if os.path.exists(project_path):
            os.remove(project_path)
        self.assertTrue(check)
    

if __name__ == '__main__':
    unittest.main()
