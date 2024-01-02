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
from PyQt5.QtTest import QTest
from PyQt5.QtCore import Qt
from pyssa.main import MainWindow
from pyssa.util import constants


class TestMainWindow(unittest.TestCase):
    """Test class."""

    project_name: str = "pyunit-test"
    project_path: str = pathlib.Path(f"{constants.DEFAULT_WORKSPACE_PATH}/{project_name}.xml")

    def setUp(self) -> None:
        """Sets up the gui for the tests."""
        self.app = QApplication(sys.argv)
        self.window = MainWindow()
        self.window.show()

    def tearDown(self) -> None:
        """Cleans all rests of the tests."""
        self.window.hide()
        del self.window
        del self.app
        if os.path.exists(self.project_path):
            os.remove(self.project_path)

    def test_create_new_project(self) -> None:
        """Tests if a new project can be created."""
        QTest.mouseClick(self.window.ui.btn_new_page, Qt.LeftButton)
        self.assertEqual(self.window.ui.lbl_page_title.text(), constants.PAGE_CREATE_NEW_PROJECT,
                         "Page names does not match!")
        self.window.ui.txt_new_project_name.setText(self.project_name)
        QTest.mouseClick(self.window.ui.btn_new_create_project, Qt.LeftButton)
        self.assertEqual(self.window.ui.lbl_page_title.text(), constants.PAGE_VIEW,
                         "Page names does not match!")
        self.assertTrue(os.path.exists(self.project_path))

    def test_create_open_project(self) -> None:
        """Tests if an existing project can be opened."""
        self.test_create_new_project()

        try:
            print(self.window.app_settings.ask_save_pymol_session)
            self.assertEqual(self.window.ui.lbl_page_title.text(), constants.PAGE_VIEW)
            QTest.mouseClick(self.window.ui.btn_close_project, Qt.LeftButton)
            self.assertEqual(self.window.ui.lbl_page_title.text(), constants.PAGE_HOME)
            QTest.mouseClick(self.window.ui.btn_open_page, Qt.LeftButton)
            self.assertEqual(self.window.ui.lbl_page_title.text(), constants.PAGE_OPEN_EXISTING_PROJECT)
            self.window.ui.txt_open_search.setText(self.project_name)
            QTest.mouseClick(self.window.ui.btn_open_open_project, Qt.LeftButton)
            self.assertEqual(self.window.ui.lbl_page_title.text(), constants.PAGE_VIEW)
            self.assertEqual(self.window.ui.lbl_current_project_name.text(), self.project_name)
        except Exception as e:
            print(e)


if __name__ == "__main__":
    unittest.main()
