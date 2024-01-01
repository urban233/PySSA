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
"""This module contains tests about the different default values for the pages."""
import pathlib
import shutil
import sys
import unittest
from PyQt5.QtWidgets import QApplication, QListWidgetItem
from pyssa.main import MainWindow
from pyssa.tests import helpers
from pyssa.io_pyssa import filesystem_io
from PyQt5.QtTest import QTest
from PyQt5.QtCore import Qt


class TestNewPageFunctionality(unittest.TestCase):
    def setUp(self):
        self.app = QApplication(sys.argv)
        self.window = MainWindow()
        self.window.show()

    def tearDown(self):
        self.window.hide()
        del self.window
        del self.app

    def test_if_project_list_is_correctly_filled(self):
        """This function test if the project list, is filled correctly."""
        helpers.change_page(self.window.ui.btn_new_page, self.window.ui.lbl_page_title, "Create new project")
        scanner = filesystem_io.WorkspaceScanner(self.window.workspace_path)
        item_list = []
        for row in range(len(self.window.ui.list_new_projects)):
            item_list.append(self.window.ui.list_new_projects.item(row).text())
        self.window.ui.list_new_projects.clear()
        scanner.scan_workspace_for_valid_projects(self.window.ui.list_new_projects)
        new_item_list = []
        for row in range(len(self.window.ui.list_new_projects)):
            new_item_list.append(self.window.ui.list_new_projects.item(row).text())
        if new_item_list == item_list:
            self.assertTrue(True)
        else:
            self.assertTrue(False)

    def test_invalid_characters_in_project_name(self):
        """This function test if the gui page behaves correctly, if an invalid character is entered."""
        helpers.change_page(self.window.ui.btn_new_page, self.window.ui.lbl_page_title, "Create new project")
        chars_to_test = [Qt.Key_Space, Qt.Key_Exclam, Qt.Key_Udiaeresis]
        for tmp_char_to_test in chars_to_test:
            # space is an invalid character
            QTest.keyPress(self.window.ui.txt_new_project_name, tmp_char_to_test, Qt.NoModifier)
            # Check if the status label changes
            label_text_to_check = "Invalid character."
            helpers.check_if_label_changed_in_500_ms(self.window.ui.lbl_new_status_project_name, label_text_to_check)
            if self.window.ui.lbl_new_status_project_name.text() != label_text_to_check:
                self.assertTrue(False)
            if self.window.ui.cb_new_add_reference.isCheckable():
                self.assertTrue(False)
            if self.window.ui.btn_new_create_project.isEnabled():
                self.assertTrue(False)
            # cleaning the line edit
            QTest.keyPress(self.window.ui.txt_new_project_name, Qt.Key_Backspace, Qt.NoModifier)
        self.assertTrue(True)

    def test_project_already_exists_in_project_name(self):
        """This function test if the gui page behaves correctly, if the project name overlaps with an existing one."""
        # <editor-fold desc="Check if underscore project exists">
        helpers.change_page(self.window.ui.btn_new_page, self.window.ui.lbl_page_title, "Create new project")
        underscore_project_exists = False
        underscore_project_rows = []
        for row in range(len(self.window.ui.list_new_projects)):
            row_content = self.window.ui.list_new_projects.item(row).text()
            if self.window.ui.list_new_projects.item(row).text().find("test") != -1:
                underscore_project_exists = True
                underscore_project_rows.append(row)
        if underscore_project_exists:
            for tmp_row in underscore_project_rows:
                self.window.ui.list_new_projects.takeItem(tmp_row)
        # </editor-fold>
        self.window.ui.list_new_projects.addItem(QListWidgetItem("test1"))

        QTest.keyPress(self.window.ui.txt_new_project_name, Qt.Key_T, Qt.NoModifier)

        label_text_to_check = "Project name already exists."
        helpers.check_if_label_changed_in_500_ms(self.window.ui.lbl_new_status_project_name, label_text_to_check)
        # Check if the status label changes
        if self.window.ui.lbl_new_status_project_name.text() != label_text_to_check:
            self.assertTrue(False)
        # Check if the add reference checkbox can be used
        if self.window.ui.cb_new_add_reference.isCheckable():
            self.assertTrue(False)
        # Check if the button to create a new project, can be clicked
        if self.window.ui.btn_new_create_project.isEnabled():
            self.assertTrue(False)
        # Check if the correct row is highlighted
        self.assertEqual(
            self.window.ui.list_new_projects.item(self.window.ui.list_new_projects.currentRow()).text(),
            "test1",
        )
        self.assertTrue(True)

    def test_valid_characters_in_project_name(self):
        """This function test if the gui page behaves correctly, if a valid character is entered."""
        # <editor-fold desc="Check if underscore project exists">
        helpers.change_page(self.window.ui.btn_new_page, self.window.ui.lbl_page_title, "Create new project")
        underscore_project_exists = False
        underscore_project_rows = []
        for row in range(len(self.window.ui.list_new_projects)):
            index = self.window.ui.list_new_projects.item(row).text().find("_")
            if self.window.ui.list_new_projects.item(row).text().find("_") != -1:
                underscore_project_exists = True
                underscore_project_rows.append(row)
        if underscore_project_exists:
            for tmp_row in underscore_project_rows:
                self.window.ui.list_new_projects.takeItem(tmp_row)
        # </editor-fold>

        # <editor-fold desc="Check if dash project exists">
        dash_project_exists = False
        dash_project_rows = []
        for row in range(len(self.window.ui.list_new_projects)):
            if self.window.ui.list_new_projects.item(row).text().find("-") != -1:
                dash_project_exists = True
                dash_project_rows.append(row)
        if dash_project_exists:
            for tmp_row in dash_project_rows:
                self.window.ui.list_new_projects.takeItem(tmp_row)
        # </editor-fold>

        chars_to_test = [Qt.Key_Q, Qt.Key_Underscore, Qt.Key_Minus]
        for tmp_char_to_test in chars_to_test:
            # space is an invalid character
            QTest.keyPress(self.window.ui.txt_new_project_name, tmp_char_to_test, Qt.NoModifier)
            # Check if the status label changes
            label_text_to_check = ""
            helpers.check_if_label_changed_in_500_ms(self.window.ui.lbl_new_status_project_name, label_text_to_check)
            if self.window.ui.lbl_new_status_project_name.text() != label_text_to_check:
                self.assertTrue(False)
            if not self.window.ui.cb_new_add_reference.isCheckable():
                self.assertTrue(False)
            if not self.window.ui.btn_new_create_project.isEnabled():
                self.assertTrue(False)
            if not self.window.ui.cb_new_add_reference.isCheckable():
                self.assertTrue(False)
            if not self.window.ui.btn_new_create_project.isEnabled():
                self.assertTrue(False)
            # cleaning the line edit
            QTest.keyPress(self.window.ui.txt_new_project_name, Qt.Key_Backspace, Qt.NoModifier)
        self.assertTrue(True)

    def test_add_reference_to_project(self):
        helpers.change_page(self.window.ui.btn_new_page, self.window.ui.lbl_page_title, "Create new project")
        QTest.keyPress(self.window.ui.txt_new_project_name, Qt.Key_Q, Qt.NoModifier)
        self.window.ui.cb_new_add_reference.toggle()
        QTest.qWait(500)
        # Check if choose_reference label is invisible
        if not self.window.ui.lbl_new_choose_reference.isVisible():
            self.assertTrue(False)
        # Check if choose_reference text field is invisible
        if not self.window.ui.txt_new_choose_reference.isVisible():
            self.assertTrue(False)
        # Check if choose_reference tool button is invisible
        if not self.window.ui.btn_new_choose_reference.isVisible():
            self.assertTrue(False)
        self.assertTrue(True)

    def test_invalid_pdb_id_input(self):
        helpers.change_page(self.window.ui.btn_new_page, self.window.ui.lbl_page_title, "Create new project")
        QTest.keyPress(self.window.ui.txt_new_project_name, Qt.Key_Q, Qt.NoModifier)
        self.window.ui.cb_new_add_reference.toggle()
        QTest.qWait(500)
        # enters the pdb id 3bmp
        QTest.keyPress(self.window.ui.txt_new_choose_reference, Qt.Key_3, Qt.NoModifier)
        QTest.keyPress(self.window.ui.txt_new_choose_reference, Qt.Key_4, Qt.NoModifier)
        QTest.keyPress(self.window.ui.txt_new_choose_reference, Qt.Key_M, Qt.NoModifier)
        QTest.keyPress(self.window.ui.txt_new_choose_reference, Qt.Key_P, Qt.NoModifier)
        QTest.qWait(500)
        if self.window.ui.lbl_new_status_choose_reference.text() != "Invalid PDB ID.":
            self.assertTrue(False)
        if self.window.ui.btn_new_create_project.isEnabled():
            self.assertTrue(False)
        self.assertTrue(True)

    def test_valid_pdb_id_input(self):
        helpers.change_page(self.window.ui.btn_new_page, self.window.ui.lbl_page_title, "Create new project")
        QTest.keyPress(self.window.ui.txt_new_project_name, Qt.Key_Q, Qt.NoModifier)
        self.window.ui.cb_new_add_reference.toggle()
        QTest.qWait(500)
        # enters the pdb id 3bmp
        QTest.keyPress(self.window.ui.txt_new_choose_reference, Qt.Key_3, Qt.NoModifier)
        QTest.keyPress(self.window.ui.txt_new_choose_reference, Qt.Key_B, Qt.NoModifier)
        QTest.keyPress(self.window.ui.txt_new_choose_reference, Qt.Key_M, Qt.NoModifier)
        QTest.keyPress(self.window.ui.txt_new_choose_reference, Qt.Key_P, Qt.NoModifier)
        QTest.qWait(500)
        if self.window.ui.lbl_new_status_choose_reference.text() != "":
            self.assertTrue(False)
        if not self.window.ui.btn_new_create_project.isEnabled():
            self.assertTrue(False)
        self.assertTrue(True)

    def test_project_creation_without_reference(self):
        helpers.change_page(self.window.ui.btn_new_page, self.window.ui.lbl_page_title, "Create new project")
        QTest.keyPress(self.window.ui.txt_new_project_name, Qt.Key_Q, Qt.NoModifier)
        QTest.keyPress(self.window.ui.txt_new_project_name, Qt.Key_A, Qt.NoModifier)
        QTest.keyPress(self.window.ui.txt_new_project_name, Qt.Key_3, Qt.NoModifier)
        QTest.mouseClick(self.window.ui.btn_new_create_project, Qt.LeftButton)
        scanner = filesystem_io.WorkspaceScanner(self.window.workspace_path)
        projects_of_workspace = scanner.scan_workspace_for_valid_projects(self.window.ui.list_new_projects)
        if "qa3" not in projects_of_workspace:
            self.assertTrue(False)
        else:
            shutil.rmtree(pathlib.Path(f"{self.window.workspace_path}/qa3"))
            self.assertTrue(True)


if __name__ == "__main__":
    unittest.main()
