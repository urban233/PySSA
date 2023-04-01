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

import sys
import unittest
from PyQt5.QtWidgets import QApplication
from pyssa.main import MainWindow
from pyssa.tests import helpers


class TestMainWindow(unittest.TestCase):
    def setUp(self):
        self.app = QApplication(sys.argv)
        self.window = MainWindow()
        self.window.show()

    def tearDown(self):
        self.window.hide()
        del self.window
        del self.app
    
    def test_clicking_new_button_opens_create_new_project_page(self):
        # Click the button
        self.window.ui.btn_new_page.animateClick()

        label_text_to_check = "Create new project"
        # Wait for the label text to update
        helpers.check_if_label_changed_in_500_ms(self.window.ui.lbl_page_title, label_text_to_check)

        # Verify that the label text has been updated
        self.assertEqual(self.window.ui.lbl_page_title.text(), label_text_to_check,
                         "Page title does not match expected value.")

    def test_clicking_open_button_opens_open_project_page(self):
        # Click the button
        self.window.ui.btn_open_page.animateClick()

        label_text_to_check = "Open existing project"
        # Wait for the label text to update
        helpers.check_if_label_changed_in_500_ms(self.window.ui.lbl_page_title, label_text_to_check)

        # Verify that the label text has been updated
        self.assertEqual(self.window.ui.lbl_page_title.text(), label_text_to_check,
                         "Page title does not match expected value.")

    def test_clicking_delete_button_opens_delete_project_page(self):
        # Click the button
        self.window.ui.btn_delete_page.animateClick()

        label_text_to_check = "Delete existing project"
        # Wait for the label text to update
        helpers.check_if_label_changed_in_500_ms(self.window.ui.lbl_page_title, label_text_to_check)

        # Verify that the label text has been updated
        self.assertEqual(self.window.ui.lbl_page_title.text(), label_text_to_check,
                         "Page title does not match expected value.")

    def test_delete_page_defaults(self):
        """This test method checks if the delete page is set up with the correct default values

        """
        # Click the button
        self.window.ui.btn_delete_page.animateClick()

        # Check if page title is correct
        label_text_to_check = "Delete existing project"
        helpers.check_if_label_changed_in_500_ms(self.window.ui.lbl_page_title, label_text_to_check)

        # Check if search status label is invisible
        if self.window.ui.lbl_delete_status_search.text() != "":
            self.assertFalse(expr=True, msg="The search status label of the delete page is visible!")
        # Check if search text field is editable
        if not self.window.ui.txt_delete_search.isEnabled():
            self.assertFalse(expr=True, msg="The search text field of the delete page is not editable!")
        # Check if search text field is empty
        if self.window.ui.txt_delete_search.text() == " ":  # the space is necessary to detect the right case
            self.assertFalse(expr=True, msg="The selected project text field of the delete page is not empty!")
        # Check if project list is not empty
        if len(self.window.ui.list_delete_projects) == 0:
            self.assertFalse(expr=True, msg="The project list of the delete page is empty!")
        # Check if selected project text field is not editable
        if not self.window.ui.txt_delete_selected_projects.isReadOnly():
            self.assertFalse(expr=True, msg="The selected project text field of the delete page is editable!")
        # Check if selected project text field is empty
        if self.window.ui.txt_delete_selected_projects.text() == " ":  # the space is necessary to detect the right case
            self.assertFalse(expr=True, msg="The selected project text field of the delete page is not empty!")
        # Check if delete button is disabled
        if self.window.ui.btn_delete_delete_project.isEnabled():
            self.assertFalse(expr=True, msg="The delete project button of the delete page is enabled!")
        self.assertFalse(expr=False)

    def test_open_page_defaults(self):
        """This test method checks if the open page is set up with the correct default values

        """
        # Click the button
        self.window.ui.btn_open_page.animateClick()

        # Check if page title is correct
        label_text_to_check = "Open existing project"
        helpers.check_if_label_changed_in_500_ms(self.window.ui.lbl_page_title, label_text_to_check)

        # Check if search status label is invisible
        if self.window.ui.lbl_open_status_search.text() != "":
            self.assertFalse(expr=True, msg="The search status label of the open page is visible!")
        # Check if search text field is editable
        if not self.window.ui.txt_open_search.isEnabled():
            self.assertFalse(expr=True, msg="The search text field of the open page is not editable!")
        # Check if search text field is empty
        if self.window.ui.txt_open_search.text() == " ":  # the space is necessary to detect the right case
            self.assertFalse(expr=True, msg="The selected project text field of the open page is not empty!")
        # Check if project list is not empty
        if len(self.window.ui.list_open_projects) == 0:
            self.assertFalse(expr=True, msg="The project list of the open page is empty!")
        # Check if selected project text field is not editable
        if not self.window.ui.txt_open_selected_project.isReadOnly():
            self.assertFalse(expr=True, msg="The selected project text field of the open page is editable!")
        # Check if selected project text field is empty
        if self.window.ui.txt_open_selected_project.text() == " ":  # the space is necessary to detect the right case
            self.assertFalse(expr=True, msg="The selected project text field of the open page is not empty!")
        # Check if open button is disabled
        if self.window.ui.btn_open_open_project.isEnabled():
            self.assertFalse(expr=True, msg="The open project button of the open page is enabled!")
        self.assertFalse(expr=False)

    def test_new_page_defaults(self):
        """This test method checks if the new page is set up with the correct default values

        """
        # Click the button
        self.window.ui.btn_new_page.animateClick()

        # Check if page title is correct
        label_text_to_check = "Create new project"
        helpers.check_if_label_changed_in_500_ms(self.window.ui.lbl_page_title, label_text_to_check)

        # Check if search status label is invisible
        if self.window.ui.lbl_new_status_project_name.text() != "":
            self.assertFalse(expr=True, msg="The project name status label of the new page is visible!")
        # Check if project list is not empty
        if len(self.window.ui.list_new_projects) == 0:
            self.assertFalse(expr=True, msg="The project list of the new page is empty!")
        # Check if the add reference protein checkbox is not enabled
        if self.window.ui.cb_new_add_reference.isCheckable():
            self.assertFalse(expr=True, msg="The the add reference protein checkbox of the new page is enabled!")
        # Check if choose_reference label is invisible
        if self.window.ui.lbl_new_choose_reference.isVisible():
            self.assertFalse(expr=True)
        # Check if choose_reference text field is invisible
        if self.window.ui.txt_new_choose_reference.isVisible():
            self.assertFalse(expr=True)
        # Check if choose_reference tool button is invisible
        if self.window.ui.btn_new_choose_reference.isVisible():
            self.assertFalse(expr=True)
        # Check if create button is disabled
        if self.window.ui.btn_new_create_project.isEnabled():
            self.assertFalse(expr=True, msg="The create project button of the new page is enabled!")
        self.assertFalse(expr=False)


if __name__ == '__main__':
    unittest.main()
