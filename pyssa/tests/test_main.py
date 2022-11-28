import os
import sys
import unittest
from PyQt5 import QtWidgets
from pyssa import main
from pyssa.gui.utilities import global_variables
from pyssa.tests import helper


class MyTestCase(unittest.TestCase):
    # basic variables for all tests
    silent = 1
    app = QtWidgets.QApplication(sys.argv)
    main_window = main.MainWindow()
    with open(os.path.join(global_variables.global_var_root_dir, "gui", "styles", "styles.css"), 'r', encoding="utf-8") as file:
        style = file.read()
        # Set the stylesheet of the application
        app.setStyleSheet(style)

    def test_open_open_page(self) -> None:
        """Tests if the open page can be opened.

        """
        self.helper = helper.Helper(self.app, self, self.silent, self.main_window)
        try:
            self.main_window.ui.btn_open_page.click()
            self.main_window.ui.txt_open_search.setText("gelis_wild")
        except:
            self.helper.set_assertTrue(False)
            return
        self.helper.set_assertTrue(True)

    def test_open_distance_plot(self) -> None:
        """Tests if the distance plot can be opened.

        Notes:
            If a FileNotFoundError is raised, it could be an invalid project.

        """
        self.helper = helper.Helper(self.app, self, self.silent, self.main_window)
        self.main_window = main.MainWindow()

        self.main_window.ui.btn_open_page.click()
        self.main_window.ui.txt_open_search.setText("gelis_wild")
        self.main_window.ui.btn_open_open_project.click()
        try:
            self.main_window.ui.btn_view_distance_plot.click()
        except:
            self.helper.set_assertTrue(False)
            return
        self.helper.set_assertTrue(True)


if __name__ == '__main__':
    unittest.main()
