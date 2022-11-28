import sys


class Helper:
    """This class contains useful functions for the test_main.py

    """
    def __init__(self, app, unittest_class, silent, main_window):
        self.app = app
        self.unittest_class = unittest_class
        self.silent = silent
        self.main_window = main_window

    def set_assertTrue(self, value: bool):
        if value is True:
            if self.silent == 0:
                try:
                    self.main_window.show()
                    sys.exit(self.app.exec_())
                except SystemExit:
                    self.unittest_class.assertTrue(True)
            else:
                self.unittest_class.assertTrue(True)
        else:
            if self.silent == 0:
                try:
                    self.main_window.show()
                    sys.exit(self.app.exec_())
                except SystemExit:
                    self.unittest_class.assertTrue(False)
            else:
                self.unittest_class.assertTrue(False)
