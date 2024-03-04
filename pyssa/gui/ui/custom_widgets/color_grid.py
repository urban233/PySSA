from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtCore import Qt


class PyMOLColorGrid(QtWidgets.QWidget):

    def __init__(self):
        super().__init__()
        grid = QtWidgets.QGridLayout()
        self.setLayout(grid)
        tmp_basic_stylesheet = """min-width: 30px; max-width: 30px; min-height: 30px; max-height: 30px;"""
        # Colors - Reds
        self.c_red = QtWidgets.QPushButton()
        self.c_tv_red = QtWidgets.QPushButton()
        self.c_salomon = QtWidgets.QPushButton()
        self.c_raspberry = QtWidgets.QPushButton()

        self.c_red.setStyleSheet(self.generate_color_stylesheet("#ff0000"))
        self.c_tv_red.setStyleSheet(self.generate_color_stylesheet("#ff3333"))
        self.c_salomon.setStyleSheet(self.generate_color_stylesheet("#ff9999"))
        self.c_raspberry.setStyleSheet(self.generate_color_stylesheet("#b24c66"))

        grid.addWidget(self.c_red, 0, 0)
        grid.addWidget(self.c_tv_red, 1, 0)
        grid.addWidget(self.c_salomon, 2, 0)
        grid.addWidget(self.c_raspberry, 3, 0)
        # Colors - Greens
        self.c_green = QtWidgets.QPushButton()
        self.c_tv_green = QtWidgets.QPushButton()
        self.c_palegreen = QtWidgets.QPushButton()
        self.c_forest = QtWidgets.QPushButton()

        self.c_green.setStyleSheet(self.generate_color_stylesheet("#00ff00"))
        self.c_tv_green.setStyleSheet(self.generate_color_stylesheet("#33ff33"))
        self.c_palegreen.setStyleSheet(self.generate_color_stylesheet("#a5e5a5"))
        self.c_forest.setStyleSheet(self.generate_color_stylesheet("#339933"))

        grid.addWidget(self.c_green, 0, 1)
        grid.addWidget(self.c_tv_green, 1, 1)
        grid.addWidget(self.c_palegreen, 2, 1)
        grid.addWidget(self.c_forest, 3, 1)
        # Colors - Blues
        self.c_blue = QtWidgets.QPushButton()
        self.c_tv_blue = QtWidgets.QPushButton()
        self.c_lightblue = QtWidgets.QPushButton()
        self.c_skyblue = QtWidgets.QPushButton()

        self.c_blue.setStyleSheet(self.generate_color_stylesheet("#0000ff"))
        self.c_tv_blue.setStyleSheet(self.generate_color_stylesheet("#4c4cff"))
        self.c_lightblue.setStyleSheet(self.generate_color_stylesheet("#bfbfff"))
        self.c_skyblue.setStyleSheet(self.generate_color_stylesheet("#337fcc"))

        grid.addWidget(self.c_blue, 0, 2)
        grid.addWidget(self.c_tv_blue, 1, 2)
        grid.addWidget(self.c_lightblue, 2, 2)
        grid.addWidget(self.c_skyblue, 3, 2)
        # Colors - Yellows
        self.c_yellow = QtWidgets.QPushButton()
        self.c_tv_yellow = QtWidgets.QPushButton()
        self.c_paleyellow = QtWidgets.QPushButton()
        self.c_sand = QtWidgets.QPushButton()

        self.c_yellow.setStyleSheet(self.generate_color_stylesheet("#ffff00"))
        self.c_tv_yellow.setStyleSheet(self.generate_color_stylesheet("#ffff33"))
        self.c_paleyellow.setStyleSheet(self.generate_color_stylesheet("#ffff7f"))
        self.c_sand.setStyleSheet(self.generate_color_stylesheet("#b78c4c"))

        grid.addWidget(self.c_yellow, 0, 3)
        grid.addWidget(self.c_tv_yellow, 1, 3)
        grid.addWidget(self.c_paleyellow, 2, 3)
        grid.addWidget(self.c_sand, 3, 3)
        # Colors - Magentas
        self.c_magenta = QtWidgets.QPushButton()
        self.c_purple = QtWidgets.QPushButton()
        self.c_pink = QtWidgets.QPushButton()
        self.c_hotpink = QtWidgets.QPushButton()

        self.c_magenta.setStyleSheet(self.generate_color_stylesheet("#ff00ff"))
        self.c_purple.setStyleSheet(self.generate_color_stylesheet("#bf00bf"))
        self.c_pink.setStyleSheet(self.generate_color_stylesheet("#ffa5d8"))
        self.c_hotpink.setStyleSheet(self.generate_color_stylesheet("#ff007f"))

        grid.addWidget(self.c_magenta, 0, 4)
        grid.addWidget(self.c_purple, 1, 4)
        grid.addWidget(self.c_pink, 2, 4)
        grid.addWidget(self.c_hotpink, 3, 4)
        # Colors - Cyan
        self.c_cyan = QtWidgets.QPushButton()
        self.c_aquamarine = QtWidgets.QPushButton()
        self.c_palecyan = QtWidgets.QPushButton()
        self.c_teal = QtWidgets.QPushButton()

        self.c_cyan.setStyleSheet(self.generate_color_stylesheet("#00ffff"))
        self.c_aquamarine.setStyleSheet(self.generate_color_stylesheet("#7fffff"))
        self.c_palecyan.setStyleSheet(self.generate_color_stylesheet("#ccffff"))
        self.c_teal.setStyleSheet(self.generate_color_stylesheet("#00bfbf"))

        grid.addWidget(self.c_cyan, 0, 5)
        grid.addWidget(self.c_aquamarine, 1, 5)
        grid.addWidget(self.c_palecyan, 2, 5)
        grid.addWidget(self.c_teal, 3, 5)
        # Colors - Orange
        self.c_orange = QtWidgets.QPushButton()
        self.c_tv_orange = QtWidgets.QPushButton()
        self.c_lightorange = QtWidgets.QPushButton()
        self.c_olive = QtWidgets.QPushButton()

        self.c_orange.setStyleSheet(self.generate_color_stylesheet("#ff7f00"))
        self.c_tv_orange.setStyleSheet(self.generate_color_stylesheet("#ff8c26"))
        self.c_lightorange.setStyleSheet(self.generate_color_stylesheet("#ffcc7f"))
        self.c_olive.setStyleSheet(self.generate_color_stylesheet("#c4b200"))

        grid.addWidget(self.c_orange, 0, 6)
        grid.addWidget(self.c_tv_orange, 1, 6)
        grid.addWidget(self.c_lightorange, 2, 6)
        grid.addWidget(self.c_olive, 3, 6)
        # Colors - Grays
        self.c_white = QtWidgets.QPushButton()
        self.c_grey_70 = QtWidgets.QPushButton()
        self.c_grey_30 = QtWidgets.QPushButton()
        self.c_black = QtWidgets.QPushButton()

        self.c_white.setStyleSheet(self.generate_color_stylesheet("#ffffff"))
        self.c_grey_70.setStyleSheet(self.generate_color_stylesheet("#b2b2b2"))
        self.c_grey_30.setStyleSheet(self.generate_color_stylesheet("#4c4c4c"))
        self.c_black.setStyleSheet(self.generate_color_stylesheet("#000000"))

        grid.addWidget(self.c_white, 0, 7)
        grid.addWidget(self.c_grey_70, 1, 7)
        grid.addWidget(self.c_grey_30, 2, 7)
        grid.addWidget(self.c_black, 3, 7)

        self.setGeometry(100, 100, self.minimumWidth(), self.minimumHeight())

    def generate_color_stylesheet(self, a_hex_color: str) -> str:
        stylesheet = """QPushButton {
                background-color: %s;
                border: none;
                border-radius: 4px;
                min-width: 20px;
                max-width: 20px;
                min-height: 20px;
                max-height: 20px;
            }
            QPushButton::hover {
                background-color: %s;
                border: solid;
                border-color: black;
                border-width: 2px;
                border-radius: 4px;
                min-width: 20px;
                max-width: 20px;
                min-height: 20px;
                max-height: 20px;
            }
        """ % (a_hex_color, a_hex_color)
        return stylesheet
