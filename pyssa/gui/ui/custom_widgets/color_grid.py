#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/zielesny/PySSA>
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
"""Module for the custom color gird widget."""
from PyQt5 import QtWidgets

from pyssa.internal.data_structures import chain


class PyMOLColorGrid(QtWidgets.QWidget):

    def __init__(self):
        super().__init__()

        self.last_clicked_color = ""
        self.last_selected_chain: "chain.Chain" = None

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

        self.set_all_tooltips()
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

    def set_all_tooltips(self) -> None:
        self.c_red.setToolTip("red")
        self.c_tv_red.setToolTip("tv_red")
        self.c_salomon.setToolTip("salmon")
        self.c_raspberry.setToolTip("raspberry")

        self.c_green.setToolTip("green")
        self.c_tv_green.setToolTip("tv_green")
        self.c_palegreen.setToolTip("palegreen")
        self.c_forest.setToolTip("forest")

        self.c_blue.setToolTip("blue")
        self.c_tv_blue.setToolTip("tv_blue")
        self.c_lightblue.setToolTip("lightblue")
        self.c_skyblue.setToolTip("skyblue")

        self.c_yellow.setToolTip("yellow")
        self.c_tv_yellow.setToolTip("tv_yellow")
        self.c_paleyellow.setToolTip("paleyellow")
        self.c_sand.setToolTip("sand")

        self.c_magenta.setToolTip("magenta")
        self.c_purple.setToolTip("purple")
        self.c_pink.setToolTip("pink")
        self.c_hotpink.setToolTip("hotpink")

        self.c_cyan.setToolTip("cyan")
        self.c_aquamarine.setToolTip("aquamarine")
        self.c_palecyan.setToolTip("palecyan")
        self.c_teal.setToolTip("teal")

        self.c_orange.setToolTip("orange")
        self.c_tv_orange.setToolTip("tv_orange")
        self.c_lightorange.setToolTip("lightorange")
        self.c_olive.setToolTip("olive")

        self.c_white.setToolTip("white")
        self.c_grey_70.setToolTip("grey70")
        self.c_grey_30.setToolTip("grey30")
        self.c_black.setToolTip("black")

    def get_all_color_buttons(self) -> dict:
        return {
            "red": self.c_red,
            "tv_red": self.c_tv_red,
            "salmon": self.c_salomon,
            "raspberry": self.c_raspberry,
            "green": self.c_green,
            "tv_green": self.c_tv_green,
            "palegreen": self.c_palegreen,
            "forest": self.c_forest,
            "blue": self.c_blue,
            "tv_blue": self.c_tv_blue,
            "lightblue": self.c_lightblue,
            "skyblue": self.c_skyblue,
            "yellow": self.c_yellow,
            "tv_yellow": self.c_tv_yellow,
            "paleyellow": self.c_paleyellow,
            "sand": self.c_sand,
            "magenta": self.c_magenta,
            "pink": self.c_pink,
            "hotpink": self.c_hotpink,
            "cyan": self.c_cyan,
            "aquamarine": self.c_aquamarine,
            "palecyan": self.c_palecyan,
            "teal": self.c_teal,
            "orange": self.c_orange,
            "tv_orange": self.c_tv_orange,
            "lightorange": self.c_lightorange,
            "olive": self.c_olive,
            "white": self.c_white,
            "grey70": self.c_grey_70,
            "grey30": self.c_grey_30,
            "black": self.c_black
        }
