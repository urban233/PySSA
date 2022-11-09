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
"""Module for all styles related functions"""


def color_button_ready(button):
    """This functions colors a button blue to signal it is ready to be pressed

    Args:
        button:
            button ui element
    """
    with open('styles/styles_start_button_ready.css', 'r') as style_sheet_file:
        button_style = style_sheet_file.read()
        # Set the stylesheet of the application
        button.setStyleSheet(button_style)


def color_button_not_ready(button):
    """This functions colors a button white to signal it is NOT ready to be pressed

    Args:
        button:
            button ui element
    """
    with open('styles/styles_start_button_not_ready.css', 'r') as style_sheet_file:
        button_style = style_sheet_file.read()
        # Set the stylesheet of the application
        button.setStyleSheet(button_style)
