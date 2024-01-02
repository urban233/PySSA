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
"""Module for all styles related functions."""
import os
import pathlib
from PyQt5 import QtWidgets
from pyssa.util import constants, global_variables


def color_button_ready(button: QtWidgets.QPushButton) -> None:
    """This functions colors a button blue to signal it is ready to be pressed.

    Args:
        button:
            button ui element
    """
    with open(
        os.path.join(global_variables.global_var_root_dir, "gui", "ui", "styles", "styles_start_button_ready.css"),
        "r",
    ) as style_sheet_file:
        button_style = style_sheet_file.read()
        # Set the stylesheet of the application
        button.setStyleSheet(button_style)


def color_button_not_ready(button: QtWidgets.QPushButton) -> None:
    """This functions colors a button white to signal it is NOT ready to be pressed.

    Args:
        button:
            button ui element
    """
    with open(
        os.path.join(global_variables.global_var_root_dir, "gui", "ui", "styles", "styles_start_button_not_ready.css"),
        "r",
    ) as style_sheet_file:
        button_style = style_sheet_file.read()
        # Set the stylesheet of the application
        button.setStyleSheet(button_style)


def set_stylesheet(self) -> None:  # noqa: ANN001
    """Sets the style sheet to the QMainWindow or a QDialog.

    Args:
        self: a QMainWindow or QDialog
    """
    with open(
        pathlib.Path(f"{constants.PLUGIN_ROOT_PATH}/pyssa/gui/ui/styles/styles.css"),
        "r",
        encoding="utf-8",
    ) as file:
        style = file.read()
        # Set the stylesheet of the application
        self.setStyleSheet(style)


def color_sidebar_buttons(
    last_button: QtWidgets.QPushButton,
    active_button: QtWidgets.QPushButton,
) -> QtWidgets.QPushButton:
    """Colors the active sidebar button and the inactive sidebar button appropriately.

    Args:
        last_button: the inactive button on the sidebar.
        active_button: the active button on the sidebar.
    """
    with open(
        pathlib.Path(f"{constants.PLUGIN_ROOT_PATH}/pyssa/gui/ui/styles/sidebar_inactive.css"),
        "r",
        encoding="utf-8",
    ) as file:
        inactive_style = file.read()
        file.close()
    last_button.setStyleSheet(inactive_style)

    with open(
        pathlib.Path(f"{constants.PLUGIN_ROOT_PATH}/pyssa/gui/ui/styles/sidebar_active.css"),
        "r",
        encoding="utf-8",
    ) as file:
        active_style = file.read()
        file.close()
    active_button.setStyleSheet(active_style)
    return active_button
