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
"""Module for all styles related functions."""
import logging
import os
import pathlib
from PyQt5 import QtWidgets
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants, global_variables, exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


def color_bottom_frame_button(button: QtWidgets.QPushButton) -> None:
    """Colors the button in the style of the bottom frame button.

    Args:
        button (QtWidgets.QPushButton): The button to apply the style to.
    
    Raises:
        exception.IllegalArgumentError: If `button` is None.
    """
    # <editor-fold desc="Checks">
    if button is None:
        logger.error("button is None.")
        raise exception.IllegalArgumentError("button is None.")
    
    # </editor-fold>
    
    with open(
        os.path.join(global_variables.global_var_root_dir, "gui", "ui", "styles", "bottom_frame_button.css"),
        "r",
    ) as style_sheet_file:
        button_style = style_sheet_file.read()
        # Set the stylesheet of the application
        button.setStyleSheet(button_style)


def color_button_not_ready(button: QtWidgets.QPushButton) -> None:
    """Colors the button in the style of the NOT ready button.

    Args:
        button (QtWidgets.QPushButton): The button to apply the style to.
    
    Raises:
        exception.IllegalArgumentError: If `button` is None.
    """
    # <editor-fold desc="Checks">
    if button is None:
        logger.error("button is None.")
        raise exception.IllegalArgumentError("button is None.")

    # </editor-fold>
    
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
    logger.info("Using the 'in-project' stylesheet.")
    with open(
        pathlib.Path(f"{constants.PLUGIN_ROOT_PATH}/pyssa/gui/ui/styles/pyssa_style.css"),
        "r",
        encoding="utf-8",
    ) as file:
        style = file.read()
        # Set the stylesheet of the application
        self.setStyleSheet(style)


def set_stylesheet_homepage(self) -> None:  # noqa: ANN001
    """Sets the style sheet to the QMainWindow or a QDialog.

    Args:
        self: a QMainWindow or QDialog
    """
    logger.info("Using the homepage stylesheet.")
    with open(
        pathlib.Path(f"{constants.PLUGIN_ROOT_PATH}/pyssa/gui/ui/styles/pyssa_style_homepage.css"),
        "r",
        encoding="utf-8",
    ) as file:
        style = file.read()
        # Set the stylesheet of the application
        self.setStyleSheet(style)
