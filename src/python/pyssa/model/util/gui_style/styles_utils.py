"""Module for all styles related functions."""
import logging
import os
import pathlib
from PyQt6 import QtWidgets
from pyssa.model.pyssa_logging import default_logging
from pyssa.model.util import exception
from pyssa.model.preference import model_definitions

__docformat__ = "google"

logger = default_logging.setup_logger(__file__)


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
  # Set the stylesheet of the application
  button.setStyleSheet(
    """
    QPushButton {
        background-color: #367AF6;
        color: #fff;
        font-family: "Segoe UI Semibold";
        font-size: 12px;
        border: solid;
        border-width: 1px;
        border-radius: 4px;
        border-color: #DCDCDC;
        padding: 2px;
        min-width: 65px;
        max-width: 65px;
        min-height: 15px;
        max-height: 15px;
        /*font-family: "Segoe UI";*/
        /*font-size: 12px;*/
        /*font: bold;*/
        /*border: none;*/
        /*border-width: 2px;*/
        /*border-radius: 4px;*/
        /*border-color: #DCDCDC;*/
        /*padding: 2px;*/
        /*min-width: 65px;*/
        /*max-width: 65px;*/
        /*min-height: 15px;*/
    }
    
    QPushButton:disabled {
        background-color: #fff;
        color: #B0B0B0;
        font-family: "Segoe UI";
        font-size: 12px;
        border: solid;
        border-width: 1px;
        border-radius: 4px;
        border-color: #DCDCDC;
        padding: 2px;
        min-width: 65px;
        max-width: 65px;
        min-height: 15px;
        /*background-color: #fff;*/
        /*color: #B0B0B0;*/
        /*font-family: "Segoe UI";*/
        /*font-size: 12px;*/
        /*font: bold;*/
        /*border: solid;*/
        /*border-width: 2px;*/
        /*border-radius: 4px;*/
        /*border-color: #fff;*/
        /*padding: 2px;*/
        /*min-width: 65px;*/
        /*max-width: 65px;*/
        /*min-height: 15px;*/
    }
    
    QPushButton::pressed {
        background-color: #204993;
        color: #fff;
        font-family: "Segoe UI";
        font-size: 12px;
        border: none;
        border-width: 2px;
        border-radius: 4px;
        border-color: #DCDCDC;
        padding: 2px;
        min-width: 65px;
        max-width: 65px;
        min-height: 15px;
    }
    """
  )


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
  button.setStyleSheet(
    """
    QPushButton {
        background-color: white;
        border-style: outset;
        border-width: 2px;
        border-radius: 6px;
        border-color: #f7f7f7;
        min-width: 6em;
        padding: 2px;
        min-width: 75px;
        max-width: 75px;
        min-height: 15px;
    }
    """
  )


def set_stylesheet(self) -> None:  # noqa: ANN001
  """Sets the style sheet to the QMainWindow or a QDialog.

  Args:
      self: a QMainWindow or QDialog
  """
  logger.info("Using the 'in-project' stylesheet.")
  with open(
          pathlib.Path(
            f"{model_definitions.ModelDefinitions.PROGRAM_SRC_PATH}/python/pyssa/model/util/gui_style/style.css"
          ),
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
          pathlib.Path(
            f"{model_definitions.ModelDefinitions.PROGRAM_SRC_PATH}/pyssa/gui/ui/styles/pyssa_style_homepage.css"
          ),
          "r",
          encoding="utf-8",
  ) as file:
    style = file.read()
    # Set the stylesheet of the application
    self.setStyleSheet(style)
