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

  with open(
          os.path.join(
            global_variables.global_var_root_dir,
            "gui",
            "ui",
            "styles",
            "bottom_frame_button.css",
          ),
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
          os.path.join(
            global_variables.global_var_root_dir,
            "gui",
            "ui",
            "styles",
            "styles_start_button_not_ready.css",
          ),
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
  #logger.info("Using the 'in-project' stylesheet.")
  with open(
          pathlib.Path(
            f"{model_definitions.ModelDefinitions.PROGRAM_SRC_PATH}/python/pydd/components/internal/core/util/gui_style/style.css"
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
