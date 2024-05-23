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
"""Module for all asynchronous functions that are related to the image creation and preview for PyMOL."""
import logging

from pyssa.controller import pymol_session_manager
from pyssa.logging_pyssa import log_handlers
from pyssa_pymol import pymol_enums

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


def preview_image(
    the_pymol_session_manager: "pymol_session_manager.PymolSessionManager",
    a_placeholder_2: int,
) -> tuple[int, str]:
  """Previews an image that would be ray-traced.

  Args:
      the_pymol_session_manager (pymol_session_manager.PymolSessionManager): The PymolSessionManager object responsible for managing the Pymol session.
      a_placeholder_2 (int): An integer placeholder that is necessary for using `LegacyTask`.

  Returns:
      A tuple containing two values:
          - The first value is 0 if operation was successful and 1 if operation failed.
          - The second value is an empty string.

  Note:
      This method is used to preview an image of the Pymol session using the specified renderer. The default renderer is 0.
      The image dimensions are set to a width and height of 2400 pixels.
  """
  # <editor-fold desc="Checks">
  if the_pymol_session_manager is None:
    logger.error("the_pymol_session_manager is None.")
    return 1, ""

  # </editor-fold>

  try:
    # TODO: the renderer should be changeable
    the_pymol_session_manager.cmd(
      pymol_enums.CommandEnum.RAY, (2400, 2400, int(0))
    )
  except Exception as e:
    logger.error(e)
  return 0, ""


def create_drawn_image(
    an_image_filepath: str,
    the_pymol_session_manager: "pymol_session_manager.PymolSessionManager",
) -> tuple[int, str]:
  """Creates a drawn image and save it as a PNG file.

  Args:
      an_image_filepath (str): The file path where the image will be saved.
      the_pymol_session_manager (pymol_session_manager.PymolSessionManager): An instance of the PymolSessionManager class responsible for managing the Pymol session.

  A tuple containing two values:
          - The first value is 0 if operation was successful and 1 if operation failed.
          - The second value is an empty string.
  """
  # <editor-fold desc="Checks">
  if an_image_filepath is None or an_image_filepath == "":
    logger.error("an_image_filepath is either None or an empty string.")
    return 1, ""
  if the_pymol_session_manager is None:
    logger.error("the_pymol_session_manager is None.")
    return 1, ""

  # </editor-fold>

  try:
    the_pymol_session_manager.cmd(
      pymol_enums.CommandEnum.DRAW, (2400, 2400, 2)
    )
    the_pymol_session_manager.cmd(
      pymol_enums.CommandEnum.PNG, (an_image_filepath, 300)
    )
  except Exception as e:
    logger.error(e)
  return 0, ""
