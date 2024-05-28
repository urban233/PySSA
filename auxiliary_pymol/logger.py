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
"""Module contains the handlers used by the loggers."""
import logging
import os
from typing import Optional

__docformat__ = "google"

from auxiliary_pymol import local_constants


def setup_logger(
    a_name: str, level: int = logging.DEBUG, log_path: str = local_constants.LOG_PATH, add_console_handler: bool = True
) -> Optional[logging.Logger]:
  """Sets up a logger with a FileHandler directing output to the specified file.

  Args:
      a_name (str): Name for the logger.
      level (int, optional): Logging level (Default: logging.INFO).
      log_path (str): The path for the log file. (Default: The default PySSA log path.)
      add_console_handler (bool): A flag for whether to add a console handler. (Default: True).
  
  Returns:
      A logger object or None if any of the arguments are None or `a_name` is an empty string.
  """
  # <editor-fold desc="Checks">
  if a_name is None or a_name == "":
    return None
  if log_path is None or log_path == "":
    return None
  if add_console_handler is None:
    return None
  
  # </editor-fold>

  # <editor-fold desc="Definitions of important paths">
  if log_path == local_constants.LOG_PATH:
    if not os.path.exists(local_constants.SETTINGS_DIR):
      os.mkdir(local_constants.SETTINGS_DIR)
    if not os.path.exists(local_constants.LOG_PATH):
      os.mkdir(local_constants.LOG_PATH)
  else:
    if not os.path.exists(log_path):
      os.mkdir(log_path)

  # </editor-fold>

  tmp_logger = logging.getLogger(a_name)
  tmp_logger.setLevel(level)

  # File handler for separate log files
  tmp_handler = logging.FileHandler(f"{log_path}\\{local_constants.LOG_FILENAME}")
  tmp_log_formatter = logging.Formatter(
      "%(asctime)s: %(name)s %(levelname)s - %(message)s"
  )
  tmp_handler.setFormatter(tmp_log_formatter)
  tmp_logger.addHandler(tmp_handler)
  
  if add_console_handler:
    # Console handler for logging to console
    tmp_console_handler = logging.StreamHandler()
    tmp_console_formatter = logging.Formatter(
        "%(levelname)s: %(message)s"
    )  # Simpler console format
    tmp_console_handler.setFormatter(tmp_console_formatter)
    tmp_logger.addHandler(tmp_console_handler)

  return tmp_logger
