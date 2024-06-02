#
# TEA - Task Event-based Async library for Python
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/urban233/TEA>
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
import datetime
import os
import pathlib
from typing import Optional


def setup_logger(
    a_name: str, a_path: pathlib.Path, level=logging.DEBUG
) -> Optional[logging.Logger]:
  """Sets up a logger with a FileHandler directing output to the specified file.

  Args:
    a_name (str): Name for the logger.
    a_path (pathlib.Path): Path where the logging file is stored.
    level (int, optional): Logging level (default: logging.INFO).
  """
  # <editor-fold desc="Checks">
  if a_name is None or a_name == "":
    return None
  if a_path is None or a_path == "":
    return None
  if not pathlib.Path.exists(a_path):
    os.mkdir(a_path)

  # </editor-fold>

  tmp_logger = logging.getLogger(a_name)
  tmp_logger.setLevel(level)

  # File handler for separate log files
  tmp_current_time = datetime.datetime.now()
  tmp_filename = f"{tmp_current_time.year}-{tmp_current_time.month:02d}-{tmp_current_time.day:02d}_{tmp_current_time.hour:02d}-{tmp_current_time.minute:02d}.log"  # noqa: E501
  tmp_filepath = str(f"{a_path}\\{tmp_filename}")
  tmp_handler = logging.FileHandler(tmp_filepath)
  tmp_log_formatter = logging.Formatter(
      "%(asctime)s: %(name)s %(levelname)s - %(message)s"
  )
  tmp_handler.setFormatter(tmp_log_formatter)
  tmp_logger.addHandler(tmp_handler)

  # Console handler for logging to console
  tmp_console_handler = logging.StreamHandler()
  tmp_console_formatter = logging.Formatter(
      "%(levelname)s: %(message)s"
  )  # Simpler console format
  tmp_console_handler.setFormatter(tmp_console_formatter)
  tmp_logger.addHandler(tmp_console_handler)

  return tmp_logger
