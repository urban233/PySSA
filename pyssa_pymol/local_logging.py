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
import datetime
import os.path
import pathlib

# current_time = datetime.datetime.now()
# filename = f"{current_time.year}-{current_time.month:02d}-{current_time.day:02d}_{current_time.hour:02d}-{current_time.minute:02d}.log"  # noqa: E501
# SETTINGS_DIR = str(pathlib.Path(f"{os.path.expanduser('~')}/.pyssa/"))
#
# log_formatter = logging.Formatter("%(asctime)s: %(name)s %(levelname)s - %(message)s")
# if not os.path.exists(f"{SETTINGS_DIR}\\user_pymol"):
#     os.mkdir(f"{SETTINGS_DIR}\\user_pymol")
# if not os.path.exists(f"{SETTINGS_DIR}\\user_pymol\\logs"):
#     os.mkdir(f"{SETTINGS_DIR}\\user_pymol\\logs")
# filepath = f"{SETTINGS_DIR}\\user_pymol\\logs\\{filename}"
# log_file_handler = logging.FileHandler(filepath)
# log_file_handler.setLevel(logging.DEBUG)
# log_file_handler.setFormatter(log_formatter)


def setup_logger(name, level=logging.DEBUG):
    """
    Sets up a logger with a FileHandler directing output to the specified file.

    Args:
        name (str): Name for the logger.
        path (str): Path to the log file.
        level (int, optional): Logging level (default: logging.INFO).
    """
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # File handler for separate log files
    current_time = datetime.datetime.now()
    filename = f"{current_time.year}-{current_time.month:02d}-{current_time.day:02d}_{current_time.hour:02d}-{current_time.minute:02d}.log"  # noqa: E501
    SETTINGS_DIR = str(pathlib.Path(f"{os.path.expanduser('~')}/.pyssa/"))
    path = f"{SETTINGS_DIR}\\user_pymol\\logs"
    filepath = f"{path}\\{filename}"
    handler = logging.FileHandler(filepath)
    log_formatter = logging.Formatter("%(asctime)s: %(name)s %(levelname)s - %(message)s")
    handler.setFormatter(log_formatter)
    logger.addHandler(handler)

    # Console handler for logging to console
    console_handler = logging.StreamHandler()
    console_formatter = logging.Formatter('%(levelname)s: %(message)s')  # Simpler console format
    console_handler.setFormatter(console_formatter)
    logger.addHandler(console_handler)

    return logger
