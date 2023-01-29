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
"""This module contains all loggers used in the pyssa project."""
import logging
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants


# adding necessary handlers
constants.PYSSA_LOGGER.addHandler(log_handlers.log_file_handler)
constants.PREDICTION_WORKER_LOGGER.addHandler(log_handlers.log_file_handler)
constants.ANALYSIS_WORKER_LOGGER.addHandler(log_handlers.log_file_handler)


def log_multiple_messages(logger, logging_level, messages: list[str]):
    """This function logs multiple messages.

    Args:
        logger:
            the actual logger from constants.py
        logging_level:
            the level of the log message
        messages:
            the messages which should get logged
    """
    for message in messages:
        logger.log(level=logging_level, msg=message)


def log_multiple_variable_values(logger, step: str, variables_and_description: list[tuple[str, vars]]):
    """This function creates a list of messages which contain the name and value of multiple variables.

    Args:
        logger:
            the actual logger from constants.py
        step:
            usage location in the code e.g. function or class + function
        variables_and_description:
            a tuple with the name and the value of the variable
    """
    for variable in variables_and_description:
        logger.log(logging.DEBUG, msg=f"Filepath: {__file__}; Step: {step}; {variable[0]}: {variable[1]}")


def log_single_variable_value(logger, step: str, variable_name: str, variable: vars):
    """This function creates a message which contains the name and value of a variable with its usage location.

    Args:
        logger:
            the actual logger from constants.py
        step:
            usage location in the code e.g. function or class + function
        variable_name:
            name of the var
        variable:
            value of the var
    """
    logger.log(logging.DEBUG, msg=f"Filepath: {__file__}; Step: {step}; {variable_name}: {variable}")
