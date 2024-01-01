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
"""Module for functions which create standardized log messages."""
from pyssa.logging_pyssa import loggers


def generate_multiple_var_value_msgs(step: str, variables_and_description: list[tuple[str, vars]]) -> list:
    """This function creates a list of messages which contain the name and value of multiple variables.

    Args:
        step:
            usage location in the code e.g. function or class + function
        variables_and_description:
            a tuple with the name and the value of the variable

    Returns:
        a list of messages which should be logged
    """
    messages = []
    for variable in variables_and_description:
        messages.append(f"Filepath: {__file__}; Step: {step}; {variable[0]}: {variable[1]}")
    return messages


def log_single_variable_value(logger, logging_level, step: str, variable_name: str, variable: vars) -> str:
    """This function creates a message which contains the name and value of a variable with its usage location.

    Args:
        logger:
            one of the application loggers:
                pyssa, analysis_worker, prediction_worker
        logging_level:
            the level of the log message
        step:
            usage location in the code e.g. function or class + function
        variable_name:
            name of the var
        variable:
            value of the var

    Returns:
        a message string which contains the name and value of a variable
    """
    message = f"Filepath: {__file__}; Step: {step}; {variable_name}: {variable}"
    loggers.prediction_worker.log(logging_level, msg=message)
