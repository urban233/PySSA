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
"""Module contains snippets that can be copy and pasted into differnt locations."""

# Logging
from pyssa.util import exception
logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"

# Basic argument check for None
if a_value_to_check is None:
  logger.error("a_value_to_check is None.")
  raise exception.IllegalArgumentError("a_value_to_check is None.")

if a_value_to_check is None or a_value_to_check == "":
  logger.error("a_value_to_check is either None or an empty string.")
  raise exception.IllegalArgumentError("a_value_to_check is either None or an empty string.")

# messages
"""
exception.IllegalArgumentError: If any of the arguments are None or if `a_filepath` is an empty string.
exception.IllegalArgumentError: If any of the arguments are None.
True: Operation successful, False: Otherwise
"""

"""Singal used to transfer data back to the previous window."""
"""Connects all UI elements to their corresponding slot functions in the class."""
"""Restores the UI."""

"""Constructor.

Args:
    the_interface_manager (interface_manager.InterfaceManager): The InterfaceManager object.

Raises:
    exception.IllegalArgumentError: If `the_interface_manager` is None.
"""
# <editor-fold desc="Checks">
if the_interface_manager is None:
  logger.error("the_interface_manager is None.")
  raise exception.IllegalArgumentError("the_interface_manager is None.")

# </editor-fold>


"""


Raises:
    exception.IllegalArgumentError: If `a_page_name` is None.
"""
# <editor-fold desc="Checks">
if a_page_name is None:
    logger.error("a_page_name is None.")
    raise exception.IllegalArgumentError("a_page_name is None.")

# </editor-fold>

"""Opens the help center and performs necessary actions based on the return value.

Args:
    return_value (tuple): The return value from opening the help center.
"""


# Scratch

"""Checks if a selection is empty."""

