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
"""Module which contains the void function."""
from pyssa.util import constants, exception


def rvoid(a_python_object) -> None:  # noqa: ANN001
  """Only logs that a return value of a function or method is not used.

  Args:
      a_python_object (object): the object to be void.

  Raises:
      exception.IllegalArgumentError: If a_python_object is None.
  """
  # <editor-fold desc="Checks">
  if a_python_object is None:
    constants.PYSSA_LOGGER.error("a_python_object is None.")
    raise exception.IllegalArgumentError("a_python_object is None.")

  # </editor-fold>

  constants.PYSSA_LOGGER.debug(
      "VOID: Return value: %s is not used.", a_python_object
  )
