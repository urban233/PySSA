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
"""Module contains the CurrentSession dataclass."""
__docformat__ = "google"

from pyssa.util import exception


# @dataclass
# class CurrentSession:
#     """Holds information about the current pymol session.
#
#     DO NOT USE THIS CLASS ANYMORE IN PRODUCTION!
#     """
#
#     type: str  # noqa: A003
#     name: str
#     session: str


class CurrentPymolSession:
  """Holds information about the current pymol session."""

  # <editor-fold desc="Class attributes">
  session_name: str
  """The name of the current pymol session."""

  object_type: str
  """The object type of the current pymol session."""

  # </editor-fold>

  def __init__(self, a_session_name: str, an_object_type) -> None:
    """Constructor.

    Args:
        a_session_name (str): The name of the session.
        an_object_type (any): The type of the object.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_session_name is None:
      raise exception.IllegalArgumentError("a_session_name is None.")
    if an_object_type is None:
      raise exception.IllegalArgumentError("an_object_type is None.")

    # </editor-fold>

    self.session_name = a_session_name
    self.object_type = an_object_type
