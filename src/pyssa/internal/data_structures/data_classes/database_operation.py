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
"""Module contains the DatabaseOperation dataclass."""
from dataclasses import dataclass
from src.pyssa.util import enums


@dataclass
class DatabaseOperation:
  """DatabaseOperation is a dataclass that represents an SQL database operation."""

  # <editor-fold desc="Class attributes">
  sql_query_type: enums.SQLQueryType
  """The type of SQL query."""

  buffered_data: tuple
  """The buffered data to be stored in the database."""
  # </editor-fold>

  def __init__(
      self, a_sql_query_type: "enums.SQLQueryType", the_buffered_data: tuple
  ) -> None:
    """Constructor.

    Args:
        a_sql_query_type (enums.SQLQueryType): The type of SQL query to perform.
        the_buffered_data (tuple): The data that is necessary for the query to be performed.

    Note:
        The buffered_data is a tuple which always starts with a 0 -> (0, your_data)
    """
    self.sql_query_type: enums.SQLQueryType = a_sql_query_type
    self.buffered_data: tuple = the_buffered_data
