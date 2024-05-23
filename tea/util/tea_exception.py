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
"""Module for custom tea exceptions."""


class ActionIsNotRunnableError(Exception):
  """Class for an action is not runnable type."""

  def __init__(self, message: str = "") -> None:
    """Constructor."""
    super().__init__("self.is_runnable is set to False. " + message)


class IllegalArgumentError(Exception):
  """Class for an illegal argument exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


# <editor-fold desc="Filesystem related">
class FileIsEmptyError(Exception):
  """Class for a file is empty exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class DirectoryNotFoundError(Exception):
  """Class for a directory does not exist exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToOpenFileError(Exception):
  """Class for an unable to open file exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToDeleteDirectoryError(Exception):
  """Class for an unable to delete a directory exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToCopyDirectoryError(Exception):
  """Class for an unable to copy a directory exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToCreateDirectoryError(Exception):
  """Class for an unable to create a directory exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToCopyFileError(Exception):
  """Class for an unable to copy a file exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


# </editor-fold>


class SubprocessExecutionError(Exception):
  """Class for a subprocess execution exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)
