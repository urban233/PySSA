#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
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
"""Module for custom pyssa exceptions."""


class IllegalArgumentError(Exception):
  """Class for an illegal argument exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class FileIsEmptyError(Exception):
  """Class for a file is empty exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class NoInternetConnectionError(Exception):
  """Class for a no internet exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class ActionIsNotRunnableError(Exception):
  """Class for an action is not runnable type."""

  def __init__(self, message: str = "") -> None:
    """Constructor."""
    super().__init__("self.is_runnable is set to False. " + message)


class NotMainThreadError(Exception):
  """Class for a not main thread exception type."""

  def __init__(self, message: str = "") -> None:
    """Constructor."""
    super().__init__(
        "Call was from a separate thread but this function must only be called from the main thread! "
        + message
    )


class PyMOLCommandFailedError(Exception):
  """Class for a pymol command failed exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)
    

class DirectoryNotFoundError(Exception):
  """Class for a directory does not exist exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class FastaFilesNotFoundError(Exception):
  """Class for an empty WSL fasta directory exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class FastaFilesNotCreatedError(Exception):
  """Class for an empty WSL fasta directory exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class PredictionEndedWithError(Exception):
  """Class for an empty WSL fasta directory exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToOpenFileError(Exception):
  """Class for an unable to open file exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToCreatePdbFileError(Exception):
  """Class for an unable to create pdb file exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToLoadPdbFileError(Exception):
  """Class for an unable to load pdb file exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToLoadProteinError(Exception):
  """Class for an unable to load a protein exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToLoadProteinPairError(Exception):
  """Class for an unable to load a protein pair exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToAddProteinPairError(Exception):
  """Class for an unable to add a protein pair exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToColorProteinPairError(Exception):
  """Class for an unable to color a protein pair exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToReinitializePymolSessionError(Exception):
  """Class for an unable to reinitialize the pymol session exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToDoAnalysisError(Exception):
  """Class for an unable to do analysis exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToTakeImageError(Exception):
  """Class for an unable to do analysis exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToSetImageError(Exception):
  """Class for an unable to do analysis exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToSafeSessionError(Exception):
  """Class for an unable to load a protein pair exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class SubprocessExecutionError(Exception):
  """Class for a subprocess execution exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class Wsl2PreparationFailedError(Exception):
  """Class for a failed wsl2 preparation exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToCreateAnalysisRunInfoError(Exception):
  """Class for an unable to create analysis run info exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToCreateProteinPairError(Exception):
  """Class for an unable to create a protein pair object exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToCreateDistanceAnalysisObjectError(Exception):
  """Class for an unable to create a distance analysis object exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToTransformDataForAnalysisError(Exception):
  """Class for an unable to transform data for analysis exception type."""

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


class UnableToSetupAnalysisError(Exception):
  """Class for an unable to set up analysis exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class UnableToFindColabfoldModelError(Exception):
  """Class for an unable to find colabfold model exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class DistanceAnalysisFailedError(Exception):
  """Class for a distance analysis failed exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class ProteinNotFoundInCurrentProjectError(Exception):
  """Class for protein not found in current project exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)


class ProteinNotFoundInPyMOLSession(Exception):
  """Class for protein not found in current pymol session exception type."""

  def __init__(self, message: str) -> None:
    """Constructor."""
    super().__init__(message)
