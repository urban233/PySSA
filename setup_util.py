import os
import pathlib
import shutil
from urllib import request


class Directory:
  """Container for common filesystem directory-based operations."""

  @staticmethod
  def exists(a_path: str | pathlib.Path) -> bool | None:
    """Checks if the given directory exists.

    Args:
      a_path: The path to the directory.

    Returns:
      A boolean indicating whether the directory exists or not or None if the argument is illegal.
    """
    # <editor-fold desc="Checks">
    if a_path is None or a_path == "":
      return False

    # </editor-fold>

    return pathlib.Path(a_path).exists()

  @staticmethod
  def create_directory(a_path: str | pathlib.Path) -> pathlib.Path | None:
    """Creates a directory and all its subdirectories.

    Args:
      a_path: The path to the directory to be deleted.

    Returns:
      The path that was created or None if the directory could not be deleted or the argument was illegal.

    Notes:
      This static method always deletes the directory if it already exists.
    """
    # <editor-fold desc="Checks">
    if a_path is None or a_path == "":
      return None

    # </editor-fold>

    tmp_path: pathlib.Path = pathlib.Path(a_path)
    if not Directory.purge(tmp_path):
      return None
    tmp_path.mkdir(parents=True)
    return tmp_path

  @staticmethod
  def copy_directory(source_directory_name: str | pathlib.Path,
                     destination_directory_name: str | pathlib.Path) -> bool:
    """Copies a directory and all its contents to a new location.

    Args:
      source_directory_name: The path to the directory to be copied.
      destination_directory_name: The path to the destination directory.

    Returns:
      A boolean indicating the success of the operation.
    """
    # <editor-fold desc="Checks">
    if source_directory_name is None or source_directory_name == "":
      return False
    if destination_directory_name is None or destination_directory_name == "":
      return False

    # </editor-fold>

    try:
      tmp_source_dir = pathlib.Path(source_directory_name)
      tmp_destination_dir = pathlib.Path(destination_directory_name)

      if not Directory.purge(tmp_destination_dir):
        return False
      if tmp_destination_dir.exists():
        tmp_destination_dir.mkdir(parents=True)

      shutil.copytree(tmp_source_dir, tmp_destination_dir)
      return True
    except Exception as e:
      return False

  @staticmethod
  def purge(a_path: str | pathlib.Path) -> bool:
    """Deletes a directory and all its contents.

    This function can also be used if the directory might not exist because
    it checks the existence before purging.

    Args:
      a_path: The path to the directory to be deleted.

    Returns:
      A boolean indicating the success of the operation.
    """
    # <editor-fold desc="Checks">
    if a_path is None or a_path == "":
      return False

    # </editor-fold>

    try:
      tmp_path = pathlib.Path(a_path)
      if tmp_path.exists():
        shutil.rmtree(tmp_path)
    except Exception as e:
      return False
    return True


class File:
  """Container for common filesystem file-based operations."""

  @staticmethod
  def copy(a_source_file_name: str | pathlib.Path, a_dest_file_name: str | pathlib.Path, overwrite: bool = False) -> bool:
    """Copies a file to a new location.

    Args:
      a_source_file_name: The path of the file to copy.
      a_dest_file_name: The destination path of the file to copy.
      overwrite: A boolean indicating whether the file should be overwritten.

    Returns:
      A boolean indicating the success of the operation.
    """
    # <editor-fold desc="Checks">
    if a_source_file_name is None or a_source_file_name == "":
      return False
    if a_dest_file_name is None or a_dest_file_name == "":
      return False
    if overwrite is None:
      return False

    # </editor-fold>

    try:
      if overwrite:
        if pathlib.Path(a_dest_file_name).exists():
          os.remove(a_dest_file_name)
      shutil.copy(a_source_file_name, pathlib.Path(a_dest_file_name))
    except Exception as e:
      return False
    return True

  @staticmethod
  def delete(a_filepath: str | pathlib.Path) -> bool:
    """Deletes a file.

    Args:
      a_filepath: The path of the file to delete.

    Returns:
      A boolean indicating the success of the operation.
    """
    # <editor-fold desc="Checks">
    if a_filepath is None or a_filepath == "":
      return False

    # </editor-fold>

    try:
      os.remove(a_filepath)
    except Exception as e:
      print(e)
      return False
    return True


def download_file(an_url: str, a_filepath: str) -> bool:
  """Downloads a single file of the given URL.

  Args:
    an_url: The URL to download.
    a_filepath: The path to the file to download.

  Returns:
    True if the download was successful, False otherwise.
  """
  try:
    request.urlretrieve(an_url, a_filepath)
    return True
  except Exception as e:
    print(e)
    return False
