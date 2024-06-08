import pathlib
import shutil


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
      else:
        print("Nothing to do. The given path does not exist.")
    except Exception as e:
      return False
    return True
