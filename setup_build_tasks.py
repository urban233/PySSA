import abc
import pathlib
import subprocess
import setup_util


class IBuildTask(abc.ABC):
  """Interface for defining build tasks used as commands in the setup.py build script."""

  @abc.abstractmethod
  def execute_task(self):
    """Defines the build task to be executed."""
    raise NotImplementedError()


class CreateWinPackageBuildTask(IBuildTask):
  """Creates a Windows package that can be used for deployment."""

  # <editor-fold desc="Class attributes">
  project_root_path = pathlib.Path(__file__).parent.absolute()
  """The root path of the project."""

  temp_path: pathlib.Path = pathlib.Path(project_root_path, "_build", "temp")
  """The temporary path for the build package."""

  base_win_package_path: pathlib.Path = pathlib.Path(project_root_path, "base_win_package")
  """The path to the base Windows package."""

  temp_base_win_package_path: pathlib.Path = pathlib.Path(temp_path, "base_win_package")
  """The temporary path to the base Windows package."""

  temp_pyssa_path: pathlib.Path = pathlib.Path(temp_base_win_package_path, "bin", "PySSA")
  """The temporary path for the PySSA package."""

  assets_path: pathlib.Path = pathlib.Path(project_root_path, "assets")
  """The path to the assets directory."""

  temp_assets_path: pathlib.Path = pathlib.Path(temp_pyssa_path, "assets")
  """The temporary path to the assets directory."""

  docs_path: pathlib.Path = pathlib.Path(project_root_path, "docs")
  """The path to the docs directory."""

  temp_docs_path: pathlib.Path = pathlib.Path(temp_pyssa_path, "docs")
  """The temporary path to the docs directory."""

  scripts_path: pathlib.Path = pathlib.Path(project_root_path, "scripts")
  """The path to the scripts directory."""

  temp_scripts_path: pathlib.Path = pathlib.Path(project_root_path, "scripts")
  """The temporary path to the scripts directory."""

  src_path: pathlib.Path = pathlib.Path(project_root_path, "src")
  """The path to the src directory."""

  temp_src_path: pathlib.Path = pathlib.Path(temp_pyssa_path, "src")
  """The temporary path to the src directory."""

  winbatch_path: pathlib.Path = pathlib.Path(project_root_path, "winbatch")
  """The path to the winbatch directory."""

  temp_winbatch_path: pathlib.Path = pathlib.Path(temp_pyssa_path, "winbatch")
  """The temporary path to the winbatch directory."""

  license_filepath: pathlib.Path = pathlib.Path(project_root_path, "LICENSE")
  """The filepath to the license file."""

  temp_license_filepath: pathlib.Path = pathlib.Path(temp_pyssa_path, "LICENSE")
  """The temporary filepath to the license file."""

  readme_filepath: pathlib.Path = pathlib.Path(project_root_path, "README.md")
  """The filepath to the README file."""

  temp_readme_filepath: pathlib.Path = pathlib.Path(temp_pyssa_path, "README.md")
  """The temporary filepath to the README file."""

  build_output_path: pathlib.Path = pathlib.Path(project_root_path, "_build", "output")
  """The path to the build output directory."""

  win_package_zip_filepath: pathlib.Path = pathlib.Path(build_output_path, "win_package.zip")
  """The path to the win_package.zip file."""

  # </editor-fold>

  def execute_task(self):
    try:
      # Clean temporary build directory
      print("Clean temporary build directory ...")
      if setup_util.Directory.exists(self.temp_path):
        setup_util.Directory.purge(self.temp_path)
      setup_util.Directory.create_directory(self.temp_path)
      # Copy base windows package from resources to tmp directory
      print("Copy base windows package from resources to temp directory ...")
      setup_util.Directory.copy_directory(self.base_win_package_path, self.temp_base_win_package_path)
      # Copy new PySSA src files from GitHub repo to tmp/offline_win_package directory
      print("Copy relevant PySSA files and directories into bin directory of base windows package ...")
      setup_util.Directory.copy_directory(self.assets_path, self.temp_assets_path)
      setup_util.Directory.copy_directory(self.docs_path, self.temp_docs_path)
      setup_util.Directory.copy_directory(self.scripts_path, self.temp_scripts_path)
      setup_util.Directory.copy_directory(self.src_path, self.temp_src_path)
      setup_util.Directory.copy_directory(self.winbatch_path, self.temp_winbatch_path)
      setup_util.File.copy(self.license_filepath, self.temp_license_filepath, overwrite=True)
      setup_util.File.copy(self.readme_filepath, self.temp_readme_filepath, overwrite=True)
      # Construct the PowerShell command
      setup_util.Directory.create_directory(self.build_output_path)
      powershell_command = f"Compress-Archive -Path {self.temp_base_win_package_path}\\* -DestinationPath {self.win_package_zip_filepath} -Force"
      # Execute the PowerShell command using subprocess.run()
      print("Compressing the base Windows package ...")
      subprocess.run(["powershell", "-Command", powershell_command], shell=True)
      print("Build process finished without errors.")
    except Exception as e:
      print(e)
    finally:
      # Clean directory
      print("Clean directory")
      setup_util.Directory.purge(self.temp_path)
