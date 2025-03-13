"""
#A* -------------------------------------------------------------------
#B* This file contains source code for running automation tasks related
#-* to the build process of the PySSA computer program
#C* Copyright 2025 by Martin Urban.
#D* -------------------------------------------------------------------
#E* It is unlawful to modify or remove this copyright notice.
#F* -------------------------------------------------------------------
#G* Please see the accompanying LICENSE file for further information.
#H* -------------------------------------------------------------------
#I* Additional authors of this source file include:
#-*
#-*
#-*
#Z* -------------------------------------------------------------------
"""
import argparse
import os
import pathlib
import subprocess
import shutil
import sys
import time
import zipfile
import setup_util


# pyproject_toml = toml.load("pyproject.toml")
# PROJECT_NAME = pyproject_toml["project"]["name"]
# PROJECT_VERSION = pyproject_toml["project"]["version"]

PROJECT_ROOT_DIR = pathlib.Path(__file__).parent

PYTHON_EXECUTABLE = sys.executable  # This gives the current Python executable
DEBUG = True


# <editor-fold desc="Automation classes">
class BuildWinPackage:
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

  temp_scripts_path: pathlib.Path = pathlib.Path(temp_pyssa_path, "scripts")
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

  version_history_filepath: pathlib.Path = pathlib.Path(project_root_path, "version_history.json")
  """The filepath to the version history file."""

  temp_version_history_filepath: pathlib.Path = pathlib.Path(temp_pyssa_path, "version_history.json")
  """The temporary filepath to the version history file."""

  build_output_path: pathlib.Path = pathlib.Path(project_root_path, "_build", "output")
  """The path to the build output directory."""

  win_package_zip_filepath: pathlib.Path = pathlib.Path(build_output_path, "win_package.zip")
  """The path to the win_package.zip file."""

  # </editor-fold>

  def build(self) -> None:
    """Builds the acutal windows package."""
    try:
      # Clean temporary build directory
      print("Clean temporary build directory ...")
      if setup_util.Directory.exists(self.temp_path):
        setup_util.Directory.purge(self.temp_path)
      setup_util.Directory.create_directory(self.temp_path)
      # Download bin.zip if not available offline
      if not pathlib.Path.exists(pathlib.Path(f"{self.base_win_package_path}\\bin.zip")):
        print("Downloading bin.zip ...")
        if not setup_util.download_file("https://w-hs.sciebo.de/s/nX5TklMAdRNnyoC/download", f"{self.base_win_package_path}\\bin.zip"):
          print("Unable to download bin.zip, build process exists.")
          return
        print("Finished downloading bin.zip.")
      # Copy base windows package from resources to tmp directory
      print("Copy base windows package from resources to temp directory ...")
      setup_util.Directory.copy_directory(self.base_win_package_path, self.temp_base_win_package_path)
      # Expand-Archive bin.zip
      powershell_command = f"Expand-Archive -Path {self.temp_base_win_package_path}\\bin.zip -DestinationPath {self.temp_base_win_package_path}\\bin"
      print("Expand-Archive bin.zip ...")
      subprocess.run(["powershell", "-Command", powershell_command], shell=True)
      setup_util.File.delete(f"{self.temp_base_win_package_path}\\bin.zip")
      # Copy new PySSA src files from GitHub repo to tmp/offline_win_package directory
      print("Copy relevant PySSA files and directories into bin directory of base windows package ...")
      setup_util.Directory.copy_directory(self.assets_path, self.temp_assets_path)
      setup_util.Directory.copy_directory(self.docs_path, self.temp_docs_path)
      setup_util.Directory.copy_directory(self.scripts_path, self.temp_scripts_path)
      setup_util.Directory.copy_directory(self.src_path, self.temp_src_path)
      setup_util.Directory.copy_directory(self.winbatch_path, self.temp_winbatch_path)
      setup_util.File.copy(self.license_filepath, self.temp_license_filepath, overwrite=True)
      setup_util.File.copy(self.readme_filepath, self.temp_readme_filepath, overwrite=True)
      # TODO: Add automated way of updating the version number project wide
      # Note: Until now there are different locations where the versions are stored
      # (1) constants.py (2) pyproject.toml (3) version_history.json => TODO: This has to change for better maintenance!
      setup_util.File.copy(self.version_history_filepath, self.temp_version_history_filepath, overwrite=True)
      setup_util.Directory.create_directory(self.build_output_path)
      setup_util.Directory.copy_directory(self.temp_base_win_package_path, self.build_output_path)
      # Construct the PowerShell command
      # powershell_command = f"Compress-Archive -Path {self.temp_base_win_package_path}\\* -DestinationPath {self.win_package_zip_filepath} -Force"
      # # Execute the PowerShell command using subprocess.run()
      # print("Compressing the base Windows package ...")
      # subprocess.run(["powershell", "-Command", powershell_command], shell=True)
      print("Build process finished without errors.")
    except Exception as e:
      print(e)
    finally:
      # Clean directory
      print("Clean directory")
      setup_util.Directory.purge(self.temp_path)


# class BuildInnoSetup:
#   """Contains the logic for building the inno setup EXE file."""
#
#   def __init__(self) -> None:
#     """Constructor."""
#     self.inno_build_path = pathlib.Path(PROJECT_ROOT_DIR / "innoBuild")
#     self.inno_sources_build_path = pathlib.Path(PROJECT_ROOT_DIR / "innoBuild/inno-sources")
#     self.inno_build_bin_path = pathlib.Path(self.inno_build_path / "inno-sources/bin")
#     self.inno_build_prerequisite_path = pathlib.Path(self.inno_build_path / "inno-sources/prerequisite")
#     self.inno_build_third_party_path = pathlib.Path(self.inno_build_path / "inno-sources/third_party")
#     self.inno_build_assets_path = pathlib.Path(self.inno_build_path / "inno-assets")
#     self.inno_setup_compiler_filepath = pathlib.Path(r"C:\Program Files (x86)\Inno Setup 6\ISCC.exe")
#     self.inno_setup_script_filepath = pathlib.Path(PROJECT_ROOT_DIR / "deployment/src/inno_setup" / "setup.iss")
#     self.build_src_path = pathlib.Path(PROJECT_ROOT_DIR / "pyssa")
#     self.build_script_filepath = pathlib.Path(PROJECT_ROOT_DIR / "pyssa", "build_win_exe.py")
#     self.build_aux_pymol_script_filepath = pathlib.Path(PROJECT_ROOT_DIR / "pyssa", "build_aux_pymol_exe.py")
#
#   def setup_build_environment(self) -> None:
#     """Sets up a temporary build environment."""
#     # <editor-fold desc="Path/Filepath definitions">
#     tmp_build_script_filepath = pathlib.Path(PROJECT_ROOT_DIR / "scripts/python", "build_win_exe.py")
#     tmp_build_aux_pymol_script_filepath = pathlib.Path(PROJECT_ROOT_DIR / "scripts/python", "build_aux_pymol_exe.py")
#     tmp_pyssa_win_build_logo_filepath = pathlib.Path(PROJECT_ROOT_DIR / "assets/convert_logo_to_ico" / "logo.ico")
#     tmp_pyssa_python_src_path = pathlib.Path(PROJECT_ROOT_DIR / "src")
#     tmp_pyssa_win_build_files_path = pathlib.Path(PROJECT_ROOT_DIR / "pyssa/build/exe.win-amd64-3.11")
#     tmp_pyssa_win_package_files_path = pathlib.Path(PROJECT_ROOT_DIR / "_build/output")
#     tmp_pyssa_win_lib_build_files_path = pathlib.Path(tmp_pyssa_win_build_files_path / "lib")
#     tmp_vc_redist_setup_filepath = pathlib.Path(PROJECT_ROOT_DIR / "third_party/microsoft" / "VC_redist.x64.exe")
#     tmp_windows_tasks_exe_filepath = pathlib.Path(PROJECT_ROOT_DIR / "deployment/offline_resources" / "WindowsTasks.exe")
#     tmp_site_packages_path = pathlib.Path(PROJECT_ROOT_DIR / ".venv/Lib/site-packages")
#     # </editor-fold>
#     """IMPORTANT
#     Use the python interpreter of the venv of pymol windows build because
#     that interpreter gets also used in the build script of the
#     pymol windows build repo!
#     """
#     # if self.inno_build_path.exists():
#     #   print("Only move PySSA sources! To rebuild the entire setup, run 'clean' beforehand!")
#     #   shutil.copytree(
#     #     tmp_pyssa_python_src_path,
#     #     pathlib.Path(self.inno_build_bin_path / "lib"), dirs_exist_ok=True
#     #   )
#     # else:
#     #   shutil.copytree(tmp_pyssa_python_src_path, self.build_src_path)
#     #   shutil.copy(tmp_build_script_filepath, self.build_script_filepath)
#     #   shutil.copy(tmp_build_aux_pymol_script_filepath, self.build_aux_pymol_script_filepath)
#     #   subprocess.run(
#     #     [PYTHON_EXECUTABLE, self.build_script_filepath],
#     #     stdout=sys.stdout, stderr=sys.stderr, text=True, cwd=self.build_src_path
#     #   )
#     # <editor-fold desc="Copy operations">
#     if self.inno_build_path.exists():
#       shutil.rmtree(self.inno_build_path)
#       self.inno_build_path.mkdir()
#     # shutil.copytree(tmp_pyssa_python_src_path, tmp_pyssa_win_lib_build_files_path, dirs_exist_ok=True)
#     # shutil.copytree(pathlib.Path(tmp_site_packages_path / "zmq"), pathlib.Path(tmp_pyssa_win_build_files_path / "lib/zmq"), dirs_exist_ok=True)
#     # shutil.copytree(pathlib.Path(tmp_site_packages_path / "pyzmq.libs"), pathlib.Path(tmp_pyssa_win_build_files_path / "lib/pyzmq.libs"), dirs_exist_ok=True)
#     #
#     # shutil.copytree(
#     #   pathlib.Path(PROJECT_ROOT_DIR / "assets"),
#     #   pathlib.Path(self.inno_build_path / "assets"),
#     #   dirs_exist_ok=True
#     # )
#     # shutil.copytree(
#     #   pathlib.Path(PROJECT_ROOT_DIR / "docs"),
#     #   pathlib.Path(self.inno_build_path / "docs"),
#     #   dirs_exist_ok=True
#     # )
#     # shutil.copytree(
#     #   pathlib.Path(PROJECT_ROOT_DIR / "winbatch"),
#     #   pathlib.Path(self.inno_build_path / "winbatch"),
#     #   dirs_exist_ok=True
#     # )
#     # shutil.copytree(tmp_pyssa_win_build_files_path, self.inno_build_bin_path, dirs_exist_ok=True)
#     # tmp_aux_pymol_exe_filepath = pathlib.Path(PROJECT_ROOT_DIR / "dist" / "AuxPyMOL.exe")
#     # if not tmp_aux_pymol_exe_filepath.exists():
#     #   subprocess.run(
#     #     [
#     #       pathlib.Path(PROJECT_ROOT_DIR / ".venv/Scripts" / "pyinstaller.exe"),
#     #       pathlib.Path(PROJECT_ROOT_DIR / "aux_pymol.spec"),
#     #     ]
#     #   )
#     # shutil.copy(
#     #   tmp_aux_pymol_exe_filepath,
#     #   pathlib.Path(tmp_pyssa_win_build_files_path / "AuxPyMOL.exe")
#     # )
#
#     shutil.copytree(tmp_pyssa_win_package_files_path, self.inno_sources_build_path, dirs_exist_ok=True)
#     self.inno_build_third_party_path.mkdir(exist_ok=True)
#     self.inno_build_prerequisite_path.mkdir(exist_ok=True)
#     shutil.copy(tmp_vc_redist_setup_filepath, pathlib.Path(self.inno_build_third_party_path / "VC_redist.x64.exe"))
#     shutil.copy(tmp_windows_tasks_exe_filepath, pathlib.Path(self.inno_build_prerequisite_path / "WindowsTasks.exe"))
#     self.inno_build_assets_path.mkdir(exist_ok=True)
#     shutil.copy(tmp_pyssa_win_build_logo_filepath, pathlib.Path(self.inno_build_assets_path / "logo.ico"))
#     # </editor-fold>
#
#   def build(self) -> None:
#     """Builds the PyMOL Windows EXE file."""
#     try:
#       tmp_start_time = time.time()
#       self.setup_build_environment()
#       subprocess.run(
#         [self.inno_setup_compiler_filepath, self.inno_setup_script_filepath],
#         stdout=sys.stdout, stderr=sys.stderr, text=True
#       )
#       tmp_end_time = time.time()
#       tmp_duration = tmp_end_time - tmp_start_time
#       print(f"The build process of the inno setup EXE took: {tmp_duration:.2f} seconds ({(tmp_duration/60):.2f} minutes).")
#     except Exception as e:
#       print(e)




class BuildInnoSetupDebug:
  """Contains the logic for building a debug version of the inno setup EXE file."""

  def __init__(self) -> None:
    """Constructor."""
    self.deployment_resources_path = pathlib.Path(PROJECT_ROOT_DIR / "deployment/resources")
    self.inno_build_path = pathlib.Path(PROJECT_ROOT_DIR / "inno-build-debug")
    self.inno_build_assets_path = pathlib.Path(self.inno_build_path / "inno-assets")
    self.inno_build_cache_path = pathlib.Path(self.inno_build_path / "inno-cache")
    self.inno_sources_build_path = pathlib.Path(self.inno_build_path / "inno-sources")
    self.inno_build_prerequisite_path = pathlib.Path(self.inno_sources_build_path / "prerequisite")
    self.inno_build_third_party_path = pathlib.Path(self.inno_sources_build_path / "third_party")
    self.inno_build_tmp_path = pathlib.Path(self.inno_sources_build_path / "tmp")
    self.inno_setup_script_path = pathlib.Path(PROJECT_ROOT_DIR / "deployment/src/inno_setup")
    self.inno_setup_script_filepath = pathlib.Path(PROJECT_ROOT_DIR / "deployment/src/inno_setup" / "setup.iss")
    self.inno_setup_compiler_filepath = pathlib.Path(r"C:\Program Files (x86)\Inno Setup 6\ISCC.exe")

  def clean_build_directory(self) -> None:
    """Cleans the inno setup build directory."""
    pass

  def setup_build_environment(self, include_wsl2_distro: bool = False, do_setup_build: bool = False) -> None:
    """Sets up a temporary build environment."""
    # <editor-fold desc="Configure for inno setup build">
    if do_setup_build:
      self.inno_build_path = pathlib.Path(PROJECT_ROOT_DIR / "inno-build-release")
      if self.inno_build_path.exists():
        shutil.rmtree(self.self.inno_build_path)
        self.inno_build_path.mkdir()
    # </editor-fold>
    # <editor-fold desc="Path/Filepath definitions">
    tmp_pyssa_win_build_logo_filepath = pathlib.Path(PROJECT_ROOT_DIR / "assets/convert_logo_to_ico" / "logo.ico")
    tmp_vc_redist_setup_filepath = pathlib.Path(PROJECT_ROOT_DIR / "third_party/microsoft" / "VC_redist.x64.exe")
    tmp_windows_tasks_exe_filepath = pathlib.Path(PROJECT_ROOT_DIR / "deployment/offline_resources" / "WindowsTasks.exe")
    # </editor-fold>
    """IMPORTANT
    Use the python interpreter of the venv of pymol windows build because
    that interpreter gets also used in the build script of the 
    pymol windows build repo!
    """
    # <editor-fold desc="Restore build directory for new build">
    if self.inno_build_assets_path.exists():
      shutil.rmtree(self.inno_build_assets_path)
    if self.inno_sources_build_path.exists():
      shutil.rmtree(self.inno_sources_build_path)
    self.inno_build_path.mkdir(exist_ok=True)
    self.inno_build_assets_path.mkdir()
    self.inno_sources_build_path.mkdir()
    self.inno_build_cache_path.mkdir(exist_ok=True)
    # </editor-fold>
    # <editor-fold desc="Get WSL2 distro from sciebo">
    pathlib.Path(self.inno_build_tmp_path).mkdir()
    if include_wsl2_distro:
      if not pathlib.Path.exists(pathlib.Path(self.inno_build_cache_path / "alma-colabfold-9-rootfs.tar")):
        print("Downloading alma-colabfold-9-rootfs.tar ...")
        if not setup_util.download_file("https://w-hs.sciebo.de/s/q5oYjcZdEzCDyEH/download", str(pathlib.Path(self.inno_build_cache_path / "alma-colabfold-9-rootfs.tar"))):
          print("Unable to download alma-colabfold-9-rootfs.tar, build process exists.")
          return
        print("Finished downloading alma-colabfold-9-rootfs.tar.")

      if not setup_util.File.copy(
              pathlib.Path(self.inno_build_cache_path / "alma-colabfold-9-rootfs.tar"),
              pathlib.Path(self.inno_build_tmp_path / "alma-colabfold-9-rootfs.tar"),
      ):
        print("Copying the alma-colabfold-9-rootfs.tar failed!")
        return
    # </editor-fold>
    # Extract bin.zip file
    # <editor-fold desc="Copy operations">
    setup_util.Directory.copy_directory(
      pathlib.Path(PROJECT_ROOT_DIR / "build/exe.win-amd64-3.11"),
      pathlib.Path(self.inno_sources_build_path)
    )
    setup_util.Directory.copy_directory(
      pathlib.Path(PROJECT_ROOT_DIR / "build/user_pymol"),
      pathlib.Path(self.inno_sources_build_path / "user_pymol")
    )
    self.inno_build_third_party_path.mkdir(exist_ok=True)
    self.inno_build_prerequisite_path.mkdir(exist_ok=True)
    shutil.copy(tmp_vc_redist_setup_filepath, pathlib.Path(self.inno_build_third_party_path / "VC_redist.x64.exe"))
    shutil.copy(tmp_windows_tasks_exe_filepath, pathlib.Path(self.inno_build_prerequisite_path / "WindowsTasks.exe"))
    self.inno_build_assets_path.mkdir(exist_ok=True)
    shutil.copy(tmp_pyssa_win_build_logo_filepath, pathlib.Path(self.inno_build_assets_path / "logo.ico"))
    # </editor-fold>

  def build(self, a_inno_script_path: pathlib.Path) -> None:
    """Builds the PyMOL Windows EXE file."""
    try:
      tmp_start_time = time.time()
      subprocess.run(
        [self.inno_setup_compiler_filepath, a_inno_script_path],
        stdout=sys.stdout, stderr=sys.stderr, text=True
      )
      tmp_end_time = time.time()
      tmp_duration = tmp_end_time - tmp_start_time
      print(f"The build process of the inno setup EXE took: {tmp_duration:.2f} seconds ({(tmp_duration/60):.2f} minutes).")
    except Exception as e:
      print(e)


class BuildLiveLinuxExec:
  """Class for building an environment to run and test the linux version."""

  def build(self) -> None:
    """Builds a live linux execution environment containing the user pymol."""
    tmp_user_pymol_path = pathlib.Path(PROJECT_ROOT_DIR / "third_party/pymol-oss/linux-build/dist/exe.linux-x86_64-3.11")
    if not tmp_user_pymol_path.exists():
      print(f"User PyMOL distribution could not be found under {tmp_user_pymol_path}!")
      exit(-1)

    subprocess.run([f"{PROJECT_ROOT_DIR}/.venv/bin/python", "setup.py", "build"])
    shutil.copytree(
      tmp_user_pymol_path,
      pathlib.Path(PROJECT_ROOT_DIR / "build/exe.linux-x86_64-3.11/user_pymol")
    )

# </editor-fold>


# <editor-fold desc="Automation functions">



def build_win_package() -> None:
  """Builds the windows package ZIP file."""
  tmp_builder = BuildWinPackage()
  tmp_builder.build()


def build_setup_exe() -> None:
  """Builds the inno setup EXE file."""
  tmp_builder = BuildInnoSetup()
  tmp_builder.setup_build_environment(include_wsl2_distro=True)
  tmp_builder.build(tmp_builder.inno_setup_script_filepath)


def build_update_src_only_exe() -> None:
  """Builds the update src only inno setup EXE file."""
  tmp_builder = BuildInnoSetup()
  tmp_builder.setup_build_environment()
  tmp_builder.build(
    pathlib.Path(tmp_builder.inno_setup_script_path / "setup_only_src.iss")
  )


def clean_build_setup_exe() -> None:
  """Cleans the inno setup build directory."""
  shutil.rmtree(pathlib.Path(PROJECT_ROOT_DIR / "innoBuild"))


def build_setup_exe_dbg():
  """Builds the setup exe directory n debug mode."""
  tmp_builder = BuildInnoSetupDebug()
  tmp_builder.setup_build_environment(include_wsl2_distro=True)
  tmp_builder.build(tmp_builder.inno_setup_script_filepath)
  # TODO: Comment build line


def build_portable() -> None:
  """Builds the portable installation ZIP file."""
  tmp_builder = BuildInnoSetup()
  tmp_builder.setup_build_environment()
  # Construct the PowerShell command
  shutil.copy(
    pathlib.Path(PROJECT_ROOT_DIR / "src/inno_setup/LICENSE.txt"),
    pathlib.Path(PROJECT_ROOT_DIR / "innoBuild/inno-sources/LICENSE.txt"),
  )
  tmp_path = pathlib.Path(PROJECT_ROOT_DIR / "innoBuild/inno-sources" / "*")
  # TODO: Make the version number variable so that it can be changed in one place
  tmp_filename = "PyMOL_Open_source_v3.1.0a0_WINx64_portable.zip"
  powershell_command = [
    "powershell",
    "-Command",
    f"Compress-Archive -Path '{tmp_path}' -DestinationPath '{tmp_filename}' -Force"
  ]
  # Run the PowerShell command
  subprocess.run(powershell_command, check=True)
  shutil.rmtree(pathlib.Path(PROJECT_ROOT_DIR / "innoBuild"))


def build_linux_live():
  """Builds the live linux exec env incl. user_pymol."""
  tmp_builder = BuildLiveLinuxExec()
  tmp_builder.build()

# </editor-fold>


def main() -> None:
  """Main function."""
  parser = argparse.ArgumentParser(description="Automation script with subcommands.")
  # <editor-fold desc="Subparsers">
  subparsers = parser.add_subparsers(dest='command')
  install_parser = subparsers.add_parser('setup-dev-env', help="Installs build dependencies.")
  install_parser.set_defaults(func=setup_dev_env)
  build_win_package_parser = subparsers.add_parser('build-win-package', help="Builds the windows package ZIP file.")
  build_win_package_parser.set_defaults(func=build_win_package)
  build_setup_exe_parser = subparsers.add_parser('build-setup-exe', help="Builds the inno setup EXE file.")
  build_setup_exe_parser.set_defaults(func=build_setup_exe)
  build_update_src_only_exe_parser = subparsers.add_parser('build-update-src-exe', help="Builds update src the inno setup EXE file.")
  build_update_src_only_exe_parser.set_defaults(func=build_update_src_only_exe)
  clean_build_setup_exe_parser = subparsers.add_parser('clean', help="Cleans the inno setup build directory.")
  clean_build_setup_exe_parser.set_defaults(func=clean_build_setup_exe)
  build_setup_exe_dbg_parser = subparsers.add_parser('build-setup-exe-dbg', help="Builds a debug mode inno setup.")
  build_setup_exe_dbg_parser.set_defaults(func=build_setup_exe_dbg)
  build_portable_parser = subparsers.add_parser('build-portable', help="Builds the portable installation zip.")
  build_portable_parser.set_defaults(func=build_portable)
  build_linux_live_parser = subparsers.add_parser('build-linux-live', help="Builds the live linux exec env incl. user_pymol.")
  build_linux_live_parser.set_defaults(func=build_linux_live)
  # </editor-fold>
  args = parser.parse_args()

  if args.command:
    args.func()


if __name__ == "__main__":
  main()
