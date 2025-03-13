import pathlib
import shutil
import subprocess
import sys
import time

from task_automator.IO import file
from task_automator.utils import web_utils

import const
import build_win_exe


class BuildInnoSetup:
  """Contains the logic for building the inno setup EXE file."""

  def __init__(self) -> None:
    """Constructor."""
    self.deployment_resources_path = pathlib.Path(const.PROJECT_ROOT_DIR / "deployment/resources")
    self.inno_build_path = pathlib.Path(const.PROJECT_ROOT_DIR / "inno-build-release")
    self.inno_build_assets_path = pathlib.Path(self.inno_build_path / "inno-assets")
    self.inno_build_cache_path = pathlib.Path(self.inno_build_path / "inno-cache")
    self.inno_sources_build_path = pathlib.Path(self.inno_build_path / "inno-sources")
    self.pyssa_program_path = pathlib.Path(self.inno_build_path / "inno-sources/bin/PySSA")
    self.inno_build_prerequisite_path = pathlib.Path(self.inno_sources_build_path / "prerequisite")
    self.inno_build_third_party_path = pathlib.Path(self.inno_sources_build_path / "third_party")
    self.inno_build_tmp_path = pathlib.Path(self.inno_sources_build_path / "tmp")
    self.inno_setup_script_path = pathlib.Path(const.PROJECT_ROOT_DIR / "deployment/src/inno_setup")
    self.inno_setup_script_filepath = pathlib.Path(const.PROJECT_ROOT_DIR / "deployment/src/inno_setup" / "setup.iss")
    self.inno_setup_update_script_filepath = pathlib.Path(const.PROJECT_ROOT_DIR / "deployment/src/inno_setup" / "setup_only_src.iss")
    self.inno_setup_compiler_filepath = pathlib.Path(r"C:\Program Files (x86)\Inno Setup 6\ISCC.exe")

  def setup_build_environment(self, include_wsl2_distro: bool = False) -> None:
    """Sets up a temporary build environment."""
    # <editor-fold desc="Path/Filepath definitions">
    tmp_pyssa_win_build_logo_filepath = pathlib.Path(const.PROJECT_ROOT_DIR / "assets/convert_logo_to_ico" / "logo.ico")
    tmp_vc_redist_setup_filepath = pathlib.Path(const.PROJECT_ROOT_DIR / "third_party/microsoft" / "VC_redist.x64.exe")
    tmp_windows_tasks_exe_filepath = pathlib.Path(const.PROJECT_ROOT_DIR / "deployment/offline_resources" / "WindowsTasks.exe")
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
    build_win_exe.build()
    # <editor-fold desc="Copy frozen PySSA Python venv">
    shutil.copytree(
      pathlib.Path(const.PROJECT_ROOT_DIR / "build/exe.win-amd64-3.11"),
      pathlib.Path(self.inno_sources_build_path),
      dirs_exist_ok=True
    )
    shutil.copytree(
      pathlib.Path(const.PROJECT_ROOT_DIR / "build/user_pymol"),
      pathlib.Path(self.inno_sources_build_path / "user_pymol"),
      dirs_exist_ok=True
    )
    # </editor-fold>
    # <editor-fold desc="Get WSL2 distro from sciebo">
    pathlib.Path(self.inno_build_tmp_path).mkdir()
    if include_wsl2_distro:
      if not pathlib.Path.exists(pathlib.Path(self.inno_build_cache_path / "alma-colabfold-9-rootfs.tar")):
        print("Downloading alma-colabfold-9-rootfs.tar ...")
        if not web_utils.download_file("https://w-hs.sciebo.de/s/q5oYjcZdEzCDyEH/download", str(pathlib.Path(self.inno_build_cache_path / "alma-colabfold-9-rootfs.tar"))):
          print("Unable to download alma-colabfold-9-rootfs.tar, build process exists.")
          return
        print("Finished downloading alma-colabfold-9-rootfs.tar.")
      else:
        print("Using cached version of alma-colabfold-9-rootfs.tar under inno-build-release/inno-cache")

      if not file.File.copy(
              pathlib.Path(self.inno_build_cache_path / "alma-colabfold-9-rootfs.tar"),
              pathlib.Path(self.inno_build_tmp_path / "alma-colabfold-9-rootfs.tar"),
      ):
        print("Copying the alma-colabfold-9-rootfs.tar failed!")
        exit(1)
    # </editor-fold>
    # <editor-fold desc="Copy operations">
    if not file.File.copy(
      pathlib.Path(self.deployment_resources_path / "setup.bat"),
      pathlib.Path(self.inno_build_tmp_path / "setup.bat"),
      overwrite=True
    ):
      print("Copying the setup.bat file failed!")
      exit(1)
    self.inno_build_third_party_path.mkdir(exist_ok=True)
    self.inno_build_prerequisite_path.mkdir(exist_ok=True)
    shutil.copy(tmp_vc_redist_setup_filepath, pathlib.Path(self.inno_build_third_party_path / "VC_redist.x64.exe"))
    shutil.copy(tmp_windows_tasks_exe_filepath, pathlib.Path(self.inno_build_prerequisite_path / "WindowsTasks.exe"))
    self.inno_build_assets_path.mkdir(exist_ok=True)
    shutil.copy(tmp_pyssa_win_build_logo_filepath, pathlib.Path(self.inno_build_assets_path / "logo.ico"))
    # </editor-fold>

  def build(self, a_inno_script_path: pathlib.Path) -> None:
    """Builds the inno setup EXE file."""
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


def build_full_setup_exe() -> None:
  """Builds the full inno setup EXE file."""
  tmp_builder = BuildInnoSetup()
  tmp_builder.setup_build_environment(include_wsl2_distro=True)
  tmp_builder.build(tmp_builder.inno_setup_script_filepath)


def build_update_setup_exe() -> None:
  """Builds the update inno setup EXE file."""
  tmp_builder = BuildInnoSetup()
  tmp_builder.setup_build_environment(include_wsl2_distro=False)
  tmp_builder.build(tmp_builder.inno_setup_update_script_filepath)
