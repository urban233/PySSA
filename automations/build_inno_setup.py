import pathlib
import shutil
import subprocess
import sys
import time
import zipfile
import tarfile

from task_automator.IO import file
from task_automator.IO import directory
from task_automator.utils import web_utils

import const
import build_win_exe


class BuildInnoSetup:
  """Contains the logic for building the inno setup EXE file."""

  def __init__(self) -> None:
    """Constructor."""
    self.deployment_resources_path = pathlib.Path(const.PROJECT_ROOT_DIR / "deployment/resources")
    self.virtual_env_path = pathlib.Path(const.PROJECT_ROOT_DIR / ".venv")
    self.original_pyssa_source_path = pathlib.Path(const.PROJECT_ROOT_DIR / "src")
    self.inno_build_path = pathlib.Path(const.PROJECT_ROOT_DIR / "ib-release")
    self.inno_build_assets_path = pathlib.Path(self.inno_build_path / "assets")
    self.inno_build_cache_path = pathlib.Path(self.inno_build_path / "cache")
    self.inno_sources_build_path = pathlib.Path(self.inno_build_path / "sources")
    self.pyssa_program_path = pathlib.Path(self.inno_build_path / "sources/bin/PySSA")
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
    # tmp_windows_tasks_exe_filepath = pathlib.Path(const.PROJECT_ROOT_DIR / "deployment/offline_resources" / "WindowsCli.exe")
    # </editor-fold>
    """IMPORTANT
    Use the python interpreter of the venv of pymol windows build because
    that interpreter gets also used in the build script of the 
    pymol windows build repo!
    """
    # <editor-fold desc="Restore build directory for new build">
    directory.Directory.purge(self.inno_build_assets_path)
    directory.Directory.purge(self.inno_sources_build_path)
    self.inno_build_path.mkdir(exist_ok=True)
    self.inno_build_assets_path.mkdir()
    self.inno_sources_build_path.mkdir()
    self.inno_build_cache_path.mkdir(exist_ok=True)
    # </editor-fold>
    pathlib.Path(self.inno_build_tmp_path).mkdir()
    # <editor-fold desc="Download cpython-3.11.14+20260211-x86_64-pc-windows-msvc-install_only.tar.gz">
    if not pathlib.Path.exists(pathlib.Path(self.inno_build_cache_path / "cpython-3.11.14+20260211-x86_64-pc-windows-msvc-install_only")):
      print("Downloading cpython-3.11.14+20260211-x86_64-pc-windows-msvc-install_only.tar.gz ...")
      tmp_ret = web_utils.download_file(
        "https://github.com/astral-sh/python-build-standalone/releases/download/20260211/cpython-3.11.14+20260211-x86_64-pc-windows-msvc-install_only.tar.gz",
        str(pathlib.Path(self.inno_build_cache_path / "cpython-3.11.14+20260211-x86_64-pc-windows-msvc-install_only.tar.gz"))
      )
      if not tmp_ret:
        print("Unable to download cpython-3.11.14+20260211-x86_64-pc-windows-msvc-install_only.tar.gz, build process exists.")
        return
      print("Finished downloading cpython-3.11.14+20260211-x86_64-pc-windows-msvc-install_only.tar.gz.")
      with tarfile.open(pathlib.Path(self.inno_build_cache_path / "cpython-3.11.14+20260211-x86_64-pc-windows-msvc-install_only.tar.gz"), "r:gz") as tar:
        tar.extractall(path=pathlib.Path(self.inno_build_cache_path / "cpython-3.11.14+20260211-x86_64-pc-windows-msvc-install_only"))
      pathlib.Path(self.inno_build_cache_path / "cpython-3.11.14+20260211-x86_64-pc-windows-msvc-install_only.tar.gz").unlink()
    # </editor-fold>

    # <editor-fold desc="Setup Python environment">
    if not directory.Directory.copy_directory(
      pathlib.Path(self.inno_build_cache_path / "cpython-3.11.14+20260211-x86_64-pc-windows-msvc-install_only"),
      pathlib.Path(self.inno_sources_build_path / "cpython-3.11.14"),
    ):
      print("Copying the cpython-3.11.14+20260211-x86_64-pc-windows-msvc-install_only directory failed!")
      exit(1)

    if not directory.Directory.copy_directory(
            pathlib.Path(self.original_pyssa_source_path),
            pathlib.Path(self.inno_sources_build_path / "cpython-3.11.14/python/Lib/site-packages/src")
    ):
      print("Copying the original source folder directory failed!")
      exit(1)

    subprocess.run(
      [
        str(pathlib.Path(self.inno_sources_build_path / "cpython-3.11.14/python/python.exe")),
        "-m", "pip", "install", "-r", str(pathlib.Path(const.PROJECT_ROOT_DIR / "requirements.txt"))
      ],
      stdout=sys.stdout, stderr=sys.stderr, text=True
    )

    subprocess.run(
      [
        str(pathlib.Path(self.inno_sources_build_path / "cpython-3.11.14/python/python.exe")),
        "-m", "pip", "install", "pymol-open-source-whl"
      ],
      stdout=sys.stdout, stderr=sys.stderr, text=True
    )
    # </editor-fold>

    # <editor-fold desc="Patching PyMOL source">
    # Patch the pymol_gl_widget.py for HighDpi support
    if not file.File.copy(pathlib.Path(self.deployment_resources_path / "pymol_gl_widget.py"),
                          pathlib.Path(self.inno_sources_build_path / "cpython-3.11.14/python/Lib/site-packages/pmg_qt/pymol_gl_widget.py")):
      print("Copying the pymol_gl_widget.py file failed!")
      exit(1)
    # Patch the invocation.py for custom options
    if not file.File.copy(pathlib.Path(self.deployment_resources_path / "invocation.py"),
                          pathlib.Path(self.inno_sources_build_path / "cpython-3.11.14/python/Lib/site-packages/pymol/invocation.py")):
      print("Copying the invocation.py file failed!")
      exit(1)
    if not file.File.copy(pathlib.Path(self.deployment_resources_path / "controlling.py"),
                          pathlib.Path(self.inno_sources_build_path / "cpython-3.11.14/python/Lib/site-packages/pymol/controlling.py")):
      print("Copying the controlling.py file failed!")
      exit(1)
    # </editor-fold>

    # <editor-fold desc="Get WSL2 distro from sciebo">
    if include_wsl2_distro:
      if not pathlib.Path.exists(pathlib.Path(self.inno_build_cache_path / "alma-colabfold-9-rootfs.tar")):
        print("Downloading alma-colabfold-9-rootfs.tar ...")
        if not web_utils.download_file("https://w-hs.sciebo.de/s/yOQ8Qo1Uvk1eaQc/download", str(pathlib.Path(self.inno_build_cache_path / "alma-colabfold-9-rootfs.tar"))):
          print("Unable to download alma-colabfold-9-rootfs.tar, build process exists.")
          exit(1)
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

    if not file.File.copy(
      pathlib.Path(self.deployment_resources_path / "start_pyssa.bat"),
      pathlib.Path(self.inno_sources_build_path / "start_pyssa.bat"),
      overwrite=True
    ):
      print("Copying the start_pyssa.bat file failed!")
      exit(1)
    if not file.File.copy(
            pathlib.Path(self.deployment_resources_path / "uninstall_helper.bat"),
            pathlib.Path(self.inno_sources_build_path / "uninstall_helper.bat"),
            overwrite=True
    ):
      print("Copying the setup.bat file failed!")
      exit(1)
    self.inno_build_third_party_path.mkdir(exist_ok=True)
    self.inno_build_prerequisite_path.mkdir(exist_ok=True)
    if not file.File.copy(tmp_vc_redist_setup_filepath, pathlib.Path(self.inno_build_third_party_path / "VC_redist.x64.exe")):
      print("Copying the VC_redist.x64.exe file failed!")
      exit(1)
    self.inno_build_assets_path.mkdir(exist_ok=True)
    if not file.File.copy(tmp_pyssa_win_build_logo_filepath, pathlib.Path(self.inno_build_assets_path / "logo.ico")):
      print("Copying the logo.ico file failed!")
      exit(1)
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
