from __future__ import annotations

import pathlib
import shutil

from cx_Freeze import Executable, setup, build_exe

PROJECT_ROOT_PATH = pathlib.Path(__file__).parent


class CustomBuildExe(build_exe):
  def run(self):
    super().run()
    print("\nRunning post-build steps...")
    build_dir: pathlib.Path = pathlib.Path(self.build_exe)
    # <editor-fold desc="Copying necessary dependencies">
    print("Copying necessary dependencies ...")
    shutil.copytree(
      pathlib.Path("./.venv/Lib/site-packages/pyzmq.libs"),
      pathlib.Path(build_dir / "lib/pyzmq.libs")
    )
    shutil.rmtree(build_dir / "lib/src")
    shutil.copytree(
      pathlib.Path(PROJECT_ROOT_PATH / "src"),
      pathlib.Path(build_dir / "lib/src"),
      dirs_exist_ok=True
    )
    # </editor-fold>
    # <editor-fold desc="Copying additional files">
    print("Copying additional files ...")
    shutil.copytree(
      pathlib.Path(PROJECT_ROOT_PATH / "assets"),
      pathlib.Path(build_dir / "assets"),
    )
    shutil.copytree(
      pathlib.Path(PROJECT_ROOT_PATH / "docs"),
      pathlib.Path(build_dir / "docs"),
    )
    shutil.copytree(
      pathlib.Path(PROJECT_ROOT_PATH / "scripts"),
      pathlib.Path(build_dir / "scripts"),
    )
    shutil.copytree(
      pathlib.Path(PROJECT_ROOT_PATH / "winbatch"),
      pathlib.Path(build_dir / "winbatch"),
    )
    shutil.copy(
      pathlib.Path(PROJECT_ROOT_PATH / "LICENSE"),
      pathlib.Path(build_dir / "LICENSE"),
    )
    shutil.copy(
      pathlib.Path(PROJECT_ROOT_PATH / "README.md"),
      pathlib.Path(build_dir / "README.md"),
    )
    shutil.copy(
      pathlib.Path(PROJECT_ROOT_PATH / ".pymolrc.py"),
      pathlib.Path(build_dir / ".pymolrc.py"),
    )
    # </editor-fold>
    print("\nFinished post-build steps.")


try:
  from cx_Freeze.hooks import get_qt_plugins_paths
except ImportError:
  get_qt_plugins_paths = None

include_files = []
if get_qt_plugins_paths:
    # Inclusion of extra plugins (since cx_Freeze 6.8b2)
    # cx_Freeze automatically imports the following plugins depending on the
    # module used, but suppose we need the following:
    include_files += get_qt_plugins_paths("PyQt5", "multimedia")

build_exe_options = {
  # exclude packages that are not really needed
  "excludes": ["tkinter"],
  "include_files": include_files,
  "packages": ["PyQt5.uic", "pymol.povray", "pymol.parser"]
}

setup(
  name="PySSA",
  version="1.0.7",
  description="Build script for PySSA",
  options={
    "build_exe": build_exe_options,
  },
  executables=[
    {"script": pathlib.Path("src/pyssa/main.py"), "target_name": "pyssa", "base": "Win32GUI"},
    {"script": pathlib.Path("src/auxiliary_pymol/main.py"), "target_name": "aux_pymol"}
  ],
  cmdclass={
    "build_exe": CustomBuildExe
  },
)
