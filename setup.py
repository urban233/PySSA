# from setuptools import setup
# from setuptools import Command
# import setup_build_tasks
#
#
# class CreateWinPackage(Command):
#     """Setup command for creating the Windows package that gets deployed."""
#     description = 'Create a custom ZIP file containing source and additional files'
#     user_options = []
#
#     def initialize_options(self):
#         """Override the default initialization options."""
#         pass
#
#     def finalize_options(self):
#         """Override the default finalization options."""
#         pass
#
#     def run(self):
#         """Override the default run method."""
#         tmp_build_task = setup_build_tasks.CreateWinPackageBuildTask()
#         tmp_build_task.execute_task()
#
#
# class MakeDocs(Command):
#     """Build command for creating sphinx documentation."""
#     description = 'Activates the .venv and runs the make.bat for building the sphinx docs.'
#     user_options = []
#
#     def initialize_options(self):
#         """Override the default initialization options."""
#         pass
#
#     def finalize_options(self):
#         """Override the default finalization options."""
#         pass
#
#     def run(self):
#         """Override the default run method."""
#         tmp_build_task = setup_build_tasks.MakeSphinxDocs()
#         tmp_build_task.execute_task()
#
#
# class MakeVenv(Command):
#     """Build command for creating sphinx documentation."""
#     description = 'Export the .venv information into a deployable format.'
#     user_options = []
#
#     def initialize_options(self):
#         """Override the default initialization options."""
#         pass
#
#     def finalize_options(self):
#         """Override the default finalization options."""
#         pass
#
#     def run(self):
#         """Override the default run method."""
#         tmp_build_task = setup_build_tasks.MakeVenvForDeployment()
#         tmp_build_task.execute_task()
#
#
# setup(
#     cmdclass={
#         'create_win_package': CreateWinPackage,
#         'make_docs': MakeDocs,
#         'make_venv': MakeVenv,
#     },
# )

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
    # shutil.copytree(
    #   pathlib.Path("./.venv/Lib/site-packages/pyzmq.libs"),
    #   pathlib.Path(build_dir / "lib/pyzmq.libs")
    # )
    # shutil.copytree(
    #   pathlib.Path("./.venv/lib64/python3.11/site-packages/pyzmq.libs"),
    #   pathlib.Path(build_dir / "lib/pyzmq.libs")
    # )
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

build_exe_options = {
  # exclude packages that are not really needed
  "excludes": ["tkinter"],
  "include_files": include_files,
  "packages": ["pymol.povray", "pymol.parser", "Xlib"] # Xlib -> linux
}

bdist_mac_options = {
  "bundle_name": "Test",
}

bdist_dmg_options = {
  "volume_label": "TEST",
}

setup(
  name="PySSA",
  version="1.0.2",
  description="Build script for PySSA",
  options={
    "build_exe": build_exe_options,
    "bdist_mac": bdist_mac_options,
    "bdist_dmg": bdist_dmg_options,
  },
  executables=[
    {"script": pathlib.Path("src/pyssa/main.py"), "target_name": "pyssa"},  # base="gui"
    {"script": pathlib.Path("src/auxiliary_pymol/main.py"), "target_name": "aux_pymol"}
  ],
  cmdclass={
    "build_exe": CustomBuildExe
  },
)
