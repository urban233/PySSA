import pathlib
import subprocess

import const


def setup_dev_env() -> None:
  """Installs the dependencies needed for building the _cmd extension module."""
  subprocess.run(["git", "clone", "https://github.com/urban233/pymol-open-source-windows-build", pathlib.Path("./vendor/pymol-open-source-windows-build")])
  subprocess.run(["powershell.exe", "pwd"], cwd=str(pathlib.Path(const.PROJECT_ROOT_DIR / 'vendor/pymol-open-source-windows-build'))
                 )
  subprocess.run(
    [
      "cmd.exe", "/c", str(pathlib.Path(r'.\setup_dev_env.bat'))
    ],
    cwd=pathlib.Path(const.PROJECT_ROOT_DIR / 'vendor/pymol-open-source-windows-build')
  )
  subprocess.run(
    [
      pathlib.Path(const.PROJECT_ROOT_DIR / "vendor/pymol-open-source-windows-build/.venv/Scripts" / "python.exe"),
      pathlib.Path(const.PROJECT_ROOT_DIR / "vendor/pymol-open-source-windows-build" / "run_automation.py"), 'setup-dev-env'
    ], cwd=pathlib.Path(const.PROJECT_ROOT_DIR / "vendor/pymol-open-source-windows-build")
  )
