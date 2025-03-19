import pathlib
import subprocess

import const


def build() -> None:
  """Freezes the PySSA application using cx_freeze."""
  subprocess.run(
    [
      const.PYTHON_EXECUTABLE,
      pathlib.Path(const.PROJECT_ROOT_DIR / "setup.py"),
      "build"
    ]
  )
