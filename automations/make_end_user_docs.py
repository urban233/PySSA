import pathlib
import subprocess

import const


def make_docs() -> None:
  """Creates the HTML files for the End-User documentation."""
  tmp_venv_activate_bat_filepath = pathlib.Path(const.PROJECT_ROOT_DIR / ".venv/Scripts" / "activate.bat")
  tmp_make_bat_filepath = pathlib.Path(const.PROJECT_ROOT_DIR / "docs" / "make.bat")

  activate_cmd = f'& "{tmp_venv_activate_bat_filepath}"'
  combined_cmd = f'{activate_cmd} ; {tmp_make_bat_filepath} html'
  subprocess.run(['powershell', '-Command', combined_cmd])
