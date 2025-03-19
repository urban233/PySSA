#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/urban233/PySSA>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""Module that is used to start PySSA."""
import os
import pathlib
import sys
import subprocess


def main():
  tmp_root_path = pathlib.Path(__file__).parent
  sys.path.append(str(tmp_root_path / "lib"))
  sys.path.append(str(tmp_root_path / "user_pymol/lib"))

  from PyQt5 import QtWidgets
  from PyQt5 import QtGui
  from PyQt5.QtCore import Qt

  # Check the session type
  # session_type = os.environ.get('XDG_SESSION_TYPE', 'x11')  # Default to 'x11' if not set
  #
  # if session_type == 'wayland':
  #   os.environ['QT_QPA_PLATFORM'] = 'wayland'
  # else:
  #   os.environ['QT_QPA_PLATFORM'] = 'xcb'

  from src.pyssa.util import constants
  from src.pyssa.gui.ui.styles import styles
  from src.pyssa.controller import main_view_controller
  from src.pyssa.controller import interface_manager

  app = QtWidgets.QApplication(sys.argv)
  # setup QSplashScreen
  pixmapi = QtGui.QPixmap(
    f"{constants.PROGRAM_BIN_ROOT_PATH}\\assets\\images\\splash_screen.png"
  )
  smaller_pixmapi = pixmapi.scaled(
    700, 700, Qt.KeepAspectRatio, Qt.SmoothTransformation
  )
  tmp_splash = QtWidgets.QSplashScreen(smaller_pixmapi)
  tmp_splash.show()
  # Begin with PySSA startup
  styles.set_stylesheet(app)
  interfaceManager = interface_manager.InterfaceManager()
  main_window = interfaceManager.get_main_view()
  main_controller = main_view_controller.MainViewController(interfaceManager)
  styles.set_stylesheet_homepage(main_window)
  main_window.show()
  tmp_splash.finish(None)
  subprocess.Popen(
    [constants.ARRANGE_WINDOWS_EXE_FILEPATH],
    creationflags=subprocess.CREATE_NO_WINDOW,
  )
  sys.exit(app.exec_())


if __name__ == "__main__":
  main()
