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
"""Main window class for the PySSA frontend application.

Authors: Martin Urban, Hannah Kullik

Version: 1.4.0
"""

import pathlib
import os

import sys

import pymol
from pymol._gui import PyMOLDesktopGUI
from pymol.Qt.utils import MainThreadCaller

from pmg_qt.pymol_gl_widget import PyMOLGLWidget
from pmg_qt import keymapping

from src.pyssa.gui import user_pymol

from src.pyssa.gui.ui.views import molecule_objects_panel

from src.pyssa.gui.ui.styles.icon_manager import IconManager
from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.qt import QtWidgets
from src.pyssa.gui.qt import QtGui

from src.pyssa.gui.ui.custom_widgets import dropdown_menu
from src.pyssa.gui.ui.custom_widgets import color_grid
from src.pyssa.gui.ui.custom_widgets import tool_window_layout
from src.pyssa.gui.ui.custom_widgets import quick_access_bar
from src.pyssa.gui.ui.custom_widgets import quick_access_bar_action

__docformat__ = "google"

Qt = QtCore.Qt


class MainWindow(QtWidgets.QMainWindow, PyMOLDesktopGUI):
  """Main window class."""

  # <editor-fold desc="Class attributes">
  viewportsignal = QtCore.Signal(int, int)
  """A signal for thread-safe viewport command."""
  dialogClosed = QtCore.pyqtSignal(tuple)
  """A signal indicating that the dialog is closed."""
  # </editor-fold>

  def __init__(self):  # noqa
    QtWidgets.QMainWindow.__init__(self)
    # <editor-fold desc="Viewer toolbar actions">
    self.viewer_toolbar_actions = {
      "open_session": quick_access_bar_action.QuickAccessBarAction(
        "Open Session", "left", 0, None, IconManager.instance().get_icon(IconManager.Icons.OPEN_IN_NEW)
      ),
      "create_scene": quick_access_bar_action.QuickAccessBarAction(
        "Create Scene", "left", 0, None, IconManager.instance().get_icon(IconManager.Icons.ADD_CIRCLE)
      ),
      "save_scene": quick_access_bar_action.QuickAccessBarAction(
        "Save Scene", "left", 0, None, IconManager.instance().get_icon(IconManager.Icons.CHANGE_CIRCLE)
      ),
      "delete_scene": quick_access_bar_action.QuickAccessBarAction(
        "Delete Scene", "left", 0, None, IconManager.instance().get_icon(IconManager.Icons.CANCEL)
      ),
      "cartoon": quick_access_bar_action.QuickAccessBarAction(
        "Cartoon", "left", 1, None, IconManager.instance().get_icon(IconManager.Icons.CARTOON_REPR)
      ),
      "sticks": quick_access_bar_action.QuickAccessBarAction(
        "Sticks", "left", 2, None, IconManager.instance().get_icon(IconManager.Icons.STICKS_REPR)
      ),
      "ribbon": quick_access_bar_action.QuickAccessBarAction(
        "Ribbon", "left", 3, None, IconManager.instance().get_icon(IconManager.Icons.RIBBON_REPR)
      ),
      "lines": quick_access_bar_action.QuickAccessBarAction(
        "Lines", "left", 4, None, IconManager.instance().get_icon(IconManager.Icons.LINES_REPR)
      ),
      "spheres": quick_access_bar_action.QuickAccessBarAction(
        "Spheres", "left", 5, None, IconManager.instance().get_icon(IconManager.Icons.SPHERES_REPR)
      ),
      "dots": quick_access_bar_action.QuickAccessBarAction(
        "Dots", "left", 6, None, IconManager.instance().get_icon(IconManager.Icons.DOTS_REPR)
      ),
      "mesh": quick_access_bar_action.QuickAccessBarAction(
        "Mesh", "left", 7, None, IconManager.instance().get_icon(IconManager.Icons.MESH_REPR)
      ),
      "surface": quick_access_bar_action.QuickAccessBarAction(
        "Surface", "left", 8, None, IconManager.instance().get_icon(IconManager.Icons.SURFACE_REPR)
      ),
      "color": quick_access_bar_action.QuickAccessBarAction(
        "Color", "left", 0, None, IconManager.instance().get_icon(IconManager.Icons.PALETTE)
      ),
      "running_jobs": quick_access_bar_action.QuickAccessBarAction(
        "Running Jobs", "left", 0, None, IconManager.instance().get_icon(IconManager.Icons.PLAY_CIRCLE)
      ),
      "notifications": quick_access_bar_action.QuickAccessBarAction(
        "Notifications", "left", 0, None, IconManager.instance().get_icon(IconManager.Icons.NOTIFICATIONS)
      ),
    }
    # </editor-fold>

    # <editor-fold desc="Representation actions">
    self.cartoon_show_hide_menu = dropdown_menu.DropDownMenu()
    self.cartoon_show_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY), "Show"
    )
    self.cartoon_hide_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY_OFF), "Hide"
    )
    self.cartoon_show_hide_menu.addAction(self.cartoon_show_action)
    self.cartoon_show_hide_menu.addAction(self.cartoon_hide_action)

    self.sticks_show_hide_menu = dropdown_menu.DropDownMenu()
    self.sticks_show_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY), "Show"
    )
    self.sticks_hide_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY_OFF), "Hide"
    )
    self.sticks_show_hide_menu.addAction(self.sticks_show_action)
    self.sticks_show_hide_menu.addAction(self.sticks_hide_action)

    self.ribbon_show_hide_menu = dropdown_menu.DropDownMenu()
    self.ribbon_show_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY), "Show"
    )
    self.ribbon_hide_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY_OFF), "Hide"
    )
    self.ribbon_show_hide_menu.addAction(self.ribbon_show_action)
    self.ribbon_show_hide_menu.addAction(self.ribbon_hide_action)

    self.lines_show_hide_menu = dropdown_menu.DropDownMenu()
    self.lines_show_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY), "Show"
    )
    self.lines_hide_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY_OFF), "Hide"
    )
    self.lines_show_hide_menu.addAction(self.lines_show_action)
    self.lines_show_hide_menu.addAction(self.lines_hide_action)

    self.spheres_show_hide_menu = dropdown_menu.DropDownMenu()
    self.spheres_show_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY), "Show"
    )
    self.spheres_hide_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY_OFF), "Hide"
    )
    self.spheres_show_hide_menu.addAction(self.spheres_show_action)
    self.spheres_show_hide_menu.addAction(self.spheres_hide_action)

    self.dots_show_hide_menu = dropdown_menu.DropDownMenu()
    self.dots_show_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY), "Show"
    )
    self.dots_hide_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY_OFF), "Hide"
    )
    self.dots_show_hide_menu.addAction(self.dots_show_action)
    self.dots_show_hide_menu.addAction(self.dots_hide_action)

    self.mesh_show_hide_menu = dropdown_menu.DropDownMenu()
    self.mesh_show_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY), "Show"
    )
    self.mesh_hide_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY_OFF), "Hide"
    )
    self.mesh_show_hide_menu.addAction(self.mesh_show_action)
    self.mesh_show_hide_menu.addAction(self.mesh_hide_action)
    self.surface_show_hide_menu = dropdown_menu.DropDownMenu()
    self.surface_show_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY), "Show"
    )
    self.surface_hide_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY_OFF), "Hide"
    )
    self.surface_show_hide_menu.addAction(self.surface_show_action)
    self.surface_show_hide_menu.addAction(self.surface_hide_action)
    # </editor-fold>
    self.color_grid = color_grid.PyMOLColorGrid()
    self.color_grid_menu = dropdown_menu.DropDownMenu()
    self.color_grid_action = QtWidgets.QWidgetAction(None)

    # <editor-fold desc="PyMOL OpenGL widget">
    # For thread-safe viewport command
    self.viewportsignal.connect(self.pymolviewport)
    self.pymolwidget = PyMOLGLWidget(self)
    cmd = self.cmd = self.pymolwidget.cmd
    self.pymolwidget.installEventFilter(self)
    # </editor-fold>
    # <editor-fold desc="Viewer + toolbar widget">
    self.viewer_toolbar = quick_access_bar.QuickAccessBar(list(self.viewer_toolbar_actions.values()), horizontal=True)
    self.viewer_widget = QtWidgets.QWidget()
    tmp_viewer_layout = QtWidgets.QVBoxLayout()
    tmp_viewer_layout.addWidget(self.viewer_toolbar)
    tmp_viewer_layout.addWidget(self.pymolwidget)
    tmp_viewer_layout.setContentsMargins(0, 0, 0, 0)
    tmp_viewer_layout.setSpacing(0)
    tmp_viewer_layout.setStretch(0, 0)  # toolbar should not stretch vertically
    tmp_viewer_layout.setStretch(1, 1)  # PyMOL widget takes all remaining space
    self.viewer_widget.setLayout(tmp_viewer_layout)
    # </editor-fold>
    """Change behaviour of PyMOL widget

    To change certain options like the display of the internal GUI, you have
    to edit the invocation.py in the site-packages/pymol package!
    Typical changes include:
      options.internal_gui = 0  # This hides the right sidebar
      # The following options are used to hide the internal PyMOL CLI
      options.internal_feedback = 0
      options.show_splash = 0

    To turn off the context menu it is necessary to edit the controlling.py
    module and change all occurrences of these:
      ('double_left','none','menu'),
      ('single_right','none', 'menu'),
    into these:
      ('double_left','none','none'),
      ('single_right','none', 'none'),
    """
    self.user_pymol = user_pymol.UserPyMOL(self.pymolwidget)
    # <editor-fold desc="Panels">
    self.left_side_panel_stacked_widget = QtWidgets.QStackedWidget()
    self.side_panel_molecule_objects = molecule_objects_panel.MoleculeObjectsPanel(self.user_pymol)
    self.right_side_panel_stacked_widget = QtWidgets.QStackedWidget()
    self.bottom_panel_stacked_widget = QtWidgets.QStackedWidget()
    # </editor-fold>
    # self.split_pane = split_pane_design.SplitPaneDesign(
    #   [self.viewer_widget],
    #   self.left_side_panel_stacked_widget,
    #   self.right_side_panel_stacked_widget,
    #   self.bottom_panel_stacked_widget,
    # )
    # self.split_pane.hide_right_side_panel()
    # self.split_pane.hide_bottom_panel()
    # # </editor-fold>
    self._setup_left_side_panels()
    # self._setup_right_side_panels()
    # self._setup_bottom_panels()
    self._setup_color_grid()
    self._central_widget = QtWidgets.QWidget()
    self._main_layout = QtWidgets.QVBoxLayout(self._central_widget)
    self._main_content_layout = QtWidgets.QHBoxLayout()
    self._main_content_layout.setContentsMargins(0, 0, 0, 0)
    self.tool_window_layout = tool_window_layout.ToolWindowLayout(
      None,
      None,
      list(self.viewer_toolbar_actions.values()),
      self.pymolwidget, self
    )
    self.pymolwidget.cmd.fetch("3bmp")
    self.tool_window_layout.set_right_panel_hidden(True)
    self.tool_window_layout.set_bottom_panel_hidden(True)
    # Register panels in ToolWindowLayout in enum order
    # Left stack
    self.tool_window_layout.add_left_panel(self.side_panel_molecule_objects)  # LeftSidePanel.PROTEIN_STRUCTURE = 0

    self._main_content_layout.addWidget(self.tool_window_layout)
    self._main_layout.addLayout(self._main_content_layout)
    self.setWindowTitle("PySSA")
    self.setCentralWidget(self._central_widget)
    self.status_bar = QtWidgets.QStatusBar()
    self.setStatusBar(self.status_bar)
    # self.menu_bar = pymol_menu_bar.PyMOLStandardMenuBar(
    #   self,
    #   self.user_pymol.get_cmd_module(),
    #   MenuBarController(self, self.user_pymol.get_cmd_module())
    # )
    self.menu_bar = QtWidgets.QMenuBar()
    menu = self.menuBar()
    # Apply themed, compact menu bar style (avoid unreliable reset via empty stylesheet)
    # menu.setStyleSheet(ThemeManager.instance().get_current_theme().get_menu_bar_style())
    self.plugins_menu_action = QtGui.QAction(QtGui.QIcon("bug.png"), "&Plugins", self)
    file_menu = menu.addMenu("&Settings")
    file_menu.addAction(self.plugins_menu_action)
    # Apply compact QMenu popup style alongside the main window background
    base_style = "QMainWindow {background-color: #eeeff0;}"
    self.setStyleSheet(base_style)

  # <editor-fold desc="Private methods">
  def _setup_color_grid(self) -> None:
    """Sets up the color grid on the ribbon bar."""
    self.color_grid_action.setDefaultWidget(self.color_grid)
    self.color_grid_menu.addAction(self.color_grid_action)

  def _setup_left_side_panels(self) -> None:
    """Sets up the left side panel."""
    self.left_side_panel_stacked_widget.addWidget(self.side_panel_molecule_objects)

  def _setup_right_side_panels(self) -> None:
    """Sets up the right side panel."""
    self.right_side_panel_stacked_widget.setContentsMargins(0, 0, 0, 0)
    self.right_side_panel_stacked_widget.addWidget(self.side_panel_command_runner)
    self.right_side_panel_stacked_widget.addWidget(self.side_panel_scripts_library)
    self.right_side_panel_stacked_widget.addWidget(self.side_panel_pymol_scenes)
    self.right_side_panel_stacked_widget.addWidget(self.side_panel_move_molecule)

  def _setup_bottom_panels(self) -> None:
    """Sets up the right side panel."""
    self.bottom_panel_stacked_widget.setContentsMargins(0, 0, 0, 0)
    # self.bottom_panel_stacked_widget.addWidget(self.base_bottom_panel)
    self.bottom_panel_stacked_widget.addWidget(self.bottom_panel_command_output)
    # self.bottom_panel_stacked_widget.addWidget(self.side_panel_pymol_scenes)
  # </editor-fold>

  # <editor-fold desc="Public methods">
  def keyPressEvent(self, ev):
    args = keymapping.keyPressEventToPyMOLButtonArgs(ev)
    if args is not None:
      self.pymolwidget.pymol.button(*args)

  def closeEvent(self, event):
    # Emit the custom signal when the window is closed
    self.dialogClosed.emit(("", event))
    self.cmd.quit()

  def pymolviewport(self, w, h):
    cw, ch = self.cmd.get_viewport()
    pw = self.pymolwidget
    scale = pw.fb_scale

    # maintain aspect ratio
    if h < 1:
      if w < 1:
        pw.pymol.reshape(int(scale * pw.width()), int(scale * pw.height()), True)
        return
      h = (w * ch) / cw
    if w < 1:
      w = (h * cw) / ch

    win_size = self.size()
    delta = QtCore.QSize(w - cw, h - ch) / scale

    # window resize
    self.resize(delta + win_size)

  def get_view(self):
    self.cmd.get_view(2, quiet=0)
    QtWidgets.QApplication.clipboard().setText(self.cmd.get_view(3))
    print(" get_view: matrix copied to clipboard.")

  def show_import_popup(self) -> None:
    """Shows the dropdown menu for the scene selection."""
    # Show the dialog momentarily to ensure its size is calculated
    self.import_popup.adjustSize()
    # Get the button's position in global coordinates
    button_pos = self.ui.btn_import_seq.mapToGlobal(
      QtCore.QPoint(0, 0)
    )  # TODO: This should be replaced by the correct ribbon bar item
    # Subtract the dialog's height to position it above the button
    dialog_height = self.import_popup.height()
    adjusted_pos = button_pos - QtCore.QPoint(0, dialog_height)
    # Move the dialog to the adjusted position
    self.import_popup.move(adjusted_pos)
    self.import_popup.show()

  def toggle_fullscreen(self, toggle=-1):
    """
    Full screen
    """
    is_fullscreen = self.windowState() == Qt.WindowFullScreen
    if toggle == -1:
      toggle = not is_fullscreen

  # </editor-fold>


def commandoverloaddecorator(func):
  name = func.__name__
  func.__doc__ = getattr(pymol.cmd, name).__doc__
  setattr(pymol.cmd, name, func)
  pymol.cmd.extend(func)
  return func


window = None


def exec_app():
  """Run PySSA as a Qt application"""
  global window
  global pymol

  # don't let exceptions stop PyMOL
  import traceback

  sys.excepthook = traceback.print_exception

  # use QT_OPENGL=desktop (auto-detection may fail on Windows)
  # if hasattr(Qt, 'AA_UseDesktopOpenGL') and pymol.IS_WINDOWS:
  if pymol.IS_WINDOWS:
    print("Handling AA_UseDesktopOpenGL")
    QtCore.QCoreApplication.setAttribute(Qt.ApplicationAttribute.AA_UseDesktopOpenGL)
    # QtCore.QCoreApplication.setAttribute(Qt.AA_UseDesktopOpenGL)

  # enable 4K scaling on Windows and Linux
  if hasattr(Qt, "AA_EnableHighDpiScaling") and not any(
          v in os.environ for v in ["QT_SCALE_FACTOR", "QT_SCREEN_SCALE_FACTORS"]
  ):
    QtCore.QCoreApplication.setAttribute(Qt.AA_EnableHighDpiScaling)

  # fix Windows taskbar icon
  if pymol.IS_WINDOWS:
    import ctypes

    ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID("PySSA")

  app = QtWidgets.QApplication(sys.argv)
  app.setWindowIcon(IconManager.instance().get_icon(IconManager.Icons.LOGO))
  window = MainWindow()
  # main_window_controller.MainWindowController(window)

  @commandoverloaddecorator
  def viewport(w=-1, h=-1, _self=None):
    window.viewportsignal.emit(int(w), int(h))

  @commandoverloaddecorator
  def full_screen(toggle=-1, _self=None):
    from pymol import viewing as v

    toggle = v.toggle_dict[v.toggle_sc.auto_err(str(toggle), "toggle")]
    window.toggle_fullscreen(toggle)

  pymol.cmd._call_in_gui_thread = MainThreadCaller()

  # Assume GUI thread, make OpenGL context current before calling func().
  def _call_with_opengl_context_gui_thread(func):
    with window.pymolwidget:
      return func()

  # Dispatch to GUI thread and make OpenGL context current before calling func().
  pymol.cmd._call_with_opengl_context = lambda func: pymol.cmd._call_in_gui_thread(
    lambda: _call_with_opengl_context_gui_thread(func)
  )

  pymol.cmd.set("internal_gui", 0)
  pymol.cmd.set("internal_feedback", 0)

  window.show()
  # window.raise_()
  #
  # # window size according to -W -H options
  # options = pymol.invocation.options
  # if options.win_xy_set:
  #     scale = window.pymolwidget.fb_scale
  #     viewport(scale * options.win_x, scale * options.win_y)
  app.exec()
