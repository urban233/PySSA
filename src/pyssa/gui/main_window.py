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

from src.pyssa.controller import main_window_controller, welcome_screen_view_controller
from src.pyssa.gui import user_pymol

from src.pyssa.gui.ui.views import pyssa_objects_panel, welcome_screen_view, help_panel

from src.pyssa.gui.ui.styles.icon_manager import IconManager
from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.qt import QtWidgets
from src.pyssa.gui.qt import QtGui

from src.pyssa.gui.ui.custom_widgets import dropdown_menu, psa_color_config, job_panel
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
    # <editor-fold desc="Create Menus">
    self.menuProject = QtWidgets.QMenu("Project", self)
    self.menuPrediction = QtWidgets.QMenu("Prediction", self)
    self.menuAnalysis = QtWidgets.QMenu("Analysis", self)
    self.menuResults = QtWidgets.QMenu("Results", self)
    self.menuImage = QtWidgets.QMenu("Image", self)
    self.menuHotspots = QtWidgets.QMenu("Hotspots", self)
    self.menuSettings = QtWidgets.QMenu("Settings", self)
    self.menuAbout = QtWidgets.QMenu("Help", self)
    # </editor-fold>

    # <editor-fold desc="Create Actions">
    # <editor-fold desc="Project Actions">
    self.action_new_project = QtGui.QAction("New", self)
    self.action_new_project.setObjectName("action_new_project")
    self.action_open_project = QtGui.QAction("Open", self)
    self.action_open_project.setObjectName("action_open_project")
    self.action_use_project = QtGui.QAction("Use", self)
    self.action_use_project.setObjectName("action_use_project")
    self.action_delete_project = QtGui.QAction("Delete", self)
    self.action_delete_project.setObjectName("action_delete_project")
    self.action_import_project = QtGui.QAction("Import", self)
    self.action_import_project.setObjectName("action_import_project")
    self.action_export_project = QtGui.QAction("Export", self)
    self.action_export_project.setObjectName("action_export_project")
    self.action_close_project = QtGui.QAction("Close", self)
    self.action_close_project.setObjectName("action_close_project")
    self.action_exit_application = QtGui.QAction("Exit Application", self)
    # </editor-fold>

    # <editor-fold desc="Prediction Actions">
    self.action_predict_monomer = QtGui.QAction("Monomer", self)
    self.action_predict_monomer.setObjectName("action_predict_monomer")
    self.action_predict_multimer = QtGui.QAction("Multimer", self)
    self.action_predict_multimer.setObjectName("action_predict_multimer")
    # </editor-fold>

    # <editor-fold desc="Analysis Actions">
    self.action_distance_analysis = QtGui.QAction("Distance", self)
    self.action_distance_analysis.setObjectName("action_distance_analysis")
    # </editor-fold>

    # <editor-fold desc="Results Actions">
    self.action_results_summary = QtGui.QAction("Summary", self)
    self.action_results_summary.setObjectName("action_results_summary")
    # </editor-fold>

    # <editor-fold desc="Image Actions">
    self.action_preview_image = QtGui.QAction("Preview", self)
    self.action_preview_image.setObjectName("action_preview_image")
    self.action_ray_tracing_image = QtGui.QAction("Ray-Tracing", self)
    self.action_ray_tracing_image.setObjectName("action_ray_tracing_image")
    self.action_simple_image = QtGui.QAction("Simple", self)
    self.action_simple_image.setObjectName("action_simple_image")
    # </editor-fold>

    # <editor-fold desc="Hotspots Actions">
    self.action_protein_regions = QtGui.QAction("Protein Regions", self)
    self.action_protein_regions.setObjectName("action_protein_regions")
    self.action_protein_regions.setCheckable(False)
    # </editor-fold>

    # <editor-fold desc="Settings Actions">
    self.action_edit_settings = QtGui.QAction("Edit", self)
    self.action_edit_settings.setObjectName("action_edit_settings")
    self.action_restore_settings = QtGui.QAction("Restore", self)
    self.action_restore_settings.setObjectName("action_restore_settings")
    # </editor-fold>

    # <editor-fold desc="About/Help Actions">
    self.action_documentation = QtGui.QAction("Documentation", self)
    self.action_documentation.setObjectName("action_documentation")
    self.action_get_demo_projects = QtGui.QAction("Get Demo Projects", self)
    self.action_get_demo_projects.setObjectName("action_get_demo_projects")
    self.action_show_log_in_explorer = QtGui.QAction("Show Logs in Explorer", self)
    self.action_show_log_in_explorer.setObjectName("action_show_log_in_explorer")
    self.action_clear_logs = QtGui.QAction("Clear Logs", self)
    self.action_clear_logs.setObjectName("action_clear_logs")
    self.action_about = QtGui.QAction("About", self)
    self.action_about.setObjectName("action_about")
    # </editor-fold>
    # </editor-fold>

    # <editor-fold desc="Set up status bar">
    self.status_bar = QtWidgets.QStatusBar()
    self.btn_update = QtWidgets.QPushButton("Update")
    self.btn_update.setStyleSheet("""
        QPushButton {
            color: #0000FF;
            text-decoration: underline;
            border: none;
            background-color: transparent;
            padding: 0px;
            margin: 0px;
        }
        QPushButton:hover {
            color: #0000CC;
        }
        QPushButton:pressed {
            color: #000088;
        }
    """)
    self.status_bar.addWidget(self.btn_update)
    self.btn_update.hide()
    self.progress_bar = QtWidgets.QProgressBar()
    self.status_bar.addWidget(self.progress_bar)
    self.progress_bar.hide()
    self.setStatusBar(self.status_bar)
    # </editor-fold>

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
      "clean": quick_access_bar_action.QuickAccessBarAction(
        "Clean", "left", 0, None, IconManager.instance().get_icon(IconManager.Icons.MOP)
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
    self.cartoon_show_action.setObjectName("popup_cartoon_show")
    self.cartoon_hide_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY_OFF), "Hide"
    )
    self.cartoon_hide_action.setObjectName("popup_cartoon_hide")
    self.cartoon_show_hide_menu.addAction(self.cartoon_show_action)
    self.cartoon_show_hide_menu.addAction(self.cartoon_hide_action)

    self.sticks_show_hide_menu = dropdown_menu.DropDownMenu()
    self.sticks_show_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY), "Show"
    )
    self.sticks_show_action.setObjectName("popup_sticks_show")
    self.sticks_hide_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY_OFF), "Hide"
    )
    self.sticks_hide_action.setObjectName("popup_sticks_hide")
    self.sticks_show_hide_menu.addAction(self.sticks_show_action)
    self.sticks_show_hide_menu.addAction(self.sticks_hide_action)

    self.ribbon_show_hide_menu = dropdown_menu.DropDownMenu()
    self.ribbon_show_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY), "Show"
    )
    self.ribbon_show_action.setObjectName("popup_ribbon_show")
    self.ribbon_hide_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY_OFF), "Hide"
    )
    self.ribbon_hide_action.setObjectName("popup_ribbon_hide")
    self.ribbon_show_hide_menu.addAction(self.ribbon_show_action)
    self.ribbon_show_hide_menu.addAction(self.ribbon_hide_action)

    self.lines_show_hide_menu = dropdown_menu.DropDownMenu()
    self.lines_show_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY), "Show"
    )
    self.lines_show_action.setObjectName("popup_lines_show")
    self.lines_hide_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY_OFF), "Hide"
    )
    self.lines_hide_action.setObjectName("popup_lines_hide")
    self.lines_show_hide_menu.addAction(self.lines_show_action)
    self.lines_show_hide_menu.addAction(self.lines_hide_action)

    self.spheres_show_hide_menu = dropdown_menu.DropDownMenu()
    self.spheres_show_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY), "Show"
    )
    self.spheres_show_action.setObjectName("popup_spheres_show")
    self.spheres_hide_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY_OFF), "Hide"
    )
    self.spheres_hide_action.setObjectName("popup_spheres_hide")
    self.spheres_show_hide_menu.addAction(self.spheres_show_action)
    self.spheres_show_hide_menu.addAction(self.spheres_hide_action)

    self.dots_show_hide_menu = dropdown_menu.DropDownMenu()
    self.dots_show_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY), "Show"
    )
    self.dots_show_action.setObjectName("popup_dots_show")
    self.dots_hide_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY_OFF), "Hide"
    )
    self.dots_hide_action.setObjectName("popup_dots_hide")
    self.dots_show_hide_menu.addAction(self.dots_show_action)
    self.dots_show_hide_menu.addAction(self.dots_hide_action)

    self.mesh_show_hide_menu = dropdown_menu.DropDownMenu()
    self.mesh_show_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY), "Show"
    )
    self.mesh_show_action.setObjectName("popup_mesh_show")
    self.mesh_hide_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY_OFF), "Hide"
    )
    self.mesh_hide_action.setObjectName("popup_mesh_hide")
    self.mesh_show_hide_menu.addAction(self.mesh_show_action)
    self.mesh_show_hide_menu.addAction(self.mesh_hide_action)
    self.surface_show_hide_menu = dropdown_menu.DropDownMenu()
    self.surface_show_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY), "Show"
    )
    self.surface_show_action.setObjectName("popup_surface_show")
    self.surface_hide_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.VISIBILITY_OFF), "Hide"
    )
    self.surface_hide_action.setObjectName("popup_surface_hide")
    self.surface_show_hide_menu.addAction(self.surface_show_action)
    self.surface_show_hide_menu.addAction(self.surface_hide_action)
    # </editor-fold>

    # <editor-fold desc="Clean">
    self.clean_solvent_organic_menu = dropdown_menu.DropDownMenu()
    self.clean_solvent_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.DELETE), "Solvent Molecules"
    )
    self.clean_solvent_action.setObjectName("popup_clean_solvent")
    self.clean_organic_action = QtGui.QAction(
      IconManager.instance().get_icon(IconManager.Icons.DELETE), "Organic Molecules"
    )
    self.clean_organic_action.setObjectName("popup_clean_organic")
    self.clean_solvent_organic_menu.addAction(self.clean_solvent_action)
    self.clean_solvent_organic_menu.addAction(self.clean_organic_action)
    # </editor-fold>

    # <editor-fold desc="Color Grid">
    self.color_grid = color_grid.PyMOLColorGrid()
    self.color_config = psa_color_config.PSAColorConfig(self.color_grid)
    self.color_config.btn_white_bg.setObjectName("popup_color_bg_white")
    self.color_config.btn_grey_bg.setObjectName("popup_color_bg_grey")
    self.color_config.btn_black_bg.setObjectName("popup_color_bg_black")
    self.color_config.btn_color_by_elements.setObjectName("popup_color_by_elements")
    self.color_grid_menu = dropdown_menu.DropDownMenu()
    self.color_grid_action = QtWidgets.QWidgetAction(None)
    # </editor-fold>

    # <editor-fold desc="Jobs">
    self.active_jobs = job_panel.JobPanel()
    self.active_jobs.setObjectName("popup_active_jobs")
    self.active_jobs_menu = dropdown_menu.DropDownMenu()
    self.active_jobs_action = QtWidgets.QWidgetAction(None)

    self.complete_jobs = job_panel.JobPanel()
    self.complete_jobs.setObjectName("popup_complete_jobs")
    self.complete_jobs_menu = dropdown_menu.DropDownMenu()
    self.complete_jobs_action = QtWidgets.QWidgetAction(None)
    # </editor-fold>

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
    self.pyssa_objects_panel = pyssa_objects_panel.PySSAObjectsPanel()
    self.pyssa_objects_panel.setObjectName("pyssa_objects_panel")
    self.help_panel = help_panel.HelpPanel()
    self.right_side_panel_stacked_widget = QtWidgets.QStackedWidget()
    self.bottom_panel_stacked_widget = QtWidgets.QStackedWidget()
    # </editor-fold>

    # # </editor-fold>
    self._setup_menu()
    self._setup_left_side_panels()
    # self._setup_right_side_panels()
    # self._setup_bottom_panels()
    self._setup_color_grid()
    self._setup_job_popups()
    self._central_widget = QtWidgets.QWidget()
    self._main_layout = QtWidgets.QVBoxLayout(self._central_widget)
    self._main_layout.setContentsMargins(4, 2, 4, 2)
    self._main_content_layout = QtWidgets.QHBoxLayout()
    self._main_content_layout.setContentsMargins(0, 0, 0, 0)
    self.tool_window_layout = tool_window_layout.ToolWindowLayout(
      None,
      None,
      list(self.viewer_toolbar_actions.values()),
      self.pymolwidget, self
    )

    self.tool_window_layout.set_right_panel_hidden(True)
    self.tool_window_layout.set_bottom_panel_hidden(True)
    # Register panels in ToolWindowLayout in enum order
    # Left stack
    self.tool_window_layout.add_left_panel(self.pyssa_objects_panel)  # LeftSidePanel.PROTEIN_STRUCTURE = 0
    self.tool_window_layout.add_right_panel(self.help_panel)

    self._main_content_layout.addWidget(self.tool_window_layout)
    self._main_layout.addLayout(self._main_content_layout)
    self.setWindowTitle("PySSA")
    self.setCentralWidget(self._central_widget)
    self._add_menu_style()
    base_style = "QMainWindow {background-color: #ebecf0;}"
    self.setStyleSheet(base_style)

  # <editor-fold desc="Private methods">
  def _setup_menu(self):
    """Builds the main menu bar, menus, and actions in the exact original order."""
    menubar = self.menuBar()

    # <editor-fold desc="Add Actions to Menus (with separators)">
    # --- Project Menu ---
    self.menuProject.addAction(self.action_new_project)
    self.menuProject.addAction(self.action_open_project)
    self.menuProject.addAction(self.action_use_project)
    self.menuProject.addAction(self.action_delete_project)
    self.menuProject.addAction(self.action_import_project)
    self.menuProject.addAction(self.action_export_project)
    self.menuProject.addSeparator()
    self.menuProject.addAction(self.action_close_project)
    self.menuProject.addSeparator()
    self.menuProject.addAction(self.action_exit_application)

    # --- Prediction Menu ---
    self.menuPrediction.addAction(self.action_predict_monomer)
    self.menuPrediction.addAction(self.action_predict_multimer)

    # --- Analysis Menu ---
    self.menuAnalysis.addAction(self.action_distance_analysis)

    # --- Results Menu ---
    self.menuResults.addAction(self.action_results_summary)

    # --- Image Menu ---
    self.menuImage.addAction(self.action_preview_image)
    self.menuImage.addAction(self.action_ray_tracing_image)
    self.menuImage.addAction(self.action_simple_image)

    # --- Hotspots Menu ---
    self.menuHotspots.addAction(self.action_protein_regions)

    # --- Settings Menu ---
    self.menuSettings.addAction(self.action_edit_settings)
    self.menuSettings.addAction(self.action_restore_settings)

    # --- About/Help Menu ---
    self.menuAbout.addAction(self.action_documentation)
    self.menuAbout.addAction(self.action_get_demo_projects)
    self.menuAbout.addSeparator()
    self.menuAbout.addAction(self.action_show_log_in_explorer)
    self.menuAbout.addAction(self.action_clear_logs)
    self.menuAbout.addSeparator()
    self.menuAbout.addAction(self.action_about)
    # </editor-fold>

    # <editor-fold desc="Add Menus to MenuBar">
    menubar.addMenu(self.menuProject)
    menubar.addMenu(self.menuPrediction)
    menubar.addMenu(self.menuAnalysis)
    menubar.addMenu(self.menuResults)
    menubar.addMenu(self.menuImage)
    menubar.addMenu(self.menuHotspots)
    menubar.addMenu(self.menuSettings)
    menubar.addMenu(self.menuAbout)
    # </editor-fold>

  def _add_menu_style(self):
    modern_light_menu_style = """
    /* Main Menu Bar */
    QMenuBar {
        background-color: #ebecf0;
        color: #000000; /* Dark gray text */
        font-family: "Segoe UI", "Helvetica Neue", sans-serif;
        font-size: 12px;
    }
    
    /* Menu Bar Items (Project, Prediction, etc.) */
    QMenuBar::item {
        background-color: transparent;
        padding: 8px 6px;
        margin: 0px 2px;
         
    }
    
    QMenuBar::item:selected {
        background-color: #f3f4f6; /* Light gray highlight */
        color: #111827; /* Nearly black text on hover */
    }
    
    QMenuBar::item:pressed {
        background-color: #e5e7eb; /* Slightly darker gray when clicked */
    }
    /* Disabled Menu Bar */
    QMenuBar:disabled {
        background-color: #ebecf0;
        color: #9ca3af; /* Gray text */
    }
    
    /* Disabled Menu Bar Items */
    QMenuBar::item:disabled {
        color: #9ca3af;              /* Medium gray text */
        background-color: transparent;
    }
    
    /* Disabled + Hover (prevents highlight when disabled) */
    QMenuBar::item:disabled:selected {
        background-color: transparent;
        color: #9ca3af;
    }
    
    /* Disabled + Pressed (safety override) */
    QMenuBar::item:disabled:pressed {
        background-color: transparent;
        color: #9ca3af;
    }
    """
    self.menuBar().setStyleSheet(modern_light_menu_style)
    tmp_menu_style = """
    /* The Dropdown Menu */
    QMenu {
        background-color: #ffffff;
        color: #374151;
        border: 1px solid #d1d5db; /* Soft gray border for shadowless definition */
        border-radius: 8px;
        padding: 0px;
        font-family: "Segoe UI", "Helvetica Neue", sans-serif;
        font-size: 12px;
    }
    
    /* Individual Dropdown Actions */
    QMenu::item {
        padding: 3px 10px 3px 24px;
        margin: 2px;
        border-radius: 5px;
        background-color: transparent;
    }
    
    QMenu::item:selected {
        background-color: #d4e2ff;
        color: #111827;
    }
    
    QMenu::item:disabled {
        color: #9ca3af; /* Faded gray for disabled items */
    }
    
    /* Horizontal Separators */
    QMenu::separator {
        height: 1px;
        background-color: #e5e7eb;
    }
    """
    self.menuProject.setStyleSheet(tmp_menu_style)
    self.menuPrediction.setStyleSheet(tmp_menu_style)
    self.menuAnalysis.setStyleSheet(tmp_menu_style)
    self.menuResults.setStyleSheet(tmp_menu_style)
    self.menuImage.setStyleSheet(tmp_menu_style)
    self.menuHotspots.setStyleSheet(tmp_menu_style)
    self.menuSettings.setStyleSheet(tmp_menu_style)
    self.menuAbout.setStyleSheet(tmp_menu_style)

  def _setup_color_grid(self) -> None:
    """Sets up the color grid on the ribbon bar."""
    self.color_grid_action.setDefaultWidget(self.color_config)
    self.color_grid_menu.addAction(self.color_grid_action)

  def _setup_job_popups(self):
    self.active_jobs_action.setDefaultWidget(self.active_jobs)
    self.active_jobs_menu.addAction(self.active_jobs_action)
    self.complete_jobs_action.setDefaultWidget(self.complete_jobs)
    self.complete_jobs_menu.addAction(self.complete_jobs_action)

  def _setup_left_side_panels(self) -> None:
    """Sets up the left side panel."""
    self.left_side_panel_stacked_widget.addWidget(self.pyssa_objects_panel)

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

  def disable_menu_bar_without_exit_application(self) -> None:
    """Disables the menu entries but not 'Exit Application'."""
    self.menuProject.setEnabled(True)
    self.action_new_project.setEnabled(False)
    self.action_open_project.setEnabled(False)
    self.action_close_project.setEnabled(False)
    self.action_use_project.setEnabled(False)
    self.action_delete_project.setEnabled(False)
    self.action_export_project.setEnabled(False)
    self.action_import_project.setEnabled(False)
    self.menuPrediction.setEnabled(False)
    self.menuAnalysis.setEnabled(False)
    self.menuResults.setEnabled(False)
    self.menuImage.setEnabled(False)
    self.menuHotspots.setEnabled(False)
    self.menuSettings.setEnabled(False)
    self.menuAbout.setEnabled(False)
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
controller = None


def exec_app():
  """Run PySSA as a Qt application"""
  global window
  global controller
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
  controller = main_window_controller.MainWindowController(window)

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
