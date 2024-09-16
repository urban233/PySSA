from PyQt6 import QtCore
from PyQt6 import QtWidgets
from PyQt6 import QtGui

from pyssa.model.preference import model_definitions
from pyssa.model.util.gui_style import icons


class RibbonBar(QtWidgets.QTabWidget):
  """Container for the entire tabbed toolbar."""

  def __init__(self):
    super().__init__()
    # <editor-fold desc="Instance attributes">
    # <editor-fold desc="Home ribbon">
    # File panel
    self.import_action = QtGui.QAction("Import", self)
    self.import_action.setIcon(icons.get_icon(model_definitions.IconsEnum.IMPORT_SEQUENCE))
    self.export_action = QtGui.QAction("Export", self)
    self.export_action.setIcon(icons.get_icon(model_definitions.IconsEnum.SAVE_SEQUENCE))
    self.delete_action = QtGui.QAction("Delete", self)
    self.delete_action.setIcon(icons.get_icon(model_definitions.IconsEnum.DELETE))
    # </editor-fold>

    # <editor-fold desc="Job ribbon">
    # <editor-fold desc="Prediction panel">
    self.monomer_action = QtGui.QAction("Monomer", self)
    self.monomer_action.setIcon(icons.get_icon(model_definitions.IconsEnum.MONOMER_PROTEIN))
    self.multimer_action = QtGui.QAction("Multimer", self)
    self.multimer_action.setIcon(icons.get_icon(model_definitions.IconsEnum.MULTIMER_PROTEIN))
    # </editor-fold>
    # <editor-fold desc="Analysis panel">
    self.distance_action = QtGui.QAction("Distance", self)
    self.distance_action.setIcon(icons.get_icon(model_definitions.IconsEnum.DISTANCE_ANALYSIS))
    # </editor-fold>
    # <editor-fold desc="Docking panel">
    self.small_molecule_action = QtGui.QAction("Small Molecule", self)
    self.small_molecule_action.setIcon(icons.get_icon(model_definitions.IconsEnum.DOCKING))
    # </editor-fold>
    # <editor-fold desc="Image panel">
    self.ray_tracing_action = QtGui.QAction("Ray-Tracing", self)
    self.ray_tracing_action.setIcon(icons.get_icon(model_definitions.IconsEnum.RAY_TRACED_IMAGE))
    self.simple_action = QtGui.QAction("Simple", self)
    self.simple_action.setIcon(icons.get_icon(model_definitions.IconsEnum.IMAGE))
    # </editor-fold>
    # <editor-fold desc="Status panel">
    self.running_action = QtGui.QAction("Running", self)
    self.running_action.setIcon(icons.get_icon(model_definitions.IconsEnum.JOBS))
    self.completed_action = QtGui.QAction("Completed", self)
    self.completed_action.setIcon(icons.get_icon(model_definitions.IconsEnum.NOTIFY))
    # </editor-fold>
    # <editor-fold desc="Results panel">
    self.summary_action = QtGui.QAction("Summary", self)
    self.summary_action.setIcon(icons.get_icon(model_definitions.IconsEnum.OPEN_SESSION))
    # </editor-fold>
    # </editor-fold>

    # <editor-fold desc="Session ribbon">
    # <editor-fold desc="Session panel">
    self.open_action = QtGui.QAction("Open", self)
    self.open_action.setIcon(icons.get_icon(model_definitions.IconsEnum.OPEN_SESSION))
    self.close_action = QtGui.QAction("Close", self)
    self.close_action.setIcon(icons.get_icon(model_definitions.IconsEnum.CLOSE))
    # </editor-fold>
    # <editor-fold desc="Scene panel">
    self.create_action = QtGui.QAction("Create", self)
    self.create_action.setIcon(icons.get_icon(model_definitions.IconsEnum.NEW))
    self.save_action = QtGui.QAction("Save", self)
    self.save_action.setIcon(icons.get_icon(model_definitions.IconsEnum.SAVE_PROTEIN))
    self.delete_action = QtGui.QAction("Delete", self)
    self.delete_action.setIcon(icons.get_icon(model_definitions.IconsEnum.DELETE))
    self.rename_action = QtGui.QAction("Rename", self)
    self.rename_action.setIcon(icons.get_icon(model_definitions.IconsEnum.EDIT))
    # </editor-fold>
    # <editor-fold desc="Representation panel">
    self.sticks_action = QtGui.QAction("Sticks", self)
    self.sticks_action.setIcon(icons.get_icon(model_definitions.IconsEnum.STICKS_REPR))
    self.cartoon_action = QtGui.QAction("Cartoon", self)
    self.cartoon_action.setIcon(icons.get_icon(model_definitions.IconsEnum.CARTOON_REPR))
    self.spheres_action = QtGui.QAction("Spheres", self)
    self.spheres_action.setIcon(icons.get_icon(model_definitions.IconsEnum.SPHERES_REPR))
    self.surface_action = QtGui.QAction("Surface", self)
    self.surface_action.setIcon(icons.get_icon(model_definitions.IconsEnum.SURFACE_REPR))
    self.more_action = QtGui.QAction("More", self)
    self.more_action.setIcon(icons.get_icon(model_definitions.IconsEnum.MORE))
    # </editor-fold>
    # <editor-fold desc="Color panel">
    self.grid_action = QtGui.QAction("Grid", self)
    self.grid_action.setIcon(icons.get_icon(model_definitions.IconsEnum.COLOR_GRID))
    self.by_elements_action = QtGui.QAction("By Elements", self)
    # TODO: icon for by elements action is missing
    # </editor-fold>
    # </editor-fold>

    # <editor-fold desc="Help ribbon">
    # <editor-fold desc="Help panel">
    self.docs_action = QtGui.QAction("Documentation", self)
    self.docs_action.setIcon(icons.get_icon(model_definitions.IconsEnum.HELP))
    self.arrange_windows_action = QtGui.QAction("Arrange\nWindows", self)
    self.arrange_windows_action.setIcon(icons.get_icon(model_definitions.IconsEnum.COLOR_GRID))
    self.restart_viewer_action = QtGui.QAction("Restart\nViewer", self)
    self.restart_viewer_action.setIcon(icons.get_icon(model_definitions.IconsEnum.UPDATE_SESSION_SCENE))
    self.demo_projects_action = QtGui.QAction("Get Demo\nProjects", self)
    self.demo_projects_action.setIcon(icons.get_icon(model_definitions.IconsEnum.IMPORT))
    self.about_action = QtGui.QAction("About", self)
    self.about_action.setIcon(icons.get_icon(model_definitions.IconsEnum.MORE))
    # </editor-fold>
    # <editor-fold desc="Logs panel">
    self.open_action = QtGui.QAction("Open", self)
    self.open_action.setIcon(icons.get_icon(model_definitions.IconsEnum.OPEN))
    self.clear_action = QtGui.QAction("Clear", self)
    self.clear_action.setIcon(icons.get_icon(model_definitions.IconsEnum.DELETE))
    # </editor-fold>
    # </editor-fold>
    # </editor-fold>

    shadow_effect = QtWidgets.QGraphicsDropShadowEffect()
    shadow_effect.setBlurRadius(8)
    shadow_effect.setOffset(3, 3)
    shadow_effect.setColor(QtGui.QColor(0, 0, 0, 15))
    self.setGraphicsEffect(shadow_effect)
    self.setStyleSheet(
      """
      QTabWidget::pane { /* The tab widget frame */
          /*border: 1px solid #DCDBE3;*/
          border: 1px solid white;
          background: white;
          margin: 5px;
          border-radius: 10px;
      }
      
      QTabWidget::tab-bar {
          left: 5px; /* move to the right by 5px */
      }
      
      QTabBar::tab {
          border: none;
          background: #f0f0f0;
          min-width: 8ex;
          padding: 3px;
          font-size: 13px;
      }
      
      QTabBar::tab:selected, QTabBar::tab:hover {
          background: #e6e6e6;
          border: 2px solid #B8B8B8;
          border-top: none;
          border-left: none;
          border-right: none;
      }
      
      QTabBar::tab:selected {
          border-color: #367AF6;
          background: #e6e6e6;
          font: bold;
      }
      
      QTabBar::tab:selected#Help {
          border-color: #367AF6;
          background: red;
          font: bold;
      }
      
      QTabBar::tab:!selected {
          margin-top: 0px
      }
      
      QToolButton {
          background-color: white;
          /*min-width: 1.5em;*/
          min-width: 24px;
          min-height: 24px;
          max-height: 64px;
          padding: 0px;
          border: none;
      }
      
      QToolButton::hover {
          background: #f5f5f5;
          color: black;
      }
      """
    )
    self.setFixedHeight(135)

    self.defaults = {
      "Project": (
        RibbonPanel(
          "Placeholder",
          [
            QtGui.QAction("New", self)
          ],
          self
        ),
      ),
      "Home": (self._create_home_panels()),
      "Job": (self._create_jobs_panels()),
      "Session": (self._create_session_panels()),
      "Help": (self._create_help_panels()),
      "Protein": (self._create_protein_panels()),
    }
    for tmp_key in self.defaults.keys():
      self.add_tab(tmp_key, self.defaults[tmp_key])
    self.setCurrentIndex(1)
    self.setTabVisible(5, False)
  
  def add_tab(self, a_tab_name, panels: tuple["RibbonPanel"]):
    """Adds a tab to the tabbed ribbon widget"""
    tmp_tab = QtWidgets.QWidget()
    tmp_layout = QtWidgets.QVBoxLayout()
    tmp_layout.setContentsMargins(3, 3, 3, 3)  # Remove margins for compactness
    # Create the toolbar for the Home tab
    tmp_toolbar = QtWidgets.QToolBar()
    tmp_toolbar.setMovable(False)  # Keep the toolbar fixed in place
    # Adding different sections to the toolbar
    for tmp_panel in panels:
      if tmp_panel.is_separator():
        tmp_toolbar.addSeparator()  # Add a separator between sections
      else:
        if tmp_panel.panel_layout == "stacked":
          tmp_toolbar.addAction(tmp_panel.create_stacked_panel())
        else:
          tmp_toolbar.addAction(tmp_panel.create_panel())
    # Add the toolbar to the layout
    tmp_layout.addWidget(tmp_toolbar)
    tmp_tab.setLayout(tmp_layout)
    self.addTab(tmp_tab, a_tab_name)

  def _create_home_panels(self):
    tmp_panels = (
      RibbonPanel(
        "File",
        [
          self.import_action,
          self.export_action,
          self.delete_action
        ],
        self
      ),
      RibbonPanel("", [], self, "separator")
    )
    return tmp_panels

  def _create_jobs_panels(self):
    tmp_panels = (
      RibbonPanel(
        "Prediction",
        [self.monomer_action, self.multimer_action],
        self, has_settings=True
      ),
      RibbonPanel("", [], self, "separator"),
      RibbonPanel(
        "Analysis",
        [self.distance_action],
        self, has_settings=True
      ),
      RibbonPanel("", [], self, "separator"),
      RibbonPanel(
        "Docking",
        [self.small_molecule_action],
        self, has_settings=True
      ),
      RibbonPanel("", [], self, "separator"),
      RibbonPanel(
        "Image",
        [self.ray_tracing_action, self.simple_action],
        self, panel_layout="stacked", has_settings=True
      ),
      RibbonPanel("", [], self, "separator"),
      RibbonPanel(
        "Status",
        [self.running_action, self.completed_action],
        self
      ),
      RibbonPanel("", [], self, "separator"),
      RibbonPanel(
        "Results",
        [self.summary_action],
        self
      ),
      RibbonPanel("", [], self, "separator")
    )
    return tmp_panels

  def _create_session_panels(self):
    tmp_panels = (
      RibbonPanel(
        "Session",
        [
          self.open_action,
          self.close_action
        ],
        self, panel_layout="stacked"
      ),
      RibbonPanel("", [], self, "separator"),
      RibbonPanel(
        "Scene",
        [
          self.create_action,
          self.save_action,
          self.delete_action,
          #self.rename_action
        ],
        self, panel_layout="stacked"
      ),
      RibbonPanel("", [], self, "separator"),
      RibbonPanel(
        "Representation",
        [
          self.sticks_action,
          self.cartoon_action,
          self.spheres_action,
          self.surface_action,
          self.more_action,
        ],
        self
      ),
      RibbonPanel("", [], self, "separator"),
      RibbonPanel(
        "Color",
        [
          self.grid_action,
          self.by_elements_action
        ],
        self
      ),
      RibbonPanel("", [], self, "separator"),
      # RibbonPanel(
      #   "Sticks",
      #   [
      #     QtGui.QAction("Show", self),
      #     QtGui.QAction("Hide", self),
      #   ],
      #   self, panel_layout="stacked"
      # ),
      # RibbonPanel("", [], self, "separator"),
      # RibbonPanel(
      #   "Cartoon",
      #   [
      #     QtGui.QAction("Show", self),
      #     QtGui.QAction("Hide", self),
      #   ],
      #   self, panel_layout="stacked"
      # ),
      # RibbonPanel("", [], self, "separator"),
      # RibbonPanel(
      #   "Ribbon",
      #   [
      #     QtGui.QAction("Show", self),
      #     QtGui.QAction("Hide", self),
      #   ],
      #   self, panel_layout="stacked"
      # ),
      # RibbonPanel("", [], self, "separator"),
      # RibbonPanel(
      #   "Lines",
      #   [
      #     QtGui.QAction("Show", self),
      #     QtGui.QAction("Hide", self),
      #   ],
      #   self, panel_layout="stacked"
      # ),
      # RibbonPanel("", [], self, "separator"),
      # RibbonPanel(
      #   "Spheres",
      #   [
      #     QtGui.QAction("Show", self),
      #     QtGui.QAction("Hide", self),
      #   ],
      #   self, panel_layout="stacked"
      # ),
      # RibbonPanel("", [], self, "separator"),
      # RibbonPanel(
      #   "Dots",
      #   [
      #     QtGui.QAction("Show", self),
      #     QtGui.QAction("Hide", self),
      #   ],
      #   self, panel_layout="stacked"
      # ),
      # RibbonPanel("", [], self, "separator"),
      # RibbonPanel(
      #   "Mesh",
      #   [
      #     QtGui.QAction("Show", self),
      #     QtGui.QAction("Hide", self),
      #   ],
      #   self, panel_layout="stacked"
      # ),
      # RibbonPanel("", [], self, "separator"),
      # RibbonPanel(
      #   "Surface",
      #   [
      #     QtGui.QAction("Show", self),
      #     QtGui.QAction("Hide", self),
      #   ],
      #   self, panel_layout="stacked"
      # ),
      # RibbonPanel("", [], self, "separator"),
    )
    return tmp_panels

  def _create_help_panels(self):
    tmp_panels = (
      RibbonPanel(
        "Help",
        [
          self.docs_action,
          self.arrange_windows_action,
          self.restart_viewer_action,
          self.demo_projects_action,
          self.about_action
        ],
        self
      ),
      RibbonPanel("", [], self, "separator"),
      RibbonPanel(
        "Logs",
        [
          self.open_action,
          self.clear_action
        ],
        self, panel_layout="stacked"
      ),
      RibbonPanel("", [], self, "separator")
    )
    return tmp_panels

  def _create_protein_panels(self):
    tmp_panels = (
      RibbonPanel(
        "Color",
        [
          self.grid_action,
          self.by_elements_action,
        ],
        self
      ),
      RibbonPanel("", [], self, "separator"),
      RibbonPanel(
        "Representation",
        [
          self.sticks_action,
          self.cartoon_action,
          self.spheres_action,
          self.surface_action,
          self.more_action,
        ],
        self
      ),
      RibbonPanel("", [], self, "separator")
    )
    return tmp_panels

class RibbonPanel(QtWidgets.QWidgetAction):

  def __init__(self, 
               a_panel_name: str, 
               actions: list, 
               a_parent: QtCore.QObject,
               panel_type: str = "content",
               panel_layout: str = "side-by-side",
               has_settings: bool = False):
    super().__init__(a_parent)
    self.parent = a_parent
    self.panel_name = a_panel_name
    self.panel_type = panel_type  # content or separator
    self.panel_layout = panel_layout  # "side-by-side" or stacked
    self.has_settings: bool = has_settings
    self.actions = actions
    
  def is_separator(self) -> bool:
    return True if self.panel_type == "separator" else False

  def create_panel(self) -> QtWidgets.QWidgetAction:
    """Create a section with a group of actions and a label beneath them."""
    tmp_panel_widget = QtWidgets.QWidget()
    tmp_panel_layout = QtWidgets.QVBoxLayout()
    tmp_panel_layout.setContentsMargins(0, 0, 0, 0)  # Remove margins for compactness
    # Horizontal layout for the actions (e.g., buttons)
    tmp_action_layout = QtWidgets.QHBoxLayout()
    for action in self.actions:
      button = QtWidgets.QToolButton()
      button.setDefaultAction(action)
      button.setToolButtonStyle(QtCore.Qt.ToolButtonStyle.ToolButtonTextUnderIcon)
      button.setIconSize(QtCore.QSize(36, 36))
      button.setSizePolicy(QtWidgets.QSizePolicy.Policy.Minimum, QtWidgets.QSizePolicy.Policy.Fixed)
      button.setFixedHeight(115)
      tmp_action_layout.addWidget(button)
    # Add the actions layout
    tmp_panel_layout.addLayout(tmp_action_layout)
    # # Add a label beneath the actions
    # tmp_bottom_layout = QtWidgets.QHBoxLayout()
    # label = QtWidgets.QLabel(self.panel_name)
    # label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
    # tmp_bottom_layout.addWidget(label)
    # if self.has_settings:
    #   settings_button = QtWidgets.QPushButton()
    #   icons.set_icon(
    #     settings_button,
    #     model_definitions.IconsEnum.RIBBON_PANEL_SETTINGS,
    #     size=QtCore.QSize(12, 12)
    #   )
    #   settings_button.setStyleSheet(
    #     """
    #     QPushButton {
    #         background-color: white;
    #         border: none;
    #         padding: 2px;
    #         min-width: 12px;
    #         max-width: 12px;
    #         min-height: 12px;
    #         max-height: 12px
    #     }
    #     QPushButton::hover {
    #         background: rgba(220, 219, 227, 0.5);
    #         color: white;
    #         color: #4B91F7;
    #         border: 2px solid #DCDBE3;
    #         border: none;
    #     }
    #     """
    #   )
    #   tmp_bottom_layout.addWidget(settings_button)
    tmp_panel_layout.addLayout(self._add_panel_name(tmp_panel_layout))
    #tmp_panel_layout.addLayout(tmp_bottom_layout)
    tmp_panel_widget.setLayout(tmp_panel_layout)
    # Embed the section widget in a QWidgetAction to add it to the toolbar
    tmp_panel_action = QtWidgets.QWidgetAction(self.parent)
    tmp_panel_action.setDefaultWidget(tmp_panel_widget)
    return tmp_panel_action

  def create_stacked_panel(self) -> QtWidgets.QWidgetAction:
    """Create a section with a group of actions and a label beneath them."""
    tmp_panel_widget = QtWidgets.QWidget()
    tmp_panel_layout = QtWidgets.QVBoxLayout()
    tmp_panel_layout.setContentsMargins(2, 0, 2, 0)  # Remove margins for compactness

    # Grid layout for the actions (e.g., buttons)
    tmp_action_layout = QtWidgets.QGridLayout()
    tmp_action_layout.setSpacing(1)  # Adjust spacing between buttons if needed
    tmp_action_layout.setContentsMargins(0, 0, 0, 0)  # Remove margins for compactness

    # Determine the number of columns needed for a maximum of 2 rows
    num_actions = len(self.actions)
    num_rows = 3
    num_columns = (num_actions + num_rows - 1) // num_rows  # Ceiling division to calculate columns

    row = 0
    column = 0
    for action in self.actions:
      button = QtWidgets.QToolButton()
      button.setDefaultAction(action)
      button.setToolButtonStyle(QtCore.Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
      button.setIconSize(QtCore.QSize(20, 20))
      button.setSizePolicy(QtWidgets.QSizePolicy.Policy.Minimum, QtWidgets.QSizePolicy.Policy.Minimum)

      tmp_action_layout.addWidget(button, row, column)
      column += 1

      # Move to the next row if we reach the maximum number of columns
      if column >= num_columns:
        column = 0
        row += 1
      # Stop adding widgets if we reach the maximum number of rows
      if row >= num_rows:
        break
    # Add the actions layout
    tmp_panel_layout.addLayout(tmp_action_layout)
    tmp_panel_layout.addLayout(self._add_panel_name(tmp_panel_layout))
    tmp_panel_widget.setLayout(tmp_panel_layout)
    # Embed the section widget in a QWidgetAction to add it to the toolbar
    tmp_panel_action = QtWidgets.QWidgetAction(self.parent)
    tmp_panel_action.setDefaultWidget(tmp_panel_widget)
    return tmp_panel_action

  def _add_panel_name(self, a_panel_layout):
    """Adds the panel name to the ribbon panel."""
    # Add a vertical spacer to push the bottom layout down
    spacer_item = QtWidgets.QSpacerItem(20, 1, QtWidgets.QSizePolicy.Policy.Minimum,
                                        QtWidgets.QSizePolicy.Policy.Expanding)
    a_panel_layout.addItem(spacer_item)
    # Add a label beneath the actions
    tmp_bottom_layout = QtWidgets.QHBoxLayout()
    label = QtWidgets.QLabel(self.panel_name)
    label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
    tmp_bottom_layout.addWidget(label)
    if self.has_settings:
      settings_button = QtWidgets.QPushButton()
      icons.set_icon(
        settings_button,
        model_definitions.IconsEnum.RIBBON_PANEL_SETTINGS,
        size=QtCore.QSize(12, 12)
      )
      settings_button.setStyleSheet(
        """
        QPushButton {
            background-color: white;
            border: none;
            padding: 2px;
            min-width: 12px;
            max-width: 12px;
            min-height: 12px;
            max-height: 12px
        }
        QPushButton::hover {
            background: rgba(220, 219, 227, 0.5);
            color: white;
            color: #4B91F7;
            border: 2px solid #DCDBE3;
            border: none;
        }
        """
      )
      tmp_bottom_layout.addWidget(settings_button)
    return tmp_bottom_layout
