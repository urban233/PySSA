from typing import Optional

from PyQt6 import QtWidgets
from PyQt6 import QtGui

from pyssa.model.preference import model_definitions
from pyssa.model.util.gui_style import icons


class InlineMenu(QtWidgets.QWidget):
  def __init__(self):
    super().__init__()
    self.setStyleSheet("""
        QMenu {
            background-color: white;
            margin: 2px; /* some spacing around the menu */ 
        }
        QMenu::item {
            padding-top: 5px;
            padding-bottom: 5px;
            padding-left: 7px;
            padding-right: 15px;
            font-size: 13px;
        }
        QMenu::item:selected {
            background: #D6E4FD;
            border-width: 2px;
            border-radius: 4px;
            border-color: white;
        }
        QMenu::icon {
            padding-left: 15px;  /* Add left padding to the icon */
        }
        QMenu::separator {
            height: 1px;
            background: #E2E2E2;
            margin-left: 0px;
            margin-right: 0px;
        }
        QLabel {
            padding-top: 5px;
            padding-bottom: 5px;
            padding-right: 10px;
            margin-left: 10px;
            font: bold;
            font-size: 12px;
            color: #242424;
        }
    """)
    # Set up the menu
    layout = QtWidgets.QVBoxLayout(self)
    self.menu = QtWidgets.QMenu()
    layout.addWidget(self.menu)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(0)

  def add_action(self, a_name, an_icon: Optional["model_definitions.IconsEnum"] = None):
    tmp_action = QtGui.QAction(a_name, self)
    if an_icon is not None:
      tmp_action.setIcon(icons.get_icon(an_icon))
    self.menu.addAction(tmp_action)

  def add_description(self, a_description: str):
    tmp_lbl_description = QtWidgets.QLabel(a_description)
    widget_action = QtWidgets.QWidgetAction(self)
    widget_action.setDefaultWidget(tmp_lbl_description)
    self.menu.addAction(widget_action)

  def add_separator(self):
    self.menu.addSeparator()

  # def __init__(self):
  #   super().__init__()
  #   self.setStyleSheet("""
  #       QDialog {
  #           background-color: white;
  #           border: 2px solid #e0e0e0;
  #           border-radius: 4px;
  #       }
  #       QMenu {
  #           background-color: white;
  #           margin: 2px; /* some spacing around the menu */
  #       }
  #       QMenu::item {
  #           padding-top: 5px;
  #           padding-bottom: 5px;
  #           padding-left: 7px;
  #           padding-right: 15px;
  #           font-size: 13px;
  #       }
  #       QMenu::item:selected {
  #           background: #D6E4FD;
  #           border-width: 2px;
  #           border-radius: 4px;
  #           border-color: white;
  #       }
  #       QMenu::icon {
  #           padding-left: 15px;  /* Add left padding to the icon */
  #       }
  #       QMenu::separator {
  #           height: 1px;
  #           background: #E2E2E2;
  #           margin-left: 0px;
  #           margin-right: 0px;
  #       }
  #       QLabel {
  #           padding-top: 5px;
  #           padding-bottom: 5px;
  #           padding-right: 10px;
  #           margin-left: 10px;
  #           font: bold;
  #           font-size: 12px;
  #           color: #242424;
  #       }
  #   """)
  #   # Set up the menu
  #   layout = QtWidgets.QVBoxLayout(self)
  #   self.menu = QtWidgets.QMenu()
  #   # Create a non-clickable descriptive text using QWidgetAction
  #   self.descriptive_text = QtWidgets.QLabel("Import Sequence From")
  #   # self.descriptive_text.setStyleSheet("padding: 0px 10px; color: black; font: bold;")
  #   widget_action = QtWidgets.QWidgetAction(self)
  #   widget_action.setDefaultWidget(self.descriptive_text)
  #   self.menu.addAction(widget_action)
  #   # Create actions
  #   self.action_copy_paste = QtGui.QAction("Copy + Paste", self)
  #   self.action_copy_paste.setIcon(icons.get_icon(model_definitions.IconsEnum.IMPORT_SEQUENCE))
  #   self.action_this_device = QtGui.QAction("This Device", self)
  #   self.action_this_device.setIcon(icons.get_icon(model_definitions.IconsEnum.IMPORT_SEQUENCE))
  #   # Add actions to the menu
  #   self.menu.addAction(self.action_copy_paste)
  #   self.menu.addAction(self.action_this_device)
  #   layout.addWidget(self.menu)
  #   layout.setContentsMargins(0, 0, 0, 0)
  #   layout.setSpacing(0)
