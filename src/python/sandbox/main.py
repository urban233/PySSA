import json
import os
import pathlib
import subprocess
import sys

from PyQt6 import QtWidgets


if __name__ == "__main__":
  script_dir = os.path.dirname(os.path.abspath(__file__))
  script_path = pathlib.Path(script_dir)
  sys.path.insert(0, str(script_path))
  sys.path.insert(1, str(script_path.parent))

  tmp_process = subprocess.Popen(r'C:\ProgramData\IBCI\PySSA\bin\PySSA\scripts\batch\start_pymol.bat')

  from sandbox.gui import simple_frame
  app = QtWidgets.QApplication(sys.argv)
  main_window = simple_frame.SimpleFrame(tmp_process)
  main_window.show()
  sys.exit(app.exec())

# from PyQt6 import QtWidgets, QtCore, QtGui
#
#
# class RibbonTab(QtWidgets.QTabWidget):
#   """Ribbon Tab with a single toolbar per tab"""
#
#   def __init__(self):
#     super().__init__()
#
#     # Home Tab
#     home_tab = QtWidgets.QWidget()
#     home_layout = QtWidgets.QVBoxLayout()
#
#     # Create the toolbar for the Home tab
#     toolbar = QtWidgets.QToolBar()
#     toolbar.setMovable(False)  # Keep the toolbar fixed in place
#
#     # Adding different sections to the toolbar
#     clipboard_section = self.create_section("Clipboard", self.create_clipboard_actions())
#     font_section = self.create_section("Font", self.create_font_actions())
#
#     toolbar.addAction(clipboard_section)
#     toolbar.addSeparator()  # Add a separator between sections
#     toolbar.addAction(font_section)
#
#     # Add the toolbar to the layout
#     home_layout.addWidget(toolbar)
#     home_tab.setLayout(home_layout)
#     self.addTab(home_tab, "Home")
#
#   def create_section(self, label_text, actions):
#     """Create a section with a group of actions and a label beneath them."""
#     section_widget = QtWidgets.QWidget()
#     section_layout = QtWidgets.QVBoxLayout()
#     section_layout.setContentsMargins(0, 0, 0, 0)  # Remove margins for compactness
#
#     # Horizontal layout for the actions (e.g., buttons)
#     action_layout = QtWidgets.QHBoxLayout()
#     for action in actions:
#       button = QtWidgets.QToolButton()
#       button.setDefaultAction(action)
#       action_layout.addWidget(button)
#
#     # Add the actions layout
#     section_layout.addLayout(action_layout)
#
#     # Add a label beneath the actions
#     label = QtWidgets.QLabel(label_text)
#     label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
#     section_layout.addWidget(label)
#
#     section_widget.setLayout(section_layout)
#
#     # Embed the section widget in a QWidgetAction to add it to the toolbar
#     section_action = QtWidgets.QWidgetAction(self)
#     section_action.setDefaultWidget(section_widget)
#
#     return section_action
#
#   def create_clipboard_actions(self):
#     """Create actions for the Clipboard section."""
#     cut_action = QtGui.QAction("Cut", self)
#     copy_action = QtGui.QAction("Copy", self)
#     paste_action = QtGui.QAction("Paste", self)
#     return [cut_action, copy_action, paste_action]
#
#   def create_font_actions(self):
#     """Create actions for the Font section."""
#     bold_action = QtGui.QAction("Bold", self)
#     italic_action = QtGui.QAction("Italic", self)
#
#     # Example custom widget (e.g., font size dropdown)
#     font_size_widget = QtWidgets.QComboBox()
#     font_size_widget.addItems(["8", "10", "12", "14", "16", "18", "20"])
#     font_size_action = QtWidgets.QWidgetAction(self)
#     font_size_action.setDefaultWidget(font_size_widget)
#
#     return [bold_action, italic_action, font_size_action]
#
#
# class MainWindow(QtWidgets.QMainWindow):
#   def __init__(self):
#     super().__init__()
#     self.setWindowTitle("Ribbon Interface Example")
#
#     ribbon = RibbonTab()
#     self.setCentralWidget(ribbon)
#
#
# if __name__ == "__main__":
#   import sys
#
#   app = QtWidgets.QApplication(sys.argv)
#   main_window = MainWindow()
#   main_window.show()
#   sys.exit(app.exec())

# from PyQt6 import QtWidgets, QtCore, QtGui
#
# class CustomTabWidget(QtWidgets.QWidget):
#   def __init__(self, icon, text):
#     super().__init__()
#     layout = QtWidgets.QVBoxLayout()
#     self.setLayout(layout)
#
#     # Create and set up the icon label
#     icon_label = QtWidgets.QLabel()
#     icon_label.setPixmap(icon.pixmap(32, 32))
#     icon_label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
#     layout.addWidget(icon_label)
#
#     # Create and set up the text label
#     text_label = QtWidgets.QLabel(text)
#     text_label.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
#     layout.addWidget(text_label)
#
#     # Optional: Add some spacing between icon and text
#     layout.setSpacing(5)
#     layout.setContentsMargins(0, 0, 0, 0)
#
# class MainWindow(QtWidgets.QMainWindow):
#   def __init__(self):
#     super().__init__()
#
#     tab_widget = QtWidgets.QTabWidget()
#     tab_widget.setTabPosition(QtWidgets.QTabWidget.TabPosition.West)  # Place tabs on the left side
#
#     # Add tabs
#     tab1 = QtWidgets.QWidget()
#     tab1_layout = QtWidgets.QVBoxLayout()
#     tab1_layout.addWidget(QtWidgets.QLabel("Content of Tab 1"))
#     tab1.setLayout(tab1_layout)
#
#     tab2 = QtWidgets.QWidget()
#     tab2_layout = QtWidgets.QVBoxLayout()
#     tab2_layout.addWidget(QtWidgets.QLabel("Content of Tab 2"))
#     tab2.setLayout(tab2_layout)
#
#     tab_widget.addTab(tab1, "")
#     tab_widget.addTab(tab2, "")
#
#     # Set custom tab widgets
#     tab_widget.setTabEnabled(0, True)
#     tab_widget.setTabEnabled(1, True)
#     tab_widget.setTabText(0, "")  # Disable default text
#     tab_widget.setTabText(1, "")  # Disable default text
#
#     # Set custom widgets as the tab labels
#     tab_widget.tabBar().setTabButton(0, QtWidgets.QTabBar.ButtonPosition.LeftSide, CustomTabWidget(QtGui.QIcon("path/to/icon1.png"), "Tab 1"))
#     tab_widget.tabBar().setTabButton(1, QtWidgets.QTabBar.ButtonPosition.LeftSide, CustomTabWidget(QtGui.QIcon("path/to/icon2.png"), "Tab 2"))
#
#     self.setCentralWidget(tab_widget)
#
# if __name__ == "__main__":
#   import sys
#   app = QtWidgets.QApplication(sys.argv)
#   main_window = MainWindow()
#   main_window.show()
#   sys.exit(app.exec())

# import sys
# import os
# from PyQt6.QtWidgets import (
#   QApplication, QMainWindow, QLineEdit, QDialog, QListWidget, QVBoxLayout, QListWidgetItem
# )
# from PyQt6.QtCore import Qt, QTimer
#
#
# class SearchDialog(QDialog):
#   def __init__(self, parent=None):
#     super().__init__(parent)
#     self.setWindowFlags(Qt.WindowType.Popup)  # Make it a pop-up dialog
#
#     self.list_widget = QListWidget()
#     layout = QVBoxLayout()
#     layout.addWidget(self.list_widget)
#     self.setLayout(layout)
#
#   def update_search_results(self, results):
#     self.list_widget.clear()
#     for result in results:
#       item = QListWidgetItem(result)
#       self.list_widget.addItem(item)
#
#   def show_below_widget(self, widget):
#     pos = widget.mapToGlobal(widget.rect().bottomLeft())
#     self.move(pos)
#     self.show()
#
#
# class MainWindow(QMainWindow):
#   def __init__(self):
#     super().__init__()
#     self.initUI()
#
#   def initUI(self):
#     self.setWindowTitle("File Search Example")
#     self.setGeometry(100, 100, 400, 100)
#
#     # Create a QLineEdit for search input
#     self.line_edit = QLineEdit(self)
#     self.line_edit.setPlaceholderText("Search files...")
#     self.line_edit.setGeometry(50, 20, 300, 30)
#     self.line_edit.textChanged.connect(self.on_text_changed)
#
#     # Create a pop-up dialog for search results
#     self.search_dialog = SearchDialog(self)
#
#     # Example directory to search
#     self.search_directory = r"C:\Users\student\.pydd\default_workspace"
#
#   def on_text_changed(self, text):
#     if text:
#       # Perform the search with a delay to avoid processing on every keystroke
#       QTimer.singleShot(300, lambda: self.search_files(text))
#     else:
#       self.search_dialog.hide()
#
#   def search_files(self, text):
#     if not text:
#       self.search_dialog.hide()
#       return
#
#     # Perform the search in the folder
#     results = []
#     for root, dirs, files in os.walk(self.search_directory):
#       for file in files:
#         if text.lower() in file.lower():
#           results.append(file)
#
#     if results:
#       self.search_dialog.update_search_results(results)
#       self.search_dialog.show_below_widget(self.line_edit)
#     else:
#       self.search_dialog.hide()
#
#
# if __name__ == '__main__':
#   app = QApplication(sys.argv)
#   window = MainWindow()
#   window.show()
#   sys.exit(app.exec())


import sys
import os
from PyQt6.QtWidgets import (
  QApplication, QMainWindow, QLineEdit, QDialog, QListWidget, QVBoxLayout, QListWidgetItem
)
from PyQt6.QtCore import Qt, QTimer


class SearchDialog(QDialog):
  def __init__(self, parent=None):
    super().__init__(parent)
    # Set the dialog to behave like a tooltip, so it doesn't steal focus
    self.setWindowFlags(Qt.WindowType.ToolTip | Qt.WindowType.FramelessWindowHint)

    self.list_widget = QListWidget()
    layout = QVBoxLayout()
    layout.addWidget(self.list_widget)
    self.setLayout(layout)

  def update_search_results(self, results):
    self.list_widget.clear()
    for result in results:
      item = QListWidgetItem(result)
      self.list_widget.addItem(item)

  def show_below_widget(self, widget):
    pos = widget.mapToGlobal(widget.rect().bottomLeft())
    self.move(pos)
    self.show()


class MainWindow(QMainWindow):
  def __init__(self):
    super().__init__()
    self.initUI()

  def initUI(self):
    self.setWindowTitle("File Search Example")
    self.setGeometry(100, 100, 400, 100)

    # Create a QLineEdit for search input
    self.line_edit = QLineEdit(self)
    self.line_edit.setPlaceholderText("Search files...")
    self.line_edit.setGeometry(50, 20, 300, 30)
    self.line_edit.textChanged.connect(self.on_text_changed)

    # Create a pop-up dialog for search results
    self.search_dialog = SearchDialog(self)

    # Example directory to search
    self.search_directory = r"C:\Users\student\.pydd\default_workspace"

  def on_text_changed(self, text):
    if text:
      # Perform the search with a delay to avoid processing on every keystroke
      QTimer.singleShot(300, lambda: self.search_files(text))
    else:
      self.search_dialog.hide()

  def search_files(self, text):
    if not text:
      self.search_dialog.hide()
      return

    # Perform the search in the folder
    results = []
    for root, dirs, files in os.walk(self.search_directory):
      for file in files:
        if text.lower() in file.lower():
          results.append(file)

    if results:
      self.search_dialog.update_search_results(results)
      self.search_dialog.show_below_widget(self.line_edit)
    else:
      self.search_dialog.hide()

    # Keep the focus on the QLineEdit
    self.line_edit.setFocus()


if __name__ == '__main__':
  app = QApplication(sys.argv)
  window = MainWindow()
  window.show()
  sys.exit(app.exec())

