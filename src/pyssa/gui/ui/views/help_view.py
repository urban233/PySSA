from PyQt5 import QtWidgets, QtGui, QtCore
from PyQt5.QtWebEngineWidgets import QWebEngineView

from src.pyssa.gui.ui.styles import styles
from src.pyssa.util import constants


class BasicBrowserView(QtWidgets.QDialog):
  def __init__(self):
    super().__init__()
    self.setWindowTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)
    self.setGeometry(100, 100, 1200, 800)

    # Layout
    layout = QtWidgets.QVBoxLayout()
    self.setLayout(layout)

    # WebEngineView
    self.browser = QWebEngineView()

    self.browser.urlChanged.connect(self.update_urlbar)
    # Navigation bar
    nav_bar = QtWidgets.QToolBar("Navigation")
    # Back button
    back_btn = QtWidgets.QAction("Back", self)
    back_btn.setStatusTip("Back to previous page")
    back_btn.triggered.connect(self.browser.back)
    nav_bar.addAction(back_btn)
    # Forward button
    forward_btn = QtWidgets.QAction("Forward", self)
    forward_btn.setStatusTip("Forward to next page")
    forward_btn.triggered.connect(self.browser.forward)
    nav_bar.addAction(forward_btn)
    # URL bar
    self.url_bar = QtWidgets.QLineEdit()
    self.url_bar.setEnabled(False)
    nav_bar.addWidget(self.url_bar)
    # Add the navigation bar to the layout
    layout.addWidget(nav_bar)
    layout.addWidget(self.browser)

    self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
    styles.set_stylesheet(self)
    self.setWindowFlags(
      self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint
    )

  def update_urlbar(self, q):
    self.url_bar.setText(q.toString())
