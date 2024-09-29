from typing import Callable, Optional
from PyQt6 import QtWidgets, QtGui, QtCore


class BigCardButton(QtWidgets.QPushButton):
  """A custom button in the style of a big card."""

  def __init__(self, a_description):
    """Constructor"""
    super().__init__()
    self.setStyleSheet(
      """
      QPushButton {
          background-color: #f0f0f0;
          color: black;
          font-family: "Segoe UI";
          font-size: 12px;
          border: none;
          padding: 2px;
      }
      QPushButton::hover {
          background-color: rgba(220, 219, 227, 0.5);
      }
      QPushButton::pressed {
          background: #d6d6d6;
          border: solid;
          border-width: 2px;
          border-radius: 4px;
          border-color: #367AF6;
      }
    """
    )
    # Set up the main layout directly on the button
    self._layout = QtWidgets.QVBoxLayout(self)
    self._layout.setContentsMargins(0, 0, 0, 0)
    self._layout.setSpacing(0)

    # Set the button to have no default padding
    self.setContentsMargins(0, 0, 0, 0)

    # Create and set up the QLabel for image
    self.lbl_image = QtWidgets.QLabel(self)
    pixmap = QtGui.QPixmap(r"C:\Users\student\Downloads\MyRastergrafik.png")  # Update with the correct path
    self.lbl_image.setPixmap(pixmap.scaled(175, 175, QtCore.Qt.AspectRatioMode.KeepAspectRatio))

    # Create and set up the QLabel for description
    self.lbl_description = QtWidgets.QLabel(a_description, self)

    # Add widgets to layout and center them
    self._layout.addWidget(self.lbl_image, alignment=QtCore.Qt.AlignmentFlag.AlignHCenter)
    self._layout.addWidget(self.lbl_description, alignment=QtCore.Qt.AlignmentFlag.AlignHCenter)

    # Ensure the layout size is updated
    self.setLayout(self._layout)
    self.setFixedSize(150, 200)


class ProjectOverviewButton(QtWidgets.QPushButton):
  """A custom button in the style of a big card."""

  button_clicked = QtCore.pyqtSignal(tuple)
  """A signal used to transfer the project name."""

  def __init__(self, a_project_name, a_date_modified, a_callable: Optional[Callable] = None):
    """Constructor

    Args:
      a_project_name: The name of the project
      a_date_modified: The date the project was last modified
      a_callable: The function to use after the button was clicked.

    Notes:
      IMPORTANT: If the button is clicked the `a_callable` function will
      receive a tuple containing the project name and the date modified.
    """
    super().__init__()
    self.setStyleSheet(
      """
      QPushButton {
          background-color: #f5f5f5;
          color: black;
          font-family: "Segoe UI";
          font-size: 12px;
          border-top: 1px solid #d1d1d1;
          border-bottom: 1px solid #d1d1d1;
          padding: 15px;
      }
      QPushButton::hover {
          background-color: rgba(220, 219, 227, 0.5);
          border: none;
      }
      QPushButton:disabled {
          background-color: #f5f5f5;
          color: #B0B0B0;
          font-family: "Segoe UI";
          font-size: 12px;
          border: solid;
          border-width: 1px;
          border-radius: 4px;
          border-color: #DCDCDC;
          padding: 2px;
      }

      QPushButton::pressed {
          background: #d6d6d6;
          color: black;
          font-family: "Segoe UI";
          font-size: 12px;
          border: solid;
          border-width: 2px;
          border-radius: 4px;
          border-color: #367AF6;
          padding: 2px;
      }
      """
    )
    # Set up the main layout directly on the button
    self._layout = QtWidgets.QHBoxLayout(self)
    self._layout.setContentsMargins(0, 0, 0, 0)
    self._layout.setSpacing(0)
    # Set the button to have no default padding
    self.setContentsMargins(0, 0, 0, 0)
    self.lbl_project_name = QtWidgets.QLabel(a_project_name, self)
    self.lbl_project_name.setStyleSheet("""margin-left: 10px;""")
    self.lbl_date_modified = QtWidgets.QLabel(a_date_modified, self)
    self.lbl_date_modified.setStyleSheet("""margin-right: 10px;""")
    # Add widgets to layout and center them
    self._layout.addWidget(self.lbl_project_name)
    self._layout.addStretch(1)
    self._layout.addWidget(self.lbl_date_modified)
    # Ensure the layout size is updated
    self.setLayout(self._layout)
    # Connect clicked signal
    if a_callable is not None:
      self.clicked.connect(self.send_button_clicked)
      self.button_clicked.connect(a_callable)

  def send_button_clicked(self) -> None:
    """Sends a button clicked signal with the project name and date modified."""
    self.button_clicked.emit((self.lbl_project_name.text(), self.lbl_date_modified.text()))


class ProjectOverviewButtonWithCheckbox(ProjectOverviewButton):
  def __init__(self, a_project_name, a_date_modified, a_callable: Optional[Callable] = None):
    super().__init__(a_project_name, a_date_modified, a_callable)
    self.cb_select = QtWidgets.QCheckBox("")
    self.cb_select.setStyleSheet("""margin-left: 10px""")
    self._layout.insertWidget(0, self.cb_select)
    self.lbl_project_name.setStyleSheet("""margin-left: 43px;""")
    if a_callable is None:
      self.connect_button_clicked_with_checkbox_selection()
    else:
      self.clicked.connect(a_callable)

  def connect_button_clicked_with_checkbox_selection(self) -> None:
    """Connects the button clicked signal with the selection of the checkbox."""
    self.clicked.connect(self._switch_checkbox_check_state)

  def _switch_checkbox_check_state(self):
    """Checks or unchecks the checkbox based on the current state."""
    if self.cb_select.isChecked():
      self.cb_select.setChecked(False)
    else:
      self.cb_select.setChecked(True)

  def is_selected(self) -> bool:
    """Checks if the checkbox is currently selected."""
    return self.cb_select.isChecked()


class ProjectSearchOverviewButton(QtWidgets.QPushButton):
  """A custom button in the style of a big card."""

  button_clicked = QtCore.pyqtSignal(tuple)
  """A signal used to transfer the project name."""

  def __init__(self, a_project_name, a_date_modified, a_callable: Optional[Callable] = None):
    """Constructor

    Args:
      a_project_name: The name of the project
      a_date_modified: The date the project was last modified
      a_callable: The function to use after the button was clicked.

    Notes:
      IMPORTANT: If the button is clicked the `a_callable` function will
      receive a tuple containing the project name and the date modified.
    """
    super().__init__()
    self.setStyleSheet(
      """
      QPushButton {
          background-color: #f5f5f5;
          color: black;
          font-family: "Segoe UI";
          font-size: 12px;
          border-top: 1px solid #d1d1d1;
          border-bottom: 1px solid #d1d1d1;
          padding: 15px;
      }
      QPushButton::hover {
          background-color: rgba(220, 219, 227, 0.5);
          border: none;
      }
      QPushButton:disabled {
          background-color: #f5f5f5;
          color: #B0B0B0;
          font-family: "Segoe UI";
          font-size: 12px;
          border: solid;
          border-width: 1px;
          border-radius: 4px;
          border-color: #DCDCDC;
          padding: 2px;
      }

      QPushButton::pressed {
          background: #d6d6d6;
          color: black;
          font-family: "Segoe UI";
          font-size: 12px;
          border: solid;
          border-width: 2px;
          border-radius: 4px;
          border-color: #367AF6;
          padding: 2px;
      }
      """
    )
    # Set up the main layout directly on the button
    self._layout = QtWidgets.QVBoxLayout(self)
    self._layout.setContentsMargins(0, 0, 0, 0)
    self._layout.setSpacing(0)
    # Set the button to have no default padding
    self.setContentsMargins(0, 0, 0, 0)
    self.lbl_project_name = QtWidgets.QLabel(a_project_name, self)
    self.lbl_project_name.setStyleSheet("""margin-left: 10px;""")
    self.lbl_date_modified = QtWidgets.QLabel(a_date_modified, self)
    self.lbl_date_modified.setStyleSheet("""margin-left: 10px;""")
    self._layout.addWidget(self.lbl_project_name)
    self._layout.addWidget(self.lbl_date_modified)
    # Ensure the layout size is updated
    self.setLayout(self._layout)
    # Connect clicked signal
    if a_callable is not None:
      self.clicked.connect(self.send_button_clicked)
      self.button_clicked.connect(a_callable)

  def send_button_clicked(self) -> None:
    """Sends a button clicked signal with the project name and date modified."""
    self.button_clicked.emit((self.lbl_project_name.text(), self.lbl_date_modified.text()))