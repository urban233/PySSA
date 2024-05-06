#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/zielesny/PySSA>
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
"""Module for custom message boxes which are designed through subclassing QDialog."""
import enum

from PyQt5.QtCore import Qt
from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5 import QtWidgets

from pyssa.gui.ui.custom_widgets import custom_line_edit
from pyssa.gui.ui.styles import styles
from pyssa.util import constants


class CustomMessageBoxIcons(enum.Enum):
    """Enumeration of custom message box icons."""
    INFORMATION = ":/icons/info_w200.png"
    WARNING = ":/icons/warning_w200.png"
    ERROR = ":/icons/error_w200.png"
    DANGEROUS = ":/icons/dangerous_w200.png"


class CustomMessageBox(QtWidgets.QDialog):
    # Define a custom signal
    dialogClosed = QtCore.pyqtSignal(bool)

    def __init__(self, parent=None) -> None:  # noqa: ANN001
        """Constructor."""
        QtWidgets.QDialog.__init__(self, parent)
        self.lbl_icon = QtWidgets.QLabel("Generic Icon Pixmap")
        self.lbl_description = QtWidgets.QLabel("Generic Message")
        self.btn_left = QtWidgets.QPushButton("Generic Left")
        self.btn_right = QtWidgets.QPushButton("Generic Right")

        # <editor-fold desc="Layouts">
        self.layout_description_msg = QtWidgets.QHBoxLayout()
        self.layout_description_msg.addWidget(self.lbl_icon)
        self.layout_description_msg.addWidget(self.lbl_description)
        self.layout_description_msg.addStretch(1)

        self.layout_buttons = QtWidgets.QHBoxLayout()  # Use QHBoxLayout for the button
        self.layout_buttons.addStretch(1)  # Add stretchable space before the button
        self.layout_buttons.addWidget(self.btn_left)
        self.layout_buttons.addWidget(self.btn_right)

        self.layout_complete = QtWidgets.QVBoxLayout()
        self.layout_complete.addLayout(self.layout_description_msg)
        self.layout_complete.addLayout(self.layout_buttons)

        self.setLayout(self.layout_complete)

        # </editor-fold>

        self.setMaximumSize(600, 80)
        self.setMinimumWidth(250)
        self.resize(450, 80)

        self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
        styles.set_stylesheet(self)
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
        self.setWindowTitle("Generic Window Title")
        self.setModal(True)

    def closeEvent(self, event):
        # Emit the custom signal when the window is closed
        self.dialogClosed.emit(False)
        event.accept()


class CustomMessageBoxDelete(CustomMessageBox):
    response: bool

    def __init__(self, a_message: str, a_window_title: str, an_icon_path: str) -> None:  # noqa: ANN001
        """Constructor."""
        super().__init__()

        self.response = False

        self.lbl_icon.setText("")
        self.lbl_icon.setPixmap(QtGui.QIcon(an_icon_path).pixmap(40, 40))
        self.lbl_description.setText(a_message)
        self.btn_left.setText("Delete")
        self.btn_left.setStyleSheet("""
            QPushButton {
                background-color: #ba1a1a; 
                color: white; 
                border: none;
            }
            QPushButton::pressed {
                background-color: #410002; 
                color: white; 
                border: none;
            }
        """)
        self.btn_right.setText("Cancel")
        self.setWindowTitle(a_window_title)

        self.btn_left.clicked.connect(self.__slot_left_button)
        self.btn_right.clicked.connect(self.__slot_right_button)

    def __slot_left_button(self):
        """Method for the Delete button."""
        self.response = True
        self.close()

    def __slot_right_button(self):
        """Method for the Cancel button."""
        self.response = False
        self.close()


class CustomMessageBoxOk(CustomMessageBox):
    response: bool

    def __init__(self, a_message: str, a_window_title: str, an_icon_path: str) -> None:  # noqa: ANN001
        """Constructor."""
        super().__init__()

        self.response = False

        self.lbl_icon.setText("")
        self.lbl_icon.setPixmap(QtGui.QIcon(an_icon_path).pixmap(40, 40))
        self.lbl_description.setText(a_message)
        self.btn_left.setText("OK")
        styles.color_bottom_frame_button(self.btn_left)
        self.btn_right.hide()
        self.setWindowTitle(a_window_title)

        self.btn_left.clicked.connect(self.__slot_left_button)

    def __slot_left_button(self):
        """Method for the Delete button."""
        self.response = True
        self.close()


class CustomMessageBoxYesNo(CustomMessageBox):
    response: bool

    def __init__(self, a_message: str, a_window_title: str, an_icon_path: str) -> None:  # noqa: ANN001
        """Constructor."""
        super().__init__()

        self.response = False

        self.lbl_icon.setText("")
        self.lbl_icon.setPixmap(QtGui.QIcon(an_icon_path).pixmap(40, 40))
        self.lbl_description.setText(a_message)
        self.btn_left.setText("Yes")
        styles.color_bottom_frame_button(self.btn_left)
        self.btn_right.setText("No")
        self.setWindowTitle(a_window_title)

        self.btn_left.clicked.connect(self.__slot_left_button)
        self.btn_right.clicked.connect(self.__slot_right_button)

    def __slot_left_button(self):
        """Method for the Delete button."""
        self.response = True
        self.close()

    def __slot_right_button(self):
        """Method for the Cancel button."""
        self.response = False
        self.close()
