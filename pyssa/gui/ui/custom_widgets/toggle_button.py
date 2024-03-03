from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5.QtCore import Qt, pyqtSignal


class ToggleButton(QtWidgets.QCheckBox):

    def __init__(
            self,
            width=40,
            bg_color="#d5d5d5",
            circle_color="#5a5a5a",
            active_circle_color="#fff",
            active_color="#367AF6"
    ):
        QtWidgets.QCheckBox.__init__(self)
        # default parameters
        self.setFixedSize(width, 19)
        self.setCursor(Qt.PointingHandCursor)

        # colors
        self._bg_color = bg_color
        self._circle_color = circle_color
        self._active_circle_color = active_circle_color
        self._active_color = active_color

    def hitButton(self, pos):
        return self.contentsRect().contains(pos)

    def paintEvent(self, e):
        # set painter
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        painter.setPen(Qt.NoPen)
        rectangle = QtCore.QRect(0, 0, self.width(), self.height())

        if not self.isChecked():
            painter.setBrush(QtGui.QColor(self._bg_color))
            painter.drawRoundedRect(0, 0, rectangle.width(), self.height(), self.height() / 2, self.height() / 2)
            painter.setBrush(QtGui.QColor(self._circle_color))
            painter.drawEllipse(3, 2, 15, 15)
        else:
            painter.setBrush(QtGui.QColor(self._active_color))
            painter.drawRoundedRect(0, 0, rectangle.width(), self.height(), self.height() / 2, self.height() / 2)
            painter.setBrush(QtGui.QColor(self._active_circle_color))
            painter.drawEllipse(self.width() - 18, 2, 15, 15)
        painter.end()


class ToggleWidget(QtWidgets.QWidget):
    """A widget that combines the toggle button with an on/off label."""

    toggleChanged = pyqtSignal(bool)
    """A custom signal that should be used instead of the stateChanged signal."""

    def __init__(self):
        super().__init__()
        self.toggle_button = ToggleButton()
        self.toggle_label = QtWidgets.QLabel("Off  ")

        self._layout = QtWidgets.QHBoxLayout()
        self._layout.addWidget(self.toggle_label)
        self._layout.addWidget(self.toggle_button)
        self.setLayout(self._layout)

        self.toggle_button.stateChanged.connect(self.switch_toggle_label_text)

    def switch_toggle_label_text(self):
        if self.toggle_button.isChecked():
            self.toggle_label.setText("On  ")
            self.toggleChanged.emit(True)
        else:
            self.toggle_label.setText("Off  ")
            self.toggleChanged.emit(False)
