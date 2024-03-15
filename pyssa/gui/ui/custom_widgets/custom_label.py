from PyQt5 import QtWidgets
from PyQt5 import QtCore


class PermanentMessageLabel(QtWidgets.QLabel):
    textChanged = QtCore.pyqtSignal(str)

    def __init__(self, parent=None):
        super().__init__(parent)

    def setText(self, text):
        super().setText(text)
        self.textChanged.emit(text)
