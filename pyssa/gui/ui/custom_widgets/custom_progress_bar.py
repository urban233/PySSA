from PyQt5 import QtWidgets
from PyQt5 import QtCore


class LoadingProgressBar(QtWidgets.QProgressBar):

    def __init__(self):
        super().__init__()
        self._setup_ui()

    def _setup_ui(self):
        self.setValue(0)
        self.setTextVisible(False)
        self._animation = QtCore.QPropertyAnimation(self, b'loading')
        self._animation.setStartValue(self.minimum())
        self._animation.setEndValue(self.maximum())
        self._animation.valueChanged.connect(self.__loading)
        self._animation.setDuration(600)
        self._animation.setEasingCurve(QtCore.QEasingCurve.Linear)
        self.setStyleSheet("""
            QProgressBar {
                 border-style: solid;
                 border-width: 2px;
                 border-radius: 4px;
                 border-color: #367af6;
                 background-color: #367af6;
                 max-height: 5px;
                 max-width: 100px;
            }
            QProgressBar::chunk {
                background-color: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 transparent, stop: 0.5 #CCCCCC, stop: 0.6 #CCCCCC, stop:1 transparent);
            }
        """)
        self._animation.start()

    def __loading(self, v):
        self.setValue(v)
        if self._animation.currentValue() == self._animation.endValue():
            self._animation.setDirection(QtCore.QAbstractAnimation.Backward)
            self.setInvertedAppearance(True)
            self._animation.start()
        elif self._animation.currentValue() == self._animation.startValue():
            self._animation.setDirection(QtCore.QAbstractAnimation.Forward)
            self.setInvertedAppearance(False)
            self._animation.start()
