from pymol import Qt

from uiForms.auto.auto_DialogWarningPredictionZip import Ui_Dialog

class DialogWarningPredictionZip(Qt.QtWidgets.QDialog):
    def __init__(self, parent=None):
        """Constructor

        Args:
            args
            kwargs
        """
        Qt.QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)

        self.setWindowTitle("Warning Prediction Exists")
