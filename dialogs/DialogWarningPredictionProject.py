from pymol import Qt

from uiForms.auto.auto_DialogWarningPredictionProject import Ui_Dialog

class DialogWarningPredictionProject(Qt.QtWidgets.QDialog):
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

        self.ui.btnOkWarning.clicked.connect(self.closeDialog)

        self.setWindowTitle("Warning Project Exists")

    # @SLOT
    def closeDialog(self):
        self.close()