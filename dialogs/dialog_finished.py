from pymol import Qt

from uiForms.auto.auto_dialog_finished import Ui_Form

class DialogFinished(Qt.QtWidgets.QDialog):
    def __init__(self, parent=None):
        """Constructor

        Args:
            args
            kwargs
        """
        Qt.QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Form()
        self.ui.setupUi(self)

        self.ui.btn_ok.clicked.connect(self.close_dialog)

        self.setWindowTitle("Analysis Finished")

    # @SLOT()
    def close_dialog(self):
        self.close()
