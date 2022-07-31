from pymol import Qt

from uiForms.auto.auto_dialog_about import Ui_Dialog


class DialogAbout(Qt.QtWidgets.QDialog):

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

        self.ui.btn_ok.clicked.connect(self.close_dialog)

        self.setWindowTitle("About")

    # @SLOT
    def close_dialog(self):
        self.close()