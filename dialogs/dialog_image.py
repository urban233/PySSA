from PyQt5.QtGui import QPixmap
from pymol import Qt

from uiForms.auto.auto_dialog_image import Ui_Dialog
#from PIL import Image


class DialogImage(Qt.QtWidgets.QDialog):

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

        # image = Image.open(
        #     "/home/matt/Documents/test_pymol/Results/images/interesting_region_98.png")
        # resized_image = image.resize((500, 500))

        # pixmap = QPixmap(resized_image)
        # pixmap.scaled(64, 64)
        # self.ui.image_label.setPixmap(pixmap)
        # self.resize(pixmap.width(), pixmap.height())
