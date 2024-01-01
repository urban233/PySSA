#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2022
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/urban233/PySSA>
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
from pymol import Qt
from PyQt5 import QtCore
from pyssa.gui.ui.forms.auto_generated.auto_dialog_image import Ui_Dialog


class DialogImage(Qt.QtWidgets.QDialog):
    def __init__(self, parent=None):
        """Constructor.

        Args:
            args
            kwargs
        """
        Qt.QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)

        # image = Image.open(
        #     "/home/matt/Documents/test_pymol/Results/images/interesting_region_98.png")
        # resized_image = image.resize((500, 500))

        # pixmap = QPixmap(resized_image)
        # pixmap.scaled(64, 64)
        # self.ui.image_label.setPixmap(pixmap)
        # self.resize(pixmap.width(), pixmap.height())
