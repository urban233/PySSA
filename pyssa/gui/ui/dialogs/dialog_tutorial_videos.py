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
""""""

import os
import pathlib
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5 import QtGui
from pyssa.gui.ui.forms.auto_generated.auto_dialog_tutorial_videos import Ui_Dialog
from pyssa.util import constants


class TutorialVideosDialog(QtWidgets.QDialog):
    def __init__(self, parent=None):
        QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        self.setGeometry(200, 200, 850, 700)
        self.setWindowTitle("PySSA Tutorials")
        self.setWindowIcon(QtGui.QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))

        # connect btn
        self.ui.btn_tutorial_videos.clicked.connect(self.load_file)

        # fill list
        self.populate_list_with_mp4_files()

        # Set window flags
        self.setWindowFlags(
            QtCore.Qt.Window |
            QtCore.Qt.WindowMaximizeButtonHint |
            QtCore.Qt.WindowCloseButtonHint
        )
        # Hide help button
        self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)

    def populate_list_with_mp4_files(self):
        # check if path exists
        if os.path.exists(constants.TUTORIAL_PATH):
            # fill list
            for filename in os.listdir(constants.TUTORIAL_PATH):
                if filename.endswith(".mp4"):
                    item = QtWidgets.QListWidgetItem(filename)
                    self.ui.list_tutorial_videos.addItem(item)

    def load_file(self):
        tmp_video_filepath = pathlib.Path(f"{constants.TUTORIAL_PATH}/{self.ui.list_tutorial_videos.currentItem().text()}")
        os.startfile(tmp_video_filepath)
