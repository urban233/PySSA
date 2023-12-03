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

# import os
# from PyQt5.QtCore import Qt
# from PyQt5.QtGui import QIcon
# from PyQt5 import QtWidgets
# from PyQt5.QtMultimediaWidgets import QVideoWidget
# from PyQt5.QtWidgets import QListWidgetItem, QFileDialog
#
# from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog, QVBoxLayout, QPushButton, QWidget
# from PyQt5.QtMultimedia import QMediaPlayer, QMediaContent
# from PyQt5.QtCore import QUrl
# import sys
#
# from pyssa.gui.ui.forms.auto_generated.auto_dialog_tutorial_videos import Ui_Dialog
# from pyssa.util import constants
#
# class TutorialVideosDialog(QtWidgets.QDialog):
#     def __init__(self, parent=None):
#         QtWidgets.QDialog.__init__(self, parent)
#         # build ui object
#         self.ui = Ui_Dialog()
#         self.ui.setupUi(self)
#         self.setWindowTitle("PySSA Tutorials")
#         self.setWindowIcon(QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))
#
#         # connect btn
#         self.ui.btn_tutorial_videos.clicked.connect(self.loadFile)
#
#         # fill list
#         self.populateListWithMP4Files()
#
#         # Set window flags
#         self.setWindowFlags(
#             Qt.Window |
#             Qt.WindowMaximizeButtonHint |
#             Qt.WindowCloseButtonHint
#         )
#         # Hide help button
#         self.setWindowFlags(self.windowFlags() ^ Qt.WindowContextHelpButtonHint)
#
#     def populateListWithMP4Files(self):
#         directory_path = constants.TUTORIAL_PATH
#
#         # check if path exists
#         if os.path.exists(directory_path):
#             # fill list
#             for filename in os.listdir(directory_path):
#                 if filename.endswith(".mp4"):
#                     item = QListWidgetItem(filename)
#                     self.ui.list_tutorial_videos.addItem(item)
#
#
#
#     def loadFile(self):
#         # show choosing file
#         options = QFileDialog.Options()
#         options |= QFileDialog.ReadOnly
#         file_dialog = QFileDialog(self)
#         file_dialog.setOptions(options)
#         file_dialog.setFileMode(QFileDialog.ExistingFiles)
#
#         # fill list with choosing files
#         if file_dialog.exec_():
#             file_names = file_dialog.selectedFiles()
#             for file_name in file_names:
#                 item = QListWidgetItem(file_name)
#                 self.ui.list_tutorial_videos.addItem(item)
#
#                 # Hier sollte der Code stehen, um den Video-Player zu aktualisieren
#                 # Ersetze 'YourVideoPlayerWidget' durch den tatsächlichen Namen deines Video-Player-Widgets
#                 self.ui.YourVideoPlayerWidget.setMedia(QMediaContent(QUrl.fromLocalFile(file_name)))
#                 self.ui.YourVideoPlayerWidget.play()  # Starte das Video
#
# class TutorialVideosDialog(QMainWindow):
#     def __init__(self):
#         super().__init__()
#
#         self.initUI()
#
#     def initUI(self):
#         self.setWindowTitle('Tutorial Videos')
#         self.setGeometry(100, 100, 800, 600)
#
#         layout = QVBoxLayout()
#
#         self.videoPlayer = QMediaPlayer(self)
#         self.videoWidget = QVideoWidget(self)
#         self.videoPlayer.setVideoOutput(self.videoWidget)
#
#         open_button = QPushButton('Open Video', self)
#         open_button.clicked.connect(self.openFile)
#         layout.addWidget(open_button)
#
#         layout.addWidget(self.videoWidget)
#
#         central_widget = QWidget()
#         central_widget.setLayout(layout)
#         self.setCentralWidget(central_widget)
#
#     def openFile(self):
#         options = QFileDialog.Options()
#         options |= QFileDialog.ReadOnly
#         file_dialog = QFileDialog(self)
#         file_dialog.setOptions(options)
#         file_dialog.setFileMode(QFileDialog.ExistingFile)
#
#         if file_dialog.exec_():
#             file_name = file_dialog.selectedFiles()[0]
#
#             media_content = QMediaContent(QUrl.fromLocalFile(file_name))
#             self.videoPlayer.setMedia(media_content)
#             self.videoPlayer.play()
#
# if __name__ == '__main__':
#     app = QApplication(sys.argv)
#     ex = TutorialVideosDialog()
#     ex.show()
#     sys.exit(app.exec_())

import os
import subprocess
from PyQt5.QtCore import Qt, QUrl
from PyQt5.QtGui import QIcon
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QListWidgetItem, QFileDialog, QVBoxLayout, QPushButton, QWidget, QMainWindow, QApplication
from PyQt5.QtMultimedia import QMediaPlayer, QMediaContent
from pyssa.gui.ui.forms.auto_generated.auto_dialog_tutorial_videos import Ui_Dialog
from pyssa.util import constants
import sys

class TutorialVideosDialog(QtWidgets.QDialog):
    def __init__(self, parent=None):
        QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)
        self.setGeometry(200, 200, 850, 700)
        self.setWindowTitle("PySSA Tutorials")
        self.setWindowIcon(QIcon(f"{constants.PLUGIN_ROOT_PATH}\\assets\\pyssa_logo.png"))

        # connect btn
        self.ui.btn_tutorial_videos.clicked.connect(self.loadFile)

        # fill list
        self.populateListWithMP4Files()

        # Set window flags
        self.setWindowFlags(
            Qt.Window |
            Qt.WindowMaximizeButtonHint |
            Qt.WindowCloseButtonHint
        )
        # Hide help button
        self.setWindowFlags(self.windowFlags() ^ Qt.WindowContextHelpButtonHint)

        # Video Player
        self.videoPlayer = QMediaPlayer(self)
        self.videoWidget = None

    def populateListWithMP4Files(self):
        directory_path = constants.TUTORIAL_PATH

        # check if path exists
        if os.path.exists(directory_path):
            # fill list
            for filename in os.listdir(directory_path):
                if filename.endswith(".mp4"):
                    item = QListWidgetItem(filename)
                    self.ui.list_tutorial_videos.addItem(item)

    # fixme: something is maybe wrong in this function, because the selected tutorial doesn't get show!
    def loadFile(self):
        # show choosing file
        options = QFileDialog.Options()
        options |= QFileDialog.ReadOnly
        file_dialog = QFileDialog(self)
        file_dialog.setOptions(options)
        file_dialog.setFileMode(QFileDialog.ExistingFiles)

        # fill list with choosing files
        if file_dialog.exec_():
            file_names = file_dialog.selectedFiles()
            for file_name in file_names:
                item = QListWidgetItem(file_name)
                self.ui.list_tutorial_videos.addItem(item)

                # start videoplayer from system
                try:
                    subprocess.Popen(['start', 'wmplayer', file_name], shell=True)
                except Exception as exc:
                    print(f"Fail to start videoplayer: {exc}")
