# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\pyssa\gui\ui\forms\dlgTutorialVideos.ui'
#
# Created by: PyQt5 UI code generator 5.12.3
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def __init__(self):
        self.fileList = None

    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(362, 341)
        self.verticalLayout = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.list_tutorial_videos = QtWidgets.QListWidget(Dialog)
        self.list_tutorial_videos.setObjectName("list_tutorial_videos")
        self.verticalLayout.addWidget(self.list_tutorial_videos)
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.label_tutorial_videos = QtWidgets.QLabel(Dialog)
        self.label_tutorial_videos.setObjectName("label_tutorial_videos")
        self.horizontalLayout_7.addWidget(self.label_tutorial_videos)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_7.addItem(spacerItem)
        self.btn_tutorial_videos = QtWidgets.QPushButton(Dialog)
        self.btn_tutorial_videos.setObjectName("btn_tutorial_videos")
        self.horizontalLayout_7.addWidget(self.btn_tutorial_videos)
        self.verticalLayout.addLayout(self.horizontalLayout_7)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.label_tutorial_videos.setText(_translate("Dialog", "Selected Tutorial"))
        self.btn_tutorial_videos.setText(_translate("Dialog", "Open"))
