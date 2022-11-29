# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'dlgNotebookManagment.ui'
#
# Created by: PyQt5 UI code generator 5.15.7
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(484, 168)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.lbl_title = QtWidgets.QLabel(Dialog)
        self.lbl_title.setAlignment(QtCore.Qt.AlignCenter)
        self.lbl_title.setObjectName("lbl_title")
        self.verticalLayout_2.addWidget(self.lbl_title)
        spacerItem = QtWidgets.QSpacerItem(20, 5, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.MinimumExpanding)
        self.verticalLayout_2.addItem(spacerItem)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.lbl_description = QtWidgets.QLabel(Dialog)
        self.lbl_description.setWordWrap(True)
        self.lbl_description.setObjectName("lbl_description")
        self.verticalLayout.addWidget(self.lbl_description)
        self.lbl_current_status = QtWidgets.QLabel(Dialog)
        self.lbl_current_status.setWordWrap(True)
        self.lbl_current_status.setObjectName("lbl_current_status")
        self.verticalLayout.addWidget(self.lbl_current_status)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        spacerItem1 = QtWidgets.QSpacerItem(20, 5, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.MinimumExpanding)
        self.verticalLayout_2.addItem(spacerItem1)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.btn_abort = QtWidgets.QPushButton(Dialog)
        self.btn_abort.setMinimumSize(QtCore.QSize(110, 0))
        self.btn_abort.setMaximumSize(QtCore.QSize(110, 100))
        self.btn_abort.setObjectName("btn_abort")
        self.horizontalLayout.addWidget(self.btn_abort)
        self.btn_check_status = QtWidgets.QPushButton(Dialog)
        self.btn_check_status.setMinimumSize(QtCore.QSize(110, 0))
        self.btn_check_status.setMaximumSize(QtCore.QSize(110, 16777215))
        self.btn_check_status.setObjectName("btn_check_status")
        self.horizontalLayout.addWidget(self.btn_check_status)
        self.btn_reconnect = QtWidgets.QPushButton(Dialog)
        self.btn_reconnect.setMinimumSize(QtCore.QSize(110, 0))
        self.btn_reconnect.setMaximumSize(QtCore.QSize(110, 16777215))
        self.btn_reconnect.setObjectName("btn_reconnect")
        self.horizontalLayout.addWidget(self.btn_reconnect)
        self.btn_disable_gpu = QtWidgets.QPushButton(Dialog)
        self.btn_disable_gpu.setMinimumSize(QtCore.QSize(110, 0))
        self.btn_disable_gpu.setMaximumSize(QtCore.QSize(110, 16777215))
        self.btn_disable_gpu.setObjectName("btn_disable_gpu")
        self.horizontalLayout.addWidget(self.btn_disable_gpu)
        spacerItem2 = QtWidgets.QSpacerItem(5, 20, QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem2)
        self.verticalLayout_2.addLayout(self.horizontalLayout)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.lbl_title.setText(_translate("Dialog", "TextLabel"))
        self.lbl_description.setText(_translate("Dialog", "Current status:"))
        self.lbl_current_status.setText(_translate("Dialog", "TextLabel"))
        self.btn_abort.setText(_translate("Dialog", "Abort"))
        self.btn_check_status.setText(_translate("Dialog", "Check status"))
        self.btn_reconnect.setText(_translate("Dialog", "Reconnect"))
        self.btn_disable_gpu.setText(_translate("Dialog", "Disable GPU"))

