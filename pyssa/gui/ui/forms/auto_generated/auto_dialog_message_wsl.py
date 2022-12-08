# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\dlgMessageWsl.ui'
#
# Created by: PyQt5 UI code generator 5.12.3
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(483, 202)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        spacerItem = QtWidgets.QSpacerItem(20, 68, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.lbl_message_wsl = QtWidgets.QLabel(Dialog)
        self.lbl_message_wsl.setObjectName("lbl_message_wsl")
        self.verticalLayout.addWidget(self.lbl_message_wsl)
        spacerItem1 = QtWidgets.QSpacerItem(20, 88, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem1)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        spacerItem2 = QtWidgets.QSpacerItem(17, 20, QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem2)
        self.btn_message_wsl_restart_later = QtWidgets.QPushButton(Dialog)
        self.btn_message_wsl_restart_later.setObjectName("btn_message_wsl_restart_later")
        self.horizontalLayout.addWidget(self.btn_message_wsl_restart_later)
        self.btn_message_wsl_restart = QtWidgets.QPushButton(Dialog)
        self.btn_message_wsl_restart.setObjectName("btn_message_wsl_restart")
        self.horizontalLayout.addWidget(self.btn_message_wsl_restart)
        self.btn_message_wsl_cancel = QtWidgets.QPushButton(Dialog)
        self.btn_message_wsl_cancel.setObjectName("btn_message_wsl_cancel")
        self.horizontalLayout.addWidget(self.btn_message_wsl_cancel)
        self.btn_message_wsl_ok = QtWidgets.QPushButton(Dialog)
        self.btn_message_wsl_ok.setObjectName("btn_message_wsl_ok")
        self.horizontalLayout.addWidget(self.btn_message_wsl_ok)
        self.verticalLayout_2.addLayout(self.horizontalLayout)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.lbl_message_wsl.setText(_translate("Dialog", "TextLabel"))
        self.btn_message_wsl_restart_later.setText(_translate("Dialog", "Restart later"))
        self.btn_message_wsl_restart.setText(_translate("Dialog", "Restart"))
        self.btn_message_wsl_cancel.setText(_translate("Dialog", "Cancel"))
        self.btn_message_wsl_ok.setText(_translate("Dialog", "Ok"))
