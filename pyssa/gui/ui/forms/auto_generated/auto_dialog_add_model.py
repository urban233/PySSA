# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'dlgAddModel.ui'
#
# Created by: PyQt5 UI code generator 5.15.7
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(400, 116)
        Dialog.setMinimumSize(QtCore.QSize(400, 0))
        Dialog.setMaximumSize(QtCore.QSize(16777215, 116))
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.label = QtWidgets.QLabel(Dialog)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.txt_model = QtWidgets.QLineEdit(Dialog)
        self.txt_model.setFrame(True)
        self.txt_model.setReadOnly(True)
        self.txt_model.setObjectName("txt_model")
        self.horizontalLayout.addWidget(self.txt_model)
        self.btn_choose_model = QtWidgets.QToolButton(Dialog)
        self.btn_choose_model.setObjectName("btn_choose_model")
        self.horizontalLayout.addWidget(self.btn_choose_model)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem)
        self.btn_cancel = QtWidgets.QPushButton(Dialog)
        self.btn_cancel.setObjectName("btn_cancel")
        self.horizontalLayout_2.addWidget(self.btn_cancel)
        self.btn_add_model = QtWidgets.QPushButton(Dialog)
        self.btn_add_model.setObjectName("btn_add_model")
        self.horizontalLayout_2.addWidget(self.btn_add_model)
        self.verticalLayout_2.addLayout(self.horizontalLayout_2)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.label.setText(_translate("Dialog", "Choose a model"))
        self.btn_choose_model.setText(_translate("Dialog", "..."))
        self.btn_cancel.setText(_translate("Dialog", "Cancel"))
        self.btn_add_model.setText(_translate("Dialog", "Add model to project"))

