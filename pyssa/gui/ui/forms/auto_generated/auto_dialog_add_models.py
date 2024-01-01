# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'dlgAddModels.ui'
#
# Created by: PyQt5 UI code generator 5.15.7
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(509, 472)
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label = QtWidgets.QLabel(Dialog)
        self.label.setObjectName("label")
        self.horizontalLayout_2.addWidget(self.label)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem)
        self.btn_remove = QtWidgets.QPushButton(Dialog)
        self.btn_remove.setMinimumSize(QtCore.QSize(100, 0))
        self.btn_remove.setObjectName("btn_remove")
        self.horizontalLayout_2.addWidget(self.btn_remove)
        self.btn_add = QtWidgets.QPushButton(Dialog)
        self.btn_add.setMinimumSize(QtCore.QSize(100, 0))
        self.btn_add.setObjectName("btn_add")
        self.horizontalLayout_2.addWidget(self.btn_add)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.list_models = QtWidgets.QListWidget(Dialog)
        self.list_models.setObjectName("list_models")
        self.verticalLayout.addWidget(self.list_models)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        self.verticalLayout_3.addLayout(self.verticalLayout_2)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem1)
        self.btn_cancel = QtWidgets.QPushButton(Dialog)
        self.btn_cancel.setMinimumSize(QtCore.QSize(100, 0))
        self.btn_cancel.setObjectName("btn_cancel")
        self.horizontalLayout.addWidget(self.btn_cancel)
        self.btn_create_projects = QtWidgets.QPushButton(Dialog)
        self.btn_create_projects.setMinimumSize(QtCore.QSize(150, 0))
        self.btn_create_projects.setObjectName("btn_create_projects")
        self.horizontalLayout.addWidget(self.btn_create_projects)
        self.verticalLayout_3.addLayout(self.horizontalLayout)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.label.setText(_translate("Dialog", "Added models"))
        self.btn_remove.setText(_translate("Dialog", "Remove"))
        self.btn_add.setText(_translate("Dialog", "Add"))
        self.btn_cancel.setText(_translate("Dialog", "Cancel"))
        self.btn_create_projects.setText(_translate("Dialog", "Create projects"))
