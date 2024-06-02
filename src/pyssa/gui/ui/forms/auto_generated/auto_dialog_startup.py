# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'dlgStartup.ui'
#
# Created by: PyQt5 UI code generator 5.15.7
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):

  def setupUi(self, Dialog):
    Dialog.setObjectName("Dialog")
    Dialog.resize(515, 161)
    self.verticalLayout_2 = QtWidgets.QVBoxLayout(Dialog)
    self.verticalLayout_2.setObjectName("verticalLayout_2")
    self.verticalLayout = QtWidgets.QVBoxLayout()
    self.verticalLayout.setObjectName("verticalLayout")
    self.lbl_title = QtWidgets.QLabel(Dialog)
    self.lbl_title.setObjectName("lbl_title")
    self.verticalLayout.addWidget(self.lbl_title)
    self.lbl_description = QtWidgets.QLabel(Dialog)
    self.lbl_description.setObjectName("lbl_description")
    self.verticalLayout.addWidget(self.lbl_description)
    self.verticalLayout_2.addLayout(self.verticalLayout)
    self.horizontalLayout = QtWidgets.QHBoxLayout()
    self.horizontalLayout.setObjectName("horizontalLayout")
    self.label_3 = QtWidgets.QLabel(Dialog)
    self.label_3.setObjectName("label_3")
    self.horizontalLayout.addWidget(self.label_3)
    self.txt_workspace = QtWidgets.QLineEdit(Dialog)
    self.txt_workspace.setReadOnly(True)
    self.txt_workspace.setObjectName("txt_workspace")
    self.horizontalLayout.addWidget(self.txt_workspace)
    self.btn_choose_workspace = QtWidgets.QPushButton(Dialog)
    self.btn_choose_workspace.setMinimumSize(QtCore.QSize(100, 0))
    self.btn_choose_workspace.setMaximumSize(QtCore.QSize(100, 100))
    self.btn_choose_workspace.setObjectName("btn_choose_workspace")
    self.horizontalLayout.addWidget(self.btn_choose_workspace)
    self.verticalLayout_2.addLayout(self.horizontalLayout)
    spacerItem = QtWidgets.QSpacerItem(
        20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
    )
    self.verticalLayout_2.addItem(spacerItem)
    self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_2.setObjectName("horizontalLayout_2")
    spacerItem1 = QtWidgets.QSpacerItem(
        40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
    )
    self.horizontalLayout_2.addItem(spacerItem1)
    self.btn_cancel = QtWidgets.QPushButton(Dialog)
    self.btn_cancel.setMinimumSize(QtCore.QSize(100, 0))
    self.btn_cancel.setMaximumSize(QtCore.QSize(100, 16777215))
    self.btn_cancel.setObjectName("btn_cancel")
    self.horizontalLayout_2.addWidget(self.btn_cancel)
    self.btn_launch = QtWidgets.QPushButton(Dialog)
    self.btn_launch.setMinimumSize(QtCore.QSize(100, 0))
    self.btn_launch.setMaximumSize(QtCore.QSize(100, 16777215))
    self.btn_launch.setObjectName("btn_launch")
    self.horizontalLayout_2.addWidget(self.btn_launch)
    self.verticalLayout_2.addLayout(self.horizontalLayout_2)

    self.retranslateUi(Dialog)
    QtCore.QMetaObject.connectSlotsByName(Dialog)

  def retranslateUi(self, Dialog):
    _translate = QtCore.QCoreApplication.translate
    Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
    self.lbl_title.setText(
        _translate("Dialog", "Select a directory as workspace")
    )
    self.lbl_description.setText(
        _translate(
            "Dialog",
            "The PySSA uses the workspace directory to store all projects",
        )
    )
    self.label_3.setText(_translate("Dialog", "Workspace"))
    self.btn_choose_workspace.setText(_translate("Dialog", "Browse ..."))
    self.btn_cancel.setText(_translate("Dialog", "Cancel"))
    self.btn_launch.setText(_translate("Dialog", "Launch"))
