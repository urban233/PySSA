# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\pyssa\gui\ui\forms\restart_user_pymol_view.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):

  def setupUi(self, Dialog):
    Dialog.setObjectName("Dialog")
    Dialog.resize(564, 204)
    self.verticalLayout_2 = QtWidgets.QVBoxLayout(Dialog)
    self.verticalLayout_2.setObjectName("verticalLayout_2")
    self.verticalLayout = QtWidgets.QVBoxLayout()
    self.verticalLayout.setObjectName("verticalLayout")
    self.lbl_logo = QtWidgets.QLabel(Dialog)
    self.lbl_logo.setObjectName("lbl_logo")
    self.verticalLayout.addWidget(self.lbl_logo)
    spacerItem = QtWidgets.QSpacerItem(
        20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
    )
    self.verticalLayout.addItem(spacerItem)
    self.horizontalLayout = QtWidgets.QHBoxLayout()
    self.horizontalLayout.setObjectName("horizontalLayout")
    self.lbl_message = QtWidgets.QLabel(Dialog)
    font = QtGui.QFont()
    font.setPointSize(14)
    self.lbl_message.setFont(font)
    self.lbl_message.setAlignment(QtCore.Qt.AlignCenter)
    self.lbl_message.setObjectName("lbl_message")
    self.horizontalLayout.addWidget(self.lbl_message)
    spacerItem1 = QtWidgets.QSpacerItem(
        40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
    )
    self.horizontalLayout.addItem(spacerItem1)
    self.btn_close_pyssa = QtWidgets.QPushButton(Dialog)
    self.btn_close_pyssa.setObjectName("btn_close_pyssa")
    self.horizontalLayout.addWidget(self.btn_close_pyssa)
    self.verticalLayout.addLayout(self.horizontalLayout)
    self.verticalLayout_2.addLayout(self.verticalLayout)

    self.retranslateUi(Dialog)
    QtCore.QMetaObject.connectSlotsByName(Dialog)

  def retranslateUi(self, Dialog):
    _translate = QtCore.QCoreApplication.translate
    Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
    self.lbl_logo.setText(_translate("Dialog", "TextLabel"))
    self.lbl_message.setText(_translate("Dialog", "TextLabel"))
    self.btn_close_pyssa.setText(_translate("Dialog", "Close PySSA"))