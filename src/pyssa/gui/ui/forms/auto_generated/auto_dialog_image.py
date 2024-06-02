# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'dlgImage.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):

  def setupUi(self, Dialog):
    Dialog.setObjectName("Dialog")
    Dialog.resize(733, 558)
    self.verticalLayout = QtWidgets.QVBoxLayout(Dialog)
    self.verticalLayout.setObjectName("verticalLayout")
    self.image_label = QtWidgets.QLabel(Dialog)
    self.image_label.setAutoFillBackground(False)
    self.image_label.setFrameShape(QtWidgets.QFrame.Box)
    self.image_label.setObjectName("image_label")
    self.verticalLayout.addWidget(self.image_label)

    self.retranslateUi(Dialog)
    QtCore.QMetaObject.connectSlotsByName(Dialog)

  def retranslateUi(self, Dialog):
    _translate = QtCore.QCoreApplication.translate
    Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
    self.image_label.setText(_translate("Dialog", "TextLabel"))
