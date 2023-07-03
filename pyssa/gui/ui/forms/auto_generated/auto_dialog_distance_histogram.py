# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\dlgDistanceHistogram.ui'
#
# Created by: PyQt5 UI code generator 5.12.3
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(565, 475)
        self.main_Layout = QtWidgets.QVBoxLayout(Dialog)
        self.main_Layout.setObjectName("main_Layout")
        self.horizontalLayout_10 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_10.setObjectName("horizontalLayout_10")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_10.addItem(spacerItem)
        self.btn_scroll_down = QtWidgets.QPushButton(Dialog)
        self.btn_scroll_down.setObjectName("btn_scroll_down")
        self.horizontalLayout_10.addWidget(self.btn_scroll_down)
        self.btn_scroll_up = QtWidgets.QPushButton(Dialog)
        self.btn_scroll_up.setObjectName("btn_scroll_up")
        self.horizontalLayout_10.addWidget(self.btn_scroll_up)
        self.main_Layout.addLayout(self.horizontalLayout_10)
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.lbl_bar_size = QtWidgets.QLabel(Dialog)
        self.lbl_bar_size.setObjectName("lbl_bar_size")
        self.horizontalLayout_7.addWidget(self.lbl_bar_size)
        self.cb_bar_size = QtWidgets.QComboBox(Dialog)
        self.cb_bar_size.setObjectName("cb_bar_size")
        self.horizontalLayout_7.addWidget(self.cb_bar_size)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_7.addItem(spacerItem1)
        self.main_Layout.addLayout(self.horizontalLayout_7)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.btn_scroll_down.setText(_translate("Dialog", "Scroll down"))
        self.btn_scroll_up.setText(_translate("Dialog", "Scroll up"))
        self.lbl_bar_size.setText(_translate("Dialog", "Bar widths:"))
