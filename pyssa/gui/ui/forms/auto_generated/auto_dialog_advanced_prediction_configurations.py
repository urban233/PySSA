# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\dlgAdvancedPredictionConfigurations.ui'
#
# Created by: PyQt5 UI code generator 5.12.3
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(390, 147)
        self.verticalLayout = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.lbl_amber = QtWidgets.QLabel(Dialog)
        self.lbl_amber.setObjectName("lbl_amber")
        self.horizontalLayout.addWidget(self.lbl_amber)
        self.cb_amber = QtWidgets.QCheckBox(Dialog)
        self.cb_amber.setText("")
        self.cb_amber.setObjectName("cb_amber")
        self.horizontalLayout.addWidget(self.cb_amber)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.lbl_template = QtWidgets.QLabel(Dialog)
        self.lbl_template.setObjectName("lbl_template")
        self.horizontalLayout_2.addWidget(self.lbl_template)
        self.combo_box_template = QtWidgets.QComboBox(Dialog)
        self.combo_box_template.setObjectName("combo_box_template")
        self.horizontalLayout_2.addWidget(self.combo_box_template)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        spacerItem = QtWidgets.QSpacerItem(20, 29, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem1)
        self.btn_cancel = QtWidgets.QPushButton(Dialog)
        self.btn_cancel.setObjectName("btn_cancel")
        self.horizontalLayout_3.addWidget(self.btn_cancel)
        self.btn_save = QtWidgets.QPushButton(Dialog)
        self.btn_save.setObjectName("btn_save")
        self.horizontalLayout_3.addWidget(self.btn_save)
        self.verticalLayout.addLayout(self.horizontalLayout_3)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.lbl_amber.setText(_translate("Dialog", "Activate amber force field"))
        self.lbl_template.setText(_translate("Dialog", "Choose template mode"))
        self.btn_cancel.setText(_translate("Dialog", "Cancel"))
        self.btn_save.setText(_translate("Dialog", "Save"))
