# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'dlgSettingsPdbPreparation.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(413, 314)
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.label = QtWidgets.QLabel(Dialog)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.txtPdbStorageDir = QtWidgets.QLineEdit(Dialog)
        self.txtPdbStorageDir.setObjectName("txt_pdb_storage_dir")
        self.horizontalLayout.addWidget(self.txtPdbStorageDir)
        self.btnPdbStorageDir = QtWidgets.QToolButton(Dialog)
        self.btnPdbStorageDir.setObjectName("btn_pdb_storage_dir")
        self.horizontalLayout.addWidget(self.btnPdbStorageDir)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.verticalLayout_3.addLayout(self.verticalLayout)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.label_2 = QtWidgets.QLabel(Dialog)
        self.label_2.setObjectName("label_2")
        self.verticalLayout_2.addWidget(self.label_2)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.txtZipStorageDir = QtWidgets.QLineEdit(Dialog)
        self.txtZipStorageDir.setObjectName("txt_zip_storage_dir")
        self.horizontalLayout_2.addWidget(self.txtZipStorageDir)
        self.btnZipStorageDir = QtWidgets.QToolButton(Dialog)
        self.btnZipStorageDir.setObjectName("btn_zip_storage_dir")
        self.horizontalLayout_2.addWidget(self.btnZipStorageDir)
        self.verticalLayout_2.addLayout(self.horizontalLayout_2)
        self.verticalLayout_3.addLayout(self.verticalLayout_2)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem)
        self.btnCancel = QtWidgets.QPushButton(Dialog)
        self.btnCancel.setObjectName("btn_cancel")
        self.horizontalLayout_3.addWidget(self.btnCancel)
        self.btnOk = QtWidgets.QPushButton(Dialog)
        self.btnOk.setObjectName("btn_ok")
        self.horizontalLayout_3.addWidget(self.btnOk)
        self.verticalLayout_3.addLayout(self.horizontalLayout_3)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.label.setText(_translate("Dialog", "Directory to store the model .pdb files:"))
        self.btnPdbStorageDir.setText(_translate("Dialog", "..."))
        self.label_2.setText(_translate("Dialog", "Directory where the prediction.zip\'s are stored:"))
        self.btnZipStorageDir.setText(_translate("Dialog", "..."))
        self.btnCancel.setText(_translate("Dialog", "Cancel"))
        self.btnOk.setText(_translate("Dialog", "OK"))

