# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\dlgRenameProtein.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):

  def setupUi(self, Dialog):
    Dialog.setObjectName("Dialog")
    Dialog.resize(508, 500)
    Dialog.setMinimumSize(QtCore.QSize(400, 0))
    Dialog.setMaximumSize(QtCore.QSize(16777215, 500))
    self.verticalLayout_3 = QtWidgets.QVBoxLayout(Dialog)
    self.verticalLayout_3.setObjectName("verticalLayout_3")
    self.verticalLayout = QtWidgets.QVBoxLayout()
    self.verticalLayout.setObjectName("verticalLayout")
    self.label = QtWidgets.QLabel(Dialog)
    self.label.setObjectName("label")
    self.verticalLayout.addWidget(self.label)
    self.txt_rename_protein = QtWidgets.QLineEdit(Dialog)
    self.txt_rename_protein.setFrame(True)
    self.txt_rename_protein.setReadOnly(False)
    self.txt_rename_protein.setObjectName("txt_rename_protein")
    self.verticalLayout.addWidget(self.txt_rename_protein)
    self.lbl_status = QtWidgets.QLabel(Dialog)
    self.lbl_status.setObjectName("lbl_status")
    self.verticalLayout.addWidget(self.lbl_status)
    self.verticalLayout_3.addLayout(self.verticalLayout)
    self.verticalLayout_2 = QtWidgets.QVBoxLayout()
    self.verticalLayout_2.setObjectName("verticalLayout_2")
    self.label_2 = QtWidgets.QLabel(Dialog)
    self.label_2.setObjectName("label_2")
    self.verticalLayout_2.addWidget(self.label_2)
    self.list_workspace_proteins = QtWidgets.QListWidget(Dialog)
    self.list_workspace_proteins.setObjectName("list_workspace_proteins")
    self.verticalLayout_2.addWidget(self.list_workspace_proteins)
    self.verticalLayout_3.addLayout(self.verticalLayout_2)
    self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_2.setObjectName("horizontalLayout_2")
    spacerItem = QtWidgets.QSpacerItem(
        40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
    )
    self.horizontalLayout_2.addItem(spacerItem)
    self.btn_cancel = QtWidgets.QPushButton(Dialog)
    self.btn_cancel.setObjectName("btn_cancel")
    self.horizontalLayout_2.addWidget(self.btn_cancel)
    self.btn_rename_protein = QtWidgets.QPushButton(Dialog)
    self.btn_rename_protein.setObjectName("btn_rename_protein")
    self.horizontalLayout_2.addWidget(self.btn_rename_protein)
    self.verticalLayout_3.addLayout(self.horizontalLayout_2)

    self.retranslateUi(Dialog)
    QtCore.QMetaObject.connectSlotsByName(Dialog)

  def retranslateUi(self, Dialog):
    _translate = QtCore.QCoreApplication.translate
    Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
    self.label.setText(_translate("Dialog", "New protein name"))
    self.lbl_status.setText(_translate("Dialog", "TextLabel"))
    self.label_2.setText(_translate("Dialog", "Proteins in current workspace"))
    self.btn_cancel.setText(_translate("Dialog", "Cancel"))
    self.btn_rename_protein.setText(_translate("Dialog", "Rename"))
