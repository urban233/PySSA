# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\pyssa\gui\ui\forms\dlgAddModel.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):

  def setupUi(self, Dialog):
    Dialog.setObjectName("Dialog")
    Dialog.resize(400, 133)
    Dialog.setMinimumSize(QtCore.QSize(400, 0))
    Dialog.setMaximumSize(QtCore.QSize(16777215, 133))
    self.verticalLayout_4 = QtWidgets.QVBoxLayout(Dialog)
    self.verticalLayout_4.setContentsMargins(0, 0, 0, 0)
    self.verticalLayout_4.setObjectName("verticalLayout_4")
    self.frame = QtWidgets.QFrame(Dialog)
    self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
    self.frame.setObjectName("frame")
    self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.frame)
    self.verticalLayout_2.setObjectName("verticalLayout_2")
    self.verticalLayout = QtWidgets.QVBoxLayout()
    self.verticalLayout.setObjectName("verticalLayout")
    self.label = QtWidgets.QLabel(self.frame)
    self.label.setObjectName("label")
    self.verticalLayout.addWidget(self.label)
    self.horizontalLayout = QtWidgets.QHBoxLayout()
    self.horizontalLayout.setObjectName("horizontalLayout")
    self.txt_add_protein = QtWidgets.QLineEdit(self.frame)
    self.txt_add_protein.setFrame(True)
    self.txt_add_protein.setReadOnly(False)
    self.txt_add_protein.setObjectName("txt_add_protein")
    self.horizontalLayout.addWidget(self.txt_add_protein)
    self.btn_choose_protein = QtWidgets.QToolButton(self.frame)
    self.btn_choose_protein.setObjectName("btn_choose_protein")
    self.horizontalLayout.addWidget(self.btn_choose_protein)
    self.verticalLayout.addLayout(self.horizontalLayout)
    self.lbl_status = QtWidgets.QLabel(self.frame)
    self.lbl_status.setObjectName("lbl_status")
    self.verticalLayout.addWidget(self.lbl_status)
    self.verticalLayout_2.addLayout(self.verticalLayout)
    self.verticalLayout_4.addWidget(self.frame)
    self.frame_bottom = QtWidgets.QFrame(Dialog)
    self.frame_bottom.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self.frame_bottom.setFrameShadow(QtWidgets.QFrame.Raised)
    self.frame_bottom.setObjectName("frame_bottom")
    self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.frame_bottom)
    self.verticalLayout_3.setObjectName("verticalLayout_3")
    self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_2.setObjectName("horizontalLayout_2")
    self.btn_help = QtWidgets.QPushButton(self.frame_bottom)
    self.btn_help.setObjectName("btn_help")
    self.horizontalLayout_2.addWidget(self.btn_help)
    spacerItem = QtWidgets.QSpacerItem(
        40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
    )
    self.horizontalLayout_2.addItem(spacerItem)
    self.btn_add_protein = QtWidgets.QPushButton(self.frame_bottom)
    self.btn_add_protein.setObjectName("btn_add_protein")
    self.horizontalLayout_2.addWidget(self.btn_add_protein)
    self.btn_cancel = QtWidgets.QPushButton(self.frame_bottom)
    self.btn_cancel.setObjectName("btn_cancel")
    self.horizontalLayout_2.addWidget(self.btn_cancel)
    self.verticalLayout_3.addLayout(self.horizontalLayout_2)
    self.verticalLayout_4.addWidget(self.frame_bottom)

    self.retranslateUi(Dialog)
    QtCore.QMetaObject.connectSlotsByName(Dialog)

  def retranslateUi(self, Dialog):
    _translate = QtCore.QCoreApplication.translate
    Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
    self.label.setText(
        _translate("Dialog", "Choose an existing .pdb file or enter a PDB ID")
    )
    self.btn_choose_protein.setText(_translate("Dialog", "..."))
    self.lbl_status.setText(_translate("Dialog", "TextLabel"))
    self.btn_help.setText(_translate("Dialog", "Help"))
    self.btn_add_protein.setText(_translate("Dialog", "Add"))
    self.btn_cancel.setText(_translate("Dialog", "Cancel"))
