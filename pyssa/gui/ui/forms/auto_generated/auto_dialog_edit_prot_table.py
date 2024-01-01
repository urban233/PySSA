# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\dlgEditProtTable.ui'
#
# Created by: PyQt5 UI code generator 5.12.3
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(508, 265)
        sizePolicy = QtWidgets.QSizePolicy(
            QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding
        )
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Dialog.sizePolicy().hasHeightForWidth())
        Dialog.setSizePolicy(sizePolicy)
        Dialog.setMinimumSize(QtCore.QSize(0, 0))
        Dialog.setMaximumSize(QtCore.QSize(1000000, 1000000))
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setSpacing(4)
        self.verticalLayout.setObjectName("verticalLayout")
        self.lbl_edit_prot_table_protein_name = QtWidgets.QLabel(Dialog)
        self.lbl_edit_prot_table_protein_name.setObjectName("lbl_edit_prot_table_protein_name")
        self.verticalLayout.addWidget(self.lbl_edit_prot_table_protein_name)
        self.txt_edit_prot_table_protein_name = QtWidgets.QLineEdit(Dialog)
        self.txt_edit_prot_table_protein_name.setObjectName("txt_edit_prot_table_protein_name")
        self.verticalLayout.addWidget(self.txt_edit_prot_table_protein_name)
        self.lbl_edit_prot_table_status_protein_name = QtWidgets.QLabel(Dialog)
        self.lbl_edit_prot_table_status_protein_name.setObjectName("lbl_edit_prot_table_status_protein_name")
        self.verticalLayout.addWidget(self.lbl_edit_prot_table_status_protein_name)
        self.verticalLayout_4.addLayout(self.verticalLayout)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setSpacing(4)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.lbl_edit_prot_table_chain = QtWidgets.QLabel(Dialog)
        self.lbl_edit_prot_table_chain.setObjectName("lbl_edit_prot_table_chain")
        self.verticalLayout_2.addWidget(self.lbl_edit_prot_table_chain)
        self.txt_edit_prot_table_chain = QtWidgets.QLineEdit(Dialog)
        self.txt_edit_prot_table_chain.setObjectName("txt_edit_prot_table_chain")
        self.verticalLayout_2.addWidget(self.txt_edit_prot_table_chain)
        self.lbl_edit_prot_table_status_chain = QtWidgets.QLabel(Dialog)
        self.lbl_edit_prot_table_status_chain.setObjectName("lbl_edit_prot_table_status_chain")
        self.verticalLayout_2.addWidget(self.lbl_edit_prot_table_status_chain)
        self.verticalLayout_4.addLayout(self.verticalLayout_2)
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setSpacing(4)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.lbl_edit_prot_table_seq = QtWidgets.QLabel(Dialog)
        self.lbl_edit_prot_table_seq.setObjectName("lbl_edit_prot_table_seq")
        self.verticalLayout_3.addWidget(self.lbl_edit_prot_table_seq)
        self.txt_edit_prot_table_seq = QtWidgets.QLineEdit(Dialog)
        self.txt_edit_prot_table_seq.setObjectName("txt_edit_prot_table_seq")
        self.verticalLayout_3.addWidget(self.txt_edit_prot_table_seq)
        self.lbl_edit_prot_table_status_seq = QtWidgets.QLabel(Dialog)
        self.lbl_edit_prot_table_status_seq.setObjectName("lbl_edit_prot_table_status_seq")
        self.verticalLayout_3.addWidget(self.lbl_edit_prot_table_status_seq)
        self.verticalLayout_4.addLayout(self.verticalLayout_3)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_4.addItem(spacerItem)
        self.btn_edit_prot_table_cancel = QtWidgets.QPushButton(Dialog)
        self.btn_edit_prot_table_cancel.setObjectName("btn_edit_prot_table_cancel")
        self.horizontalLayout_4.addWidget(self.btn_edit_prot_table_cancel)
        self.btn_edit_prot_table_save = QtWidgets.QPushButton(Dialog)
        self.btn_edit_prot_table_save.setObjectName("btn_edit_prot_table_save")
        self.horizontalLayout_4.addWidget(self.btn_edit_prot_table_save)
        self.verticalLayout_4.addLayout(self.horizontalLayout_4)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.lbl_edit_prot_table_protein_name.setText(_translate("Dialog", "Protein name"))
        self.lbl_edit_prot_table_status_protein_name.setText(_translate("Dialog", "TextLabel"))
        self.lbl_edit_prot_table_chain.setText(_translate("Dialog", "Chain"))
        self.lbl_edit_prot_table_status_chain.setText(_translate("Dialog", "TextLabel"))
        self.lbl_edit_prot_table_seq.setText(_translate("Dialog", "Sequence"))
        self.lbl_edit_prot_table_status_seq.setText(_translate("Dialog", "TextLabel"))
        self.btn_edit_prot_table_cancel.setText(_translate("Dialog", "Cancel"))
        self.btn_edit_prot_table_save.setText(_translate("Dialog", "Save"))
