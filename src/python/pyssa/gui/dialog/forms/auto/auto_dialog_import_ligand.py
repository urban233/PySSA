# Form implementation generated from reading ui file '.\src\python\pydd\components\internal\import_ligand\gui\forms\dialog_import_ligand.ui'
#
# Created by: PyQt6 UI code generator 6.4.2
#
# WARNING: Any manual changes made to this file will be lost when pyuic6 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt6 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(400, 133)
        Dialog.setMinimumSize(QtCore.QSize(400, 0))
        Dialog.setMaximumSize(QtCore.QSize(16777215, 133))
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.frame = QtWidgets.QFrame(parent=Dialog)
        self.frame.setFrameShape(QtWidgets.QFrame.Shape.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Shadow.Raised)
        self.frame.setObjectName("frame")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.frame)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.label = QtWidgets.QLabel(parent=self.frame)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.txt_import_ligand = QtWidgets.QLineEdit(parent=self.frame)
        self.txt_import_ligand.setFrame(True)
        self.txt_import_ligand.setReadOnly(False)
        self.txt_import_ligand.setObjectName("txt_import_ligand")
        self.horizontalLayout.addWidget(self.txt_import_ligand)
        self.btn_choose_ligand_file = QtWidgets.QToolButton(parent=self.frame)
        self.btn_choose_ligand_file.setObjectName("btn_choose_ligand_file")
        self.horizontalLayout.addWidget(self.btn_choose_ligand_file)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.lbl_status = QtWidgets.QLabel(parent=self.frame)
        self.lbl_status.setObjectName("lbl_status")
        self.verticalLayout.addWidget(self.lbl_status)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        self.verticalLayout_4.addWidget(self.frame)
        self.frame_bottom = QtWidgets.QFrame(parent=Dialog)
        self.frame_bottom.setFrameShape(QtWidgets.QFrame.Shape.StyledPanel)
        self.frame_bottom.setFrameShadow(QtWidgets.QFrame.Shadow.Raised)
        self.frame_bottom.setObjectName("frame_bottom")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.frame_bottom)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.btn_help = QtWidgets.QPushButton(parent=self.frame_bottom)
        self.btn_help.setObjectName("btn_help")
        self.horizontalLayout_2.addWidget(self.btn_help)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem)
        self.btn_import_ligand = QtWidgets.QPushButton(parent=self.frame_bottom)
        self.btn_import_ligand.setObjectName("btn_import_ligand")
        self.horizontalLayout_2.addWidget(self.btn_import_ligand)
        self.btn_cancel = QtWidgets.QPushButton(parent=self.frame_bottom)
        self.btn_cancel.setObjectName("btn_cancel")
        self.horizontalLayout_2.addWidget(self.btn_cancel)
        self.verticalLayout_3.addLayout(self.horizontalLayout_2)
        self.verticalLayout_4.addWidget(self.frame_bottom)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.label.setText(_translate("Dialog", "Choose an existing .sdf file"))
        self.btn_choose_ligand_file.setText(_translate("Dialog", "..."))
        self.lbl_status.setText(_translate("Dialog", "TextLabel"))
        self.btn_help.setText(_translate("Dialog", "Help"))
        self.btn_import_ligand.setText(_translate("Dialog", "Import"))
        self.btn_cancel.setText(_translate("Dialog", "Cancel"))
