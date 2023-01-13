# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\dlgSettingsGlobal.ui'
#
# Created by: PyQt5 UI code generator 5.15.7
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(383, 431)
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.label_3 = QtWidgets.QLabel(Dialog)
        self.label_3.setObjectName("label_3")
        self.verticalLayout_3.addWidget(self.label_3)
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.txt_workspace_dir = QtWidgets.QLineEdit(Dialog)
        self.txt_workspace_dir.setObjectName("txt_workspace_dir")
        self.horizontalLayout_4.addWidget(self.txt_workspace_dir)
        self.btn_workspace_dir = QtWidgets.QToolButton(Dialog)
        self.btn_workspace_dir.setObjectName("btn_workspace_dir")
        self.horizontalLayout_4.addWidget(self.btn_workspace_dir)
        self.verticalLayout_3.addLayout(self.horizontalLayout_4)
        self.verticalLayout_5.addLayout(self.verticalLayout_3)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.label_2 = QtWidgets.QLabel(Dialog)
        self.label_2.setObjectName("label_2")
        self.verticalLayout_2.addWidget(self.label_2)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.txt_zip_storage_dir = QtWidgets.QLineEdit(Dialog)
        self.txt_zip_storage_dir.setObjectName("txt_zip_storage_dir")
        self.horizontalLayout_2.addWidget(self.txt_zip_storage_dir)
        self.btn_zip_storage_dir = QtWidgets.QToolButton(Dialog)
        self.btn_zip_storage_dir.setObjectName("btn_zip_storage_dir")
        self.horizontalLayout_2.addWidget(self.btn_zip_storage_dir)
        self.verticalLayout_2.addLayout(self.horizontalLayout_2)
        self.verticalLayout_5.addLayout(self.verticalLayout_2)
        self.verticalLayout_4 = QtWidgets.QVBoxLayout()
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.line = QtWidgets.QFrame(Dialog)
        self.line.setFrameShape(QtWidgets.QFrame.HLine)
        self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line.setObjectName("line")
        self.verticalLayout_4.addWidget(self.line)
        self.label_6 = QtWidgets.QLabel(Dialog)
        self.label_6.setObjectName("label_6")
        self.verticalLayout_4.addWidget(self.label_6)
        self.verticalLayout_5.addLayout(self.verticalLayout_4)
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.label_4 = QtWidgets.QLabel(Dialog)
        self.label_4.setObjectName("label_4")
        self.horizontalLayout_5.addWidget(self.label_4)
        self.spb_cycles = QtWidgets.QSpinBox(Dialog)
        self.spb_cycles.setMaximumSize(QtCore.QSize(80, 16777215))
        self.spb_cycles.setObjectName("spb_cycles")
        self.horizontalLayout_5.addWidget(self.spb_cycles)
        self.verticalLayout_5.addLayout(self.horizontalLayout_5)
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.label_5 = QtWidgets.QLabel(Dialog)
        self.label_5.setObjectName("label_5")
        self.horizontalLayout_6.addWidget(self.label_5)
        self.dspb_cutoff = QtWidgets.QDoubleSpinBox(Dialog)
        self.dspb_cutoff.setMaximumSize(QtCore.QSize(80, 16777215))
        self.dspb_cutoff.setObjectName("dspb_cutoff")
        self.horizontalLayout_6.addWidget(self.dspb_cutoff)
        self.verticalLayout_5.addLayout(self.horizontalLayout_6)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.lbl_wsl2 = QtWidgets.QLabel(Dialog)
        self.lbl_wsl2.setObjectName("lbl_wsl2")
        self.horizontalLayout.addWidget(self.lbl_wsl2)
        self.btn_install_wsl2 = QtWidgets.QPushButton(Dialog)
        self.btn_install_wsl2.setMaximumSize(QtCore.QSize(121, 16777215))
        self.btn_install_wsl2.setObjectName("btn_install_wsl2")
        self.horizontalLayout.addWidget(self.btn_install_wsl2)
        self.verticalLayout_5.addLayout(self.horizontalLayout)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.lbl_local_prediction = QtWidgets.QLabel(Dialog)
        self.lbl_local_prediction.setObjectName("lbl_local_prediction")
        self.horizontalLayout_3.addWidget(self.lbl_local_prediction)
        self.btn_install_local_prediction = QtWidgets.QPushButton(Dialog)
        self.btn_install_local_prediction.setMaximumSize(QtCore.QSize(121, 16777215))
        self.btn_install_local_prediction.setObjectName("btn_install_local_prediction")
        self.horizontalLayout_3.addWidget(self.btn_install_local_prediction)
        self.verticalLayout_5.addLayout(self.horizontalLayout_3)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.line_wsl_config = QtWidgets.QFrame(Dialog)
        self.line_wsl_config.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_wsl_config.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_wsl_config.setObjectName("line_wsl_config")
        self.verticalLayout.addWidget(self.line_wsl_config)
        self.lbl_wsl_config = QtWidgets.QLabel(Dialog)
        self.lbl_wsl_config.setObjectName("lbl_wsl_config")
        self.verticalLayout.addWidget(self.lbl_wsl_config)
        self.verticalLayout_5.addLayout(self.verticalLayout)
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.lbl_wsl_username = QtWidgets.QLabel(Dialog)
        self.lbl_wsl_username.setObjectName("lbl_wsl_username")
        self.horizontalLayout_7.addWidget(self.lbl_wsl_username)
        self.box_wsl_username = QtWidgets.QComboBox(Dialog)
        self.box_wsl_username.setObjectName("box_wsl_username")
        self.horizontalLayout_7.addWidget(self.box_wsl_username)
        self.verticalLayout_5.addLayout(self.horizontalLayout_7)
        spacerItem = QtWidgets.QSpacerItem(20, 48, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.MinimumExpanding)
        self.verticalLayout_5.addItem(spacerItem)
        self.horizontalLayout_8 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_8.setObjectName("horizontalLayout_8")
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_8.addItem(spacerItem1)
        self.btn_cancel = QtWidgets.QPushButton(Dialog)
        self.btn_cancel.setObjectName("btn_cancel")
        self.horizontalLayout_8.addWidget(self.btn_cancel)
        self.btn_ok = QtWidgets.QPushButton(Dialog)
        self.btn_ok.setObjectName("btn_ok")
        self.horizontalLayout_8.addWidget(self.btn_ok)
        self.verticalLayout_5.addLayout(self.horizontalLayout_8)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.label_3.setText(_translate("Dialog", "Current Workspace"))
        self.btn_workspace_dir.setText(_translate("Dialog", "..."))
        self.label_2.setText(_translate("Dialog", "Directory where the predictions get downloaded"))
        self.btn_zip_storage_dir.setText(_translate("Dialog", "..."))
        self.label_6.setText(_translate("Dialog", "Parameters for structure alignment:"))
        self.label_4.setText(_translate("Dialog", "Cycles"))
        self.label_5.setText(_translate("Dialog", "Cutoff in Å"))
        self.lbl_wsl2.setText(_translate("Dialog", "WSL2"))
        self.btn_install_wsl2.setText(_translate("Dialog", "Install"))
        self.lbl_local_prediction.setText(_translate("Dialog", "Local Prediction"))
        self.btn_install_local_prediction.setText(_translate("Dialog", "Install"))
        self.lbl_wsl_config.setText(_translate("Dialog", "WSL2 Configuration:"))
        self.lbl_wsl_username.setText(_translate("Dialog", "Username"))
        self.btn_cancel.setText(_translate("Dialog", "Cancel"))
        self.btn_ok.setText(_translate("Dialog", "OK"))
