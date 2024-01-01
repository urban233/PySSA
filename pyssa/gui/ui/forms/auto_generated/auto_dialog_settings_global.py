# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\dlgSettingsGlobal.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(383, 431)
        Dialog.setContextMenuPolicy(QtCore.Qt.NoContextMenu)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem)
        self.btn_info = QtWidgets.QPushButton(Dialog)
        self.btn_info.setObjectName("btn_info")
        self.horizontalLayout_2.addWidget(self.btn_info)
        self.verticalLayout_2.addLayout(self.horizontalLayout_2)
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
        self.verticalLayout_2.addLayout(self.verticalLayout_3)
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
        self.verticalLayout_2.addLayout(self.verticalLayout_4)
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.label_4 = QtWidgets.QLabel(Dialog)
        self.label_4.setObjectName("label_4")
        self.horizontalLayout_5.addWidget(self.label_4)
        self.spb_cycles = QtWidgets.QSpinBox(Dialog)
        self.spb_cycles.setMaximumSize(QtCore.QSize(80, 16777215))
        self.spb_cycles.setObjectName("spb_cycles")
        self.horizontalLayout_5.addWidget(self.spb_cycles)
        self.verticalLayout_2.addLayout(self.horizontalLayout_5)
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.label_5 = QtWidgets.QLabel(Dialog)
        self.label_5.setObjectName("label_5")
        self.horizontalLayout_6.addWidget(self.label_5)
        self.dspb_cutoff = QtWidgets.QDoubleSpinBox(Dialog)
        self.dspb_cutoff.setMaximumSize(QtCore.QSize(80, 16777215))
        self.dspb_cutoff.setObjectName("dspb_cutoff")
        self.horizontalLayout_6.addWidget(self.dspb_cutoff)
        self.verticalLayout_2.addLayout(self.horizontalLayout_6)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.line_2 = QtWidgets.QFrame(Dialog)
        self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
        self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.line_2.setObjectName("line_2")
        self.verticalLayout.addWidget(self.line_2)
        self.label = QtWidgets.QLabel(Dialog)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.lbl_ask_to_save_pymol_session = QtWidgets.QLabel(Dialog)
        self.lbl_ask_to_save_pymol_session.setObjectName("lbl_ask_to_save_pymol_session")
        self.horizontalLayout_7.addWidget(self.lbl_ask_to_save_pymol_session)
        self.cb_ask_to_save_pymol_session = QtWidgets.QCheckBox(Dialog)
        self.cb_ask_to_save_pymol_session.setText("")
        self.cb_ask_to_save_pymol_session.setObjectName("cb_ask_to_save_pymol_session")
        self.horizontalLayout_7.addWidget(self.cb_ask_to_save_pymol_session)
        self.verticalLayout_2.addLayout(self.horizontalLayout_7)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.lbl_color_vision_mode = QtWidgets.QLabel(Dialog)
        self.lbl_color_vision_mode.setObjectName("lbl_color_vision_mode")
        self.horizontalLayout.addWidget(self.lbl_color_vision_mode)
        self.cb_color_vision_mode = QtWidgets.QComboBox(Dialog)
        self.cb_color_vision_mode.setObjectName("cb_color_vision_mode")
        self.horizontalLayout.addWidget(self.cb_color_vision_mode)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        spacerItem1 = QtWidgets.QSpacerItem(
            20, 48, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.MinimumExpanding
        )
        self.verticalLayout_2.addItem(spacerItem1)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        spacerItem2 = QtWidgets.QSpacerItem(148, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem2)
        self.btn_cancel = QtWidgets.QPushButton(Dialog)
        self.btn_cancel.setObjectName("btn_cancel")
        self.horizontalLayout_3.addWidget(self.btn_cancel)
        self.btn_ok = QtWidgets.QPushButton(Dialog)
        self.btn_ok.setObjectName("btn_ok")
        self.horizontalLayout_3.addWidget(self.btn_ok)
        self.verticalLayout_2.addLayout(self.horizontalLayout_3)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.btn_info.setText(_translate("Dialog", "Info"))
        self.label_3.setText(_translate("Dialog", "Current Workspace"))
        self.btn_workspace_dir.setText(_translate("Dialog", "..."))
        self.label_6.setText(_translate("Dialog", "Parameters for structure alignment"))
        self.label_4.setText(_translate("Dialog", "Cycles"))
        self.label_5.setText(_translate("Dialog", "Cutoff in Å"))
        self.label.setText(_translate("Dialog", "Parameters for PyMOL session"))
        self.lbl_ask_to_save_pymol_session.setText(_translate("Dialog", "Ask to save pymol session"))
        self.lbl_color_vision_mode.setText(_translate("Dialog", "Color vision mode"))
        self.btn_cancel.setText(_translate("Dialog", "Cancel"))
        self.btn_ok.setText(_translate("Dialog", "OK"))
