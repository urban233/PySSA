# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\pyssa\gui\ui\forms\advanced_prediction_configurations_view.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):

  def setupUi(self, Dialog):
    Dialog.setObjectName("Dialog")
    Dialog.resize(390, 308)
    self.verticalLayout_3 = QtWidgets.QVBoxLayout(Dialog)
    self.verticalLayout_3.setContentsMargins(0, 0, 0, 0)
    self.verticalLayout_3.setObjectName("verticalLayout_3")
    self.frame_2 = QtWidgets.QFrame(Dialog)
    self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
    self.frame_2.setObjectName("frame_2")
    self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.frame_2)
    self.verticalLayout_2.setObjectName("verticalLayout_2")
    self.horizontalLayout = QtWidgets.QHBoxLayout()
    self.horizontalLayout.setObjectName("horizontalLayout")
    self.lbl_amber = QtWidgets.QLabel(self.frame_2)
    self.lbl_amber.setObjectName("lbl_amber")
    self.horizontalLayout.addWidget(self.lbl_amber)
    self.cb_amber = QtWidgets.QCheckBox(self.frame_2)
    self.cb_amber.setText("")
    self.cb_amber.setObjectName("cb_amber")
    self.horizontalLayout.addWidget(self.cb_amber)
    self.verticalLayout_2.addLayout(self.horizontalLayout)
    self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_2.setObjectName("horizontalLayout_2")
    self.lbl_template = QtWidgets.QLabel(self.frame_2)
    self.lbl_template.setObjectName("lbl_template")
    self.horizontalLayout_2.addWidget(self.lbl_template)
    self.combo_box_template = QtWidgets.QComboBox(self.frame_2)
    self.combo_box_template.setObjectName("combo_box_template")
    self.horizontalLayout_2.addWidget(self.combo_box_template)
    self.verticalLayout_2.addLayout(self.horizontalLayout_2)
    self.verticalLayout_3.addWidget(self.frame_2)
    spacerItem = QtWidgets.QSpacerItem(
        20, 29, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
    )
    self.verticalLayout_3.addItem(spacerItem)
    self.frame_bottom = QtWidgets.QFrame(Dialog)
    self.frame_bottom.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self.frame_bottom.setFrameShadow(QtWidgets.QFrame.Raised)
    self.frame_bottom.setObjectName("frame_bottom")
    self.verticalLayout = QtWidgets.QVBoxLayout(self.frame_bottom)
    self.verticalLayout.setObjectName("verticalLayout")
    self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_3.setObjectName("horizontalLayout_3")
    self.btn_help = QtWidgets.QPushButton(self.frame_bottom)
    self.btn_help.setObjectName("btn_help")
    self.horizontalLayout_3.addWidget(self.btn_help)
    spacerItem1 = QtWidgets.QSpacerItem(
        40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
    )
    self.horizontalLayout_3.addItem(spacerItem1)
    self.btn_ok = QtWidgets.QPushButton(self.frame_bottom)
    self.btn_ok.setObjectName("btn_ok")
    self.horizontalLayout_3.addWidget(self.btn_ok)
    self.btn_cancel = QtWidgets.QPushButton(self.frame_bottom)
    self.btn_cancel.setObjectName("btn_cancel")
    self.horizontalLayout_3.addWidget(self.btn_cancel)
    self.verticalLayout.addLayout(self.horizontalLayout_3)
    self.verticalLayout_3.addWidget(self.frame_bottom)

    self.retranslateUi(Dialog)
    QtCore.QMetaObject.connectSlotsByName(Dialog)

  def retranslateUi(self, Dialog):
    _translate = QtCore.QCoreApplication.translate
    Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
    self.lbl_amber.setText(_translate("Dialog", "Activate Amber force field"))
    self.lbl_template.setText(_translate("Dialog", "Choose template mode"))
    self.btn_help.setText(_translate("Dialog", "Help"))
    self.btn_ok.setText(_translate("Dialog", "OK"))
    self.btn_cancel.setText(_translate("Dialog", "Cancel"))
