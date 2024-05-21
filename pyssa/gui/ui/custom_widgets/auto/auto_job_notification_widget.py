# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\pyssa\gui\ui\custom_widgets\forms\job_notification_widget.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Form(object):

  def setupUi(self, Form):
    Form.setObjectName("Form")
    Form.resize(492, 45)
    self.verticalLayout_2 = QtWidgets.QVBoxLayout(Form)
    self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
    self.verticalLayout_2.setObjectName("verticalLayout_2")
    self.frame = QtWidgets.QFrame(Form)
    self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
    self.frame.setObjectName("frame")
    self.verticalLayout = QtWidgets.QVBoxLayout(self.frame)
    self.verticalLayout.setObjectName("verticalLayout")
    self.horizontalLayout = QtWidgets.QHBoxLayout()
    self.horizontalLayout.setObjectName("horizontalLayout")
    self.lbl_icon = QtWidgets.QLabel(self.frame)
    self.lbl_icon.setObjectName("lbl_icon")
    self.horizontalLayout.addWidget(self.lbl_icon)
    self.lbl_job_description = QtWidgets.QLabel(self.frame)
    self.lbl_job_description.setFrameShape(QtWidgets.QFrame.NoFrame)
    self.lbl_job_description.setObjectName("lbl_job_description")
    self.horizontalLayout.addWidget(self.lbl_job_description)
    spacerItem = QtWidgets.QSpacerItem(
        40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
    )
    self.horizontalLayout.addItem(spacerItem)
    self.btn_refresh = QtWidgets.QPushButton(self.frame)
    self.btn_refresh.setObjectName("btn_refresh")
    self.horizontalLayout.addWidget(self.btn_refresh)
    self.btn_open = QtWidgets.QPushButton(self.frame)
    self.btn_open.setObjectName("btn_open")
    self.horizontalLayout.addWidget(self.btn_open)
    self.btn_open_image = QtWidgets.QPushButton(self.frame)
    self.btn_open_image.setObjectName("btn_open_image")
    self.horizontalLayout.addWidget(self.btn_open_image)
    self.btn_clear = QtWidgets.QPushButton(self.frame)
    self.btn_clear.setObjectName("btn_clear")
    self.horizontalLayout.addWidget(self.btn_clear)
    self.verticalLayout.addLayout(self.horizontalLayout)
    self.verticalLayout_2.addWidget(self.frame)

    self.retranslateUi(Form)
    QtCore.QMetaObject.connectSlotsByName(Form)

  def retranslateUi(self, Form):
    _translate = QtCore.QCoreApplication.translate
    Form.setWindowTitle(_translate("Form", "Form"))
    self.lbl_icon.setText(_translate("Form", "Icon"))
    self.lbl_job_description.setText(_translate("Form", "A job description"))
    self.btn_refresh.setText(_translate("Form", "Refresh"))
    self.btn_open.setText(_translate("Form", "Open"))
    self.btn_open_image.setText(_translate("Form", "Open Image"))
    self.btn_clear.setText(_translate("Form", "Clear"))
