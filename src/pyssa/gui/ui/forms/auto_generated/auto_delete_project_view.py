# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\pyssa\gui\ui\forms\delete_project_view.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):

  def setupUi(self, Dialog):
    Dialog.setObjectName("Dialog")
    Dialog.resize(459, 605)
    self.verticalLayout_5 = QtWidgets.QVBoxLayout(Dialog)
    self.verticalLayout_5.setContentsMargins(0, 0, 0, 0)
    self.verticalLayout_5.setObjectName("verticalLayout_5")
    self.frame = QtWidgets.QFrame(Dialog)
    self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
    self.frame.setObjectName("frame")
    self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.frame)
    self.verticalLayout_3.setObjectName("verticalLayout_3")
    self.verticalLayout = QtWidgets.QVBoxLayout()
    self.verticalLayout.setObjectName("verticalLayout")
    self.label_31 = QtWidgets.QLabel(self.frame)
    self.label_31.setObjectName("label_31")
    self.verticalLayout.addWidget(self.label_31)
    self.txt_delete_search = QtWidgets.QLineEdit(self.frame)
    self.txt_delete_search.setObjectName("txt_delete_search")
    self.verticalLayout.addWidget(self.txt_delete_search)
    self.lbl_delete_status_search = QtWidgets.QLabel(self.frame)
    self.lbl_delete_status_search.setObjectName("lbl_delete_status_search")
    self.verticalLayout.addWidget(self.lbl_delete_status_search)
    self.verticalLayout_3.addLayout(self.verticalLayout)
    self.verticalLayout_2 = QtWidgets.QVBoxLayout()
    self.verticalLayout_2.setObjectName("verticalLayout_2")
    self.label_32 = QtWidgets.QLabel(self.frame)
    self.label_32.setObjectName("label_32")
    self.verticalLayout_2.addWidget(self.label_32)
    self.list_delete_projects_view = QtWidgets.QListView(self.frame)
    self.list_delete_projects_view.setObjectName("list_delete_projects_view")
    self.verticalLayout_2.addWidget(self.list_delete_projects_view)
    self.verticalLayout_3.addLayout(self.verticalLayout_2)
    self.verticalLayout_34 = QtWidgets.QVBoxLayout()
    self.verticalLayout_34.setSpacing(2)
    self.verticalLayout_34.setObjectName("verticalLayout_34")
    self.label_33 = QtWidgets.QLabel(self.frame)
    self.label_33.setObjectName("label_33")
    self.verticalLayout_34.addWidget(self.label_33)
    self.txt_delete_selected_projects = QtWidgets.QLineEdit(self.frame)
    self.txt_delete_selected_projects.setReadOnly(True)
    self.txt_delete_selected_projects.setObjectName(
        "txt_delete_selected_projects"
    )
    self.verticalLayout_34.addWidget(self.txt_delete_selected_projects)
    self.verticalLayout_3.addLayout(self.verticalLayout_34)
    self.verticalLayout_5.addWidget(self.frame)
    self.frame_bottom = QtWidgets.QFrame(Dialog)
    self.frame_bottom.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self.frame_bottom.setFrameShadow(QtWidgets.QFrame.Raised)
    self.frame_bottom.setObjectName("frame_bottom")
    self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.frame_bottom)
    self.verticalLayout_4.setObjectName("verticalLayout_4")
    self.horizontalLayout = QtWidgets.QHBoxLayout()
    self.horizontalLayout.setObjectName("horizontalLayout")
    self.btn_help = QtWidgets.QPushButton(self.frame_bottom)
    self.btn_help.setObjectName("btn_help")
    self.horizontalLayout.addWidget(self.btn_help)
    spacerItem = QtWidgets.QSpacerItem(
        40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
    )
    self.horizontalLayout.addItem(spacerItem)
    self.btn_delete_delete_project = QtWidgets.QPushButton(self.frame_bottom)
    self.btn_delete_delete_project.setMinimumSize(QtCore.QSize(131, 0))
    self.btn_delete_delete_project.setObjectName("btn_delete_delete_project")
    self.horizontalLayout.addWidget(self.btn_delete_delete_project)
    self.btn_cancel = QtWidgets.QPushButton(self.frame_bottom)
    self.btn_cancel.setObjectName("btn_cancel")
    self.horizontalLayout.addWidget(self.btn_cancel)
    self.verticalLayout_4.addLayout(self.horizontalLayout)
    self.verticalLayout_5.addWidget(self.frame_bottom)

    self.retranslateUi(Dialog)
    QtCore.QMetaObject.connectSlotsByName(Dialog)

  def retranslateUi(self, Dialog):
    _translate = QtCore.QCoreApplication.translate
    Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
    self.label_31.setText(_translate("Dialog", "Search"))
    self.lbl_delete_status_search.setText(_translate("Dialog", "TextLabel"))
    self.label_32.setText(_translate("Dialog", "Existing projects"))
    self.label_33.setText(_translate("Dialog", "Selected project"))
    self.btn_help.setText(_translate("Dialog", "Help"))
    self.btn_delete_delete_project.setText(_translate("Dialog", "Delete"))
    self.btn_cancel.setText(_translate("Dialog", "Cancel"))
