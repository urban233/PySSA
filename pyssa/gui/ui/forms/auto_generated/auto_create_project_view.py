# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\pyssa\gui\ui\forms\create_project_view.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):

  def setupUi(self, Dialog):
    Dialog.setObjectName("Dialog")
    Dialog.resize(414, 648)
    self.verticalLayout_4 = QtWidgets.QVBoxLayout(Dialog)
    self.verticalLayout_4.setContentsMargins(0, 0, 0, 0)
    self.verticalLayout_4.setObjectName("verticalLayout_4")
    self.frame = QtWidgets.QFrame(Dialog)
    self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
    self.frame.setObjectName("frame")
    self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.frame)
    self.verticalLayout_2.setObjectName("verticalLayout_2")
    self.verticalLayout_26 = QtWidgets.QVBoxLayout()
    self.verticalLayout_26.setSpacing(2)
    self.verticalLayout_26.setObjectName("verticalLayout_26")
    self.label_17 = QtWidgets.QLabel(self.frame)
    self.label_17.setObjectName("label_17")
    self.verticalLayout_26.addWidget(self.label_17)
    self.txt_new_project_name = QtWidgets.QLineEdit(self.frame)
    self.txt_new_project_name.setObjectName("txt_new_project_name")
    self.verticalLayout_26.addWidget(self.txt_new_project_name)
    self.lbl_new_status_project_name = QtWidgets.QLabel(self.frame)
    self.lbl_new_status_project_name.setAlignment(
        QtCore.Qt.AlignLeading | QtCore.Qt.AlignLeft | QtCore.Qt.AlignTop
    )
    self.lbl_new_status_project_name.setIndent(-1)
    self.lbl_new_status_project_name.setObjectName(
        "lbl_new_status_project_name"
    )
    self.verticalLayout_26.addWidget(self.lbl_new_status_project_name)
    self.verticalLayout_2.addLayout(self.verticalLayout_26)
    self.verticalLayout = QtWidgets.QVBoxLayout()
    self.verticalLayout.setObjectName("verticalLayout")
    self.label_27 = QtWidgets.QLabel(self.frame)
    self.label_27.setObjectName("label_27")
    self.verticalLayout.addWidget(self.label_27)
    self.list_create_projects_view = QtWidgets.QListView(self.frame)
    self.list_create_projects_view.setObjectName("list_create_projects_view")
    self.verticalLayout.addWidget(self.list_create_projects_view)
    self.cb_new_add_reference = QtWidgets.QCheckBox(self.frame)
    self.cb_new_add_reference.setCheckable(True)
    self.cb_new_add_reference.setObjectName("cb_new_add_reference")
    self.verticalLayout.addWidget(self.cb_new_add_reference)
    self.verticalLayout_2.addLayout(self.verticalLayout)
    self.verticalLayout_36 = QtWidgets.QVBoxLayout()
    self.verticalLayout_36.setSpacing(2)
    self.verticalLayout_36.setObjectName("verticalLayout_36")
    self.lbl_new_choose_reference = QtWidgets.QLabel(self.frame)
    self.lbl_new_choose_reference.setObjectName("lbl_new_choose_reference")
    self.verticalLayout_36.addWidget(self.lbl_new_choose_reference)
    self.horizontalLayout_26 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_26.setObjectName("horizontalLayout_26")
    self.txt_new_choose_reference = QtWidgets.QLineEdit(self.frame)
    self.txt_new_choose_reference.setObjectName("txt_new_choose_reference")
    self.horizontalLayout_26.addWidget(self.txt_new_choose_reference)
    self.btn_new_choose_reference = QtWidgets.QToolButton(self.frame)
    self.btn_new_choose_reference.setObjectName("btn_new_choose_reference")
    self.horizontalLayout_26.addWidget(self.btn_new_choose_reference)
    self.verticalLayout_36.addLayout(self.horizontalLayout_26)
    self.lbl_new_status_choose_reference = QtWidgets.QLabel(self.frame)
    self.lbl_new_status_choose_reference.setObjectName(
        "lbl_new_status_choose_reference"
    )
    self.verticalLayout_36.addWidget(self.lbl_new_status_choose_reference)
    self.verticalLayout_2.addLayout(self.verticalLayout_36)
    self.verticalLayout_4.addWidget(self.frame)
    self.frame_bottom = QtWidgets.QFrame(Dialog)
    self.frame_bottom.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self.frame_bottom.setFrameShadow(QtWidgets.QFrame.Raised)
    self.frame_bottom.setObjectName("frame_bottom")
    self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.frame_bottom)
    self.verticalLayout_3.setObjectName("verticalLayout_3")
    self.horizontalLayout = QtWidgets.QHBoxLayout()
    self.horizontalLayout.setObjectName("horizontalLayout")
    self.btn_help = QtWidgets.QPushButton(self.frame_bottom)
    self.btn_help.setObjectName("btn_help")
    self.horizontalLayout.addWidget(self.btn_help)
    spacerItem = QtWidgets.QSpacerItem(
        148, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
    )
    self.horizontalLayout.addItem(spacerItem)
    self.btn_new_create_project = QtWidgets.QPushButton(self.frame_bottom)
    self.btn_new_create_project.setMinimumSize(QtCore.QSize(131, 0))
    self.btn_new_create_project.setObjectName("btn_new_create_project")
    self.horizontalLayout.addWidget(self.btn_new_create_project)
    self.btn_cancel = QtWidgets.QPushButton(self.frame_bottom)
    self.btn_cancel.setObjectName("btn_cancel")
    self.horizontalLayout.addWidget(self.btn_cancel)
    self.verticalLayout_3.addLayout(self.horizontalLayout)
    self.verticalLayout_4.addWidget(self.frame_bottom)

    self.retranslateUi(Dialog)
    QtCore.QMetaObject.connectSlotsByName(Dialog)

  def retranslateUi(self, Dialog):
    _translate = QtCore.QCoreApplication.translate
    Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
    self.label_17.setText(_translate("Dialog", "Enter new project name"))
    self.lbl_new_status_project_name.setText(_translate("Dialog", "TextLabel"))
    self.label_27.setText(_translate("Dialog", "Existing projects"))
    self.cb_new_add_reference.setText(
        _translate("Dialog", "Import existing protein")
    )
    self.lbl_new_choose_reference.setText(
        _translate("Dialog", "Import an existing .pdb file or enter a PDB ID")
    )
    self.btn_new_choose_reference.setText(_translate("Dialog", "..."))
    self.lbl_new_status_choose_reference.setText(
        _translate("Dialog", "TextLabel")
    )
    self.btn_help.setText(_translate("Dialog", "Help"))
    self.btn_new_create_project.setText(_translate("Dialog", "Create"))
    self.btn_cancel.setText(_translate("Dialog", "Cancel"))
