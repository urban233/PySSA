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
        Dialog.resize(592, 1157)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout_26 = QtWidgets.QVBoxLayout()
        self.verticalLayout_26.setSpacing(2)
        self.verticalLayout_26.setObjectName("verticalLayout_26")
        self.label_17 = QtWidgets.QLabel(Dialog)
        self.label_17.setObjectName("label_17")
        self.verticalLayout_26.addWidget(self.label_17)
        self.txt_new_project_name = QtWidgets.QLineEdit(Dialog)
        self.txt_new_project_name.setObjectName("txt_new_project_name")
        self.verticalLayout_26.addWidget(self.txt_new_project_name)
        self.lbl_new_status_project_name = QtWidgets.QLabel(Dialog)
        self.lbl_new_status_project_name.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.lbl_new_status_project_name.setIndent(-1)
        self.lbl_new_status_project_name.setObjectName("lbl_new_status_project_name")
        self.verticalLayout_26.addWidget(self.lbl_new_status_project_name)
        self.verticalLayout_2.addLayout(self.verticalLayout_26)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.label_27 = QtWidgets.QLabel(Dialog)
        self.label_27.setObjectName("label_27")
        self.verticalLayout.addWidget(self.label_27)
        self.list_new_projects = QtWidgets.QListWidget(Dialog)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.list_new_projects.sizePolicy().hasHeightForWidth())
        self.list_new_projects.setSizePolicy(sizePolicy)
        self.list_new_projects.setMinimumSize(QtCore.QSize(0, 0))
        self.list_new_projects.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.list_new_projects.setObjectName("list_new_projects")
        self.verticalLayout.addWidget(self.list_new_projects)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        self.cb_new_add_reference = QtWidgets.QCheckBox(Dialog)
        self.cb_new_add_reference.setCheckable(True)
        self.cb_new_add_reference.setObjectName("cb_new_add_reference")
        self.verticalLayout_2.addWidget(self.cb_new_add_reference)
        self.verticalLayout_36 = QtWidgets.QVBoxLayout()
        self.verticalLayout_36.setSpacing(2)
        self.verticalLayout_36.setObjectName("verticalLayout_36")
        self.lbl_new_choose_reference = QtWidgets.QLabel(Dialog)
        self.lbl_new_choose_reference.setObjectName("lbl_new_choose_reference")
        self.verticalLayout_36.addWidget(self.lbl_new_choose_reference)
        self.horizontalLayout_26 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_26.setObjectName("horizontalLayout_26")
        self.txt_new_choose_reference = QtWidgets.QLineEdit(Dialog)
        self.txt_new_choose_reference.setObjectName("txt_new_choose_reference")
        self.horizontalLayout_26.addWidget(self.txt_new_choose_reference)
        self.btn_new_choose_reference = QtWidgets.QToolButton(Dialog)
        self.btn_new_choose_reference.setObjectName("btn_new_choose_reference")
        self.horizontalLayout_26.addWidget(self.btn_new_choose_reference)
        self.verticalLayout_36.addLayout(self.horizontalLayout_26)
        self.lbl_new_status_choose_reference = QtWidgets.QLabel(Dialog)
        self.lbl_new_status_choose_reference.setObjectName("lbl_new_status_choose_reference")
        self.verticalLayout_36.addWidget(self.lbl_new_status_choose_reference)
        self.verticalLayout_2.addLayout(self.verticalLayout_36)
        self.horizontalLayout_23 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_23.setObjectName("horizontalLayout_23")
        spacerItem = QtWidgets.QSpacerItem(198, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_23.addItem(spacerItem)
        self.btn_new_create_project = QtWidgets.QPushButton(Dialog)
        self.btn_new_create_project.setMinimumSize(QtCore.QSize(131, 0))
        self.btn_new_create_project.setObjectName("btn_new_create_project")
        self.horizontalLayout_23.addWidget(self.btn_new_create_project)
        self.btn_cancel = QtWidgets.QPushButton(Dialog)
        self.btn_cancel.setObjectName("btn_cancel")
        self.horizontalLayout_23.addWidget(self.btn_cancel)
        self.verticalLayout_2.addLayout(self.horizontalLayout_23)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.label_17.setText(_translate("Dialog", "Enter new project name"))
        self.lbl_new_status_project_name.setText(_translate("Dialog", "TextLabel"))
        self.label_27.setText(_translate("Dialog", "Existing projects"))
        self.cb_new_add_reference.setText(_translate("Dialog", "Add existing protein"))
        self.lbl_new_choose_reference.setText(_translate("Dialog", "Choose an existing .pdb file or enter a PDB ID"))
        self.btn_new_choose_reference.setText(_translate("Dialog", "..."))
        self.lbl_new_status_choose_reference.setText(_translate("Dialog", "TextLabel"))
        self.btn_new_create_project.setText(_translate("Dialog", "Create"))
        self.btn_cancel.setText(_translate("Dialog", "Cancel"))
