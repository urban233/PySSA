# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\pyssa\gui\ui\forms\open_project_view.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(480, 523)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout_28 = QtWidgets.QVBoxLayout()
        self.verticalLayout_28.setSpacing(2)
        self.verticalLayout_28.setObjectName("verticalLayout_28")
        self.label_28 = QtWidgets.QLabel(Dialog)
        self.label_28.setObjectName("label_28")
        self.verticalLayout_28.addWidget(self.label_28)
        self.txt_open_search = QtWidgets.QLineEdit(Dialog)
        self.txt_open_search.setObjectName("txt_open_search")
        self.verticalLayout_28.addWidget(self.txt_open_search)
        self.lbl_open_status_search = QtWidgets.QLabel(Dialog)
        self.lbl_open_status_search.setObjectName("lbl_open_status_search")
        self.verticalLayout_28.addWidget(self.lbl_open_status_search)
        self.verticalLayout_2.addLayout(self.verticalLayout_28)
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.label = QtWidgets.QLabel(Dialog)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.projects_list_view = QtWidgets.QListView(Dialog)
        self.projects_list_view.setObjectName("projects_list_view")
        self.verticalLayout.addWidget(self.projects_list_view)
        self.verticalLayout_2.addLayout(self.verticalLayout)
        self.verticalLayout_30 = QtWidgets.QVBoxLayout()
        self.verticalLayout_30.setSpacing(2)
        self.verticalLayout_30.setObjectName("verticalLayout_30")
        self.label_30 = QtWidgets.QLabel(Dialog)
        self.label_30.setObjectName("label_30")
        self.verticalLayout_30.addWidget(self.label_30)
        self.txt_open_selected_project = QtWidgets.QLineEdit(Dialog)
        self.txt_open_selected_project.setReadOnly(True)
        self.txt_open_selected_project.setObjectName("txt_open_selected_project")
        self.verticalLayout_30.addWidget(self.txt_open_selected_project)
        self.verticalLayout_2.addLayout(self.verticalLayout_30)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.btn_open_project = QtWidgets.QPushButton(Dialog)
        self.btn_open_project.setObjectName("btn_open_project")
        self.horizontalLayout.addWidget(self.btn_open_project)
        self.verticalLayout_2.addLayout(self.horizontalLayout)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.label_28.setText(_translate("Dialog", "Search"))
        self.lbl_open_status_search.setText(_translate("Dialog", "TextLabel"))
        self.label.setText(_translate("Dialog", "Existing projects"))
        self.label_30.setText(_translate("Dialog", "Selected project"))
        self.btn_open_project.setText(_translate("Dialog", "Open"))