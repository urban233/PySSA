# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'prediction_and_analysis_widget.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(513, 578)
        self.verticalLayout = QtWidgets.QVBoxLayout(Form)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label_18 = QtWidgets.QLabel(Form)
        self.label_18.setObjectName("label_18")
        self.verticalLayout.addWidget(self.label_18)
        self.txt_prediction_project_name = QtWidgets.QLineEdit(Form)
        self.txt_prediction_project_name.setObjectName("txt_prediction_project_name")
        self.verticalLayout.addWidget(self.txt_prediction_project_name)
        self.label_15 = QtWidgets.QLabel(Form)
        self.label_15.setObjectName("label_15")
        self.verticalLayout.addWidget(self.label_15)
        self.horizontalLayout_17 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_17.setObjectName("horizontalLayout_17")
        self.txt_prediction_load_reference = QtWidgets.QLineEdit(Form)
        self.txt_prediction_load_reference.setObjectName("txt_prediction_load_reference")
        self.horizontalLayout_17.addWidget(self.txt_prediction_load_reference)
        self.btn_prediction_load_reference = QtWidgets.QToolButton(Form)
        self.btn_prediction_load_reference.setObjectName("btn_prediction_load_reference")
        self.horizontalLayout_17.addWidget(self.btn_prediction_load_reference)
        self.verticalLayout.addLayout(self.horizontalLayout_17)
        self.cb_prediction_chain_info = QtWidgets.QCheckBox(Form)
        self.cb_prediction_chain_info.setObjectName("cb_prediction_chain_info")
        self.verticalLayout.addWidget(self.cb_prediction_chain_info)
        self.verticalLayout_14 = QtWidgets.QVBoxLayout()
        self.verticalLayout_14.setObjectName("verticalLayout_14")
        self.label_16 = QtWidgets.QLabel(Form)
        self.label_16.setObjectName("label_16")
        self.verticalLayout_14.addWidget(self.label_16)
        self.txt_prediction_chain_ref = QtWidgets.QLineEdit(Form)
        self.txt_prediction_chain_ref.setObjectName("txt_prediction_chain_ref")
        self.verticalLayout_14.addWidget(self.txt_prediction_chain_ref)
        self.label_17 = QtWidgets.QLabel(Form)
        self.label_17.setObjectName("label_17")
        self.verticalLayout_14.addWidget(self.label_17)
        self.txt_prediction_chain_model = QtWidgets.QLineEdit(Form)
        self.txt_prediction_chain_model.setObjectName("txt_prediction_chain_model")
        self.verticalLayout_14.addWidget(self.txt_prediction_chain_model)
        self.verticalLayout.addLayout(self.verticalLayout_14)
        spacerItem = QtWidgets.QSpacerItem(20, 235, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.horizontalLayout_18 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_18.setObjectName("horizontalLayout_18")
        self.btn_prediction_cancel = QtWidgets.QPushButton(Form)
        self.btn_prediction_cancel.setObjectName("btn_prediction_cancel")
        self.horizontalLayout_18.addWidget(self.btn_prediction_cancel)
        self.btn_prediction_start = QtWidgets.QPushButton(Form)
        self.btn_prediction_start.setObjectName("btn_prediction_start")
        self.horizontalLayout_18.addWidget(self.btn_prediction_start)
        self.verticalLayout.addLayout(self.horizontalLayout_18)
        self.progress_bar_prediction = QtWidgets.QProgressBar(Form)
        self.progress_bar_prediction.setProperty("value", 24)
        self.progress_bar_prediction.setObjectName("progress_bar_prediction")
        self.verticalLayout.addWidget(self.progress_bar_prediction)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Form"))
        self.label_18.setText(_translate("Form", "Enter a project name:"))
        self.label_15.setText(_translate("Form", "Choose a .pdb file or enter a PDB ID:"))
        self.btn_prediction_load_reference.setText(_translate("Form", "..."))
        self.cb_prediction_chain_info.setText(_translate("Form", "Add chain information"))
        self.label_16.setText(_translate("Form", "Chains in reference (input style: G,H)"))
        self.label_17.setText(_translate("Form", "Chains in models (input style: A,B)"))
        self.btn_prediction_cancel.setText(_translate("Form", "Cancel"))
        self.btn_prediction_start.setText(_translate("Form", "Predict with Colab Notebook"))

