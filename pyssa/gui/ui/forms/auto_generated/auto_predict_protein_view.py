# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\pyssa\gui\ui\forms\predict_protein_view.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(460, 815)
        self.verticalLayout = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout.setContentsMargins(0, 9, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.tab_widget = QtWidgets.QTabWidget(Dialog)
        self.tab_widget.setObjectName("tab_widget")
        self.tab_prediction = QtWidgets.QWidget()
        self.tab_prediction.setObjectName("tab_prediction")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.tab_prediction)
        self.verticalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.frame = QtWidgets.QFrame(self.tab_prediction)
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.frame)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.verticalLayout_123 = QtWidgets.QVBoxLayout()
        self.verticalLayout_123.setObjectName("verticalLayout_123")
        self.verticalLayout_124 = QtWidgets.QVBoxLayout()
        self.verticalLayout_124.setObjectName("verticalLayout_124")
        self.lbl_proteins_to_predict = QtWidgets.QLabel(self.frame)
        self.lbl_proteins_to_predict.setObjectName("lbl_proteins_to_predict")
        self.verticalLayout_124.addWidget(self.lbl_proteins_to_predict)
        self.table_proteins_to_predict = QtWidgets.QTableWidget(self.frame)
        self.table_proteins_to_predict.setObjectName("table_proteins_to_predict")
        self.table_proteins_to_predict.setColumnCount(2)
        self.table_proteins_to_predict.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.table_proteins_to_predict.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.table_proteins_to_predict.setHorizontalHeaderItem(1, item)
        self.verticalLayout_124.addWidget(self.table_proteins_to_predict)
        self.verticalLayout_123.addLayout(self.verticalLayout_124)
        self.horizontalLayout_108 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_108.setObjectName("horizontalLayout_108")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_108.addItem(spacerItem)
        self.btn_prediction_remove = QtWidgets.QPushButton(self.frame)
        self.btn_prediction_remove.setObjectName("btn_prediction_remove")
        self.horizontalLayout_108.addWidget(self.btn_prediction_remove)
        self.verticalLayout_123.addLayout(self.horizontalLayout_108)
        self.verticalLayout_3.addLayout(self.verticalLayout_123)
        self.horizontalLayout_112 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_112.setObjectName("horizontalLayout_112")
        self.lbl_advanced_config = QtWidgets.QLabel(self.frame)
        self.lbl_advanced_config.setObjectName("lbl_advanced_config")
        self.horizontalLayout_112.addWidget(self.lbl_advanced_config)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_112.addItem(spacerItem1)
        self.btn_edit_advanced_config = QtWidgets.QPushButton(self.frame)
        self.btn_edit_advanced_config.setObjectName("btn_edit_advanced_config")
        self.horizontalLayout_112.addWidget(self.btn_edit_advanced_config)
        self.verticalLayout_3.addLayout(self.horizontalLayout_112)
        self.checkbox_add_analysis = QtWidgets.QCheckBox(self.frame)
        self.checkbox_add_analysis.setObjectName("checkbox_add_analysis")
        self.verticalLayout_3.addWidget(self.checkbox_add_analysis)
        self.verticalLayout_4.addWidget(self.frame)
        spacerItem2 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.MinimumExpanding)
        self.verticalLayout_4.addItem(spacerItem2)
        self.frame_bottom = QtWidgets.QFrame(self.tab_prediction)
        self.frame_bottom.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_bottom.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_bottom.setObjectName("frame_bottom")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.frame_bottom)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.btn_help = QtWidgets.QPushButton(self.frame_bottom)
        self.btn_help.setObjectName("btn_help")
        self.horizontalLayout.addWidget(self.btn_help)
        spacerItem3 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem3)
        self.lbl_go_to_analysis_setup = QtWidgets.QLabel(self.frame_bottom)
        self.lbl_go_to_analysis_setup.setObjectName("lbl_go_to_analysis_setup")
        self.horizontalLayout.addWidget(self.lbl_go_to_analysis_setup)
        self.btn_go_to_analysis_setup = QtWidgets.QPushButton(self.frame_bottom)
        self.btn_go_to_analysis_setup.setObjectName("btn_go_to_analysis_setup")
        self.horizontalLayout.addWidget(self.btn_go_to_analysis_setup)
        self.btn_cancel = QtWidgets.QPushButton(self.frame_bottom)
        self.btn_cancel.setObjectName("btn_cancel")
        self.horizontalLayout.addWidget(self.btn_cancel)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        self.verticalLayout_4.addWidget(self.frame_bottom)
        self.tab_widget.addTab(self.tab_prediction, "")
        self.tab_analysis = QtWidgets.QWidget()
        self.tab_analysis.setObjectName("tab_analysis")
        self.verticalLayout_7 = QtWidgets.QVBoxLayout(self.tab_analysis)
        self.verticalLayout_7.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.frame_3 = QtWidgets.QFrame(self.tab_analysis)
        self.frame_3.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_3.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_3.setObjectName("frame_3")
        self.verticalLayout_6 = QtWidgets.QVBoxLayout(self.frame_3)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.verticalLayout_131 = QtWidgets.QVBoxLayout()
        self.verticalLayout_131.setObjectName("verticalLayout_131")
        self.lbl_analysis_overview = QtWidgets.QLabel(self.frame_3)
        self.lbl_analysis_overview.setObjectName("lbl_analysis_overview")
        self.verticalLayout_131.addWidget(self.lbl_analysis_overview)
        self.list_analysis_overview = QtWidgets.QListWidget(self.frame_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.list_analysis_overview.sizePolicy().hasHeightForWidth())
        self.list_analysis_overview.setSizePolicy(sizePolicy)
        self.list_analysis_overview.setMinimumSize(QtCore.QSize(0, 0))
        self.list_analysis_overview.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.list_analysis_overview.setObjectName("list_analysis_overview")
        self.verticalLayout_131.addWidget(self.list_analysis_overview)
        self.horizontalLayout_116 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_116.setObjectName("horizontalLayout_116")
        spacerItem4 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_116.addItem(spacerItem4)
        self.btn_analysis_remove = QtWidgets.QPushButton(self.frame_3)
        self.btn_analysis_remove.setObjectName("btn_analysis_remove")
        self.horizontalLayout_116.addWidget(self.btn_analysis_remove)
        self.btn_analysis_add = QtWidgets.QPushButton(self.frame_3)
        self.btn_analysis_add.setObjectName("btn_analysis_add")
        self.horizontalLayout_116.addWidget(self.btn_analysis_add)
        self.verticalLayout_131.addLayout(self.horizontalLayout_116)
        self.verticalLayout_6.addLayout(self.verticalLayout_131)
        self.verticalLayout_7.addWidget(self.frame_3)
        spacerItem5 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.MinimumExpanding)
        self.verticalLayout_7.addItem(spacerItem5)
        self.frame_bottom_2 = QtWidgets.QFrame(self.tab_analysis)
        self.frame_bottom_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_bottom_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_bottom_2.setObjectName("frame_bottom_2")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.frame_bottom_2)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.btn_help_2 = QtWidgets.QPushButton(self.frame_bottom_2)
        self.btn_help_2.setObjectName("btn_help_2")
        self.horizontalLayout_2.addWidget(self.btn_help_2)
        spacerItem6 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem6)
        self.btn_analysis_back_4 = QtWidgets.QPushButton(self.frame_bottom_2)
        self.btn_analysis_back_4.setObjectName("btn_analysis_back_4")
        self.horizontalLayout_2.addWidget(self.btn_analysis_back_4)
        self.btn_start_prediction_analysis = QtWidgets.QPushButton(self.frame_bottom_2)
        self.btn_start_prediction_analysis.setObjectName("btn_start_prediction_analysis")
        self.horizontalLayout_2.addWidget(self.btn_start_prediction_analysis)
        self.verticalLayout_5.addLayout(self.horizontalLayout_2)
        self.verticalLayout_7.addWidget(self.frame_bottom_2)
        self.tab_widget.addTab(self.tab_analysis, "")
        self.verticalLayout.addWidget(self.tab_widget)

        self.retranslateUi(Dialog)
        self.tab_widget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.lbl_proteins_to_predict.setText(_translate("Dialog", "Proteins to predict"))
        item = self.table_proteins_to_predict.horizontalHeaderItem(0)
        item.setText(_translate("Dialog", "Chain"))
        item = self.table_proteins_to_predict.horizontalHeaderItem(1)
        item.setText(_translate("Dialog", "Sequence"))
        self.btn_prediction_remove.setText(_translate("Dialog", "Remove"))
        self.lbl_advanced_config.setText(_translate("Dialog", "Advanced configurations"))
        self.btn_edit_advanced_config.setText(_translate("Dialog", "Edit"))
        self.checkbox_add_analysis.setText(_translate("Dialog", "Add analysis"))
        self.btn_help.setText(_translate("Dialog", "Help"))
        self.lbl_go_to_analysis_setup.setText(_translate("Dialog", "To analysis setup"))
        self.btn_go_to_analysis_setup.setText(_translate("Dialog", "Go"))
        self.btn_cancel.setText(_translate("Dialog", "Cancel"))
        self.tab_widget.setTabText(self.tab_widget.indexOf(self.tab_prediction), _translate("Dialog", "Prediction"))
        self.lbl_analysis_overview.setText(_translate("Dialog", "Protein structure analysis\'"))
        self.btn_analysis_remove.setText(_translate("Dialog", "Remove"))
        self.btn_analysis_add.setText(_translate("Dialog", "Add"))
        self.btn_help_2.setText(_translate("Dialog", "Help"))
        self.btn_analysis_back_4.setText(_translate("Dialog", "Back"))
        self.btn_start_prediction_analysis.setText(_translate("Dialog", "Start"))
        self.tab_widget.setTabText(self.tab_widget.indexOf(self.tab_analysis), _translate("Dialog", "Analysis"))
