# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\pyssa\gui\ui\forms\predict_monomer_view.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(441, 866)
        self.verticalLayout = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.tabWidget = QtWidgets.QTabWidget(Dialog)
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.tab)
        self.verticalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.frame = QtWidgets.QFrame(self.tab)
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.frame)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.verticalLayout_104 = QtWidgets.QVBoxLayout()
        self.verticalLayout_104.setObjectName("verticalLayout_104")
        self.verticalLayout_105 = QtWidgets.QVBoxLayout()
        self.verticalLayout_105.setObjectName("verticalLayout_105")
        self.lbl_pred_analysis_mono_prot_to_predict = QtWidgets.QLabel(self.frame)
        self.lbl_pred_analysis_mono_prot_to_predict.setObjectName("lbl_pred_analysis_mono_prot_to_predict")
        self.verticalLayout_105.addWidget(self.lbl_pred_analysis_mono_prot_to_predict)
        self.table_pred_analysis_mono_prot_to_predict = QtWidgets.QTableWidget(self.frame)
        self.table_pred_analysis_mono_prot_to_predict.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.table_pred_analysis_mono_prot_to_predict.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustIgnored)
        self.table_pred_analysis_mono_prot_to_predict.setObjectName("table_pred_analysis_mono_prot_to_predict")
        self.table_pred_analysis_mono_prot_to_predict.setColumnCount(2)
        self.table_pred_analysis_mono_prot_to_predict.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.table_pred_analysis_mono_prot_to_predict.setHorizontalHeaderItem(0, item)
        item = QtWidgets.QTableWidgetItem()
        self.table_pred_analysis_mono_prot_to_predict.setHorizontalHeaderItem(1, item)
        self.verticalLayout_105.addWidget(self.table_pred_analysis_mono_prot_to_predict)
        self.verticalLayout_104.addLayout(self.verticalLayout_105)
        self.horizontalLayout_87 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_87.setObjectName("horizontalLayout_87")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_87.addItem(spacerItem)
        self.btn_pred_analysis_mono_seq_to_predict_remove = QtWidgets.QPushButton(self.frame)
        self.btn_pred_analysis_mono_seq_to_predict_remove.setObjectName("btn_pred_analysis_mono_seq_to_predict_remove")
        self.horizontalLayout_87.addWidget(self.btn_pred_analysis_mono_seq_to_predict_remove)
        self.btn_pred_analysis_mono_seq_to_predict = QtWidgets.QPushButton(self.frame)
        self.btn_pred_analysis_mono_seq_to_predict.setMinimumSize(QtCore.QSize(90, 0))
        self.btn_pred_analysis_mono_seq_to_predict.setMaximumSize(QtCore.QSize(131, 16777215))
        self.btn_pred_analysis_mono_seq_to_predict.setObjectName("btn_pred_analysis_mono_seq_to_predict")
        self.horizontalLayout_87.addWidget(self.btn_pred_analysis_mono_seq_to_predict)
        self.verticalLayout_104.addLayout(self.horizontalLayout_87)
        self.verticalLayout_2.addLayout(self.verticalLayout_104)
        self.verticalLayout_106 = QtWidgets.QVBoxLayout()
        self.verticalLayout_106.setObjectName("verticalLayout_106")
        self.lbl_pred_analysis_mono_prot_name = QtWidgets.QLabel(self.frame)
        self.lbl_pred_analysis_mono_prot_name.setObjectName("lbl_pred_analysis_mono_prot_name")
        self.verticalLayout_106.addWidget(self.lbl_pred_analysis_mono_prot_name)
        self.txt_pred_analysis_mono_prot_name = QtWidgets.QLineEdit(self.frame)
        self.txt_pred_analysis_mono_prot_name.setObjectName("txt_pred_analysis_mono_prot_name")
        self.verticalLayout_106.addWidget(self.txt_pred_analysis_mono_prot_name)
        self.lbl_pred_analysis_mono_prot_name_status = QtWidgets.QLabel(self.frame)
        self.lbl_pred_analysis_mono_prot_name_status.setObjectName("lbl_pred_analysis_mono_prot_name_status")
        self.verticalLayout_106.addWidget(self.lbl_pred_analysis_mono_prot_name_status)
        self.verticalLayout_2.addLayout(self.verticalLayout_106)
        self.horizontalLayout_90 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_90.setObjectName("horizontalLayout_90")
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_90.addItem(spacerItem1)
        self.btn_pred_analysis_mono_back = QtWidgets.QPushButton(self.frame)
        self.btn_pred_analysis_mono_back.setObjectName("btn_pred_analysis_mono_back")
        self.horizontalLayout_90.addWidget(self.btn_pred_analysis_mono_back)
        self.btn_pred_analysis_mono_next = QtWidgets.QPushButton(self.frame)
        self.btn_pred_analysis_mono_next.setObjectName("btn_pred_analysis_mono_next")
        self.horizontalLayout_90.addWidget(self.btn_pred_analysis_mono_next)
        self.verticalLayout_2.addLayout(self.horizontalLayout_90)
        self.verticalLayout_107 = QtWidgets.QVBoxLayout()
        self.verticalLayout_107.setObjectName("verticalLayout_107")
        self.lbl_pred_analysis_mono_seq_name = QtWidgets.QLabel(self.frame)
        self.lbl_pred_analysis_mono_seq_name.setObjectName("lbl_pred_analysis_mono_seq_name")
        self.verticalLayout_107.addWidget(self.lbl_pred_analysis_mono_seq_name)
        self.txt_pred_analysis_mono_seq_name = QtWidgets.QTextEdit(self.frame)
        self.txt_pred_analysis_mono_seq_name.setObjectName("txt_pred_analysis_mono_seq_name")
        self.verticalLayout_107.addWidget(self.txt_pred_analysis_mono_seq_name)
        self.lbl_pred_analysis_mono_seq_name_status = QtWidgets.QLabel(self.frame)
        self.lbl_pred_analysis_mono_seq_name_status.setObjectName("lbl_pred_analysis_mono_seq_name_status")
        self.verticalLayout_107.addWidget(self.lbl_pred_analysis_mono_seq_name_status)
        self.verticalLayout_2.addLayout(self.verticalLayout_107)
        self.horizontalLayout_89 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_89.setObjectName("horizontalLayout_89")
        spacerItem2 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_89.addItem(spacerItem2)
        self.btn_pred_analysis_mono_back_2 = QtWidgets.QPushButton(self.frame)
        self.btn_pred_analysis_mono_back_2.setObjectName("btn_pred_analysis_mono_back_2")
        self.horizontalLayout_89.addWidget(self.btn_pred_analysis_mono_back_2)
        self.btn_pred_analysis_mono_add_protein = QtWidgets.QPushButton(self.frame)
        self.btn_pred_analysis_mono_add_protein.setObjectName("btn_pred_analysis_mono_add_protein")
        self.horizontalLayout_89.addWidget(self.btn_pred_analysis_mono_add_protein)
        self.verticalLayout_2.addLayout(self.horizontalLayout_89)
        self.horizontalLayout_88 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_88.setObjectName("horizontalLayout_88")
        self.lbl_pred_mono_advanced_config_2 = QtWidgets.QLabel(self.frame)
        self.lbl_pred_mono_advanced_config_2.setObjectName("lbl_pred_mono_advanced_config_2")
        self.horizontalLayout_88.addWidget(self.lbl_pred_mono_advanced_config_2)
        spacerItem3 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_88.addItem(spacerItem3)
        self.btn_pred_mono_advanced_config_2 = QtWidgets.QPushButton(self.frame)
        self.btn_pred_mono_advanced_config_2.setMinimumSize(QtCore.QSize(90, 0))
        self.btn_pred_mono_advanced_config_2.setMaximumSize(QtCore.QSize(90, 16777215))
        self.btn_pred_mono_advanced_config_2.setObjectName("btn_pred_mono_advanced_config_2")
        self.horizontalLayout_88.addWidget(self.btn_pred_mono_advanced_config_2)
        self.verticalLayout_2.addLayout(self.horizontalLayout_88)
        self.checkbox_add_analysis = QtWidgets.QCheckBox(self.frame)
        self.checkbox_add_analysis.setObjectName("checkbox_add_analysis")
        self.verticalLayout_2.addWidget(self.checkbox_add_analysis)
        self.verticalLayout_4.addWidget(self.frame)
        spacerItem4 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.MinimumExpanding)
        self.verticalLayout_4.addItem(spacerItem4)
        self.frame_bottom = QtWidgets.QFrame(self.tab)
        self.frame_bottom.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_bottom.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_bottom.setObjectName("frame_bottom")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.frame_bottom)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.horizontalLayout_91 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_91.setObjectName("horizontalLayout_91")
        self.btn_help = QtWidgets.QPushButton(self.frame_bottom)
        self.btn_help.setObjectName("btn_help")
        self.horizontalLayout_91.addWidget(self.btn_help)
        spacerItem5 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_91.addItem(spacerItem5)
        self.lbl_pred_analysis_mono_to_analysis_setup = QtWidgets.QLabel(self.frame_bottom)
        self.lbl_pred_analysis_mono_to_analysis_setup.setObjectName("lbl_pred_analysis_mono_to_analysis_setup")
        self.horizontalLayout_91.addWidget(self.lbl_pred_analysis_mono_to_analysis_setup)
        self.btn_pred_analysis_mono_go_analysis_setup = QtWidgets.QPushButton(self.frame_bottom)
        self.btn_pred_analysis_mono_go_analysis_setup.setMinimumSize(QtCore.QSize(90, 0))
        self.btn_pred_analysis_mono_go_analysis_setup.setObjectName("btn_pred_analysis_mono_go_analysis_setup")
        self.horizontalLayout_91.addWidget(self.btn_pred_analysis_mono_go_analysis_setup)
        self.verticalLayout_3.addLayout(self.horizontalLayout_91)
        self.verticalLayout_4.addWidget(self.frame_bottom)
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.verticalLayout_7 = QtWidgets.QVBoxLayout(self.tab_2)
        self.verticalLayout_7.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.frame_2 = QtWidgets.QFrame(self.tab_2)
        self.frame_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_2.setObjectName("frame_2")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.frame_2)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.verticalLayout_113 = QtWidgets.QVBoxLayout()
        self.verticalLayout_113.setObjectName("verticalLayout_113")
        self.lbl_pred_analysis_mono_overview = QtWidgets.QLabel(self.frame_2)
        self.lbl_pred_analysis_mono_overview.setObjectName("lbl_pred_analysis_mono_overview")
        self.verticalLayout_113.addWidget(self.lbl_pred_analysis_mono_overview)
        self.list_pred_analysis_mono_overview = QtWidgets.QListWidget(self.frame_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.list_pred_analysis_mono_overview.sizePolicy().hasHeightForWidth())
        self.list_pred_analysis_mono_overview.setSizePolicy(sizePolicy)
        self.list_pred_analysis_mono_overview.setMinimumSize(QtCore.QSize(0, 0))
        self.list_pred_analysis_mono_overview.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.list_pred_analysis_mono_overview.setObjectName("list_pred_analysis_mono_overview")
        self.verticalLayout_113.addWidget(self.list_pred_analysis_mono_overview)
        self.horizontalLayout_98 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_98.setObjectName("horizontalLayout_98")
        spacerItem6 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_98.addItem(spacerItem6)
        self.btn_pred_analysis_mono_remove = QtWidgets.QPushButton(self.frame_2)
        self.btn_pred_analysis_mono_remove.setObjectName("btn_pred_analysis_mono_remove")
        self.horizontalLayout_98.addWidget(self.btn_pred_analysis_mono_remove)
        self.btn_pred_analysis_mono_add = QtWidgets.QPushButton(self.frame_2)
        self.btn_pred_analysis_mono_add.setObjectName("btn_pred_analysis_mono_add")
        self.horizontalLayout_98.addWidget(self.btn_pred_analysis_mono_add)
        self.verticalLayout_113.addLayout(self.horizontalLayout_98)
        self.verticalLayout_5.addLayout(self.verticalLayout_113)
        self.horizontalLayout_94 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_94.setObjectName("horizontalLayout_94")
        self.verticalLayout_111 = QtWidgets.QVBoxLayout()
        self.verticalLayout_111.setObjectName("verticalLayout_111")
        self.lbl_pred_analysis_mono_prot_struct_1 = QtWidgets.QLabel(self.frame_2)
        self.lbl_pred_analysis_mono_prot_struct_1.setAlignment(QtCore.Qt.AlignCenter)
        self.lbl_pred_analysis_mono_prot_struct_1.setObjectName("lbl_pred_analysis_mono_prot_struct_1")
        self.verticalLayout_111.addWidget(self.lbl_pred_analysis_mono_prot_struct_1)
        self.box_pred_analysis_mono_prot_struct_1 = QtWidgets.QComboBox(self.frame_2)
        self.box_pred_analysis_mono_prot_struct_1.setObjectName("box_pred_analysis_mono_prot_struct_1")
        self.verticalLayout_111.addWidget(self.box_pred_analysis_mono_prot_struct_1)
        self.horizontalLayout_94.addLayout(self.verticalLayout_111)
        self.lbl_analysis_batch_vs_2 = QtWidgets.QLabel(self.frame_2)
        self.lbl_analysis_batch_vs_2.setAlignment(QtCore.Qt.AlignCenter)
        self.lbl_analysis_batch_vs_2.setObjectName("lbl_analysis_batch_vs_2")
        self.horizontalLayout_94.addWidget(self.lbl_analysis_batch_vs_2)
        self.verticalLayout_112 = QtWidgets.QVBoxLayout()
        self.verticalLayout_112.setObjectName("verticalLayout_112")
        self.lbl_pred_analysis_mono_prot_struct_2 = QtWidgets.QLabel(self.frame_2)
        self.lbl_pred_analysis_mono_prot_struct_2.setAlignment(QtCore.Qt.AlignCenter)
        self.lbl_pred_analysis_mono_prot_struct_2.setObjectName("lbl_pred_analysis_mono_prot_struct_2")
        self.verticalLayout_112.addWidget(self.lbl_pred_analysis_mono_prot_struct_2)
        self.box_pred_analysis_mono_prot_struct_2 = QtWidgets.QComboBox(self.frame_2)
        self.box_pred_analysis_mono_prot_struct_2.setObjectName("box_pred_analysis_mono_prot_struct_2")
        self.verticalLayout_112.addWidget(self.box_pred_analysis_mono_prot_struct_2)
        self.horizontalLayout_94.addLayout(self.verticalLayout_112)
        self.verticalLayout_5.addLayout(self.horizontalLayout_94)
        self.horizontalLayout_97 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_97.setObjectName("horizontalLayout_97")
        spacerItem7 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_97.addItem(spacerItem7)
        self.btn_pred_analysis_mono_back_3 = QtWidgets.QPushButton(self.frame_2)
        self.btn_pred_analysis_mono_back_3.setObjectName("btn_pred_analysis_mono_back_3")
        self.horizontalLayout_97.addWidget(self.btn_pred_analysis_mono_back_3)
        self.btn_pred_analysis_mono_next_2 = QtWidgets.QPushButton(self.frame_2)
        self.btn_pred_analysis_mono_next_2.setObjectName("btn_pred_analysis_mono_next_2")
        self.horizontalLayout_97.addWidget(self.btn_pred_analysis_mono_next_2)
        self.verticalLayout_5.addLayout(self.horizontalLayout_97)
        self.verticalLayout_110 = QtWidgets.QVBoxLayout()
        self.verticalLayout_110.setObjectName("verticalLayout_110")
        self.lbl_pred_analysis_mono_ref_chains = QtWidgets.QLabel(self.frame_2)
        self.lbl_pred_analysis_mono_ref_chains.setObjectName("lbl_pred_analysis_mono_ref_chains")
        self.verticalLayout_110.addWidget(self.lbl_pred_analysis_mono_ref_chains)
        self.list_pred_analysis_mono_ref_chains = QtWidgets.QListWidget(self.frame_2)
        self.list_pred_analysis_mono_ref_chains.setMinimumSize(QtCore.QSize(0, 0))
        self.list_pred_analysis_mono_ref_chains.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.list_pred_analysis_mono_ref_chains.setObjectName("list_pred_analysis_mono_ref_chains")
        self.verticalLayout_110.addWidget(self.list_pred_analysis_mono_ref_chains)
        self.verticalLayout_5.addLayout(self.verticalLayout_110)
        self.horizontalLayout_96 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_96.setObjectName("horizontalLayout_96")
        spacerItem8 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_96.addItem(spacerItem8)
        self.btn_pred_analysis_mono_back_4 = QtWidgets.QPushButton(self.frame_2)
        self.btn_pred_analysis_mono_back_4.setObjectName("btn_pred_analysis_mono_back_4")
        self.horizontalLayout_96.addWidget(self.btn_pred_analysis_mono_back_4)
        self.btn_pred_analysis_mono_next_3 = QtWidgets.QPushButton(self.frame_2)
        self.btn_pred_analysis_mono_next_3.setObjectName("btn_pred_analysis_mono_next_3")
        self.horizontalLayout_96.addWidget(self.btn_pred_analysis_mono_next_3)
        self.verticalLayout_5.addLayout(self.horizontalLayout_96)
        self.verticalLayout_109 = QtWidgets.QVBoxLayout()
        self.verticalLayout_109.setObjectName("verticalLayout_109")
        self.lbl_pred_analysis_mono_model_chains = QtWidgets.QLabel(self.frame_2)
        self.lbl_pred_analysis_mono_model_chains.setObjectName("lbl_pred_analysis_mono_model_chains")
        self.verticalLayout_109.addWidget(self.lbl_pred_analysis_mono_model_chains)
        self.list_pred_analysis_mono_model_chains = QtWidgets.QListWidget(self.frame_2)
        self.list_pred_analysis_mono_model_chains.setMinimumSize(QtCore.QSize(0, 0))
        self.list_pred_analysis_mono_model_chains.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.list_pred_analysis_mono_model_chains.setObjectName("list_pred_analysis_mono_model_chains")
        self.verticalLayout_109.addWidget(self.list_pred_analysis_mono_model_chains)
        self.verticalLayout_5.addLayout(self.verticalLayout_109)
        self.horizontalLayout_92 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_92.setObjectName("horizontalLayout_92")
        spacerItem9 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_92.addItem(spacerItem9)
        self.btn_pred_analysis_mono_back_5 = QtWidgets.QPushButton(self.frame_2)
        self.btn_pred_analysis_mono_back_5.setObjectName("btn_pred_analysis_mono_back_5")
        self.horizontalLayout_92.addWidget(self.btn_pred_analysis_mono_back_5)
        self.btn_pred_analysis_mono_next_4 = QtWidgets.QPushButton(self.frame_2)
        self.btn_pred_analysis_mono_next_4.setObjectName("btn_pred_analysis_mono_next_4")
        self.horizontalLayout_92.addWidget(self.btn_pred_analysis_mono_next_4)
        self.verticalLayout_5.addLayout(self.horizontalLayout_92)
        self.horizontalLayout_95 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_95.setObjectName("horizontalLayout_95")
        self.cb_pred_analysis_mono_images = QtWidgets.QCheckBox(self.frame_2)
        self.cb_pred_analysis_mono_images.setObjectName("cb_pred_analysis_mono_images")
        self.horizontalLayout_95.addWidget(self.cb_pred_analysis_mono_images)
        self.verticalLayout_5.addLayout(self.horizontalLayout_95)
        self.verticalLayout_7.addWidget(self.frame_2)
        spacerItem10 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.MinimumExpanding)
        self.verticalLayout_7.addItem(spacerItem10)
        self.frame_bottom_2 = QtWidgets.QFrame(self.tab_2)
        self.frame_bottom_2.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame_bottom_2.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_bottom_2.setObjectName("frame_bottom_2")
        self.verticalLayout_6 = QtWidgets.QVBoxLayout(self.frame_bottom_2)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.horizontalLayout_93 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_93.setObjectName("horizontalLayout_93")
        self.btn_help_2 = QtWidgets.QPushButton(self.frame_bottom_2)
        self.btn_help_2.setObjectName("btn_help_2")
        self.horizontalLayout_93.addWidget(self.btn_help_2)
        spacerItem11 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_93.addItem(spacerItem11)
        self.btn_pred_analysis_mono_back_pred_setup = QtWidgets.QPushButton(self.frame_bottom_2)
        self.btn_pred_analysis_mono_back_pred_setup.setObjectName("btn_pred_analysis_mono_back_pred_setup")
        self.horizontalLayout_93.addWidget(self.btn_pred_analysis_mono_back_pred_setup)
        self.btn_pred_analysis_mono_start = QtWidgets.QPushButton(self.frame_bottom_2)
        self.btn_pred_analysis_mono_start.setObjectName("btn_pred_analysis_mono_start")
        self.horizontalLayout_93.addWidget(self.btn_pred_analysis_mono_start)
        self.verticalLayout_6.addLayout(self.horizontalLayout_93)
        self.verticalLayout_7.addWidget(self.frame_bottom_2)
        self.tabWidget.addTab(self.tab_2, "")
        self.verticalLayout.addWidget(self.tabWidget)

        self.retranslateUi(Dialog)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.lbl_pred_analysis_mono_prot_to_predict.setText(_translate("Dialog", "Proteins to predict"))
        self.table_pred_analysis_mono_prot_to_predict.setSortingEnabled(False)
        item = self.table_pred_analysis_mono_prot_to_predict.horizontalHeaderItem(0)
        item.setText(_translate("Dialog", "Chain"))
        item = self.table_pred_analysis_mono_prot_to_predict.horizontalHeaderItem(1)
        item.setText(_translate("Dialog", "Sequence"))
        self.btn_pred_analysis_mono_seq_to_predict_remove.setText(_translate("Dialog", "Remove"))
        self.btn_pred_analysis_mono_seq_to_predict.setText(_translate("Dialog", "Add"))
        self.lbl_pred_analysis_mono_prot_name.setText(_translate("Dialog", "Enter new protein name"))
        self.lbl_pred_analysis_mono_prot_name_status.setText(_translate("Dialog", "TextLabel"))
        self.btn_pred_analysis_mono_back.setText(_translate("Dialog", "Back"))
        self.btn_pred_analysis_mono_next.setText(_translate("Dialog", "Next"))
        self.lbl_pred_analysis_mono_seq_name.setText(_translate("Dialog", "Enter new protein sequence"))
        self.lbl_pred_analysis_mono_seq_name_status.setText(_translate("Dialog", "TextLabel"))
        self.btn_pred_analysis_mono_back_2.setText(_translate("Dialog", "Back"))
        self.btn_pred_analysis_mono_add_protein.setText(_translate("Dialog", "Add"))
        self.lbl_pred_mono_advanced_config_2.setText(_translate("Dialog", "Advanced configurations"))
        self.btn_pred_mono_advanced_config_2.setText(_translate("Dialog", "Edit"))
        self.checkbox_add_analysis.setText(_translate("Dialog", "Add analysis"))
        self.btn_help.setText(_translate("Dialog", "Help"))
        self.lbl_pred_analysis_mono_to_analysis_setup.setText(_translate("Dialog", "To analysis setup"))
        self.btn_pred_analysis_mono_go_analysis_setup.setText(_translate("Dialog", "Go"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("Dialog", "Prediction"))
        self.lbl_pred_analysis_mono_overview.setText(_translate("Dialog", "Protein structure analysis\'"))
        self.btn_pred_analysis_mono_remove.setText(_translate("Dialog", "Remove"))
        self.btn_pred_analysis_mono_add.setText(_translate("Dialog", "Add"))
        self.lbl_pred_analysis_mono_prot_struct_1.setText(_translate("Dialog", "Protein structure 1"))
        self.lbl_analysis_batch_vs_2.setText(_translate("Dialog", "vs."))
        self.lbl_pred_analysis_mono_prot_struct_2.setText(_translate("Dialog", "Protein structure 2"))
        self.btn_pred_analysis_mono_back_3.setText(_translate("Dialog", "Back"))
        self.btn_pred_analysis_mono_next_2.setText(_translate("Dialog", "Next"))
        self.lbl_pred_analysis_mono_ref_chains.setText(_translate("Dialog", "Select chains in protein structure 1"))
        self.btn_pred_analysis_mono_back_4.setText(_translate("Dialog", "Back"))
        self.btn_pred_analysis_mono_next_3.setText(_translate("Dialog", "Next"))
        self.lbl_pred_analysis_mono_model_chains.setText(_translate("Dialog", "Select chains in protein structure 2"))
        self.btn_pred_analysis_mono_back_5.setText(_translate("Dialog", "Back"))
        self.btn_pred_analysis_mono_next_4.setText(_translate("Dialog", "Next"))
        self.cb_pred_analysis_mono_images.setText(_translate("Dialog", "Create ray-traced images (slow)"))
        self.btn_help_2.setText(_translate("Dialog", "Help"))
        self.btn_pred_analysis_mono_back_pred_setup.setText(_translate("Dialog", "Back"))
        self.btn_pred_analysis_mono_start.setText(_translate("Dialog", "Start"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("Dialog", "Analysis"))
