# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\pyssa\gui\ui\forms\results_view.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):

  def setupUi(self, Dialog):
    Dialog.setObjectName("Dialog")
    Dialog.resize(375, 449)
    self.verticalLayout_4 = QtWidgets.QVBoxLayout(Dialog)
    self.verticalLayout_4.setContentsMargins(0, 0, 0, 0)
    self.verticalLayout_4.setObjectName("verticalLayout_4")
    self.frame = QtWidgets.QFrame(Dialog)
    self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
    self.frame.setObjectName("frame")
    self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.frame)
    self.verticalLayout_2.setObjectName("verticalLayout_2")
    self.verticalLayout = QtWidgets.QVBoxLayout()
    self.verticalLayout.setObjectName("verticalLayout")
    self.label = QtWidgets.QLabel(self.frame)
    self.label.setObjectName("label")
    self.verticalLayout.addWidget(self.label)
    self.table_widget_results = QtWidgets.QTableWidget(self.frame)
    self.table_widget_results.setObjectName("table_widget_results")
    self.table_widget_results.setColumnCount(0)
    self.table_widget_results.setRowCount(0)
    self.verticalLayout.addWidget(self.table_widget_results)
    self.verticalLayout_2.addLayout(self.verticalLayout)
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
        40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
    )
    self.horizontalLayout.addItem(spacerItem)
    self.btn_view_plots = QtWidgets.QPushButton(self.frame_bottom)
    self.btn_view_plots.setObjectName("btn_view_plots")
    self.horizontalLayout.addWidget(self.btn_view_plots)
    self.btn_export_data = QtWidgets.QPushButton(self.frame_bottom)
    self.btn_export_data.setObjectName("btn_export_data")
    self.horizontalLayout.addWidget(self.btn_export_data)
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
    self.label.setText(_translate("Dialog", "Results"))
    self.btn_help.setText(_translate("Dialog", "Help"))
    self.btn_view_plots.setText(_translate("Dialog", "View Plots"))
    self.btn_export_data.setText(_translate("Dialog", "Export Data"))
    self.btn_cancel.setText(_translate("Dialog", "Cancel"))