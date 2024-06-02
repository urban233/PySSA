# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\dlgDistancePlot.ui'
#
# Created by: PyQt5 UI code generator 5.12.3
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):

  def setupUi(self, Dialog):
    Dialog.setObjectName("Dialog")
    Dialog.resize(565, 650)
    self.main_Layout = QtWidgets.QVBoxLayout(Dialog)
    self.main_Layout.setObjectName("main_Layout")
    self.horizontalLayout_8 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_8.setObjectName("horizontalLayout_8")
    self.label_6 = QtWidgets.QLabel(Dialog)
    self.label_6.setObjectName("label_6")
    self.horizontalLayout_8.addWidget(self.label_6)
    self.dsp_distance_plot_from_range = QtWidgets.QDoubleSpinBox(Dialog)
    self.dsp_distance_plot_from_range.setMinimum(0.0)
    self.dsp_distance_plot_from_range.setObjectName(
        "dsp_distance_plot_from_range"
    )
    self.horizontalLayout_8.addWidget(self.dsp_distance_plot_from_range)
    self.label_7 = QtWidgets.QLabel(Dialog)
    self.label_7.setObjectName("label_7")
    self.horizontalLayout_8.addWidget(self.label_7)
    self.dsp_distance_plot_to_range = QtWidgets.QDoubleSpinBox(Dialog)
    self.dsp_distance_plot_to_range.setMinimum(0.0)
    self.dsp_distance_plot_to_range.setObjectName("dsp_distance_plot_to_range")
    self.horizontalLayout_8.addWidget(self.dsp_distance_plot_to_range)
    self.main_Layout.addLayout(self.horizontalLayout_8)
    self.horizontalLayout_6 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_6.setObjectName("horizontalLayout_6")
    self.label_2 = QtWidgets.QLabel(Dialog)
    self.label_2.setObjectName("label_2")
    self.horizontalLayout_6.addWidget(self.label_2)
    self.sp_distance_plot_from = QtWidgets.QSpinBox(Dialog)
    self.sp_distance_plot_from.setMaximum(9999)
    self.sp_distance_plot_from.setObjectName("sp_distance_plot_from")
    self.horizontalLayout_6.addWidget(self.sp_distance_plot_from)
    self.label_5 = QtWidgets.QLabel(Dialog)
    self.label_5.setObjectName("label_5")
    self.horizontalLayout_6.addWidget(self.label_5)
    self.sp_distance_plot_to = QtWidgets.QSpinBox(Dialog)
    self.sp_distance_plot_to.setMaximum(9999)
    self.sp_distance_plot_to.setObjectName("sp_distance_plot_to")
    self.horizontalLayout_6.addWidget(self.sp_distance_plot_to)
    self.main_Layout.addLayout(self.horizontalLayout_6)
    self.cb_sync_with_pymol = QtWidgets.QCheckBox(Dialog)
    self.cb_sync_with_pymol.setObjectName("cb_sync_with_pymol")
    self.main_Layout.addWidget(self.cb_sync_with_pymol)
    self.line = QtWidgets.QFrame(Dialog)
    self.line.setFrameShape(QtWidgets.QFrame.HLine)
    self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
    self.line.setObjectName("line")
    self.main_Layout.addWidget(self.line)
    self.label = QtWidgets.QLabel(Dialog)
    self.label.setText("")
    self.label.setObjectName("label")
    self.main_Layout.addWidget(self.label)
    self.horizontalLayout_10 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_10.setObjectName("horizontalLayout_10")
    self.btn_distance_plot_reset = QtWidgets.QPushButton(Dialog)
    self.btn_distance_plot_reset.setObjectName("btn_distance_plot_reset")
    self.horizontalLayout_10.addWidget(self.btn_distance_plot_reset)
    spacerItem = QtWidgets.QSpacerItem(
        40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
    )
    self.horizontalLayout_10.addItem(spacerItem)
    self.btn_distance_plot_update = QtWidgets.QPushButton(Dialog)
    self.btn_distance_plot_update.setObjectName("btn_distance_plot_update")
    self.horizontalLayout_10.addWidget(self.btn_distance_plot_update)
    self.main_Layout.addLayout(self.horizontalLayout_10)
    self.horizontalLayout_7 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_7.setObjectName("horizontalLayout_7")
    self.cb_turn_on_grid = QtWidgets.QCheckBox(Dialog)
    self.cb_turn_on_grid.setObjectName("cb_turn_on_grid")
    self.horizontalLayout_7.addWidget(self.cb_turn_on_grid)
    spacerItem1 = QtWidgets.QSpacerItem(
        40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
    )
    self.horizontalLayout_7.addItem(spacerItem1)
    self.btn_distance_plot_save = QtWidgets.QPushButton(Dialog)
    self.btn_distance_plot_save.setObjectName("btn_distance_plot_save")
    self.horizontalLayout_7.addWidget(self.btn_distance_plot_save)
    self.main_Layout.addLayout(self.horizontalLayout_7)

    self.retranslateUi(Dialog)
    QtCore.QMetaObject.connectSlotsByName(Dialog)

  def retranslateUi(self, Dialog):
    _translate = QtCore.QCoreApplication.translate
    Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
    self.label_6.setText(_translate("Dialog", "From Å"))
    self.label_7.setText(_translate("Dialog", "to"))
    self.label_2.setText(_translate("Dialog", "From residue pair"))
    self.label_5.setText(_translate("Dialog", "to"))
    self.cb_sync_with_pymol.setText(_translate("Dialog", "Sync with PyMOL"))
    self.btn_distance_plot_reset.setText(_translate("Dialog", "Reset Plot"))
    self.btn_distance_plot_update.setText(_translate("Dialog", "Update Plot"))
    self.cb_turn_on_grid.setText(_translate("Dialog", "Turn on grid"))
    self.btn_distance_plot_save.setText(_translate("Dialog", "Save as image"))
