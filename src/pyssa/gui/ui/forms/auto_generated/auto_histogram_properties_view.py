# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\pyssa\gui\ui\forms\histogram_properties_view.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):

  def setupUi(self, Dialog):
    Dialog.setObjectName("Dialog")
    Dialog.resize(306, 337)
    self.verticalLayout_7 = QtWidgets.QVBoxLayout(Dialog)
    self.verticalLayout_7.setContentsMargins(0, 0, 0, 0)
    self.verticalLayout_7.setObjectName("verticalLayout_7")
    self.frame = QtWidgets.QFrame(Dialog)
    self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
    self.frame.setObjectName("frame")
    self.horizontalLayout_4 = QtWidgets.QHBoxLayout(self.frame)
    self.horizontalLayout_4.setObjectName("horizontalLayout_4")
    self.tab_widget = QtWidgets.QTabWidget(self.frame)
    self.tab_widget.setObjectName("tab_widget")
    self.tab = QtWidgets.QWidget()
    self.tab.setObjectName("tab")
    self.verticalLayout_6 = QtWidgets.QVBoxLayout(self.tab)
    self.verticalLayout_6.setContentsMargins(9, -1, -1, -1)
    self.verticalLayout_6.setObjectName("verticalLayout_6")
    self.frame_x_axis = QtWidgets.QFrame(self.tab)
    self.frame_x_axis.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self.frame_x_axis.setFrameShadow(QtWidgets.QFrame.Raised)
    self.frame_x_axis.setObjectName("frame_x_axis")
    self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.frame_x_axis)
    self.verticalLayout_3.setObjectName("verticalLayout_3")
    self.verticalLayout_2 = QtWidgets.QVBoxLayout()
    self.verticalLayout_2.setObjectName("verticalLayout_2")
    self.label = QtWidgets.QLabel(self.frame_x_axis)
    self.label.setObjectName("label")
    self.verticalLayout_2.addWidget(self.label)
    self.line_x_axis = QtWidgets.QFrame(self.frame_x_axis)
    self.line_x_axis.setFrameShadow(QtWidgets.QFrame.Plain)
    self.line_x_axis.setFrameShape(QtWidgets.QFrame.HLine)
    self.line_x_axis.setObjectName("line_x_axis")
    self.verticalLayout_2.addWidget(self.line_x_axis)
    self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_2.setObjectName("horizontalLayout_2")
    self.label_4 = QtWidgets.QLabel(self.frame_x_axis)
    self.label_4.setObjectName("label_4")
    self.horizontalLayout_2.addWidget(self.label_4)
    spacerItem = QtWidgets.QSpacerItem(
        40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
    )
    self.horizontalLayout_2.addItem(spacerItem)
    self.le_units_x_axis = QtWidgets.QLineEdit(self.frame_x_axis)
    self.le_units_x_axis.setObjectName("le_units_x_axis")
    self.horizontalLayout_2.addWidget(self.le_units_x_axis)
    self.verticalLayout_2.addLayout(self.horizontalLayout_2)
    self.verticalLayout_3.addLayout(self.verticalLayout_2)
    self.verticalLayout_6.addWidget(self.frame_x_axis)
    self.frame_y_axis = QtWidgets.QFrame(self.tab)
    self.frame_y_axis.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self.frame_y_axis.setFrameShadow(QtWidgets.QFrame.Raised)
    self.frame_y_axis.setObjectName("frame_y_axis")
    self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.frame_y_axis)
    self.verticalLayout_5.setObjectName("verticalLayout_5")
    self.verticalLayout_4 = QtWidgets.QVBoxLayout()
    self.verticalLayout_4.setObjectName("verticalLayout_4")
    self.label_2 = QtWidgets.QLabel(self.frame_y_axis)
    self.label_2.setObjectName("label_2")
    self.verticalLayout_4.addWidget(self.label_2)
    self.line_y_axis = QtWidgets.QFrame(self.frame_y_axis)
    self.line_y_axis.setContextMenuPolicy(QtCore.Qt.DefaultContextMenu)
    self.line_y_axis.setFrameShadow(QtWidgets.QFrame.Plain)
    self.line_y_axis.setFrameShape(QtWidgets.QFrame.HLine)
    self.line_y_axis.setObjectName("line_y_axis")
    self.verticalLayout_4.addWidget(self.line_y_axis)
    self.horizontalLayout = QtWidgets.QHBoxLayout()
    self.horizontalLayout.setObjectName("horizontalLayout")
    self.label_3 = QtWidgets.QLabel(self.frame_y_axis)
    self.label_3.setObjectName("label_3")
    self.horizontalLayout.addWidget(self.label_3)
    spacerItem1 = QtWidgets.QSpacerItem(
        40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
    )
    self.horizontalLayout.addItem(spacerItem1)
    self.cb_distance_interval = QtWidgets.QComboBox(self.frame_y_axis)
    self.cb_distance_interval.setObjectName("cb_distance_interval")
    self.horizontalLayout.addWidget(self.cb_distance_interval)
    self.verticalLayout_4.addLayout(self.horizontalLayout)
    self.verticalLayout_5.addLayout(self.verticalLayout_4)
    self.verticalLayout_6.addWidget(self.frame_y_axis)
    spacerItem2 = QtWidgets.QSpacerItem(
        20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
    )
    self.verticalLayout_6.addItem(spacerItem2)
    self.tab_widget.addTab(self.tab, "")
    self.horizontalLayout_4.addWidget(self.tab_widget)
    self.verticalLayout_7.addWidget(self.frame)
    self.frame_bottom = QtWidgets.QFrame(Dialog)
    self.frame_bottom.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self.frame_bottom.setFrameShadow(QtWidgets.QFrame.Raised)
    self.frame_bottom.setObjectName("frame_bottom")
    self.verticalLayout = QtWidgets.QVBoxLayout(self.frame_bottom)
    self.verticalLayout.setObjectName("verticalLayout")
    self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_3.setObjectName("horizontalLayout_3")
    self.btn_help = QtWidgets.QPushButton(self.frame_bottom)
    self.btn_help.setObjectName("btn_help")
    self.horizontalLayout_3.addWidget(self.btn_help)
    spacerItem3 = QtWidgets.QSpacerItem(
        40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
    )
    self.horizontalLayout_3.addItem(spacerItem3)
    self.btn_save = QtWidgets.QPushButton(self.frame_bottom)
    self.btn_save.setObjectName("btn_save")
    self.horizontalLayout_3.addWidget(self.btn_save)
    self.btn_cancel = QtWidgets.QPushButton(self.frame_bottom)
    self.btn_cancel.setObjectName("btn_cancel")
    self.horizontalLayout_3.addWidget(self.btn_cancel)
    self.verticalLayout.addLayout(self.horizontalLayout_3)
    self.verticalLayout_7.addWidget(self.frame_bottom)

    self.retranslateUi(Dialog)
    self.tab_widget.setCurrentIndex(0)
    QtCore.QMetaObject.connectSlotsByName(Dialog)

  def retranslateUi(self, Dialog):
    _translate = QtCore.QCoreApplication.translate
    Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
    self.label.setText(_translate("Dialog", "X-Axis"))
    self.label_4.setText(_translate("Dialog", "Units"))
    self.label_2.setText(_translate("Dialog", "Y-Axis"))
    self.label_3.setText(_translate("Dialog", "Distance Interval"))
    self.tab_widget.setTabText(
        self.tab_widget.indexOf(self.tab), _translate("Dialog", "Plot")
    )
    self.btn_help.setText(_translate("Dialog", "Help"))
    self.btn_save.setText(_translate("Dialog", "Save"))
    self.btn_cancel.setText(_translate("Dialog", "Cancel"))