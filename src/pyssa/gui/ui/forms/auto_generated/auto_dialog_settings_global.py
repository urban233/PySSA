# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '.\pyssa\gui\ui\forms\dlgSettingsGlobal.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):

  def setupUi(self, Dialog):
    Dialog.setObjectName("Dialog")
    Dialog.resize(518, 611)
    Dialog.setContextMenuPolicy(QtCore.Qt.NoContextMenu)
    self.verticalLayout_9 = QtWidgets.QVBoxLayout(Dialog)
    self.verticalLayout_9.setContentsMargins(0, 0, 0, 0)
    self.verticalLayout_9.setObjectName("verticalLayout_9")
    self.frame = QtWidgets.QFrame(Dialog)
    self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
    self.frame.setObjectName("frame")
    self.verticalLayout_8 = QtWidgets.QVBoxLayout(self.frame)
    self.verticalLayout_8.setObjectName("verticalLayout_8")
    self.frame_settings_general = QtWidgets.QFrame(self.frame)
    self.frame_settings_general.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self.frame_settings_general.setFrameShadow(QtWidgets.QFrame.Raised)
    self.frame_settings_general.setObjectName("frame_settings_general")
    self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.frame_settings_general)
    self.verticalLayout_5.setObjectName("verticalLayout_5")
    self.verticalLayout_3 = QtWidgets.QVBoxLayout()
    self.verticalLayout_3.setObjectName("verticalLayout_3")
    self.label_3 = QtWidgets.QLabel(self.frame_settings_general)
    self.label_3.setObjectName("label_3")
    self.verticalLayout_3.addWidget(self.label_3)
    self.horizontalLayout_4 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_4.setObjectName("horizontalLayout_4")
    self.txt_workspace_dir = QtWidgets.QLineEdit(self.frame_settings_general)
    self.txt_workspace_dir.setObjectName("txt_workspace_dir")
    self.horizontalLayout_4.addWidget(self.txt_workspace_dir)
    self.btn_workspace_dir = QtWidgets.QToolButton(self.frame_settings_general)
    self.btn_workspace_dir.setObjectName("btn_workspace_dir")
    self.horizontalLayout_4.addWidget(self.btn_workspace_dir)
    self.verticalLayout_3.addLayout(self.horizontalLayout_4)
    self.verticalLayout_5.addLayout(self.verticalLayout_3)
    self.check_box_start_help = QtWidgets.QCheckBox(self.frame_settings_general)
    self.check_box_start_help.setObjectName("check_box_start_help")
    self.verticalLayout_5.addWidget(self.check_box_start_help)
    self.verticalLayout_8.addWidget(self.frame_settings_general)
    self.frame_settings_distance_parameters = QtWidgets.QFrame(self.frame)
    self.frame_settings_distance_parameters.setFrameShape(
        QtWidgets.QFrame.StyledPanel
    )
    self.frame_settings_distance_parameters.setFrameShadow(
        QtWidgets.QFrame.Raised
    )
    self.frame_settings_distance_parameters.setObjectName(
        "frame_settings_distance_parameters"
    )
    self.verticalLayout_6 = QtWidgets.QVBoxLayout(
        self.frame_settings_distance_parameters
    )
    self.verticalLayout_6.setObjectName("verticalLayout_6")
    self.verticalLayout_4 = QtWidgets.QVBoxLayout()
    self.verticalLayout_4.setObjectName("verticalLayout_4")
    self.label_6 = QtWidgets.QLabel(self.frame_settings_distance_parameters)
    self.label_6.setObjectName("label_6")
    self.verticalLayout_4.addWidget(self.label_6)
    self.verticalLayout_6.addLayout(self.verticalLayout_4)
    self.line = QtWidgets.QFrame(self.frame_settings_distance_parameters)
    self.line.setFrameShape(QtWidgets.QFrame.HLine)
    self.line.setFrameShadow(QtWidgets.QFrame.Sunken)
    self.line.setObjectName("line")
    self.verticalLayout_6.addWidget(self.line)
    self.horizontalLayout_5 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_5.setObjectName("horizontalLayout_5")
    self.label_4 = QtWidgets.QLabel(self.frame_settings_distance_parameters)
    self.label_4.setObjectName("label_4")
    self.horizontalLayout_5.addWidget(self.label_4)
    self.spb_cycles = QtWidgets.QSpinBox(
        self.frame_settings_distance_parameters
    )
    self.spb_cycles.setMaximumSize(QtCore.QSize(60, 16777215))
    self.spb_cycles.setObjectName("spb_cycles")
    self.horizontalLayout_5.addWidget(self.spb_cycles)
    self.verticalLayout_6.addLayout(self.horizontalLayout_5)
    self.horizontalLayout_6 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_6.setObjectName("horizontalLayout_6")
    self.label_5 = QtWidgets.QLabel(self.frame_settings_distance_parameters)
    self.label_5.setObjectName("label_5")
    self.horizontalLayout_6.addWidget(self.label_5)
    self.dspb_cutoff = QtWidgets.QDoubleSpinBox(
        self.frame_settings_distance_parameters
    )
    self.dspb_cutoff.setMaximumSize(QtCore.QSize(60, 16777215))
    self.dspb_cutoff.setObjectName("dspb_cutoff")
    self.horizontalLayout_6.addWidget(self.dspb_cutoff)
    self.verticalLayout_6.addLayout(self.horizontalLayout_6)
    self.verticalLayout_8.addWidget(self.frame_settings_distance_parameters)
    self.frame_settings_image_parameters = QtWidgets.QFrame(self.frame)
    self.frame_settings_image_parameters.setFrameShape(
        QtWidgets.QFrame.StyledPanel
    )
    self.frame_settings_image_parameters.setFrameShadow(QtWidgets.QFrame.Raised)
    self.frame_settings_image_parameters.setObjectName(
        "frame_settings_image_parameters"
    )
    self.verticalLayout_7 = QtWidgets.QVBoxLayout(
        self.frame_settings_image_parameters
    )
    self.verticalLayout_7.setObjectName("verticalLayout_7")
    self.verticalLayout = QtWidgets.QVBoxLayout()
    self.verticalLayout.setObjectName("verticalLayout")
    self.label = QtWidgets.QLabel(self.frame_settings_image_parameters)
    self.label.setObjectName("label")
    self.verticalLayout.addWidget(self.label)
    self.line_2 = QtWidgets.QFrame(self.frame_settings_image_parameters)
    self.line_2.setFrameShape(QtWidgets.QFrame.HLine)
    self.line_2.setFrameShadow(QtWidgets.QFrame.Sunken)
    self.line_2.setObjectName("line_2")
    self.verticalLayout.addWidget(self.line_2)
    self.verticalLayout_7.addLayout(self.verticalLayout)
    self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_2.setObjectName("horizontalLayout_2")
    self.label_24 = QtWidgets.QLabel(self.frame_settings_image_parameters)
    self.label_24.setObjectName("label_24")
    self.horizontalLayout_2.addWidget(self.label_24)
    self.box_bg_color = QtWidgets.QComboBox(
        self.frame_settings_image_parameters
    )
    self.box_bg_color.setObjectName("box_bg_color")
    self.horizontalLayout_2.addWidget(self.box_bg_color)
    self.verticalLayout_7.addLayout(self.horizontalLayout_2)
    self.horizontalLayout_7 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_7.setObjectName("horizontalLayout_7")
    self.label_8 = QtWidgets.QLabel(self.frame_settings_image_parameters)
    self.label_8.setObjectName("label_8")
    self.horizontalLayout_7.addWidget(self.label_8)
    self.box_renderer = QtWidgets.QComboBox(
        self.frame_settings_image_parameters
    )
    self.box_renderer.setObjectName("box_renderer")
    self.horizontalLayout_7.addWidget(self.box_renderer)
    self.verticalLayout_7.addLayout(self.horizontalLayout_7)
    self.horizontalLayout_8 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_8.setObjectName("horizontalLayout_8")
    self.label_10 = QtWidgets.QLabel(self.frame_settings_image_parameters)
    self.label_10.setObjectName("label_10")
    self.horizontalLayout_8.addWidget(self.label_10)
    self.box_ray_trace_mode = QtWidgets.QComboBox(
        self.frame_settings_image_parameters
    )
    self.box_ray_trace_mode.setObjectName("box_ray_trace_mode")
    self.horizontalLayout_8.addWidget(self.box_ray_trace_mode)
    self.verticalLayout_7.addLayout(self.horizontalLayout_8)
    self.horizontalLayout_9 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_9.setObjectName("horizontalLayout_9")
    self.label_14 = QtWidgets.QLabel(self.frame_settings_image_parameters)
    self.label_14.setObjectName("label_14")
    self.horizontalLayout_9.addWidget(self.label_14)
    self.box_ray_texture = QtWidgets.QComboBox(
        self.frame_settings_image_parameters
    )
    self.box_ray_texture.setObjectName("box_ray_texture")
    self.horizontalLayout_9.addWidget(self.box_ray_texture)
    self.verticalLayout_7.addLayout(self.horizontalLayout_9)
    self.verticalLayout_8.addWidget(self.frame_settings_image_parameters)
    self.horizontalLayout = QtWidgets.QHBoxLayout()
    self.horizontalLayout.setObjectName("horizontalLayout")
    self.lbl_color_vision_mode = QtWidgets.QLabel(self.frame)
    self.lbl_color_vision_mode.setObjectName("lbl_color_vision_mode")
    self.horizontalLayout.addWidget(self.lbl_color_vision_mode)
    self.cb_color_vision_mode = QtWidgets.QComboBox(self.frame)
    self.cb_color_vision_mode.setObjectName("cb_color_vision_mode")
    self.horizontalLayout.addWidget(self.cb_color_vision_mode)
    self.verticalLayout_8.addLayout(self.horizontalLayout)
    self.verticalLayout_9.addWidget(self.frame)
    spacerItem = QtWidgets.QSpacerItem(
        20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding
    )
    self.verticalLayout_9.addItem(spacerItem)
    self.frame_bottom = QtWidgets.QFrame(Dialog)
    self.frame_bottom.setFrameShape(QtWidgets.QFrame.StyledPanel)
    self.frame_bottom.setFrameShadow(QtWidgets.QFrame.Raised)
    self.frame_bottom.setObjectName("frame_bottom")
    self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.frame_bottom)
    self.verticalLayout_2.setObjectName("verticalLayout_2")
    self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
    self.horizontalLayout_3.setObjectName("horizontalLayout_3")
    self.btn_help = QtWidgets.QPushButton(self.frame_bottom)
    self.btn_help.setObjectName("btn_help")
    self.horizontalLayout_3.addWidget(self.btn_help)
    spacerItem1 = QtWidgets.QSpacerItem(
        148, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum
    )
    self.horizontalLayout_3.addItem(spacerItem1)
    self.btn_ok = QtWidgets.QPushButton(self.frame_bottom)
    self.btn_ok.setObjectName("btn_ok")
    self.horizontalLayout_3.addWidget(self.btn_ok)
    self.verticalLayout_2.addLayout(self.horizontalLayout_3)
    self.verticalLayout_9.addWidget(self.frame_bottom)

    self.retranslateUi(Dialog)
    QtCore.QMetaObject.connectSlotsByName(Dialog)

  def retranslateUi(self, Dialog):
    _translate = QtCore.QCoreApplication.translate
    Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
    self.label_3.setText(_translate("Dialog", "Current Workspace"))
    self.btn_workspace_dir.setText(_translate("Dialog", "..."))
    self.check_box_start_help.setText(
        _translate("Dialog", "Start PySSA Help Center Automatically At Startup")
    )
    self.label_6.setText(_translate("Dialog", "Distance Analysis Parameters"))
    self.label_4.setText(_translate("Dialog", "Cycles"))
    self.label_5.setText(_translate("Dialog", "Cutoff in Å"))
    self.label.setText(_translate("Dialog", "Image Parameters"))
    self.label_24.setText(_translate("Dialog", "Background Color"))
    self.label_8.setText(_translate("Dialog", "Renderer"))
    self.label_10.setText(_translate("Dialog", "Ray-Trace-Mode"))
    self.label_14.setText(_translate("Dialog", "Ray Texture"))
    self.lbl_color_vision_mode.setText(
        _translate("Dialog", "Color vision mode")
    )
    self.btn_help.setText(_translate("Dialog", "Help"))
    self.btn_ok.setText(_translate("Dialog", "OK"))