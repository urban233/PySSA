#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2022
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/urban233/PySSA>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""Module for the job entry widget."""
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore
from PyQt5.QtCore import Qt
from pyssa.gui.ui import icon_resources  # this import is used for the icons! DO NOT DELETE THIS
from pyssa.gui.ui.custom_widgets.auto import auto_job_entry_widget
from pyssa.gui.ui.custom_widgets.auto import auto_job_notification_widget
from pyssa.internal.data_structures.data_classes import job_summary
from pyssa.internal.thread.async_pyssa import custom_signals
from pyssa.util import enums


class JobEntryWidget(QtWidgets.QWidget):
    """A widget for the job overview."""
    def __init__(self, a_job_description, a_job_base_information_object: "job_summary.JobBaseInformation"):
        super().__init__()
        self.ui = auto_job_entry_widget.Ui_Form()
        self.ui.setupUi(self)
        self.job_base_information: "job_summary.JobBaseInformation" = a_job_base_information_object

        self.ui.lbl_job_description.setText(f"{a_job_description} ({self.job_base_information.project_name})")

        tmp_cancel_job_icon = QtGui.QIcon(QtGui.QPixmap(":icons/do_not_disturb_on_w200.svg"))
        self.ui.btn_cancel_job.setIcon(tmp_cancel_job_icon)
        self.ui.btn_cancel_job.setText("")
        self.ui.btn_cancel_job.setIconSize(tmp_cancel_job_icon.actualSize(QtCore.QSize(24, 24)))

        self.ui.lbl_job_description.setStyleSheet("""
        QLabel {
            color: black;
            border-style: none;
        }
        """)
        self.ui.progress_bar_job.setStyleSheet("""
            QProgressBar {
                border-style: solid;
                border-width: 2px;
                border-radius: 4px;
                border-color: #367af6;
                background-color: #efefef;
                max-height: 5px;
                max-width: 100px;
            }
            QProgressBar::chunk {
                background-color: #367af6;
                width: 10px;
            }
        """)
        self.ui.btn_cancel_job.setStyleSheet("""
        QPushButton {
            background-color: rgba(220, 219, 227, 0.01);
            border: none;
            border-radius: 14px;
            min-width: 24px;
            max-width: 24px;
            min-height: 24px;
            max-height: 24px;
        }
        QPushButton::hover {
            background-color: rgba(220, 219, 227, 0.5);
            border: none;
            min-width: 24px;
            max-width: 24px;
            min-height: 24px;
            max-height: 24px;
        }
        """)
        self.ui.frame.setStyleSheet("""
        QFrame {
            border-style: solid;
            border-width: 1px;
            border-radius: 4px;
            border-color: #DCDBE3;
            background-color: #f2f2f2;
            margin: 5px;
        }
        """)

        self.ui.progress_bar_job.setFormat("")


class JobNotificationWidget(QtWidgets.QWidget):
    def __init__(self,
                 a_description,
                 a_job_base_information_object: "job_summary.JobBaseInformation",
                 job_is_from_current_project: bool,
                 a_refresh_after_job_finished_signal: "custom_signals.RefreshAfterJobFinishedSignal"):
        super().__init__()
        self.ui = auto_job_notification_widget.Ui_Form()
        self.ui.setupUi(self)
        
        self.job_base_information: "job_summary.JobBaseInformation" = a_job_base_information_object
        self._refresh_after_job_finished_signal: "custom_signals.RefreshAfterJobFinishedSignal" = a_refresh_after_job_finished_signal

        # <editor-fold desc="Widget setup">
        self.ui.lbl_job_description.setText(f"{a_description} ({self.job_base_information.project_name})")
        if self.job_base_information.job_progress == enums.JobProgress.FINISHED:
            pixmap = QtGui.QPixmap(":icons/info_w200.svg")
            if job_is_from_current_project:
                self.ui.btn_refresh.show()
                self.ui.btn_open.hide()
                self.ui.btn_open_image.hide()
                self.ui.btn_clear.hide()
            else:
                self.ui.btn_refresh.hide()
                self.ui.btn_open.show()
                self.ui.btn_open_image.hide()
                self.ui.btn_clear.hide()
            if a_job_base_information_object.job_type == enums.JobType.RAY_TRACING:
                self.ui.btn_refresh.hide()
                self.ui.btn_open.hide()
                self.ui.btn_clear.hide()
                self.ui.btn_open_image.show()
                self.ui.btn_open_image.setText("Show")
        elif self.job_base_information.job_progress == enums.JobProgress.FAILED:
            pixmap = QtGui.QPixmap(":icons/error_w200.svg")
            self.ui.btn_refresh.hide()
            self.ui.btn_open.hide()
            self.ui.btn_open_image.hide()
            self.ui.btn_clear.show()
        else:
            pixmap = QtGui.QPixmap(":icons/warning_w200.svg")
            self.ui.btn_refresh.hide()
            self.ui.btn_open.hide()
            self.ui.btn_open_image.hide()
            self.ui.btn_clear.show()

        # Set the scaled pixmap to the QLabel
        scaled_pixmap = pixmap.scaled(26, 26,
                                      aspectRatioMode=Qt.KeepAspectRatio,
                                      transformMode=Qt.SmoothTransformation)
        self.ui.lbl_icon.setPixmap(scaled_pixmap)
        self.ui.lbl_icon.setAlignment(Qt.AlignCenter)

        self.ui.lbl_job_description.setStyleSheet("""
            QLabel {
                color: black;
                background-color: #f2f2f2;
                border-style: none;
            }
        """)
        self.ui.lbl_icon.setStyleSheet("""
            QLabel {
                color: black;
                background-color: #f2f2f2;
                border-style: none;
            }
        """)
        self.ui.frame.setStyleSheet("""
                QFrame#frame {
                    border-style: solid;
                    border-width: 1px;
                    border-radius: 4px;
                    border-color: #DCDBE3;
                    background-color: #f2f2f2;
                    margin: 5px;
                }
                """)
        # </editor-fold>

        # Connect signals
        self.ui.btn_open.clicked.connect(self._send_open_other_project_request)
        self.ui.btn_open_image.clicked.connect(self._send_update_main_view_request)
        self.ui.btn_refresh.clicked.connect(self._send_update_main_view_request)
        self.ui.btn_clear.clicked.connect(self._send_update_main_view_request)

    def _send_open_other_project_request(self):
        self._refresh_after_job_finished_signal.emit_signal(False, self.job_base_information, self)

    def _send_update_main_view_request(self):
        self._refresh_after_job_finished_signal.emit_signal(True, self.job_base_information, self)
