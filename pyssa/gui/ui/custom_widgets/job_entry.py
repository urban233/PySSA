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
from pyssa.gui.ui import icon_resources  # this import is used for the icons! DO NOT DELETE THIS
from pyssa.internal.data_structures.data_classes import job_summary
from pyssa.internal.thread.async_pyssa import custom_signals
from pyssa.util import enums


class JobEntryWidget(QtWidgets.QWidget):
    """A widget for the job overview."""
    def __init__(self, a_job_description, a_job_base_information_object: "job_summary.JobBaseInformation"):
        super().__init__()
        self.job_base_information: "job_summary.JobBaseInformation" = a_job_base_information_object

        # <editor-fold desc="Widget setup">
        self.lbl_job_description = QtWidgets.QLabel(f"{a_job_description} ({self.job_base_information.project_name})")
        self.progress_bar_job = QtWidgets.QProgressBar()
        self.btn_cancel_job = QtWidgets.QPushButton('Cancel')
        tmp_cancel_job_icon = QtGui.QIcon(QtGui.QPixmap(":icons/cancel_w200.svg"))
        tmp_cancel_job_icon.addPixmap(QtGui.QPixmap(":icons/cancel_disabled_w200.svg"),
                                              mode=QtGui.QIcon.Mode.Disabled)
        self.btn_cancel_job.setIcon(tmp_cancel_job_icon)
        self.btn_cancel_job.setText("")
        self.btn_cancel_job.setIconSize(tmp_cancel_job_icon.actualSize(QtCore.QSize(24, 24)))
        self.progress_bar_job.setStyleSheet("""
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
        self.btn_cancel_job.setStyleSheet("""
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

        self.main_layout = QtWidgets.QHBoxLayout()
        self.main_layout.addWidget(self.lbl_job_description)
        self.main_layout.addWidget(self.progress_bar_job)
        self.main_layout.addWidget(self.btn_cancel_job)
        self.setLayout(self.main_layout)
        # </editor-fold>


class JobNotificationWidget(QtWidgets.QWidget):
    def __init__(self,
                 a_job_base_information_object: "job_summary.JobBaseInformation",
                 job_is_from_current_project: bool,
                 a_refresh_after_job_finished_signal: "custom_signals.RefreshAfterJobFinishedSignal"):
        super().__init__()
        self.job_base_information: "job_summary.JobBaseInformation" = a_job_base_information_object
        self._refresh_after_job_finished_signal: "custom_signals.RefreshAfterJobFinishedSignal" = a_refresh_after_job_finished_signal

        # <editor-fold desc="Widget setup">
        if self.job_base_information.job_type == enums.JobType.PREDICTION:
            self.lbl_job_description = QtWidgets.QLabel(f"A structure prediction job finished. ({self.job_base_information.project_name})")
        elif self.job_base_information.job_type == enums.JobType.DISTANCE_ANALYSIS:
            self.lbl_job_description = QtWidgets.QLabel(f"A distance analysis job finished. ({self.job_base_information.project_name})")
        elif self.job_base_information.job_type == enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS:
            self.lbl_job_description = QtWidgets.QLabel(f"A ColabFold prediction + distance analysis job finished. ({self.job_base_information.project_name})")
        elif self.job_base_information.job_type == enums.JobType.RAY_TRACING:
            self.lbl_job_description = QtWidgets.QLabel(f"Create ray-tracing image job finished. ({self.job_base_information.project_name})")
        else:
            self.lbl_job_description = QtWidgets.QLabel(f"Job finished. ({self.job_base_information.project_name})")

        self.btn_refresh = QtWidgets.QPushButton('Refresh')
        self.btn_open = QtWidgets.QPushButton('Open')

        self.main_layout = QtWidgets.QHBoxLayout()
        self.main_layout.addWidget(self.lbl_job_description)
        self.main_layout.addWidget(self.btn_refresh)
        self.main_layout.addWidget(self.btn_open)
        self.setLayout(self.main_layout)
        if job_is_from_current_project:
            self.btn_refresh.show()
            self.btn_open.hide()
        else:
            self.btn_refresh.hide()
            self.btn_open.show()
        # </editor-fold>

        # Connect signals
        self.btn_open.clicked.connect(self._send_open_other_project_request)
        self.btn_refresh.clicked.connect(self._send_update_main_view_request)

    def _send_open_other_project_request(self):
        self._refresh_after_job_finished_signal.emit_signal(False, self.job_base_information, self)

    def _send_update_main_view_request(self):
        self._refresh_after_job_finished_signal.emit_signal(True, self.job_base_information, self)
