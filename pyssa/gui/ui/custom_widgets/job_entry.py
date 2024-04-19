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


class JobEntry(QtWidgets.QWidget):
    """A widget for the job overview."""
    def __init__(self, a_job_description, a_project_name: str):
        super().__init__()
        self.lbl_job_description = QtWidgets.QLabel(f"{a_job_description} ({a_project_name})")
        self.progress_bar_job = QtWidgets.QProgressBar()
        self.btn_cancel_job = QtWidgets.QPushButton('Cancel')
        # TODO: add icon to cancel button
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
