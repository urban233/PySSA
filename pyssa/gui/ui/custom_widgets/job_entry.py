from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore


class JobEntry(QtWidgets.QWidget):
    """A widget for the job overview."""
    def __init__(self, a_job_description, a_project_name: str):
        super().__init__()
        self.lbl_job_description = QtWidgets.QLabel(f"{a_job_description} ({a_project_name})")
        self.progress_bar_job = QtWidgets.QProgressBar()
        self.btn_cancel_job = QtWidgets.QPushButton('Cancel')
        # TODO: add icon to cancel button
        delete_protein_session_icon = QtGui.QIcon(QtGui.QPixmap(":icons/cancel_w200.svg"))
        self.btn_cancel_job.setIcon(delete_protein_session_icon)
        self.btn_cancel_job.setText("")
        self.btn_cancel_job.setIconSize(delete_protein_session_icon.actualSize(QtCore.QSize(30, 30)))

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
            border-radius: 4px;
            min-width: 36px;
            max-width: 36px;
            min-height: 36px;
            max-height: 36px;
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
