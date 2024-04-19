import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QStatusBar, QPushButton, QDialog, QVBoxLayout, QProgressBar
from PyQt5.QtCore import QTimer


class JobProgressPopup(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Job Progress")
        self.layout = QVBoxLayout()
        self.progressbars = {}
        self.setup_ui()

    def setup_ui(self):
        self.setLayout(self.layout)

    def add_job(self, job_name):
        progressbar = QProgressBar()
        self.progressbars[job_name] = progressbar
        self.layout.addWidget(progressbar)

    def update_progress(self, job_name, progress):
        if job_name in self.progressbars:
            self.progressbars[job_name].setValue(progress)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Job Progress Example")
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        self.init_ui()

    def init_ui(self):
        self.popup = None

        button = QPushButton("Show Progress")
        button.clicked.connect(self.show_progress_popup)
        self.statusBar.addPermanentWidget(button)

        # Simulating background jobs
        self.jobs_progress = {
            "Job 1": 0,
            "Job 2": 0,
            "Job 3": 0
        }
        self.timer = QTimer()
        self.timer.timeout.connect(self.update_jobs_progress)
        self.timer.start(1000)

    def show_progress_popup(self):
        if not self.popup:
            self.popup = JobProgressPopup(self)
            for job_name in self.jobs_progress:
                self.popup.add_job(job_name)
            self.popup.show()

    def update_jobs_progress(self):
        for job_name in self.jobs_progress:
            self.jobs_progress[job_name] += 10
            if self.jobs_progress[job_name] > 100:
                self.jobs_progress[job_name] = 0
            if self.popup:
                self.popup.update_progress(job_name, self.jobs_progress[job_name])


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
