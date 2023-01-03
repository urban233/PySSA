import datetime
from pymol import Qt
from PyQt5.QtWidgets import (QHBoxLayout, QVBoxLayout,
                             QPushButton, QLabel)
from PyQt5 import QtCore
from PyQt5.QtCore import QUrl
from PyQt5.QtWebEngineWidgets import QWebEngineView, QWebEnginePage, QWebEngineDownloadItem
from pyssa.gui.utilities import constants
from pyssa.gui.utilities import gui_utils
from pyssa.gui.ui.dialogs import dialog_notebook_managment


class WebInterface(Qt.QtWidgets.QDialog):
    google_login_page = "https://accounts.google.com/v3/signin/identifier?dsh=S487182181%3A1668865110328278&continue=https%3A%2F%2Fcolab.research.google.com%2Fgithub%2Fsokrypton%2FColabFold%2Fblob%2Fmain%2FAlphaFold2.ipynb&ec=GAZAqQM&passive=true&flowName=GlifWebSignIn&flowEntry=ServiceLogin&ifkv=ARgdvAu0jNh8qAP4PLZtLHO1saEbt8gL7tdnPxJNzq-nQ8HD1od95scPNdlXBiSB5aGTVZWwLQ-f-A"
    colab_login_page = "https://accounts.google.com/ServiceLogin?passive=true&continue=https%3A%2F%2Fcolab.research.google.com%2Fgithub%2Fsokrypton%2FColabFold%2Fblob%2Fmain%2FAlphaFold2.ipynb&ec=GAZAqQM"
    login_page = "https://accounts.google.com/v3/signin/identifier?dsh=S-301988812%3A1672417305836195&continue=https%3A%2F%2Fwww.google.com%2F&ec=GAZAmgQ&hl=de&passive=true&flowName=GlifWebSignIn&flowEntry=ServiceLogin&ifkv=AeAAQh6nGazOofUFZXEEKfU4PB58c6ZxzdbUX2maYjsI4mCm7wqSvzZ3SCIux6lnC0hiRjVaT4p97w"
    url_without_pre_login = QUrl(
        "https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb")
    url_with_pre_login = QUrl(
        'https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2.ipynb?pli=1')
    account_page = QUrl("https://myaccount.google.com/?hl=de&utm_source=OGB&utm_medium=act&pli=1")

    abort_msg_box = None
    """Exit codes:
    0: Prediction successful
    1: Prediction failed
    2: init code
    3: gpu is not available
    4: colab session needs reconnect
    """
    exit_code = 2
    # the values need to be in single quotes!
    prediction_params = {
        "seq": 'AFRGQEAFRGQEAFRGQEAFRGQEAFRGQE',
        "job_name": '',
        "amber": 'false',
        "templates": ('pdb70', '1'),
    }

    def __init__(self, parent=None):
        """Constructor

        Args:
            args
            kwargs
        """
        Qt.QtWidgets.QDialog.__init__(self, parent)
        self.dialog: dialog_notebook_managment.DialogNotebookManagment
        self.status = ""
        self.current_time = ""
        # QWebEngineView is used to display the content from the QWebEnginePage
        self.browser = QWebEngineView()
        self.web_page = QWebEnginePage()
        self.browser.setPage(self.web_page)
        # additional gui elements to check the prediction status
        self.label = QLabel("The status of the prediction will be displayed here.")
        self.button = QPushButton("Check status")
        self.button_abort = QPushButton("Abort")
        self.button.setMinimumWidth(120)
        self.button_abort.setMinimumWidth(120)
        # basic layout
        outer_layout = QVBoxLayout()
        label_layout = QHBoxLayout()
        button_layout = QHBoxLayout()
        outer_layout.addWidget(self.browser)
        label_layout.addWidget(self.label)
        button_layout.addStretch()
        button_layout.addWidget(self.button)
        button_layout.addWidget(self.button_abort)
        outer_layout.addLayout(label_layout)
        outer_layout.addLayout(button_layout)
        self.setLayout(outer_layout)

        self.button.hide()
        self.label.hide()
        # Connections
        self.browser.page().profile().downloadRequested.connect(self.on_downloadRequested)
        self.button.clicked.connect(self.control_status)
        self.button_abort.clicked.connect(self.disconnect_active_runtime)
        # window settings
        self.setMinimumSize(575, 150)
        self.setMaximumSize(575, 150)
        self.setWindowTitle("Sign-In Google")
        self.show()

    @QtCore.pyqtSlot("QWebEngineDownloadItem*")
    def on_downloadRequested(self, download):
        """This function is called if the notebook wants to download the results

        Args:
            download:
                necessary argument which holds information about the requested download
        """
        path = f"{constants.SCRATCH_DIR}/{constants.NOTEBOOK_RESULTS_ZIP_NAME}.zip"
        download.setPath(path)
        download.accept()
        download.finished.connect(self.finished_prediction)
        self.exit_code = 0

    # ----- JavaScript & Python Hierarchy

    def open_colab_notebook(self):
        """This function opens the Google login page

        """
        url = QUrl(self.url_without_pre_login)
        self.web_page.setUrl(url)
        self.web_page.loadFinished.connect(self.check_authentication_status)

    def check_authentication_status(self):
        self.web_page.runJavaScript("""
                    if (document.querySelectorAll('a')[2].className === "gb_7 gb_8 gb_de gb_dd") {
                        document.title = "login needed";
                    } else if (document.querySelectorAll('a')[2].className === "gb_d gb_Ra gb_l") {
                        document.title = "no login needed";
                    } 
                    """, self.process_authentication_status)

    def process_authentication_status(self, value) -> None:
        if value == "login needed":
            # Login is needed to proceed!
            #TODO: finish this branch
            self.web_page.setUrl(QUrl(self.colab_login_page))
        elif value == "no login needed":
            # Login already happened.
            self.run_colab_notebook_js()
        else:
            print("Unexpected case happened.")

    def run_colab_notebook_js(self):
        self.browser.hide()
        self.label.show()
        self.button.show()
        self.current_time = datetime.datetime.now().strftime("%H:%M:%S")
        #self.web_page.runJavaScript("colab.global.notebook.kernel.attemptAutoconnect()")
        self.web_page.runJavaScript(
            f"document.querySelectorAll('paper-input').item(0).setAttribute('value', '{self.prediction_params['seq']}')")
        self.web_page.runJavaScript(
            "document.querySelectorAll('paper-input').item(0).dispatchEvent(new Event('change'))"
        )
        # sets the job name
        self.web_page.runJavaScript(
            f"document.querySelectorAll('paper-input').item(1).setAttribute('value', '{self.prediction_params['job_name']}')")
        self.web_page.runJavaScript(
            "document.querySelectorAll('paper-input').item(1).dispatchEvent(new Event('change'))"
        )
        # checks the amber force field checkbox
        self.web_page.runJavaScript(
            f"document.querySelectorAll('input[type=checkbox]').item(0).checked = '{self.prediction_params['amber']}'"
        )
        self.web_page.runJavaScript(
            "document.querySelectorAll('input[type=checkbox]').item(0).dispatchEvent(new Event('change'))"
        )
        # sets the templates option to pdb70
        self.web_page.runJavaScript(
            f"document.querySelectorAll('select').item(0).selectedIndex = {self.prediction_params['templates'][1]}"
        )
        self.web_page.runJavaScript(
            "document.querySelectorAll('select').item(0).dispatchEvent(new Event('change'))"
        )
        # runs all cells
        self.web_page.runJavaScript("colab.global.notebook.runAll()")
        # confirm to "run anyway" -> reason: notebook from GitHub
        self.web_page.runJavaScript("document.getElementsByTagName('paper-button').item(2).click()")

    def finished_prediction(self):
        """This function does steps to finish the prediction with the google colab notebook

        """
        gui_utils.warning_prediction_is_finished(self)

    def check_status(self) -> int:
        """This function checks the state of the Google colab kernel and returns the exit code

        Returns:
            exit code: defines in which state the kernel is
        """
        self.web_page.runJavaScript("""
                    if (colab.global.notebook.kernel.getState() === "kernel idle") {
                        document.title = "error";
                    } else if (colab.global.notebook.kernel.getState() === "allocating") {
                        document.title = "allocating";
                    } else if (colab.global.notebook.kernel.getState() === "busy") {
                        document.title = "normal";
                    } else if (colab.global.notebook.kernel.getState() === "connect") {
                        document.title = "need_to_connect";
                    }
                    """)
        if self.web_page.title() == "error":
            self.exit_code = 1
            self.close()
        elif self.web_page.title() == "allocating":
            self.exit_code = 3
        elif self.web_page.title() == "normal":
            self.exit_code = 2
        elif self.web_page.title() == "need_to_connect":
            self.exit_code = 4
        return self.exit_code

    def control_status(self):
        exit_code = self.check_status()
        current_time = datetime.datetime.now().strftime("%H:%M:%S")
        if exit_code == 2:
            self.label.setText(f"Prediction runs normal. Prediction started at: {self.current_time} and last checked: {current_time}")
        elif exit_code == 3:
            self.label.setText(f"Google Colab resources need to be allocated. Prediction started at: {self.current_time} and last checked: {current_time}")
        elif exit_code == 4:
            self.label.setText(f"You have to reconnect the session! Prediction started at: {self.current_time} and last checked: {current_time}")

    def set_protein_sequence(self, sequence):
        self.prediction_params.update(seq=sequence)

    def set_job_name(self, jobname):
        self.prediction_params.update(job_name=jobname)

    def toggle_notebook_view(self):
        if self.browser.isVisible():
            self.browser.hide()
            self.setMaximumSize(575, 150)
        else:
            self.browser.show()
            self.setMaximumSize(1080, 900)

    def show_interface(self):
        """This function is called from outside the class to initiate the opening of the interface

        """
        self.open_colab_notebook()
        self.exec_()

    def show_account_page(self):
        self.web_page.setUrl(self.account_page)
        self.exec_()

    def disconnect_active_runtime(self):
        self.web_page.runJavaScript("colab.global.notebook.kernel.disconnect()")
        self.close()
