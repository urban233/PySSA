import time
from pymol import Qt
from PyQt5.QtWidgets import (QApplication, QMainWindow,
                             QHBoxLayout, QVBoxLayout,
                             QPushButton, QWidget)
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtCore import QUrl
from PyQt5.QtWebEngineWidgets import QWebEngineView, QWebEnginePage, QWebEngineDownloadItem
from pyssa.gui.utilities import constants
from pyssa.gui.utilities import gui_utils
from pyssa.gui.ui.dialogs import dialog_notebook_managment


class WebInterface(Qt.QtWidgets.QDialog):
    google_login_page = "https://accounts.google.com/v3/signin/identifier?dsh=S487182181%3A1668865110328278&continue=https%3A%2F%2Fcolab.research.google.com%2Fgithub%2Fsokrypton%2FColabFold%2Fblob%2Fmain%2FAlphaFold2.ipynb&ec=GAZAqQM&passive=true&flowName=GlifWebSignIn&flowEntry=ServiceLogin&ifkv=ARgdvAu0jNh8qAP4PLZtLHO1saEbt8gL7tdnPxJNzq-nQ8HD1od95scPNdlXBiSB5aGTVZWwLQ-f-A"
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
        # QWebEngineView is used to display the content from the QWebEnginePage
        self.browser = QWebEngineView()
        self.web_page = QWebEnginePage()
        self.browser.setPage(self.web_page)
        # basic layout
        horizontal_layout = QHBoxLayout(self)
        horizontal_layout.addWidget(self.browser)
        # Connections
        self.browser.page().profile().downloadRequested.connect(self.on_downloadRequested)
        # window settings
        self.setMinimumSize(450, 450)
        self.setWindowTitle("Sign-In Google")
        self.show()

    def open_login_page(self):
        """This function opens the Google login page

        """
        url = QUrl(self.google_login_page)
        self.web_page.setUrl(url)
        self.web_page.loadFinished.connect(self.run_js_code)

    def run_js_code(self):
        """This function calls the javascript which is needed to automate the colabfold
        google colab notebook

        """
        self._page_loaded()
        # sets the sequence
        test = f"document.querySelectorAll('select').item(0).selectedIndex = '{self.prediction_params['templates'][1]}'"
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
        if self.web_page.url() == self.url_without_pre_login or self.web_page.url() == self.url_with_pre_login:
            self.setWindowTitle("Prediction started ...")
            #self.hide()
            # sets up a dialog window which is used to control the behaviour of the QWebEnginePage
            self.dialog = dialog_notebook_managment.DialogNotebookManagment()
            self.dialog.web_page = self.web_page
            self.dialog.check_status = self.check_status
            self.dialog.interface = self
            self.dialog.exec_()
            self.close()

    @QtCore.pyqtSlot("QWebEngineDownloadItem*")
    def on_downloadRequested(self, download):
        """This function is called if the notebook wants to download the results

        Args:
            download:
                necessary argument which holds information about the requested download
        """
        path = f"{constants.SCRATCH_DIR}/prediction.zip"
        download.setPath(path)
        download.accept()
        download.finished.connect(self.finished_prediction)
        self.exit_code = 0

    def finished_prediction(self):
        """This function does steps to finish the prediction with the google colab notebook

        """
        self.dialog.close()
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

    def _page_loaded(self):
        """This function is a helper function to indicate that a web page is loaded

        """
        print("Page loaded")
        print(self.web_page.url())

    def show_interface(self):
        """This function is called from outside the class to initiate opening of the interface

        """
        self.open_login_page()
        self.exec_()

    def show_account_page(self):
        self.web_page.setUrl(self.account_page)
        self.exec_()

    def get_exit_code(self):
        """This function simply gets the exit code from the class

        """
        return self.exit_code

    def set_protein_sequence(self, sequence):
        self.prediction_params.update(seq=sequence)

    def set_job_name(self, jobname):
        self.prediction_params.update(job_name=jobname)
