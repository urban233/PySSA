from PyQt5.QtWidgets import QApplication, QDialog, QVBoxLayout, QPushButton, QHBoxLayout, QWidget
from PyQt5.QtWebEngineWidgets import QWebEngineView
import sys


class HelpDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.initUI()

    def initUI(self):
        self.setGeometry(100, 100, 700, 800)
        self.setWindowTitle('My PyQt Help System')
        central_layout = QVBoxLayout(self)
        self.browser = QWebEngineView()
        central_layout.addWidget(self.browser)

    def loadHelpContent(self, content_path):
        with open(content_path, 'r', encoding='utf-8') as file:
            html_content = file.read()

        self.browser.setHtml(html_content)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    helpDialog = HelpDialog()
    helpDialog.loadHelpContent(r'C:\Users\martin\github_repos\PySSA\docs\internal_help\html\home.html')

    # You can show the dialog using exec_() or show()
    # helpDialog.exec_()  # If you want a modal dialog
    helpDialog.show()

    sys.exit(app.exec_())


