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
import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QTextBrowser

from pyssa.util import constants


class HTMLViewer(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("HTML Viewer")
        self.setGeometry(100, 100, 800, 600)

        self.text_browser = QTextBrowser(self)
        self.text_browser.setGeometry(10, 10, 780, 580)

        # Load the static HTML content
        html_content = """
        <html>
            <head>
                <style>
                    h1 {
                        color: blue;
                    }
                    p {
                        font-size: 18px;
                    }
                </style>
            </head>
            <body>
                <h1>Hello, World!</h1>
                <p>This is a static HTML page styled with CSS.</p>
            </body>
        </html>
        """
        # Load an existing HTML file
        with open(f'{constants.DOCS_HTML}', 'r', encoding='utf-8') as file:
            html_content = file.read()
        self.text_browser.setHtml(html_content)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = HTMLViewer()
    window.show()
    sys.exit(app.exec_())
