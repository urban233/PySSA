#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
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
"""Module for the base side panel class in PySSA.

This module defines the BaseSidePanel class, which provides a customizable
side panel widget for the PySSA main window. The panel includes a header,
close button, and a content frame for adding custom widgets.

Authors: Martin Urban, Hannah Kullik

Version: 1.4.0
"""
from src.pyssa.gui.ui.styles.icon_manager import IconManager
from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.qt import QtWidgets


__docformat__ = "google"


class BaseSidePanel(QtWidgets.QWidget):
    """Class to create a custom side panel for the main window.

    Attributes:
        panelClosed: Signal emitted when the panel is closed.
        panelOpened: Signal emitted when the panel is opened.
        lbl_header: Label displaying the panel title.
        btn_close: Button to close the panel.
        content_frame: Frame to store widgets for the panel.
        layout_header: Layout for the panel header.
        layout_content_frame: Layout for the content frame.
        global_layout: Layout for the entire panel.
    """

    # <editor-fold desc="Class attributes">
    panelClosed = QtCore.pyqtSignal()
    panelOpened = QtCore.pyqtSignal()
    # </editor-fold>

    # <editor-fold desc="Constructor">
    def __init__(self, a_title: str) -> None:
        """Initializes a new BaseSidePanel instance.

        Args:
            a_title: The title of the panel.

        Raises:
            ValueError: If `a_title` is None.
        """
        # <editor-fold desc="Checks">
        if a_title is None:
            raise ValueError("a_title cannot be None.")
        # </editor-fold>
        super().__init__()
        # <editor-fold desc="Instance attributes">
        self.lbl_header: QtWidgets.QLabel = QtWidgets.QLabel(a_title)
        self.btn_close: QtWidgets.QPushButton = QtWidgets.QPushButton()
        self.content_frame: QtWidgets.QFrame = QtWidgets.QFrame()
        self.layout_header: QtWidgets.QHBoxLayout = QtWidgets.QHBoxLayout()
        self.layout_content_frame: QtWidgets.QVBoxLayout = QtWidgets.QVBoxLayout()
        self.global_layout: QtWidgets.QVBoxLayout = QtWidgets.QVBoxLayout()
        # </editor-fold>
        self._setup_panel_ui()
        self._connect_all_signals()

    # </editor-fold>

    # <editor-fold desc="Public methods">
    def add_global_stretch(self) -> None:
        """Adds a stretch to the bottom of the global layout."""
        self.global_layout.addStretch()

    def hide_panel(self) -> None:
        """Emits the panelClosed signal to indicate the panel should be closed."""
        self.panelClosed.emit()

    def show_panel(self) -> None:
        """Emits the panelOpened signal to indicate the panel should be shown."""
        self.panelOpened.emit()

    # </editor-fold>

    # <editor-fold desc="Private methods">
    def _connect_all_signals(self) -> None:
        """Connects all signals with their appropriate slot method."""
        self.btn_close.clicked.connect(self.hide_panel)

    def _setup_style(self) -> None:
        """Sets up any stylesheets needed for the panel."""
        self.btn_close.setStyleSheet(
            """
            QPushButton {
                background-color: rgba(220, 219, 227, 0.01);
                border: none;
                border-radius: 0.3em;
                min-width: 2.1em;
                max-width: 2.1em;
                min-height: 2.1em;
                max-height: 2.1em;
            }
            QPushButton::hover {
                background-color: rgba(220, 219, 227, 0.5);
                border: none;
                min-width: 2.1em;
                max-width: 2.1em;
                min-height: 2.1em;
                max-height: 2.1em;
            }
            """
        )
        self.lbl_header.setStyleSheet(
            """
            QLabel {
                font-size: 12pt;
            }
            """
        )

    def _setup_header(self) -> None:
        """Sets up the panel header layout and widgets."""
        self.layout_header.setContentsMargins(0, 0, 0, 0)
        self.layout_header.addWidget(self.lbl_header)
        self.layout_header.addStretch()
        self.layout_header.addWidget(self.btn_close)

    def _setup_content_frame(self) -> None:
        """Sets up the main content frame and its layout."""
        self.content_frame.setContentsMargins(0, 0, 0, 0)
        self.layout_content_frame.setContentsMargins(0, 0, 0, 0)
        self.content_frame.setLayout(self.layout_content_frame)

    def _setup_panel_ui(self) -> None:
        """Sets up the complete panel UI, including header, content, and style."""
        self._setup_header()
        self._setup_content_frame()
        self.global_layout.setContentsMargins(0, 0, 0, 0)
        self.global_layout.addLayout(self.layout_header)
        self.global_layout.addWidget(self.content_frame)
        self.setLayout(self.global_layout)
        self.setContentsMargins(0, 0, 0, 0)
        self._setup_style()
        IconManager.instance().set_icon(
            self.btn_close, IconManager.Icons.CLOSE, the_size=QtCore.QSize(24, 24)
        )

    # </editor-fold>
