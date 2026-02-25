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
"""Module for the HelpPanel class.

This module defines the HelpPanel, which displays available help topics
organized by category in a tree view. Users can click on any help topic
to open the corresponding help dialog.

Authors: Martin Urban, Hannah Kullik

Version: 1.4.0
"""
from src.pyssa.gui.ui.styles.icon_manager import IconManager
from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.qt import QtWidgets
from src.pyssa.gui.qt import QtGui
from src.pyssa.gui.ui.views import base_side_panel
from src.pyssa.gui.ui.custom_widgets import quick_access_bar_action, quick_access_bar, dropdown_menu

__docformat__ = "google"


class HelpPanel(base_side_panel.BaseSidePanel):
    """Panel for displaying available help topics.

    This panel provides a tree view of all registered help topics organized
    by category. Users can click on any topic to open the help dialog.
    """

    helpDialogRequested = QtCore.pyqtSignal(str)
    """Signal emitted when a help dialog is requested. Contains the dialog ID."""

    # <editor-fold desc="Constructor">
    def __init__(self) -> None:
        """Initializes the HelpPanel.

        Sets up the UI components, including the tree view for help topics.
        """
        super().__init__("Help")
        # <editor-fold desc="Instance attributes">
        self.container_widget: QtWidgets.QWidget = QtWidgets.QWidget()
        self.help_text_browser: QtWidgets.QTextBrowser = QtWidgets.QTextBrowser()
        # </editor-fold>
        self._setup_panel_extra_ui()
        self.setMinimumWidth(250)
    # </editor-fold>

    # <editor-fold desc="Public methods">
    def populate_help_topics(self, categories: dict) -> None:
        """Populates the tree view with help topics organized by category.

        Args:
            categories: Dictionary mapping category names to lists of HelpDialog objects
        """
        self.help_tree_view.clear()

        for category, dialogs in sorted(categories.items()):
            category_item = QtWidgets.QTreeWidgetItem([category])
            font = QtGui.QFont("Arial", 10)
            font.setBold(True)
            category_item.setFont(0, font)
            self.help_tree_view.addTopLevelItem(category_item)

            for dialog in sorted(dialogs, key=lambda d: d.title):
                dialog_item = QtWidgets.QTreeWidgetItem([dialog.title])
                dialog_item.setData(0, QtCore.Qt.UserRole, dialog.id)
                if dialog.description:
                    dialog_item.setToolTip(0, dialog.description)
                category_item.addChild(dialog_item)

            category_item.setExpanded(True)

    def set_help_panel_header(self, header_text: str) -> None:
        """Sets the header text for the help panel."""
        self.lbl_header.setText(f"Help \u2014 {header_text}")

    def filter_help_topics(self, search_text: str) -> None:
        """Filters the tree view based on the search text.

        Args:
            search_text: Text to filter the help topics
        """
        search_lower = search_text.lower()

        for i in range(self.help_tree_view.topLevelItemCount()):
            category_item = self.help_tree_view.topLevelItem(i)
            category_visible = False

            for j in range(category_item.childCount()):
                child_item = category_item.child(j)
                text = child_item.text(0).lower()
                tooltip = child_item.toolTip(0).lower()

                if search_lower in text or search_lower in tooltip:
                    child_item.setHidden(False)
                    category_visible = True
                else:
                    child_item.setHidden(True)

            category_item.setHidden(not category_visible)
    # </editor-fold>

    # <editor-fold desc="Private methods">
    def _setup_tree_view(self) -> None:
        """Configures the help tree view."""
        self.help_tree_view.setHeaderLabel("Help Topics")
        self.help_tree_view.setAlternatingRowColors(True)
        self.help_tree_view.itemClicked.connect(self._on_item_clicked)

    def _setup_search_box(self) -> None:
        """Configures the search box."""
        self.search_line_edit.setPlaceholderText("Search help topics...")
        self.search_line_edit.textChanged.connect(self.filter_help_topics)
        self.search_line_edit.setClearButtonEnabled(True)

    def _on_item_clicked(self, item: QtWidgets.QTreeWidgetItem, column: int) -> None:
        """Handles tree item click events.

        Args:
            item: The clicked tree widget item
            column: The column that was clicked
        """
        dialog_id = item.data(0, QtCore.Qt.UserRole)
        if dialog_id:
            self.helpDialogRequested.emit(dialog_id)

    def _setup_extra_styles(self) -> None:
        """Applies stylesheets to panel widgets."""
        self.help_tree_view.setStyleSheet(
            """
            QTreeWidget {
                border-style: solid;
                border-width: 0.15em;
                border-radius: 0.3em;
                border-color: #DCDBE3;
                background-color: white;
            }
            QTreeWidget::item:hover {
                background-color: #E8F4F8;
            }
            QTreeWidget::item:selected {
                background-color: #0078D7;
                color: white;
            }
            """
        )

    def _setup_container_widget(self) -> None:
        """Initializes and configures the container widget."""
        tmp_container_widget_layout: QtWidgets.QVBoxLayout = QtWidgets.QVBoxLayout()
        tmp_container_widget_layout.setContentsMargins(5, 5, 5, 5)
        tmp_container_widget_layout.setSpacing(5)

        # Only show text browser for hover help
        self.help_text_browser.setStyleSheet("""
            QTextBrowser {
                border: 1px solid #CCC;
                border-radius: 3px;
                padding: 10px;
                background-color: white;
                font-size: 11pt;
            }
        """)
        tmp_container_widget_layout.addWidget(self.help_text_browser)

        self.container_widget.setLayout(tmp_container_widget_layout)

    def _setup_panel_extra_ui(self) -> None:
        """Sets up additional UI elements specific to this panel."""
        self._setup_container_widget()
        self.layout_content_frame.addWidget(self.container_widget)

    # </editor-fold>
