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
"""Module for the MoleculeObjectsPanel class.

This module defines the MoleculeObjectsPanel, which is responsible for
displaying all molecule objects associated with a project in the PySSA
application. The panel provides a tree view and controls for expanding and
collapsing the structure hierarchy.

Authors: Martin Urban, Hannah Kullik

Version: 1.4.0

TODO: Work in progress, this panel is not yet functional.
"""
from src.pyssa.gui.ui.styles.icon_manager import IconManager
from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.qt import QtWidgets
from src.pyssa.gui.qt import QtGui
from src.pyssa.gui.ui.views import base_side_panel
from src.pyssa.gui.ui.custom_widgets import quick_access_bar_action, quick_access_bar, dropdown_menu

__docformat__ = "google"


class PySSAObjectsPanel(base_side_panel.BaseSidePanel):
    """Panel for displaying all molecule objects of a project.

    This panel provides a tree view for molecule object objects and
    includes controls for expanding and collapsing the entire tree.
    """

    # <editor-fold desc="Constructor">
    def __init__(self) -> None:
        """Initializes the MoleculeObjectsPanel.

        Sets up the UI components, including the tree view and expand/collapse
        controls, and applies style and icons.
        """
        super().__init__("PySSA Objects")
        # <editor-fold desc="Instance attributes">
        self.container_widget: QtWidgets.QWidget = QtWidgets.QWidget()
        self.expand_all: QtWidgets.QPushButton = QtWidgets.QPushButton()
        self.collapse_all: QtWidgets.QPushButton = QtWidgets.QPushButton()
        self.collapse_expand_layout: QtWidgets.QHBoxLayout = QtWidgets.QHBoxLayout()
        self.tree_view: QtWidgets.QTreeView = QtWidgets.QTreeView()
        # </editor-fold>
        self._setup_panel_extra_ui()
    # </editor-fold>

    def get_toolbar(self) -> "quick_access_bar.QuickAccessBar":
        return self._quick_access_bar

    # <editor-fold desc="Private methods">
    def _setup_extra_styles(self) -> None:
        """Applies stylesheets to panel widgets.

        Sets the style for expand/collapse buttons and the tree view
        according to the current theme.
        """
        self.tree_view.setStyleSheet(
            """
            border-style: solid;
            border-width: 0.15em;
            border-radius: 0.3em;
            border-color: #DCDBE3;
            background-color: white;
            """
        )

    def _set_icons(self) -> None:
        """Assigns icons to the expand and collapse buttons.

        Uses the IconManager singleton to set appropriate icons for
        the expand_all and collapse_all buttons.
        """
        IconManager.instance().set_icon(
            self.expand_all,
            IconManager.Icons.EXPAND_ALL,
            the_size=QtCore.QSize(16, 16),
        )
        IconManager.instance().set_icon(
            self.collapse_all,
            IconManager.Icons.COLLAPSE_ALL,
            the_size=QtCore.QSize(16, 16),
        )

    def _setup_top_toolbar(self):
        self.expand_all = quick_access_bar_action.QuickAccessBarAction(
            "Expand All", "left", 0, None, IconManager.instance().get_icon(IconManager.Icons.EXPAND_ALL)
        )
        self.collapse_all = quick_access_bar_action.QuickAccessBarAction(
            "Collapse All", "left", 0, None, IconManager.instance().get_icon(IconManager.Icons.COLLAPSE_ALL)
        )
        self.import_file_action = quick_access_bar_action.QuickAccessBarAction(
            "Import File", "left", 0, None, IconManager.instance().get_icon(IconManager.Icons.UPLOAD_FILE)
        )
        self.add_sequence_action = quick_access_bar_action.QuickAccessBarAction(
            "Add Sequence", "left", 0, None, IconManager.instance().get_icon(IconManager.Icons.NOTE_ADD)
        )
        self.export_file_action = quick_access_bar_action.QuickAccessBarAction(
            "Export File", "left", 1, None, IconManager.instance().get_icon(IconManager.Icons.FILE_SAVE)
        )
        self.delete_object_action = quick_access_bar_action.QuickAccessBarAction(
            "Delete Object", "left", 1, None, IconManager.instance().get_icon(IconManager.Icons.SCAN_DELETE)
        )

        # Horizontal quick access bar above the editor
        self._quick_access_bar = quick_access_bar.QuickAccessBar([
            # self.expand_all, self.collapse_all,
            self.import_file_action, self.add_sequence_action, self.export_file_action, self.delete_object_action
        ], horizontal=True, button_size=(24, 24))

    def _setup_import_popup(self):
        self.import_seq_prot_menu = dropdown_menu.DropDownMenu()
        self.import_seq_action = QtGui.QAction("Sequence")
        self.import_prot_action = QtGui.QAction("Protein")
        self.import_seq_prot_menu.addAction(self.import_seq_action)
        self.import_seq_prot_menu.addAction(self.import_prot_action)

    def _setup_expand_collapse_header(self) -> None:
        """Configures the header layout for expand/collapse controls.

        Adds the expand_all and collapse_all buttons to the layout and
        ensures proper spacing and alignment.
        """
        self.collapse_expand_layout.setContentsMargins(5, 0, 0, 0)
        self.collapse_expand_layout.addWidget(self.expand_all)
        self.collapse_expand_layout.addWidget(self.collapse_all)
        self.collapse_expand_layout.addStretch()

    def _setup_tree_view(self) -> None:
        """Configures the tree view widget.

        Sets margins, hides the header, enables extended selection, and
        disables editing for the tree view.
        """
        self.tree_view.setContentsMargins(0, 0, 0, 0)
        self.tree_view.setHeaderHidden(True)
        self.tree_view.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        self.tree_view.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)

    def _setup_container_widget(self) -> None:
        """Initializes and configures the container widget.

        Creates a vertical layout, adds the expand/collapse controls and
        the tree view, and sets the layout for the container widget.
        """
        tmp_container_widget_layout: QtWidgets.QVBoxLayout = QtWidgets.QVBoxLayout()
        tmp_container_widget_layout.setContentsMargins(0, 0, 0, 0)
        tmp_container_widget_layout.addWidget(self._quick_access_bar)
        tmp_container_widget_layout.addLayout(self.collapse_expand_layout)
        tmp_container_widget_layout.addWidget(self.tree_view)
        tmp_container_widget_layout.setSpacing(0)
        self.container_widget.setLayout(tmp_container_widget_layout)

    def _setup_panel_extra_ui(self) -> None:
        """Sets up additional UI elements specific to this panel.

        Calls setup methods for the expand/collapse header, tree view,
        icons, style, and container widget, and adds the container to
        the main layout.
        """
        self._setup_top_toolbar()
        self._setup_import_popup()
        # self._setup_expand_collapse_header()
        self._setup_tree_view()
        # self._set_icons()
        self._setup_extra_styles()
        self._setup_container_widget()
        # This panel should always be shown to make it easier for the end-user
        self.btn_close.hide()
        self.layout_content_frame.addWidget(self.container_widget)

    # </editor-fold>
