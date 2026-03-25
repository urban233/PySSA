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


class ProjectOverviewPanel(base_side_panel.BaseSidePanel):
    """Panel for displaying all relevant project information."""

    # <editor-fold desc="Constructor">
    def __init__(self) -> None:
        """Constructor."""
        super().__init__("")
        # <editor-fold desc="Instance attributes">
        self.container_widget: QtWidgets.QWidget = QtWidgets.QWidget()
        self.lbl_project_name: QtWidgets.QLabel = QtWidgets.QLabel("Project Name: ")
        self.lbl_session_name: QtWidgets.QLabel = QtWidgets.QLabel("Session Name: ")
        self.lbl_scene_name: QtWidgets.QLabel = QtWidgets.QLabel("Scene Name: ")
        # </editor-fold>
        self._setup_panel_extra_ui()
        self.global_layout.removeItem(self.layout_header)
    # </editor-fold>

    def set_project_name(self, project_name: str) -> None:
        self.lbl_project_name.setText(f"Project Name: {project_name}")

    def set_session_name(self, session_name: str) -> None:
        self.lbl_session_name.setText(f"Session Name: {session_name}")

    def set_scene_name(self, scene_name: str) -> None:
        self.lbl_scene_name.setText(f"Scene Name: {scene_name}")

    # <editor-fold desc="Private methods">
    def _setup_container_widget(self) -> None:
        """Initializes and configures the container widget.

        Creates a vertical layout, adds the expand/collapse controls and
        the tree view, and sets the layout for the container widget.
        """
        tmp_container_widget_layout: QtWidgets.QVBoxLayout = QtWidgets.QVBoxLayout()
        tmp_container_widget_layout.setContentsMargins(0, 0, 0, 0)
        tmp_container_widget_layout.addWidget(self.lbl_project_name)
        tmp_container_widget_layout.addWidget(self.lbl_session_name)
        tmp_container_widget_layout.addWidget(self.lbl_scene_name)
        tmp_container_widget_layout.setSpacing(0)
        self.container_widget.setLayout(tmp_container_widget_layout)

    def _setup_panel_extra_ui(self) -> None:
        """Sets up additional UI elements specific to this panel.

        Calls setup methods for the expand/collapse header, tree view,
        icons, style, and container widget, and adds the container to
        the main layout.
        """
        self._setup_container_widget()
        # This panel should always be shown to make it easier for the end-user
        self.btn_close.hide()
        self.layout_content_frame.addWidget(self.container_widget)

    # </editor-fold>
