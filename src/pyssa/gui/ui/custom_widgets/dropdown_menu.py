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
"""DropDownMenu widget for persistent menu display in the PySSA frontend.

This module provides the DropDownMenu class, a custom QMenu that remains open
after clicking on checkable actions. This behavior is useful for menus containing
multiple checkable options where users need to select multiple items without
the menu closing after each selection.

The DropDownMenu automatically applies the current theme's dropdown menu styling
and provides enhanced user interaction for multi-selection scenarios.

Authors: Martin Urban, Hannah Kullik

Version: 1.4.0
"""

from typing import Optional

from src.pyssa.gui.qt import QtWidgets, QtCore

__docformat__ = "google"


class DropDownMenu(QtWidgets.QMenu):
    """A QMenu that remains open after clicking on checkable actions.

    This custom menu widget extends QMenu to provide persistent display behavior
    when users interact with checkable actions. Unlike the standard QMenu which
    closes after any action click, this menu only closes when non-checkable
    actions are clicked or when the user clicks outside the menu.

    The menu automatically applies the current theme's dropdown menu styling
    for consistent appearance with the rest of the application.

    Attributes:
        Inherits all attributes from QtWidgets.QMenu.
    """

    # <editor-fold desc="Constructor">
    def __init__(self, a_parent: Optional[QtWidgets.QWidget] = None) -> None:
        """Initialize the DropDownMenu with persistent behavior.

        Creates a new DropDownMenu that will remain open when checkable actions
        are clicked. The menu is automatically styled using the current theme's
        dropdown menu style.

        Args:
            a_parent: The parent widget for this menu. Can be None for a
                     top-level menu.

        Note:
            The menu will automatically apply the current theme's styling
            upon initialization.
        """
        super().__init__(a_parent)
       
    # </editor-fold>

    # <editor-fold desc="Public methods">
    def mouseReleaseEvent(self, a_event: QtCore.QEvent) -> None:
        """Handle mouse release events with persistent menu behavior.

        Overrides the default QMenu behavior to prevent the menu from closing
        when checkable actions are clicked. For checkable actions, the method
        toggles the action's checked state and prevents menu closure. For
        non-checkable actions, the default behavior is preserved.

        Args:
            a_event: The mouse release event containing position and button information.
                    Must not be None.

        Raises:
            ValueError: If a_event is None (via psa_comm_api validation).

        Note:
            This method enables multi-selection workflows by keeping the menu
            open for checkable actions while maintaining normal behavior for
            regular actions.
        """
        tmp_action = self.actionAt(a_event.pos())

        if tmp_action and tmp_action.isCheckable():
            group = tmp_action.actionGroup()
            if group is not None and group.isExclusive():
                # For exclusive groups, never allow deselecting the currently selected action.
                # If the clicked action is not checked, check it (Qt will uncheck others automatically).
                if not tmp_action.isChecked():
                    tmp_action.setChecked(True)
                # Keep menu open
                return
            else:
                # Non-exclusive: toggle normally and keep the menu open
                tmp_action.setChecked(not tmp_action.isChecked())
                return

        super().mouseReleaseEvent(a_event)

    # </editor-fold>
