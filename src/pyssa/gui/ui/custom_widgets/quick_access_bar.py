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
"""QuickAccessBar widget for rapid action access in the PySSA frontend.

This module provides the QuickAccessBar class, a customizable toolbar widget
that displays frequently used actions as icon-only buttons. The quick access bar
is designed for efficient user workflows by providing immediate access to
commonly used functionality without navigating through menus.

The QuickAccessBar automatically manages button creation, styling, and layout
for a collection of QuickAccessBarAction objects. Each action is represented
as a compact tool button with consistent sizing and theming.

Authors: Martin Urban, Hannah Kullik

Version: 1.4.0
"""

from typing import Optional, Dict

from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.qt import QtWidgets
from src.pyssa.gui.ui.custom_widgets import quick_access_bar_action

__docformat__ = "google"


class QuickAccessBar(QtWidgets.QWidget):
    """A customizable toolbar widget for rapid access to frequently used actions.

    The QuickAccessBar provides a compact, icon-only interface for commonly used
    functionality. It automatically creates and manages tool buttons for a collection
    of QuickAccessBarAction objects, applying consistent styling and layout.

    The widget uses a box layout (vertical by default, horizontal optional) with
    fixed-size buttons and automatically applies the current theme's quick access
    bar styling. The layout is split into two groups (start and end) separated by
    a central stretch, leaving free space in the middle. Items can be added to the
    top/left (start) or bottom/right (end) side as needed.

    Attributes:
        quick_access_bar_actions: List of QuickAccessBarAction objects that define
                                 the actions displayed in this toolbar.
        item_to_button: Dictionary mapping QuickAccessBarAction objects to their
                       corresponding QToolButton widgets for easy retrieval.
    """

    # <editor-fold desc="Constructor">
    def __init__(
            self,
            the_quick_access_bar_actions: list["quick_access_bar_action.QuickAccessBarAction"],
            horizontal: bool = False,
            button_size: tuple[int, int] = (30, 30)
    ) -> None:
        """Initialize the QuickAccessBar with the specified actions.

        Creates a new QuickAccessBar widget and populates it with tool buttons
        for each provided QuickAccessBarAction. The widget is automatically
        styled and laid out according to the current theme.

        Args:
            the_quick_access_bar_actions: List of QuickAccessBarAction objects
                                       that define the actions to display.
                                       Must not be empty for meaningful usage.

        Note:
            The widget automatically applies the current theme's quick access
            bar styling and creates a mapping between actions and buttons for
            easy retrieval.
        """
        super().__init__()

        # <editor-fold desc="Instance attributes">
        self._button_size: tuple[int, int] = button_size
        self._is_horizontal: bool = horizontal
        if horizontal:
            self._layout: QtWidgets.QHBoxLayout = QtWidgets.QHBoxLayout()
            self._start_layout: QtWidgets.QHBoxLayout = QtWidgets.QHBoxLayout()
            self._end_layout: QtWidgets.QHBoxLayout = QtWidgets.QHBoxLayout()
        else:
            self._layout: QtWidgets.QVBoxLayout = QtWidgets.QVBoxLayout()
            self._start_layout: QtWidgets.QVBoxLayout = QtWidgets.QVBoxLayout()
            self._end_layout: QtWidgets.QVBoxLayout = QtWidgets.QVBoxLayout()
        self.quick_access_bar_actions: list["quick_access_bar_action.QuickAccessBarAction"] = (
            the_quick_access_bar_actions
        )
        self.item_to_button: Dict[
            "quick_access_bar_action.QuickAccessBarAction", QtWidgets.QToolButton
        ] = {}
        # </editor-fold>

        self._init_widget()

        # </editor-fold>

    # <editor-fold desc="Public methods">
    def get_tool_button_for_action(
        self, a_action: "quick_access_bar_action.QuickAccessBarAction"
    ) -> Optional[QtWidgets.QToolButton]:
        """Retrieve the tool button associated with a specific QuickAccessBarAction.

        Looks up the QToolButton widget that corresponds to the given
        QuickAccessBarAction in the internal mapping dictionary. This is useful
        for programmatically accessing buttons to modify their properties or
        connect additional signals.

        Args:
            a_action: The QuickAccessBarAction object for which to retrieve
                     the corresponding tool button.

        Returns:
            The QToolButton widget associated with the specified action,
            or None if no button could be found for the given action.

        Note:
            The mapping is established during widget initialization when
            buttons are created for each action.

        Example:
            button = quick_bar.get_tool_button_for_action(save_action)
            if button:
                button.setEnabled(False)  # Disable the save button
        """
        return self.item_to_button.get(a_action, None)

    def add_action(self, a_quick_access_bar_action, position: str = "top"):
        """Add a single action button to the bar.

        By default, items are added to the top (for vertical bars) or left (for horizontal bars).
        You can also add items to the opposite side (bottom/right). A stretch between both sides
        keeps free space in the middle.

        Args:
            a_quick_access_bar_action: The QuickAccessBarAction to add.
            position: Where to add the item: one of "top", "bottom", "start", "end",
                      "left", "right". Defaults to "top" (aka start/left).
        """
        # Normalize position
        pos = (position or "top").strip().lower()
        start_aliases = {"top", "start", "left"}
        end_aliases = {"bottom", "end", "right"}
        if pos in start_aliases:
            target_layout = self._start_layout
        elif pos in end_aliases:
            target_layout = self._end_layout
        else:
            # Fallback to start if invalid input
            target_layout = self._start_layout

        tmp_button = QtWidgets.QToolButton()
        tmp_button.setDefaultAction(a_quick_access_bar_action.get_action())
        tmp_button.setToolButtonStyle(QtCore.Qt.ToolButtonStyle.ToolButtonIconOnly)
        tmp_button.setIconSize(QtCore.QSize(self._button_size[0], self._button_size[1]))
        tmp_button.setFixedSize(self._button_size[0], self._button_size[1])
        target_layout.addWidget(tmp_button)
        self.item_to_button[a_quick_access_bar_action] = tmp_button

    def add_top_action(self, a_quick_access_bar_action):
        """Convenience: add action to the top (or left for horizontal bars)."""
        self.add_action(a_quick_access_bar_action, position="top")

    def add_bottom_action(self, a_quick_access_bar_action):
        """Convenience: add action to the bottom (or right for horizontal bars)."""
        self.add_action(a_quick_access_bar_action, position="bottom")

    # </editor-fold>

    # <editor-fold desc="Private methods">
    def _init_widget(self) -> None:
        """Initialize the widget with actions, layout, and styling.

        Sets up the complete widget by adding all actions as tool buttons,
        configuring the layout with zero margins, and applying the current
        theme's quick access bar styling. This method is called automatically
        during widget construction.

        Note:
            This method coordinates the widget setup process by calling
            helper methods in the correct order and applying final styling.
        """
        # Configure margins/spacings
        self._layout.setContentsMargins(0, 0, 0, 0)
        self._start_layout.setContentsMargins(0, 0, 0, 0)
        self._end_layout.setContentsMargins(0, 0, 0, 0)

        # Assemble the main layout: start layout | stretch | end layout
        self._layout.addLayout(self._start_layout)
        self._layout.addStretch(1)
        self._layout.addLayout(self._end_layout)
        self.setLayout(self._layout)

        # Populate default actions on the start side
        self._add_actions()

        # Ensure the bar only grows in the intended direction
        if self._is_horizontal:
            self.setSizePolicy(
                QtWidgets.QSizePolicy.Policy.Expanding,
                QtWidgets.QSizePolicy.Policy.Fixed,
            )
        else:
            self.setSizePolicy(
                QtWidgets.QSizePolicy.Policy.Fixed,
                QtWidgets.QSizePolicy.Policy.Expanding,
            )
        self.setStyleSheet(
            """
            QToolButton {
                font-size: 7.9pt;
                background-color: white;
                padding: 0.15em;
                border: none;
                border-radius: 0.375em;
            }
            
            QToolButton::hover {
                background: #f5f5f5;
                color: black;
                border-radius: 0.375em;
            }
        """
        )

    def _add_actions(self) -> None:
        """Create and add tool buttons for all QuickAccessBarAction objects.

        Iterates through the list of QuickAccessBarAction objects and creates
        a corresponding QToolButton for each one. Each button is configured
        with icon-only display and the configured fixed size, and is added to both
        the start layout (top/left) and the action-to-button mapping dictionary.

        Note:
            The main layout contains a central stretch separating start and end groups,
            so the middle area remains free. By default, initial actions populate the
            start group.
        """
        for tmp_quick_access_bar_action in self.quick_access_bar_actions:
            tmp_button = QtWidgets.QToolButton()
            tmp_button.setDefaultAction(tmp_quick_access_bar_action.get_action())
            tmp_button.setToolButtonStyle(QtCore.Qt.ToolButtonStyle.ToolButtonIconOnly)
            # Smaller icon sizes tend to look better in this case
            tmp_button.setIconSize(QtCore.QSize(self._button_size[0], self._button_size[1]))
            tmp_button.setFixedSize(self._button_size[0], self._button_size[1])
            self._start_layout.addWidget(tmp_button)
            self.item_to_button[tmp_quick_access_bar_action] = tmp_button

    # </editor-fold>
