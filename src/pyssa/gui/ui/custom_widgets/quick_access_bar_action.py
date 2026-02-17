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
"""QuickAccessBarAction wrapper for QAction objects in quick access toolbars.

This module provides the QuickAccessBarAction class, which serves as a wrapper
around QAction objects specifically designed for use in quick access toolbars.
It extends the basic QAction functionality with positioning information that
helps organize actions within different toolbar locations.

The QuickAccessBarAction encapsulates both the action itself and metadata about
its intended placement, including which toolbar section it belongs to and its
position within that section. This enables flexible toolbar organization and
consistent action management across the application.

Authors: Martin Urban, Hannah Kullik

Version: 1.4.0
"""

from typing import Union

from src.pyssa.gui.qt import QtGui, QtWidgets

__docformat__ = "google"


class QuickAccessBarAction:
    """A wrapper class for QAction objects used in quick access toolbars.

    This class encapsulates a QAction along with positioning metadata that
    defines where the action should be placed within quick access toolbars.
    It provides a structured way to organize toolbar actions with both
    sectional and positional information.

    The class maintains the underlying QAction while adding toolbar-specific
    positioning data that can be used by toolbar managers to arrange actions
    in a consistent and organized manner.

    Attributes:
        _action: The underlying QAction object that defines the action's
                behavior, icon, text, and parent widget.
        _bar_position: String identifier for the toolbar section where this
                      action should be placed (e.g., "left", "right", "bottom").
        _in_bar_position: Integer position within the specified toolbar section,
                         determining the order of actions (0, 1, 2, ...).
    """

    # <editor-fold desc="Constructor">
    def __init__(
        self,
        a_item_name: str,
        a_bar_position: str,
        an_in_bar_position: int,
        a_parent: Union[QtWidgets.QWidget, None],
        an_icon_filepath: QtGui.QIcon,
    ) -> None:
        """Initialize a QuickAccessBarAction with positioning metadata.

        Creates a new QuickAccessBarAction that wraps a QAction with additional
        positioning information for use in quick access toolbars. The action
        is created with the specified icon, text, and parent widget.

        Args:
            a_item_name: The display text for the action. This text may be
                        shown in tooltips or context menus.
            a_bar_position: String identifier for the toolbar section where
                           this action should be placed (e.g., "left", "right",
                           "bottom"). Used by toolbar managers for organization.
            an_in_bar_position: Integer position within the specified toolbar
                               section, determining the display order. Lower
                               numbers appear first (0, 1, 2, ...).
            a_parent: The parent widget for the QAction. Can be None for
                     actions without a specific parent.
            an_icon_filepath: The QIcon to display for this action in toolbars
                             and menus.

        Note:
            The positioning information is used by toolbar managers to organize
            actions but does not directly affect the QAction's behavior.
        """
        # <editor-fold desc="Instance attributes">
        self._action: QtGui.QAction = QtGui.QAction(
            icon=an_icon_filepath, text=a_item_name, parent=a_parent
        )
        self._bar_position: str = a_bar_position
        self._in_bar_position: int = an_in_bar_position
        # </editor-fold>

    # </editor-fold>

    # <editor-fold desc="Public methods">
    def get_action(self) -> QtGui.QAction:
        """Retrieve the underlying QAction object.

        Returns the wrapped QAction that can be used directly with Qt widgets
        that accept QAction objects, such as menus, toolbars, and buttons.

        Returns:
            The QAction object containing the action's icon, text, behavior,
            and parent widget information.
        """
        return self._action

    # </editor-fold>
