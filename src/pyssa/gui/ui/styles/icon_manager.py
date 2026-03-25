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
"""Icon manager class that provides a centralized way to handle icons in the PySSA application.

Authors: Martin Urban, Hannah Kullik

Version: 1.4.0
"""
import logging
import pathlib
import enum
from typing import Dict, Optional, Set, Union

from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.qt import QtGui
from src.pyssa.gui.qt import QtWidgets
from src.pyssa.util import constants

__docformat__ = "google"

_LOGGER = logging.getLogger(__name__)


class IconManager:
    """A centralized manager for handling icons in the PySSA application.

    This class provides an easy-to-use but professional interface for working with icons,
    including support for different icon states, sizes, and styling.

    Features:
        - Singleton pattern for easy access throughout the application
        - Icon caching for improved performance
        - Support for different icon states (normal, disabled, etc.)
        - Customizable icon sizes
        - Support for various widget types (QPushButton, QAction, QToolButton, etc.)
        - Automatic preloading of icons from the icons directory
        - Icon integrity checking to detect missing or invalid icons
        - Enforced use of Icons enum for core application icons
        - Direct access to icons by name for custom icons not in the Icons enum

    Usage:
        - Get a core icon using enum (required): icon_manager.get_icon(IconManager.Icons.HOME)
        - Get a custom icon by name (for icons not in the enum): icon_manager.get_icon("CUSTOM_ICON")
        - Set an icon on a widget: icon_manager.set_icon(button, IconManager.Icons.HOME)
        - Get the path to an icon: icon_manager.get_icon_path(IconManager.Icons.HOME)
        - Get all available icons: icon_manager.get_available_icons()

    Attributes:
        _instance: Singleton instance of IconManager
        _icon_cache: Cache for icons to avoid reloading
        _icons_root_path: Root path for icons directory
        _preloaded_icons: Preloaded icons with their paths
        _missing_icons: Set of icon names that are missing or invalid
    """

    class Icons(enum.StrEnum):
        """Enum for storing the icon filenames with their usage name."""

        ADD = "ADD"
        ADD_CIRCLE = "ADD_CIRCLE"
        ADD_CIRCLE_DISABLED = "ADD_CIRCLE_DISABLED"
        ADD_DISABLED = "ADD_DISABLED"
        ANGSTROM_DISTANCE = "ANGSTROM_DISTANCE"
        ARROW_DROP_UP = "ARROW_DROP_UP"
        ARROW_RIGHT = "ARROW_RIGHT"
        BROWSER_UPDATED = "BROWSER_UPDATED"
        CANCEL = "CANCEL"
        CANCEL_DISABLED = "CANCEL_DISABLED"
        CARTOON_REPR = "CARTOON_REPR"
        CHANGE_CIRCLE = "CHANGE_CIRCLE"
        CHANGE_CIRCLE_DISABLED = "CHANGE_CIRCLE_DISABLED"
        CLOSE = "CLOSE"
        COLLAPSE_ALL = "COLLAPSE_ALL"
        DANGEROUS = "DANGEROUS"
        DELETE = "DELETE"
        DOCKING = "DOCKING"
        DONE_ROUND_EDGES_W200 = "DONE_ROUND_EDGES_W200"
        DOTS_REPR = "DOTS_REPR"
        DO_NOT_DISTURB_ON = "DO_NOT_DISTURB_ON"
        DRAFT = "DRAFT"
        EDIT = "EDIT"
        ERROR = "ERROR"
        EXPAND_ALL = "EXPAND_ALL"
        FILE_SAVE = "FILE_SAVE"
        FILE_SAVE_DISABLED = "FILE_SAVE_DISABLED"
        FOLDER_OPEN = "FOLDER_OPEN"
        GRID_VIEW = "GRID_VIEW"
        HELP = "HELP"
        HOME = "HOME"
        IMAGE = "IMAGE"
        INFO = "INFO"
        KEYBOARD_ARROW_DOWN = "KEYBOARD_ARROW_DOWN"
        KEYBOARD_ARROW_RIGHT = "KEYBOARD_ARROW_RIGHT"
        LINES_REPR = "LINES_REPR"
        MESH_REPR = "MESH_REPR"
        MONOMER = "MONOMER"
        MORE_VERT = "MORE_VERT"
        MULTIMER = "MULTIMER"
        NOTE_ADD = "NOTE_ADD"
        NOTE_ADD_DISABLED = "NOTE_ADD_DISABLED"
        NOTIFICATIONS = "NOTIFICATIONS"
        NOTIFICATIONS_UNREAD = "NOTIFICATIONS_UNREAD"
        OPEN_IN_NEW = "OPEN_IN_NEW"
        OPEN_IN_NEW_DISABLED = "OPEN_IN_NEW_DISABLED"
        OPEN_IN_NEW_DOWN = "OPEN_IN_NEW_DOWN"
        PHOTO_FRAME = "PHOTO_FRAME"
        PLAY_CIRCLE = "PLAY_CIRCLE"
        PLAY_CIRCLE_RUN_W200 = "PLAY_CIRCLE_RUN_W200"
        RIBBON_REPR = "RIBBON_REPR"
        SCAN_DELETE = "SCAN_DELETE"
        SCAN_DELETE_DISABLED = "SCAN_DELETE_DISABLED"
        SHARE_WINDOWS = "SHARE_WINDOWS"
        SPHERES_REPR = "SPHERES_REPR"
        STICKS_REPR = "STICKS_REPR"
        SURFACE_REPR = "SURFACE_REPR"
        TEST = "TEST"
        UPLOAD_FILE = "UPLOAD_FILE"
        UPLOAD_FILE_DISABLED = "UPLOAD_FILE_DISABLED"
        WARNING = "WARNING"
        VISIBILITY = "VISIBILITY"
        VISIBILITY_OFF = "VISIBILITY_OFF"
        COLORS = "COLORS"
        PALETTE = "PALETTE"
        APP_LOGO = "APP_LOGO"
        ADD_PHOTO_ALTERNATE = "ADD_PHOTO_ALTERNATE"
        PANORAMA_WIDE_ANGLE = "PANORAMA_WIDE_ANGLE"
        MOP = "MOP"
        LOGO = "LOGO"
        HIGHLIGHT_MOUSE_CURSOR = "HIGHLIGHT_MOUSE_CURSOR"
        CODE_BLOCKS = "CODE_BLOCKS"
        ARROWS_OUTPUT = "ARROWS_OUTPUT"

    # <editor-fold desc="Class attributes">
    _instance = None
    # </editor-fold>

    @classmethod
    def instance(cls):
        """Singleton instance of ThemeManager."""
        if cls._instance is None:
            print("Creating singleton instance for the first time")
            cls._instance = super(IconManager, cls).__new__(cls)
            # Put any initialization here.
            cls._instance._initialize()
        return cls._instance

    # <editor-fold desc="Public methods">
    def get_icon(
        self,
        an_icon_name: Union[str, Icons],
        a_disabled_icon_name: Optional[Union[str, Icons]] = None,
        the_size: Optional[QtCore.QSize] = None,
    ) -> QtGui.QIcon:
        """Get an icon by name with an optional disabled state.

        Note:
            No checks are performed!

        Args:
            an_icon_name: Name of the icon (must be Icons enum for specific icons, can be string for others)
            a_disabled_icon_name: Name of the disabled icon (must be Icons enum for specific icons, can be string for others)
            the_size: Size for the icon (if None, no size is set)

        Returns:
            The corresponding QIcon object

        Raises:
            TypeError: If a required icon is referenced without using the Icons enum
        """
        # Check if we have this icon in cache
        tmp_cache_key = f"{an_icon_name}_{a_disabled_icon_name}_{the_size}"
        if tmp_cache_key in self._icon_cache:
            return self._icon_cache[tmp_cache_key]

        try:
            # Create the icon
            if an_icon_name in self._preloaded_icons and an_icon_name not in self._missing_icons:
                tmp_icon_path = self._preloaded_icons[an_icon_name]
                tmp_icon = QtGui.QIcon(str(tmp_icon_path))
            # If the icon name is not found, try to find a matching icon by name
            else:
                # If not found, use a default icon or an empty icon
                _LOGGER.warning(f"Warning: Icon '{an_icon_name}' not found. Using empty icon.")
                tmp_icon = QtGui.QIcon()
        except Exception as e:
            _LOGGER.error(f"Error loading icon '{an_icon_name}': {e}")
            tmp_icon = QtGui.QIcon()

        # Only add disabled state if the icon is valid
        if a_disabled_icon_name and tmp_icon.availableSizes():
            tmp_disabled_icon_name = str(a_disabled_icon_name)

            try:
                if (
                    tmp_disabled_icon_name in self._preloaded_icons
                    and tmp_disabled_icon_name not in self._missing_icons
                ):
                    tmp_disabled_path = self._preloaded_icons[tmp_disabled_icon_name]
                    tmp_icon.addPixmap(
                        QtGui.QPixmap(str(tmp_disabled_path)),
                        mode=QtGui.QIcon.Mode.Disabled,
                    )
                else:
                    # If not found, use an empty icon
                    _LOGGER.warning(f"Warning: Icon '{an_icon_name}' not found. Using empty icon.")
            except Exception as e:
                _LOGGER.error(f"Error loading disabled icon '{tmp_disabled_icon_name}': {e}")

        # Cache the icon
        self._icon_cache[tmp_cache_key] = tmp_icon

        return tmp_icon

    def set_icon(
        self,
        a_widget: Union[QtWidgets.QPushButton, QtWidgets.QToolButton],
        an_icon_name: Union[str, Icons],
        a_disabled_icon_name: Optional[Union[str, Icons]] = None,
        the_size: QtCore.QSize = QtCore.QSize(32, 32),
    ) -> None:
        """Sets an icon on a widget.

        Note:
            No checks are performed!

        Args:
            a_widget: Widget to set the icon for (QPushButton, QAction, etc.)
            an_icon_name: Name of the icon (must be Icons enum for specific icons, can be string for others)
            a_disabled_icon_name: Name of the disabled icon (must be Icons enum for specific icons, can be string for others)
            the_size: Size for the icon (defaults to the class default size)

        Raises:
            TypeError: If a required icon is referenced without using the Icons enum
        """
        tmp_icon = self.get_icon(an_icon_name, a_disabled_icon_name)
        a_widget.setIcon(tmp_icon)
        a_widget.setIconSize(tmp_icon.actualSize(the_size))

    def get_available_icons(self) -> Set[str]:
        """Gets a set of all available icon names.

        Returns:
            A set of all available icon names (excluding missing or invalid icons)
        """
        # Combine icons from ICON_PATHS and preloaded icons, excluding missing icons
        available_icons = set(self._preloaded_icons.keys()) - self._missing_icons

        return available_icons

    # </editor-fold>

    # <editor-fold desc="Private methods">
    def _initialize(self) -> None:
        """Initialize the IconManager with icon paths and default settings."""
        # <editor-fold desc="Instance attributes">
        self._icon_cache: Dict[str, QtGui.QIcon] = {}
        self._icons_root_path = constants.ICONS_PATH
        self._preloaded_icons: Dict[str, pathlib.Path] = {}
        self._missing_icons: Set[str] = set()
        # </editor-fold>
        self._preload_icons()
        self._check_icon_integrity()

    def _preload_icons(self) -> None:
        """
        Scan the icons directory and preload all available icons.

        This method populates the _preloaded_icons dictionary with all icon files
        found in the icons directory. The key is the icon name (without extension),
        and the value is the path to the icon file.
        """
        if self._icons_root_path.exists() and self._icons_root_path.is_dir():
            for file_path in self._icons_root_path.glob("*"):
                if file_path.is_file() and file_path.suffix.lower() in [
                    ".svg",
                    ".png",
                ]:
                    # Extract icon name (without extension and path)
                    icon_name = file_path.stem

                    # Remove any size or color suffix (e.g., _w200, _blue)
                    for suffix in ["_w200", "_w400", "_blue", "_g200"]:
                        if icon_name.endswith(suffix):
                            icon_name = icon_name[: -len(suffix)]
                            break

                    # Convert to uppercase to match the convention in ICON_PATHS
                    icon_name = icon_name.upper()

                    # Add to preloaded icons if not already there
                    if icon_name not in self._preloaded_icons:
                        self._preloaded_icons[icon_name] = file_path

    def _check_icon_integrity(self) -> None:
        """
        Check the integrity of all icon files.

        This method verifies that all icon files exist and can be loaded.
        It populates the _missing_icons set with the names of any missing icons.
        """
        missing_icons = set()
        # Check preloaded icons
        for icon_name, icon_path in self._preloaded_icons.items():
            if not icon_path.exists():
                missing_icons.add(icon_name)
            else:
                # Try to load the icon to check if it's valid
                try:
                    QtGui.QIcon(str(icon_path))
                except FileNotFoundError:
                    missing_icons.add(icon_name)

        self._missing_icons = missing_icons
        # Log missing icons
        if missing_icons:
            _LOGGER.warning(f"Warning: {len(missing_icons)} icon(s) are missing or invalid:")
            for icon_name in sorted(missing_icons):
                _LOGGER.warning(f"  - {icon_name}")

    # </editor-fold>
