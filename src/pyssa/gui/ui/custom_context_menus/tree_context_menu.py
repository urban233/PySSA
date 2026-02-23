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
"""Module for the unified tree view context menu."""
import logging
from typing import Callable

from src.pyssa.gui.qt import QtWidgets, QtCore
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.model import selection_snapshot

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"

# Valid section names that correspond to distinct selection contexts.
_VALID_SECTIONS = frozenset({
    "sequence",
    "standalone_protein",
    "protein_pair",
    "protein_pair_child",
})


class TreeContextMenu:
    """A unified, declarative context menu for the PySSA Objects Panel tree view.

    Items are grouped into named *sections* that map directly to selection
    contexts derived from a ``SelectionSnapshot``.  Each section is shown
    only when the current snapshot matches that context.

    Sections:
        ``sequence``: Visible when one or more sequences are selected.
        ``standalone_protein``: Visible when standalone proteins are selected.
        ``protein_pair``: Visible when protein pairs are selected.
        ``protein_pair_child``: Visible when proteins inside a pair are selected.

    Usage::

        menu = TreeContextMenu()
        menu.register_action(
            section="sequence",
            key="rename_sequence",
            label="Rename Sequence",
            callback=self.__slot_rename_sequence,
        )
        # In refresh_ui:
        menu.configure(snapshot)
        # On customContextMenuRequested:
        menu.show_at(tree_view.mapToGlobal(pos))
    """

    def __init__(self) -> None:
        """Constructor."""
        self._menu: QtWidgets.QMenu = QtWidgets.QMenu()
        # Maps section name -> ordered list of QAction objects in that section.
        self._sections: dict[str, list[QtWidgets.QAction]] = {
            section: [] for section in _VALID_SECTIONS
        }
        # Maps action key -> QAction for quick lookup.
        self._actions: dict[str, QtWidgets.QAction] = {}
        # Tracks which section each key belongs to, so configure() can target it.
        self._action_section: dict[str, str] = {}
        # Separator QActions injected between sections (one per section boundary).
        self._section_separators: dict[str, QtWidgets.QAction] = {}

    def register_action(
        self,
        section: str,
        key: str,
        label: str,
        callback: Callable,
    ) -> None:
        """Register a new context menu action under a named section.

        The action is appended to the underlying ``QMenu`` in section order,
        with separators automatically inserted between sections that contain
        items.  Actions do not reorder after registration, so call this method
        for all items during initialisation (typically inside
        ``_register_tree_context_menu_actions``).

        Args:
            section: The selection context this action belongs to.
                Must be one of ``'sequence'``, ``'standalone_protein'``,
                ``'protein_pair'``, or ``'protein_pair_child'``.
            key: A unique string identifier for this action.
                Used internally for lookup; never shown to the user.
            label: The human-readable text displayed in the menu.
            callback: The callable invoked when the action is triggered.

        Raises:
            ValueError: If ``section`` is not a recognised section name.
            ValueError: If ``key`` is already registered.
        """
        if section not in _VALID_SECTIONS:
            raise ValueError(
                f"Unknown section '{section}'. Valid sections are: {sorted(_VALID_SECTIONS)}"
            )
        if key in self._actions:
            raise ValueError(
                f"Action key '{key}' is already registered. Keys must be unique."
            )

        # Insert a separator before the first item of each new section, but
        # only if a previous section already has items (to avoid a leading sep).
        section_list = self._sections[section]
        if not section_list and any(self._sections[s] for s in _VALID_SECTIONS if s != section):
            separator = self._menu.addSeparator()
            self._section_separators[f"_sep_{section}"] = separator

        action = QtWidgets.QAction(label)
        action.triggered.connect(callback)
        self._menu.addAction(action)

        section_list.append(action)
        self._actions[key] = action
        self._action_section[key] = section

    def configure(self, snapshot: "selection_snapshot.SelectionSnapshot | None") -> None:
        """Show and enable only the actions that match the current selection context.

        This method is called by ``MainWindowController.refresh_ui`` so that the
        menu state is always consistent with the rest of the application UI.
        When ``snapshot`` is ``None`` (no project open or nothing selected) all
        actions are hidden and the menu will appear empty.

        Args:
            snapshot: The current immutable selection snapshot, or ``None``.
        """
        # Determine which sections are active for this snapshot.
        active_sections: set[str] = set()
        if snapshot is not None:
            if snapshot.raw_sequences:
                active_sections.add("sequence")
            if snapshot.raw_standalone_proteins:
                active_sections.add("standalone_protein")
            if snapshot.raw_protein_pairs:
                active_sections.add("protein_pair")
            if snapshot.raw_protein_pair_children:
                active_sections.add("protein_pair_child")

        # Apply visibility to every action and its separator.
        for section, actions in self._sections.items():
            is_visible = section in active_sections
            for action in actions:
                action.setVisible(is_visible)
            sep_key = f"_sep_{section}"
            if sep_key in self._section_separators:
                self._section_separators[sep_key].setVisible(is_visible)

    def show_at(self, global_pos: QtCore.QPoint) -> None:
        """Display the context menu at the given global screen position.

        This method should be called from the ``customContextMenuRequested``
        slot after mapping the local widget position to global screen coordinates
        via ``tree_view.mapToGlobal(pos)``.

        The menu only appears if at least one action is currently visible.

        Args:
            global_pos: Screen-space position of the right-click event.
        """
        has_visible = any(
            action.isVisible()
            for actions in self._sections.values()
            for action in actions
        )
        if has_visible:
            self._menu.exec(global_pos)
