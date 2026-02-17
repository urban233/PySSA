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
"""Tool window layout for managing the layout of the main window."""
from typing import Optional

from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.qt import QtWidgets
from src.pyssa.gui.ui.custom_widgets import quick_access_bar


class ToolWindowLayout(QtWidgets.QWidget):
  """Modern IDE-style layout with collapsible left, right, and bottom panels.

  Mimics IntelliJ-style tool windows using nested QSplitters and QStackedWidgets.
  The central area typically hosts a main viewer widget, surrounded by collapsible
  tool panels on the left, right, and bottom sides.
  """

  def __init__(
          self,
          left_quick_access_bar_items,
          right_quick_access_bar_items,
          viewer_quick_access_bar_items,
          molecular_viewer: QtWidgets.QWidget,
          parent: Optional[QtWidgets.QWidget] = None
  ) -> None:
    super().__init__(parent)

    # --- Main components ---
    self.viewer_panel = molecular_viewer
    self.viewer_toolbar = quick_access_bar.QuickAccessBar(
      viewer_quick_access_bar_items, horizontal=True
    )
    # TODO: PySSA does not need these toolbars right now
    # self.left_toolbar = quick_access_bar.QuickAccessBar(left_quick_access_bar_items)
    # self.right_toolbar = quick_access_bar.QuickAccessBar(right_quick_access_bar_items)

    # --- QStackedWidgets for tool panels ---
    self.left_stack = QtWidgets.QStackedWidget()
    self.right_stack = QtWidgets.QStackedWidget()
    self.bottom_stack = QtWidgets.QStackedWidget()

    # --- QFrames for stylable containers ---
    self.left_frame = QtWidgets.QFrame()
    self.right_frame = QtWidgets.QFrame()
    self.bottom_frame = QtWidgets.QFrame()
    self.viewer_frame = QtWidgets.QFrame()

    # --- Splitters ---
    self.right_splitter = QtWidgets.QSplitter()
    self.left_splitter = QtWidgets.QSplitter()
    self.bottom_splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)

    # --- Panel state ---
    self._last_left_size = 200
    self._last_right_size = 250
    self._last_bottom_size = 200
    self._left_hidden = False
    self._right_hidden = False
    self._bottom_hidden = False

    self._init_ui()

  # -------------------------------------------------------------------------
  # UI setup
  # -------------------------------------------------------------------------
  def _init_ui(self) -> None:
    """Initialize UI layout with stylable frames."""
    # Configure stacked widget frames
    self._setup_stacked_panel_frame(self.left_frame, self.left_stack)
    self._setup_stacked_panel_frame(self.right_frame, self.right_stack)
    self._setup_stacked_panel_frame(self.bottom_frame, self.bottom_stack)

    # Viewer container (toolbar + viewer)
    content_panel = QtWidgets.QWidget()
    content_layout = QtWidgets.QVBoxLayout(content_panel)
    content_layout.setContentsMargins(0, 0, 0, 0)
    content_layout.addWidget(self.viewer_toolbar)
    tmp_viewer_frame_layout = QtWidgets.QVBoxLayout()
    tmp_viewer_frame_layout.addWidget(self.viewer_panel)
    tmp_viewer_frame_layout.setContentsMargins(4, 4, 4, 4)
    self.viewer_frame.setLayout(tmp_viewer_frame_layout)
    content_layout.addWidget(self.viewer_frame)
    self._apply_viewer_background("black")

    # Splitter hierarchy
    self.right_splitter.addWidget(content_panel)
    self.right_splitter.addWidget(self.right_frame)
    self.right_splitter.setSizes([800, self._last_right_size])
    self.right_splitter.setChildrenCollapsible(False)

    self.left_splitter.addWidget(self.left_frame)
    self.left_splitter.addWidget(self.right_splitter)
    self.left_splitter.setSizes([self._last_left_size, 800])
    self.left_splitter.setChildrenCollapsible(False)

    self.bottom_splitter.addWidget(self.left_splitter)
    self.bottom_splitter.addWidget(self.bottom_frame)
    self.bottom_splitter.setSizes([600, self._last_bottom_size])
    self.bottom_splitter.setChildrenCollapsible(False)

    # Toolbars + splitter
    main_panel = QtWidgets.QWidget()
    main_layout = QtWidgets.QHBoxLayout(main_panel)
    main_layout.setContentsMargins(0, 0, 0, 0)
    main_layout.addWidget(self.bottom_splitter)
    # TODO: PySSA does not need these toolbars right now
    # main_layout.addWidget(self.left_toolbar)
    # main_layout.addWidget(self.right_toolbar)

    # Root layout
    root_layout = QtWidgets.QVBoxLayout(self)
    root_layout.setContentsMargins(0, 0, 0, 0)
    root_layout.addWidget(main_panel)

    # Deferred divider setup
    QtCore.QTimer.singleShot(0, lambda: self.left_splitter.setSizes([self._last_left_size, 800]))

    # Default styling (you can override externally)
    self._apply_default_styles()

  @staticmethod
  def _setup_stacked_panel_frame(frame: QtWidgets.QFrame, stacked_widget: QtWidgets.QStackedWidget) -> None:
    """Embed a stacked widget into a frame for styling."""
    frame_layout = QtWidgets.QVBoxLayout(frame)
    frame_layout.setContentsMargins(2, 2, 2, 2)
    frame_layout.addWidget(stacked_widget)
    # TODO: Fancy stuff, is only available in JBioMOL
    # tmp_shadow_effect = QtWidgets.QGraphicsDropShadowEffect()
    # tmp_shadow_effect.setBlurRadius(8)
    # tmp_shadow_effect.setOffset(3, 3)
    # tmp_shadow_effect.setColor(QtGui.QColor(0, 0, 0, 15))
    # frame.setGraphicsEffect(tmp_shadow_effect)

  def _apply_default_styles(self) -> None:
    """Apply a simple border style to demonstrate panel styling."""
    self.left_frame.setStyleSheet(
      """
      QFrame {
            border: 0.075em solid white;
            background: white;
            border-radius: 0.75em;
        }
      """
    )
    self.right_frame.setStyleSheet(
      """
      QFrame {
            border: 0.075em solid white;
            background: white;
            border-radius: 0.75em;
        }
      """
    )
    self.bottom_frame.setStyleSheet(
      """
      QFrame {
            border: 0.075em solid white;
            background: white;
            border-radius: 0.75em;
        }
      """
    )

  def _apply_viewer_background(self, bg_color: str):
    """Applies the background color of PyMOL to the underlying QFrame.

    Args:
      bg_color: The background color for PyMOL as hex or color word (e.g. black, white)
    """
    tmp_stylesheet = """
    QFrame {
            border: 0.075em solid %s;
            background: %s;
            border-radius: 0.75em;
        }
      """ % (
      bg_color,
      bg_color,
    )
    self.viewer_frame.setStyleSheet(tmp_stylesheet)

  # -------------------------------------------------------------------------
  # Public methods for adding content
  # -------------------------------------------------------------------------
  def add_left_panel(self, panel: QtWidgets.QWidget) -> None:
    self.left_stack.addWidget(panel)

  def add_right_panel(self, panel: QtWidgets.QWidget) -> None:
    self.right_stack.addWidget(panel)

  def add_bottom_panel(self, panel: QtWidgets.QWidget) -> None:
    self.bottom_stack.addWidget(panel)

  # -------------------------------------------------------------------------
  # Visibility management
  # -------------------------------------------------------------------------
  def toggle_left_panel(self) -> None:
    self.set_left_panel_hidden(not self._left_hidden)

  def toggle_right_panel(self) -> None:
    self.set_right_panel_hidden(not self._right_hidden)

  def toggle_bottom_panel(self) -> None:
    self.set_bottom_panel_hidden(not self._bottom_hidden)

  def set_left_panel_hidden(self, hidden: bool) -> None:
    if hidden == self._left_hidden:
      return
    if hidden:
      # Save current size only when panel is visible
      current_sizes = self.left_splitter.sizes()
      if current_sizes[0] > 0:
        self._last_left_size = current_sizes[0]
      self.left_frame.hide()
    else:
      self.left_frame.show()
      # Calculate the size for the right section based on current total width
      total_width = sum(self.left_splitter.sizes())
      right_size = max(total_width - self._last_left_size, 100)
      QtCore.QTimer.singleShot(0, lambda: self.left_splitter.setSizes([self._last_left_size, right_size]))
    self._left_hidden = hidden

  def set_right_panel_hidden(self, hidden: bool) -> None:
    if hidden == self._right_hidden:
      return
    if hidden:
      # Save current size only when panel is visible
      current_sizes = self.right_splitter.sizes()
      if current_sizes[1] > 0:
        self._last_right_size = current_sizes[1]
      self.right_frame.hide()
    else:
      self.right_frame.show()
      # Calculate the size for the left section based on current total width
      total_width = sum(self.right_splitter.sizes())
      left_size = max(total_width - self._last_right_size, 100)
      QtCore.QTimer.singleShot(0, lambda: self.right_splitter.setSizes([left_size, self._last_right_size]))
    self._right_hidden = hidden

  def set_bottom_panel_hidden(self, hidden: bool) -> None:
    if hidden == self._bottom_hidden:
      return
    if hidden:
      # Save current size only when panel is visible
      current_sizes = self.bottom_splitter.sizes()
      if current_sizes[1] > 0:
        self._last_bottom_size = current_sizes[1]
      self.bottom_frame.hide()
    else:
      # Ensure the currently selected bottom page is visible
      try:
        current_widget = self.bottom_stack.currentWidget()
        if current_widget is not None:
          current_widget.show()
      except Exception:
        pass
      self.bottom_frame.show()
      # Calculate the size for the top section based on current total height
      total_height = sum(self.bottom_splitter.sizes())
      top_size = max(total_height - self._last_bottom_size, 100)
      QtCore.QTimer.singleShot(0, lambda: self.bottom_splitter.setSizes([top_size, self._last_bottom_size]))
    self._bottom_hidden = hidden

  # -------------------------------------------------------------------------
  # Accessors
  # -------------------------------------------------------------------------
  @property
  def is_left_panel_hidden(self) -> bool:
    return self._left_hidden

  @property
  def is_right_panel_hidden(self) -> bool:
    return self._right_hidden

  @property
  def is_bottom_panel_hidden(self) -> bool:
    return self._bottom_hidden
