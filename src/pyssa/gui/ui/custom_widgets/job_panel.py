"""Dropdown panel that displays active and completed jobs in two table views.

The panel is designed to be shown as a popup from a ``QToolButton`` in the
Quick Access Bar toolbar above the PyMOL viewer.  It contains two
``QTableView`` widgets stacked vertically — one for active (queued/running)
jobs and one for completed (finished/failed/cancelled) jobs.
"""
from __future__ import annotations

from typing import Optional

from src.pyssa.gui.qt import QtWidgets, QtCore


class JobPanel(QtWidgets.QWidget):
  """Popup panel showing active and completed jobs in table views.

  Args:
      model: The ``JobModel`` that backs both table views.
      parent: Optional parent widget.
  """

  def __init__(
      self,
      parent: Optional[QtWidgets.QWidget] = None,
  ) -> None:
    super().__init__(parent)
    self._lbl_jobs = QtWidgets.QLabel("Active Jobs")
    self._table_view = self._create_table_view()
    self._build_ui()
    self.setMinimumSize(450, 250)

  # ------------------------------------------------------------------
  # UI construction
  # ------------------------------------------------------------------

  def _build_ui(self) -> None:
    """Assemble the panel layout with two labelled table views."""
    layout = QtWidgets.QVBoxLayout(self)
    layout.setContentsMargins(8, 8, 8, 8)
    layout.setSpacing(6)

    self._lbl_jobs.setStyleSheet("font-weight: bold; font-size: 12px;")
    layout.addWidget(self._lbl_jobs)
    layout.addWidget(self._table_view, stretch=1)

    self.setLayout(layout)
    self.setStyleSheet(
      """
      JobPanel {
        background-color: #ffffff;
        border: 1px solid #DCDBE3;
        border-radius: 6px;
      }
      QTableView {
        border: 1px solid #e0e0e0;
        border-radius: 4px;
        background-color: #fafafa;
        selection-background-color: #e3f0ff;
        gridline-color: #eeeeee;
        font-size: 11px;
      }
      QTableView::item {
        padding: 4px 6px;
      }
      QHeaderView::section {
        background-color: #f5f5f5;
        border: none;
        border-bottom: 1px solid #e0e0e0;
        padding: 4px 6px;
        font-weight: bold;
        font-size: 11px;
      }
      QLabel {
        color: #333333;
      }
      """
    )

  def show_active_jobs(self, an_active_proxy):
    self._lbl_jobs.setText("Active Jobs")
    self._table_view.setModel(an_active_proxy)
    self._table_view.resizeColumnsToContents()

  def show_completed_jobs(self, an_completed_proxy):
    self._lbl_jobs.setText("Completed Jobs")
    self._table_view.setModel(an_completed_proxy)
    self._table_view.resizeColumnsToContents()

  @staticmethod
  def _create_table_view(
      proxy: QtCore.QSortFilterProxyModel | None = None,
  ) -> QtWidgets.QTableView:
    """Create and configure a read-only ``QTableView``.

    Args:
        proxy: The proxy model to display.

    Returns:
        A configured ``QTableView`` instance.
    """
    view = QtWidgets.QTableView()
    view.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
    view.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
    view.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
    view.verticalHeader().setVisible(False)
    view.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.ResizeToContents)
    view.horizontalHeader().setStretchLastSection(True)
    view.setAlternatingRowColors(True)
    view.setShowGrid(False)
    return view

  # ------------------------------------------------------------------
  # Public helpers
  # ------------------------------------------------------------------

  def show_below(self, widget: QtWidgets.QWidget) -> None:
    """Position the panel below *widget* and show it.

    Args:
        widget: The toolbar button (or any widget) used as the anchor.
    """
    pos = widget.mapToGlobal(
      QtCore.QPoint(0, widget.height()),
    )
    # Align right edge of panel with right edge of anchor
    pos.setX(pos.x() + widget.width() - self.width())
    self.move(pos)
    self.show()
