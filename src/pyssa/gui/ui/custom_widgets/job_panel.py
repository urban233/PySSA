"""Dropdown panel that displays active and completed jobs in two table views.

The panel is designed to be shown as a popup from a ``QToolButton`` in the
Quick Access Bar toolbar above the PyMOL viewer.  It contains two
``QTableView`` widgets stacked vertically — one for active (queued/running)
jobs and one for completed (finished/failed/cancelled) jobs.
"""
from __future__ import annotations

from typing import Optional

from src.pyssa.gui.qt import QtWidgets, QtCore, Qt
from src.pyssa.model.job_model import (
  JobModel,
  ActiveJobsProxyModel,
  CompletedJobsProxyModel,
)


class JobPanel(QtWidgets.QWidget):
  """Popup panel showing active and completed jobs in table views.

  Args:
      model: The ``JobModel`` that backs both table views.
      parent: Optional parent widget.
  """

  def __init__(
      self,
      model: JobModel,
      parent: Optional[QtWidgets.QWidget] = None,
  ) -> None:
    super().__init__(parent, Qt.Popup | Qt.FramelessWindowHint)
    if model is None:
      raise ValueError("model must not be None")

    self._model = model
    self._active_proxy = ActiveJobsProxyModel(self)
    self._active_proxy.setSourceModel(model)
    self._completed_proxy = CompletedJobsProxyModel(self)
    self._completed_proxy.setSourceModel(model)

    self._build_ui()
    self.setMinimumSize(520, 320)

  # ------------------------------------------------------------------
  # UI construction
  # ------------------------------------------------------------------

  def _build_ui(self) -> None:
    """Assemble the panel layout with two labelled table views."""
    layout = QtWidgets.QVBoxLayout(self)
    layout.setContentsMargins(8, 8, 8, 8)
    layout.setSpacing(6)

    active_label = QtWidgets.QLabel("Active Jobs")
    active_label.setStyleSheet("font-weight: bold; font-size: 12px;")
    layout.addWidget(active_label)

    self._active_view = self._create_table_view(self._active_proxy)
    layout.addWidget(self._active_view, stretch=1)

    completed_label = QtWidgets.QLabel("Completed Jobs")
    completed_label.setStyleSheet("font-weight: bold; font-size: 12px;")
    layout.addWidget(completed_label)

    self._completed_view = self._create_table_view(self._completed_proxy)
    layout.addWidget(self._completed_view, stretch=1)

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

  @staticmethod
  def _create_table_view(
      proxy: QtCore.QSortFilterProxyModel,
  ) -> QtWidgets.QTableView:
    """Create and configure a read-only ``QTableView``.

    Args:
        proxy: The proxy model to display.

    Returns:
        A configured ``QTableView`` instance.
    """
    view = QtWidgets.QTableView()
    view.setModel(proxy)
    view.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
    view.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
    view.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
    view.verticalHeader().setVisible(False)
    view.horizontalHeader().setStretchLastSection(True)
    view.setAlternatingRowColors(True)
    view.setShowGrid(False)
    return view

  # ------------------------------------------------------------------
  # Public helpers
  # ------------------------------------------------------------------

  @property
  def active_view(self) -> QtWidgets.QTableView:
    """The table view displaying queued and running jobs."""
    return self._active_view

  @property
  def completed_view(self) -> QtWidgets.QTableView:
    """The table view displaying finished, failed, and cancelled jobs."""
    return self._completed_view

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
