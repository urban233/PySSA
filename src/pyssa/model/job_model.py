"""Table model for tracking all jobs in the application.

``JobModel`` is a ``QAbstractTableModel`` that acts as the single source of
truth for every job's lifecycle.  Two ``QSortFilterProxyModel`` instances
(active and completed) consume this model so that separate ``QTableView``
widgets can display running/queued and finished/failed jobs respectively.
"""
from __future__ import annotations

import logging
from typing import Any, Optional

from src.pyssa.gui.qt import QtCore, pyqtSignal, Qt
from src.pyssa.internal.data_structures.data_classes import job_descriptor
from src.pyssa.util import enums

logger = logging.getLogger(__name__)

_COLUMN_HEADERS: list[str] = ["Job", "Name", "Project", "Status"]
_COL_JOB = 0
_COL_NAME = 1
_COL_PROJECT = 2
_COL_STATUS = 3
_COLUMN_COUNT = len(_COLUMN_HEADERS)

_JOB_TYPE_LABELS: dict[enums.JobType, str] = {
  enums.JobType.PREDICTION: "Prediction",
  enums.JobType.DISTANCE_ANALYSIS: "Distance Analysis",
  enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS: "Prediction + Dist. Analysis",
  enums.JobType.RAY_TRACING: "Ray Tracing",
  enums.JobType.SIMPLE_IMAGE: "Simple Image",
  enums.JobType.GENERAL_PURPOSE: "General Purpose",
}


class _JobEntry:
  """Mutable row-level storage used only inside ``JobModel``.

  Attributes:
      descriptor: The immutable job metadata.
      status: The current lifecycle status (mutable).
      result: The return value of the job's worker function, set on success.
  """

  __slots__ = ("descriptor", "status", "result")

  def __init__(
      self,
      descriptor: job_descriptor.JobDescriptor,
      status: enums.JobStatus = enums.JobStatus.QUEUED,
  ) -> None:
    self.descriptor = descriptor
    self.status = status
    self.result: Any = None


class JobModel(QtCore.QAbstractTableModel):
  """Table model tracking every job submitted to the scheduler.

  Signals:
      job_finished: Emitted when a job transitions to ``FINISHED``.
                    Carries the ``JobDescriptor`` and the result object.
      job_failed: Emitted when a job transitions to ``FAILED``.
                  Carries the ``JobDescriptor`` and the exception.
  """

  job_finished = pyqtSignal(object, object)
  job_failed = pyqtSignal(object, object)

  def __init__(self, parent: Optional[QtCore.QObject] = None) -> None:
    """Initialise an empty job model.

    Args:
        parent: Optional Qt parent object.
    """
    super().__init__(parent)
    self._entries: list[_JobEntry] = []

  # ------------------------------------------------------------------
  # QAbstractTableModel interface
  # ------------------------------------------------------------------

  def rowCount(
      self, parent: QtCore.QModelIndex = QtCore.QModelIndex(),
  ) -> int:
    """Return the number of tracked jobs.

    Args:
        parent: Must be invalid for a flat table model.

    Returns:
        The number of rows in the model.
    """
    if parent.isValid():
      return 0
    return len(self._entries)

  def columnCount(
      self, parent: QtCore.QModelIndex = QtCore.QModelIndex(),
  ) -> int:
    """Return the fixed column count (4).

    Args:
        parent: Must be invalid for a flat table model.

    Returns:
        The number of columns.
    """
    if parent.isValid():
      return 0
    return _COLUMN_COUNT

  def data(
      self,
      index: QtCore.QModelIndex,
      role: int = Qt.DisplayRole,
  ) -> Any:
    """Return cell data for the given *index* and *role*.

    Args:
        index: The model index to query.
        role: The data role (only ``Qt.DisplayRole`` is supported).

    Returns:
        The string value for the cell, or ``None`` if the role is not
        ``DisplayRole``.
    """
    if not index.isValid():
      return None
    if role != Qt.DisplayRole:
      return None

    entry = self._entries[index.row()]
    col = index.column()
    if col == _COL_JOB:
      return _JOB_TYPE_LABELS.get(
        entry.descriptor.job_type,
        entry.descriptor.job_type.value,
      )
    if col == _COL_NAME:
      return entry.descriptor.display_name
    if col == _COL_PROJECT:
      return entry.descriptor.project_name
    if col == _COL_STATUS:
      return entry.status.value
    return None

  def headerData(
      self,
      section: int,
      orientation: Qt.Orientation,
      role: int = Qt.DisplayRole,
  ) -> Any:
    """Return column header labels.

    Args:
        section: The column index.
        orientation: Must be ``Qt.Horizontal`` for column headers.
        role: Only ``Qt.DisplayRole`` is supported.

    Returns:
        The header string for the given section, or ``None``.
    """
    if orientation == Qt.Horizontal and role == Qt.DisplayRole:
      if 0 <= section < _COLUMN_COUNT:
        return _COLUMN_HEADERS[section]
    return None

  # ------------------------------------------------------------------
  # Public mutation API
  # ------------------------------------------------------------------

  def add_job(self, descriptor: job_descriptor.JobDescriptor) -> int:
    """Append a new job to the model.

    Args:
        descriptor: The immutable job metadata.

    Returns:
        The row index of the newly added job.

    Raises:
        ValueError: If *descriptor* is ``None``.
    """
    if descriptor is None:
      raise ValueError("descriptor must not be None")

    row = len(self._entries)
    self.beginInsertRows(QtCore.QModelIndex(), row, row)
    self._entries.append(_JobEntry(descriptor))
    self.endInsertRows()
    logger.debug(
      "Job added at row %d: type=%s, project=%s",
      row, descriptor.job_type.value, descriptor.project_name,
    )
    return row

  def update_status(
      self,
      row: int,
      status: enums.JobStatus,
      result: Any = None,
  ) -> None:
    """Transition a job's status and notify attached views.

    When the status changes to ``FINISHED``, ``job_finished`` is emitted.
    When the status changes to ``FAILED``, ``job_failed`` is emitted.

    Args:
        row: The row index of the job.
        status: The new status.
        result: Optional result or exception to store alongside the status.

    Raises:
        IndexError: If *row* is out of range.
    """
    if row < 0 or row >= len(self._entries):
      raise IndexError(f"Row {row} is out of range (0..{len(self._entries) - 1})")

    entry = self._entries[row]
    entry.status = status
    entry.result = result

    first_index = self.index(row, 0)
    last_index = self.index(row, _COLUMN_COUNT - 1)
    self.dataChanged.emit(first_index, last_index, [Qt.DisplayRole])

    if status == enums.JobStatus.FINISHED:
      self.job_finished.emit(entry.descriptor, result)
    elif status == enums.JobStatus.FAILED:
      self.job_failed.emit(entry.descriptor, result)

  def remove_job(self, row: int) -> None:
    """Remove a job from the model.

    Args:
        row: The row index of the job to remove.

    Raises:
        IndexError: If *row* is out of range.
    """
    if row < 0 or row >= len(self._entries):
      raise IndexError(f"Row {row} is out of range (0..{len(self._entries) - 1})")

    self.beginRemoveRows(QtCore.QModelIndex(), row, row)
    del self._entries[row]
    self.endRemoveRows()

  def remove_completed_jobs(self) -> None:
    completed_statuses = {
      enums.JobStatus.FINISHED,
      enums.JobStatus.FAILED,
      enums.JobStatus.CANCELLED,
    }

    rows_to_remove = [
      i for i, entry in enumerate(self._entries)
      if entry.status in completed_statuses
    ]

    for row in reversed(rows_to_remove):
      self.remove_job(row)

  def get_descriptor(self, row: int) -> job_descriptor.JobDescriptor:
    """Return the ``JobDescriptor`` for a given row.

    Args:
        row: The row index.

    Returns:
        The descriptor for the job at the given row.

    Raises:
        IndexError: If *row* is out of range.
    """
    if row < 0 or row >= len(self._entries):
      raise IndexError(f"Row {row} is out of range (0..{len(self._entries) - 1})")
    return self._entries[row].descriptor

  def get_status(self, row: int) -> enums.JobStatus:
    """Return the current status for a given row.

    Args:
        row: The row index.

    Returns:
        The current ``JobStatus`` of the job.

    Raises:
        IndexError: If *row* is out of range.
    """
    if row < 0 or row >= len(self._entries):
      raise IndexError(f"Row {row} is out of range (0..{len(self._entries) - 1})")
    return self._entries[row].status

  def find_row_for_descriptor(
      self, descriptor: job_descriptor.JobDescriptor,
  ) -> int:
    """Find the row index of a descriptor by identity comparison.

    Args:
        descriptor: The descriptor to search for.

    Returns:
        The row index, or ``-1`` if not found.
    """
    for i, entry in enumerate(self._entries):
      if entry.descriptor is descriptor:
        return i
    return -1


class ActiveJobsProxyModel(QtCore.QSortFilterProxyModel):
  """Proxy that accepts only ``QUEUED`` and ``RUNNING`` jobs."""

  _ACTIVE_STATUSES = frozenset({
    enums.JobStatus.QUEUED.value,
    enums.JobStatus.RUNNING.value,
  })

  def filterAcceptsRow(
      self, source_row: int, source_parent: QtCore.QModelIndex,
  ) -> bool:
    """Accept rows whose status column value is in the active set.

    Args:
        source_row: Row in the source model.
        source_parent: Parent index (unused for flat models).

    Returns:
        ``True`` if the row's status is ``QUEUED`` or ``RUNNING``.
    """
    index = self.sourceModel().index(source_row, _COL_STATUS, source_parent)
    status_text = self.sourceModel().data(index, Qt.DisplayRole)
    return status_text in self._ACTIVE_STATUSES


class CompletedJobsProxyModel(QtCore.QSortFilterProxyModel):
  """Proxy that accepts only ``FINISHED``, ``FAILED``, and ``CANCELLED`` jobs."""

  _COMPLETED_STATUSES = frozenset({
    enums.JobStatus.FINISHED.value,
    enums.JobStatus.FAILED.value,
    enums.JobStatus.CANCELLED.value,
  })

  def filterAcceptsRow(
      self, source_row: int, source_parent: QtCore.QModelIndex,
  ) -> bool:
    """Accept rows whose status column value is in the completed set.

    Args:
        source_row: Row in the source model.
        source_parent: Parent index (unused for flat models).

    Returns:
        ``True`` if the row's status is ``FINISHED``, ``FAILED``, or
        ``CANCELLED``.
    """
    index = self.sourceModel().index(source_row, _COL_STATUS, source_parent)
    status_text = self.sourceModel().data(index, Qt.DisplayRole)
    return status_text in self._COMPLETED_STATUSES
