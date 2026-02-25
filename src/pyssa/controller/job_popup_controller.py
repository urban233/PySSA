import os
import pathlib

from src.pyssa.gui.qt import QtWidgets, QtCore
from src.pyssa.gui.ui.custom_widgets.job_panel import JobPanel
from src.pyssa.model.job_model import JobModel, ActiveJobsProxyModel, CompletedJobsProxyModel, _COL_STATUS
from src.pyssa.util import enums


class JobPopupController:
  """Controller for the job popup widget."""

  def __init__(self, view: JobPanel, model: JobModel, filter_active_jobs: bool = True):
    self._view = view
    self._model = model
    self._is_active_mode = filter_active_jobs
    self._active_proxy = ActiveJobsProxyModel(self._view)
    self._active_proxy.setSourceModel(model)
    self._active_proxy.setFilterKeyColumn(_COL_STATUS)
    self._active_proxy.setDynamicSortFilter(True)
    self._completed_proxy = CompletedJobsProxyModel(self._view)
    self._completed_proxy.setSourceModel(model)
    self._completed_proxy.setFilterKeyColumn(_COL_STATUS)
    self._completed_proxy.setDynamicSortFilter(True)
    self._view.get_table_view().setContextMenuPolicy(QtCore.Qt.ContextMenuPolicy.CustomContextMenu)

    if self._is_active_mode:
      self._view.button_container_widget.hide()
    self._connect_all_signals_with_their_slots()

  def _connect_all_signals_with_their_slots(self):
    if not self._is_active_mode:
      self._completed_proxy.rowsInserted.connect(self._update_clear_history_action)
      self._completed_proxy.rowsRemoved.connect(self._update_clear_history_action)
      self._view.get_table_view().clicked.connect(self._check_row_for_job_result)
      self._view.btn_clear_history.clicked.connect(self._clear_history)
      self._view.btn_open_result.clicked.connect(self._open_selected_job_result)

  def show_popup(self, a_widget: QtWidgets.QWidget | None = None):
    if self._is_active_mode:
        self.show_active_jobs()
    else:
        self.show_completed_jobs()
    self._view.show_below(a_widget)

  def show_active_jobs(self):
    self._is_active_mode = True
    self._view.show_active_jobs(self._active_proxy)

  def show_completed_jobs(self):
    self._is_active_mode = False
    self._view.show_completed_jobs(self._completed_proxy)

  def _get_job_descriptor_of_selected_row(self):
    table = self._view.get_table_view()
    selection_model = table.selectionModel()
    selected_indexes = selection_model.selectedRows()  # returns a list
    # Take the single selected row
    proxy_index = selected_indexes[0]
    # Map from proxy to source model
    source_index = self._completed_proxy.mapToSource(proxy_index)
    row = source_index.row()
    # Check if job is finished
    return self._model.get_descriptor(row)

  def _update_clear_history_action(self):
    """Enable the Clear History action only if completed proxy has at least one row."""
    row_count = self._completed_proxy.rowCount()
    self._view.btn_clear_history.setEnabled(row_count > 0)

  def _check_row_for_job_result(self):
    tmp_descriptor = self._get_job_descriptor_of_selected_row()
    match tmp_descriptor.job_type:
      case enums.JobType.RAY_TRACING:
        self._view.btn_open_result.setEnabled(True)
      case enums.JobType.SIMPLE_IMAGE:
        self._view.btn_open_result.setEnabled(True)
      case _:
        self._view.btn_open_result.setEnabled(False)

  def _clear_history(self):
    if not self._is_active_mode:
      self._model.remove_completed_jobs()
    self._view.btn_open_result.setEnabled(False)

  def _open_selected_job_result(self):
    """Open the result of the currently selected finished job."""
    tmp_descriptor = self._get_job_descriptor_of_selected_row()
    match tmp_descriptor.job_type:
      case enums.JobType.RAY_TRACING:
        os.startfile(str(pathlib.Path(tmp_descriptor.run_args[0])))
      case enums.JobType.SIMPLE_IMAGE:
        os.startfile(str(pathlib.Path(tmp_descriptor.run_args[0])))
      case _:
        self._view.btn_open_result.setEnabled(False)
