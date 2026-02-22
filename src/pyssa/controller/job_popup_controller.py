from src.pyssa.gui.qt import QtWidgets
from src.pyssa.gui.ui.custom_widgets.job_panel import JobPanel
from src.pyssa.model.job_model import JobModel, ActiveJobsProxyModel, CompletedJobsProxyModel


class JobPopupController:
  """Controller for the job popup widget."""

  def __init__(self, view: JobPanel, model: JobModel):
    self._view = view
    self._model = model
    self._active_proxy = ActiveJobsProxyModel(self._view)
    self._active_proxy.setSourceModel(model)
    self._completed_proxy = CompletedJobsProxyModel(self._view)
    self._completed_proxy.setSourceModel(model)

  def show_popup(self, a_widget: QtWidgets.QWidget | None = None):
    self._view.show_below(a_widget)

  def show_active_jobs(self):
    self._view.show_active_jobs(self._active_proxy)

  def show_completed_jobs(self):
    self._view.show_completed_jobs(self._completed_proxy)
