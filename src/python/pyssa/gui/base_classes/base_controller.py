from PyQt6 import QtWidgets, QtCore


class BaseController(QtCore.QObject):
  """Base class for all controller classes."""

  def get_dialog(self) -> QtWidgets.QDialog:
    """Gets the dialog of the controller."""
    raise NotImplementedError

  def restore_ui(self) -> None:
    """Restores the UI to default values."""
    raise NotImplementedError

  def set_dialog_close_as_canceled(self) -> None:
    """Sets the was_canceled flag to true."""
    raise NotImplementedError
