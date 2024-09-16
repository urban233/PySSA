from PyQt6 import QtCore
from PyQt6 import QtWidgets


class BaseDialog(QtWidgets.QDialog):
  """Base class for all dialogs."""

  dialogClosed = QtCore.pyqtSignal()
  """A signal indicating that the dialog is closed."""

  def setup_ui(self) -> None:
    """Sets up the initial ui."""
    raise NotImplementedError()
