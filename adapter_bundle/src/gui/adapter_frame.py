from PyQt6 import QtWidgets, QtGui, QtCore
from PyQt6.QtGui import QAction
from PyQt6.QtWidgets import QMenuBar, QMenu
from chimerax.ui import MainToolWindow

from pyssa.gui.dialog.forms.auto import auto_adapter_frame
from . import adapter_frame_controller


class AdapterFrame(QtWidgets.QDialog):
  """Adapter form to use in ChimeraX as tool."""

  def __init__(self, the_tool_instance, the_tool_window: MainToolWindow, a_session):
    super().__init__()
    self.ui = auto_adapter_frame.Ui_Dialog()
    self.ui.setupUi(the_tool_window.ui_area)
    self.adapter_frame_controller = adapter_frame_controller.AdapterFrameController(the_tool_instance, self, a_session)
    self.adapter_frame_controller.start_checking()
