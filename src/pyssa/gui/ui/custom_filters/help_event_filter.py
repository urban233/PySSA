from src.pyssa.gui.qt import QtCore
import logging

logger = logging.getLogger(__name__)


class HelpEventFilter(QtCore.QObject):
  def __init__(self, help_browser, help_map, help_panel=None):
    super().__init__()
    self.help_browser = help_browser
    self.help_map = help_map
    self.help_panel = help_panel

  def eventFilter(self, obj, event):
    if event.type() == QtCore.QEvent.Type.Enter:
      object_name = obj.objectName()
      logger.debug(f"Enter event on object: {object_name}")
      help_text = self.help_map.get(object_name, "No help available")
      logger.debug(f"Setting help text: {help_text[:50]}...")
      self.help_browser.setHtml(help_text)

    elif event.type() == QtCore.QEvent.Type.Leave:
      logger.debug(f"Leave event on object: {obj.objectName()}")
      self.help_browser.clear()

    return super().eventFilter(obj, event)
