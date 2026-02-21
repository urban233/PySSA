from src.pyssa.gui.qt import QtCore


class HelpEventFilter(QtCore.QObject):
  def __init__(self, help_browser, help_map):
    super().__init__()
    self.help_browser = help_browser
    self.help_map = help_map

  def eventFilter(self, obj, event):
    if event.type() == QtCore.QEvent.Type.Enter:
      help_text = self.help_map.get(obj.objectName(), "No help available")
      self.help_browser.setHtml(help_text)

    elif event.type() == QtCore.QEvent.Type.Leave:
      self.help_browser.clear()

    return super().eventFilter(obj, event)
