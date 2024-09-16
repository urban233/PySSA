import pathlib

from PyQt6 import QtWidgets
from PyQt6 import QtCore
from PyQt6 import QtGui


class MoleculeNaming(QtWidgets.QWidget):
  """A custom widget class that is used for displaying a molecule that can then be named."""

  def __init__(self, an_image_filepath: pathlib.Path, a_name: str, parent=None):
    """Constructor."""
    super().__init__(parent)
    self.layout = QtWidgets.QVBoxLayout(self)
    self.line_edit = QtWidgets.QLineEdit(self)
    self.line_edit.setText(a_name)
    self.image_label = QtWidgets.QLabel(self)
    pixmap = QtGui.QPixmap(str(an_image_filepath))
    self.image_label.setPixmap(pixmap)
    self.layout.addWidget(self.line_edit)
    self.layout.addWidget(self.image_label)
    self.setLayout(self.layout)
