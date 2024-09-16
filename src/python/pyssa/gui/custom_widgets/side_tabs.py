from typing import Optional

from PyQt6 import QtWidgets
from PyQt6 import QtCore

from pyssa.model.preference import model_definitions
from pyssa.model.util.gui_style import icons


class SideTabs(QtWidgets.QWidget):
  def __init__(self, text, icon: Optional["model_definitions.IconsEnum"] = None):
    super().__init__()
    layout = QtWidgets.QHBoxLayout()
    self.setLayout(layout)
    spacer_item = QtWidgets.QSpacerItem(
      20, 1,
      QtWidgets.QSizePolicy.Policy.Maximum,
      QtWidgets.QSizePolicy.Policy.Minimum
    )
    # Create and set up the icon label
    icon_label = QtWidgets.QLabel()
    if icon is not None:
      icon_label.setPixmap(icons.get_icon(icon).pixmap(28, 28))
      icon_label.setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft)
    layout.addWidget(icon_label)
    # Create and set up the text label
    text_label = QtWidgets.QLabel(text)
    text_label.setAlignment(QtCore.Qt.AlignmentFlag.AlignLeft)
    layout.addWidget(text_label)

    # Optional: Add some spacing between icon and text
    #layout.setSpacing(5)
    layout.setContentsMargins(0, 0, 0, 0)
