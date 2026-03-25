from src.pyssa.gui.ui.custom_widgets import color_grid
from src.pyssa.gui.qt import QtWidgets

class PSAColorConfig(QtWidgets.QWidget):

  def __init__(self, a_color_grid: "color_grid.PyMOLColorGrid"):
    super().__init__()
    self._color_grid = a_color_grid
    self._main_layout = QtWidgets.QVBoxLayout()

    self._lbl_bg = QtWidgets.QLabel("Background Color")
    self.btn_white_bg = QtWidgets.QPushButton()
    self.btn_grey_bg = QtWidgets.QPushButton()
    self.btn_black_bg = QtWidgets.QPushButton()
    self._bg_layout = QtWidgets.QHBoxLayout()

    self._lbl_color_by_elements = QtWidgets.QLabel("Atoms by Element")
    self.btn_color_by_elements = QtWidgets.QPushButton("Color")
    self._color_by_elements_layout = QtWidgets.QHBoxLayout()

    self._init_ui()
    self._set_all_tooltips()

  def _init_ui(self):
    self._color_by_elements_layout.setContentsMargins(8, 0, 8, 2)
    self._color_by_elements_layout.addWidget(self._lbl_color_by_elements)
    self._color_by_elements_layout.addStretch()
    self._color_by_elements_layout.addWidget(self.btn_color_by_elements)

    self.btn_white_bg.setStyleSheet(self._generate_color_stylesheet("#ffffff"))
    self.btn_grey_bg.setStyleSheet(self._generate_color_stylesheet("#666666"))
    self.btn_black_bg.setStyleSheet(self._generate_color_stylesheet("#000000"))
    self._bg_layout.setContentsMargins(8, 2, 8, 6)
    self._bg_layout.addWidget(self._lbl_bg)
    self._bg_layout.addStretch()
    self._bg_layout.addWidget(self.btn_white_bg)
    self._bg_layout.addWidget(self.btn_grey_bg)
    self._bg_layout.addWidget(self.btn_black_bg)

    self._main_layout.setContentsMargins(0, 0, 0, 0)
    self._main_layout.addWidget(self._color_grid)
    self._main_layout.addLayout(self._color_by_elements_layout)
    self._main_layout.addLayout(self._bg_layout)
    self.setLayout(self._main_layout)

  def _set_all_tooltips(self) -> None:
    self.btn_white_bg.setToolTip("white")
    self.btn_grey_bg.setToolTip("grey40")
    self.btn_black_bg.setToolTip("black")

  def _generate_color_stylesheet(self, a_hex_color: str) -> str:
    """Generates a color stylesheet for a QPushButton with the specified background color.

    Args:
        a_hex_color: A string representing a hexadecimal color code.

    Returns:
        A string representing the generated color stylesheet.

    Raises:
        exception.IllegalArgumentError: If a_hex_color is either None or an empty string.
    """
    stylesheet = """QPushButton {
                background-color: %s;
                border: solid;
                border-width: 1px;
                border-radius: 4px;
                min-width: 20px;
                max-width: 20px;
                min-height: 20px;
                max-height: 20px;
            }
            QPushButton::hover {
                background-color: %s;
                border: solid;
                border-color: black;
                border-width: 2px;
                border-radius: 4px;
                min-width: 20px;
                max-width: 20px;
                min-height: 20px;
                max-height: 20px;
            }
        """ % (
      a_hex_color,
      a_hex_color,
    )
    return stylesheet  # noqa: RET504
