#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/urban233/PySSA>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""Module for the About dialog."""
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore
from src.pyssa.gui.ui.forms.auto_generated.auto_dialog_about import Ui_Dialog
from src.pyssa.gui.ui.styles import styles
from src.pyssa.util import constants


class DialogAbout(QtWidgets.QDialog):
  """Class representing an About dialog."""

  def __init__(self, parent=None) -> None:  # noqa: ANN001
    """Constructor.

    Args:
        parent: The parent.
    """
    QtWidgets.QDialog.__init__(self, parent)
    # build ui object
    self.ui = Ui_Dialog()
    self.ui.setupUi(self)

    self.ui.btn_ok.clicked.connect(self.close_dialog)
    original_pixmap = QtGui.QPixmap(
        f"{constants.PROGRAM_BIN_ROOT_PATH}\\assets\\images\\pyssa_logo.png"
    )
    scaled_pixmap = original_pixmap.scaled(150, 150)
    self.ui.lbl_pyssa_logo.setPixmap(scaled_pixmap)
    styles.set_stylesheet(self)
    self.ui.label_2.setText(f"Version: {constants.VERSION_NUMBER}")
    self.ui.label_2.setStyleSheet("font-weight: bold;")
    self.ui.label.setStyleSheet("font-size: 19px")

    self._fill_table_view()

    self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
    self.setWindowTitle("About")
    self.setWindowFlags(
        self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint
    )

  def _fill_table_view(self) -> None:
    """Fill the table with the packages."""
    tmp_table_model = QtGui.QStandardItemModel()
    tmp_table_model.setHorizontalHeaderLabels(["Name", "Version", "License"])
    # Sample data
    data = [
        ["Biopython", "1.78", "BSD 3-Clause License"],
        ["Material Design 3 Icons", "4.0.0", "Apache License 2.0"],
        ["Matplotlib", "3.8.0", "Python Software Foundation License (PSF)"],
        ["Numpy", "1.26.2", "BSD 3-Clause 'New' or 'Revised' License"],
        ["PyMOL Open-Source", "3.0.0", "BSD-like license"],
        ["Pandas", "2.1.1", "BSD 3-Clause 'New' or 'Revised' License"],
        ["PyQt5", "5.15.10", "GNU General Public License (GPL)"],
        ["SQLite", "3.41.2", "Public Domain"],
    ]

    # Populate the model with data
    for row, row_data in enumerate(data):
      for col, value in enumerate(row_data):
        item = QtGui.QStandardItem(
            str(value)
        )  # or QStandardItem(str(value)) for QStandardItemModel
        tmp_table_model.setItem(row, col, item)
    self.ui.tableView.setModel(tmp_table_model)
    self.ui.tableView.resizeColumnsToContents()
    self.ui.tableView.setEditTriggers(QtWidgets.QTableView.NoEditTriggers)

  # @SLOT
  def close_dialog(self) -> None:
    """Closes the dialog."""
    self.close()
