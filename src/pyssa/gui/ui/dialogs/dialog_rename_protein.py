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
"""Module for the renaming protein dialog."""
import pathlib

from PyQt5 import QtCore
from PyQt5 import QtGui
from PyQt5 import QtWidgets
from PyQt5.QtCore import pyqtSignal

from src.pyssa.gui.ui.custom_widgets import custom_line_edit
from src.pyssa.gui.ui.forms.auto_generated import auto_dialog_rename_protein
from src.pyssa.gui.ui.styles import styles
from src.pyssa.util import constants, workspace_util, input_validator

global_var_rename_protein = ("", False)


class DialogRenameProtein(QtWidgets.QDialog):
  """Class for renaming protein dialog."""

  # TODO: If this class should be used a new implementation is needed

  # """
  # A pyqtsignal that is used to hand over the protein structure information.
  # """
  # return_value = pyqtSignal(tuple)
  #
  # def __init__(self, the_workspace_path: pathlib.Path, parent=None) -> None:  # noqa: ANN001
  #     """Constructor.
  #
  #     Args:
  #         the_workspace_path: the path to the workspace.
  #         parent: the parent.
  #     """
  #     QtWidgets.QDialog.__init__(self, parent)
  #     # build ui object
  #     self.ui = auto_dialog_rename_protein.Ui_Dialog()
  #     self.ui.setupUi(self)
  #     self.ui.btn_rename_protein.setEnabled(False)
  #     self.ui.btn_cancel.clicked.connect(self.close_dialog)
  #     self.ui.btn_rename_protein.clicked.connect(self.rename_protein)
  #     self.ui.txt_rename_protein.textChanged.connect(self.validate_protein_name)
  #     self.ui.lbl_status.setText("")
  #
  #     tmp_protein_infos = workspace_util.scan_workspace_for_non_duplicate_proteins(the_workspace_path)
  #     if len(tmp_protein_infos) > 0:
  #         self.ui.list_workspace_proteins.clear()
  #         for tmp_info in tmp_protein_infos:
  #             self.ui.list_workspace_proteins.addItem(QtWidgets.QListWidgetItem(tmp_info.name))
  #     self.setWindowTitle("Add an existing protein to the current project")
  #     self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
  #     styles.set_stylesheet(self)
  #     self.ui.list_workspace_proteins.setStyleSheet(
  #         "border-style: solid;"
  #         "border-width: 2px;"
  #         "border-radius: 8px;"
  #         "border-color: #DCDBE3;"
  #         "background: #ffffff;",
  #     )
  #     self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
  #     self._temporary_ui_redesign()
  #
  # def _temporary_ui_redesign(self):
  #     self.ui.list_workspace_proteins.hide()
  #     self.ui.label_2.hide()
  #     self.ui.txt_rename_protein = custom_line_edit.CustomLineEdit()
  #
  # def rename_protein(self) -> None:
  #     """Renames the current selected protein."""
  #     self.return_value.emit((self.ui.txt_rename_protein.text(), True))
  #     self.close()
  #
  # def close_dialog(self) -> None:
  #     """Closes dialog."""
  #     self.return_value.emit((self.ui.txt_rename_protein.text(), False))
  #     self.close()
  #
  # def validate_protein_name(self) -> None:
  #     """Validates the protein name."""
  #     input_validator.InputValidator.validate_search_input(
  #         self.ui.list_workspace_proteins,
  #         self.ui.txt_rename_protein,
  #         self.ui.lbl_status,
  #     )
  #     input_validator.InputValidator.validate_protein_name(
  #         self.ui.txt_rename_protein,
  #         self.ui.lbl_status,
  #         self.ui.btn_rename_protein,
  #     )
