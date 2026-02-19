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
"""PySSAObjectsPanelController class for the PySSA frontend application.

Authors: Martin Urban, Hannah Kullik

Version: 1.4.0
"""
from pyssa.controller import interface_manager
from pyssa.gui.ui.views import pyssa_objects_panel
from src.pyssa.gui import user_pymol
from src.pyssa.gui.qt import QtCore, QtWidgets


class PySSAObjectsPanelController:
  """Controller class for the PySSAObjectsPanel."""

  def __init__(
          self,
          a_pyssa_objects_panel: "pyssa_objects_panel.PySSAObjectsPanel",
          an_interface_manager: "interface_manager.InterfaceManager",
          a_user_pymol: "user_pymol.UserPyMOL"
  ):
    """Constructor."""
    # <editor-fold desc="Instance attributes">
    self._panel = a_pyssa_objects_panel
    self._interface_manager = an_interface_manager
    self._user_pymol: "user_pymol.UserPyMOL" = a_user_pymol
    self._model = self._interface_manager.get_pyssa_objects_model()
    # </editor-fold>
    self._set_model()
    self._connect_signals()

  def _set_model(self) -> None:
    self._panel.tree_view.setModel(self._model)

  def _connect_signals(self) -> None:
    self._panel.import_file_action.get_action().triggered.connect(
      self.__slot_import_file
    )
    self._panel.export_file_action.get_action().triggered.connect(
      self.__slot_export_file
    )

    # self.expand_all.clicked.connect(self._panel.tree_view.expandAll)
    # self.collapse_all.clicked.connect(self._panel.tree_view.collapseAll)

  def __slot_import_file(self):
    file_path, _ = QtWidgets.QFileDialog.getOpenFileName(
      self._panel,
      "Open PDB/mmCIF File",
      "",
      "PDB/mmcif Files (*.pdb *.mmcif *.cif)"
    )

    if file_path:
      print(file_path)
      self._user_pymol.get_cmd_module().load(file_path)
      self._model.add_protein(self._user_pymol.get_cmd_module().get_model())

  def __slot_export_file(self):
    # Determine selected top-level molecule name from tree view
    selection_model = self._panel.tree_view.selectionModel()
    if selection_model is None:
      return
    selected_indexes = selection_model.selectedIndexes()
    if not selected_indexes:
      QtWidgets.QMessageBox.information(
        self._panel,
        "Export",
        "Please select a molecule (or any of its children) in the list to export.",
      )
      return

    # Use the first selected index and climb to the top-level item to get the molecule name
    index = selected_indexes[0]
    model = self._panel.tree_view.model()
    while index.isValid() and model.parent(index).isValid():
      index = model.parent(index)

    molecule_name = index.data(QtCore.Qt.ItemDataRole.DisplayRole) if index.isValid() else None
    if not molecule_name:
      QtWidgets.QMessageBox.warning(self._panel, "Export", "Could not determine the molecule name for export.")
      return

    # Ask for destination file path. Provide a default filename based on the molecule name.
    default_filename = f"{molecule_name}.pdb"
    file_path, selected_filter = QtWidgets.QFileDialog.getSaveFileName(
      self._panel,
      "Save Molecule As",
      default_filename,
      "PDB File (*.pdb);;mmCIF File (*.cif *.mmcif);;All Files (*.*)",
    )

    if not file_path:
      return

    # If no extension was provided, infer from the chosen filter.
    if "." not in file_path.split("\\")[-1]:
      if "mmCIF" in selected_filter:
        file_path = file_path + ".cif"
      else:
        file_path = file_path + ".pdb"

    try:
      # Export the selected molecule by name from PyMOL
      self._user_pymol.get_cmd_module().save(file_path, molecule_name)
    except Exception as e:
      QtWidgets.QMessageBox.critical(
        self._panel,
        "Export Failed",
        f"Failed to export '{molecule_name}': {e}",
      )
