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
import logging

from src.pyssa.gui.ui.views import pyssa_objects_panel
from src.pyssa.gui import user_pymol, app_state
from src.pyssa.gui.qt import QtCore, QtWidgets
from src.pyssa.logging_pyssa import log_levels, log_handlers
from src.pyssa.model.selection_snapshot import SelectionSnapshot

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class PySSAObjectsPanelController(QtCore.QObject):
  """Controller class for the PySSAObjectsPanel."""
  
  # Custom signal emitted whenever the debounced selection changes
  selectionSnapshotUpdated = QtCore.Signal(object)

  def __init__(
          self,
          the_app_state: "app_state.AppState",
          a_pyssa_objects_panel: "pyssa_objects_panel.PySSAObjectsPanel",
          a_user_pymol: "user_pymol.UserPyMOL"
  ):
    """Constructor."""
    super().__init__()
    # <editor-fold desc="Instance attributes">
    self._app_state = the_app_state
    self._panel = a_pyssa_objects_panel
    self._user_pymol: "user_pymol.UserPyMOL" = a_user_pymol
    self._model = the_app_state.pyssa_objects_model
    
    # Setup debounce timer for selection changes
    self._selection_debounce_timer = QtCore.QTimer(self)
    self._selection_debounce_timer.setSingleShot(True)
    self._selection_debounce_timer.setInterval(50)  # 50ms debounce
    self._selection_debounce_timer.timeout.connect(self._evaluate_tree_selection)
    
    # </editor-fold>
    self._set_model()
    self._connect_signals()

  def _set_model(self) -> None:
    self._panel.tree_view.setModel(self._model)
    self._connect_selection_signal()

  def _connect_selection_signal(self) -> None:
    """Connect (or reconnect) the selection changed signal to the tree view's selection model."""
    selection_model = self._panel.tree_view.selectionModel()
    if selection_model:
        # Disconnect any existing connections to avoid duplicates
        try:
            selection_model.selectionChanged.disconnect(self._on_tree_selection_changed)
        except (RuntimeError, TypeError):
            # No previous connection exists, which is fine
            pass
        # Now connect
        selection_model.selectionChanged.connect(self._on_tree_selection_changed)
        logger.info("Selection model signal connected successfully")
    else:
        logger.warning("Selection model is None, cannot connect signal")

  def _connect_signals(self) -> None:
    self._panel.import_file_action.get_action().triggered.connect(
      self.__slot_display_import_popup
    )
    self._panel.import_prot_action.triggered.connect(self.__slot_import_protein)
    self._panel.import_seq_action.triggered.connect(self.__slot_import_sequence)
    self._panel.export_file_action.get_action().triggered.connect(
      self.__slot_export_file
    )

    # self.expand_all.clicked.connect(self._panel.tree_view.expandAll)
    # self.collapse_all.clicked.connect(self._panel.tree_view.collapseAll)

  def suppress_selection_signal(self) -> None:
    """Temporarily disconnect the tree's selectionChanged signal and stop the debounce timer.

    Call this before programmatically modifying the tree selection from the
    reverse-sync path (PyMOL → tree) to prevent re-entrant snapshot updates.
    Must be paired with a subsequent call to :meth:`restore_selection_signal`.
    """
    self._selection_debounce_timer.stop()
    selection_model = self._panel.tree_view.selectionModel()
    if selection_model:
      try:
        selection_model.selectionChanged.disconnect(self._on_tree_selection_changed)
      except (RuntimeError, TypeError):
        pass

  def restore_selection_signal(self) -> None:
    """Reconnect the tree's selectionChanged signal after a suppressed update.

    This is the counterpart to :meth:`suppress_selection_signal`.  It
    reconnects the debounced selection handling so normal tree interactions
    resume producing snapshots.
    """
    self._connect_selection_signal()

  def _on_tree_selection_changed(self) -> None:
    """Triggered on every selection change, resets the debounce timer."""
    self._selection_debounce_timer.start()

  def _evaluate_tree_selection(self) -> None:
    """Triggered by the debounce timer to safely resolve selection and emit snapshot."""
    selection_model = self._panel.tree_view.selectionModel()
    if not selection_model:
        return
        
    indexes = selection_model.selectedIndexes()
    snapshot = self._model.resolve_selection(indexes)
    self.selectionSnapshotUpdated.emit(snapshot)

  def __slot_display_import_popup(self):
    tmp_button = self._panel.get_toolbar().get_tool_button_for_action(
      self._panel.import_file_action
    )
    if tmp_button is not None and tmp_button.isVisible():
      # Position just below the button
      self._panel.import_seq_prot_menu.exec(
        tmp_button.mapToGlobal(tmp_button.rect().bottomLeft())
      )

  def __slot_import_sequence(self):
    from src.pyssa.controller.import_sequence_view_controller import ImportSequenceViewController
    from src.pyssa.internal.thread.thread_api import thread_runtime
    
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Import sequence' button on the 'Sequence Tab' was clicked.",
    )
    self._external_controller = ImportSequenceViewController(
      the_app_state=self._app_state,
      a_parent=self._panel.window()
    )
    self._external_controller.restore_ui()
    self._external_controller.get_view().show()
    
  def __slot_import_protein(self):
    from src.pyssa.controller.add_protein_view_controller import AddProteinViewController
    
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Import protein' button was clicked.",
    )
    self._external_controller = AddProteinViewController(
      the_app_state=self._app_state,
      a_parent=self._panel.window()
    )
    self._external_controller.restore_ui()
    self._external_controller.get_view().show()

  def __slot_export_file(self):
    pass
    # fixme: Legacy implementation
    # Determine selected top-level molecule name from tree view
    # selection_model = self._panel.tree_view.selectionModel()
    # if selection_model is None:
    #   return
    # selected_indexes = selection_model.selectedIndexes()
    # if not selected_indexes:
    #   QtWidgets.QMessageBox.information(
    #     self._panel,
    #     "Export",
    #     "Please select a molecule (or any of its children) in the list to export.",
    #   )
    #   return
    #
    # # Use the first selected index and climb to the top-level item to get the molecule name
    # index = selected_indexes[0]
    # model = self._panel.tree_view.model()
    # while index.isValid() and model.parent(index).isValid():
    #   index = model.parent(index)
    #
    # molecule_name = index.data(QtCore.Qt.ItemDataRole.DisplayRole) if index.isValid() else None
    # if not molecule_name:
    #   QtWidgets.QMessageBox.warning(self._panel, "Export", "Could not determine the molecule name for export.")
    #   return
    #
    # # Ask for destination file path. Provide a default filename based on the molecule name.
    # default_filename = f"{molecule_name}.pdb"
    # file_path, selected_filter = QtWidgets.QFileDialog.getSaveFileName(
    #   self._panel,
    #   "Save Molecule As",
    #   default_filename,
    #   "PDB File (*.pdb);;mmCIF File (*.cif *.mmcif);;All Files (*.*)",
    # )
    #
    # if not file_path:
    #   return
    #
    # # If no extension was provided, infer from the chosen filter.
    # if "." not in file_path.split("\\")[-1]:
    #   if "mmCIF" in selected_filter:
    #     file_path = file_path + ".cif"
    #   else:
    #     file_path = file_path + ".pdb"
    #
    # try:
    #   # Export the selected molecule by name from PyMOL
    #   self._user_pymol.get_cmd_module().save(file_path, molecule_name)
    # except Exception as e:
    #   QtWidgets.QMessageBox.critical(
    #     self._panel,
    #     "Export Failed",
    #     f"Failed to export '{molecule_name}': {e}",
    #   )
