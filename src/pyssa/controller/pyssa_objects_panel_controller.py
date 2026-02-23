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
import pathlib

from Bio import SeqIO

from src.pyssa.gui.ui.views import pyssa_objects_panel
from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.gui import user_pymol, app_state
from src.pyssa.gui.qt import QtCore, QtWidgets
from src.pyssa.io_pyssa import bio_data
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

    # Holds the most recently resolved selection snapshot so that slot methods
    # can access the current selection without re-querying the Qt model.
    self._current_snapshot: "SelectionSnapshot | None" = None
    
    # Setup debounce timer for selection changes
    self._selection_debounce_timer = QtCore.QTimer(self)
    self._selection_debounce_timer.setSingleShot(True)
    self._selection_debounce_timer.setInterval(50)  # 50ms debounce
    self._selection_debounce_timer.timeout.connect(self._evaluate_tree_selection)
    
    # </editor-fold>
    self._set_model()
    self._connect_all_signals_with_their_slots()

  @property
  def _model(self) -> "psa_objects_model.PSAObjectsModel":
      return self._app_state.pyssa_objects_model

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

  def _connect_all_signals_with_their_slots(self) -> None:
    self._panel.expand_all.get_action().triggered.connect(
      self._panel.tree_view.expandAll
    )
    self._panel.collapse_all.get_action().triggered.connect(
      self._panel.tree_view.collapseAll
    )
    self._panel.import_file_action.get_action().triggered.connect(
      self.__slot_display_import_popup
    )
    self._panel.import_prot_action.triggered.connect(self.__slot_import_protein)
    self._panel.import_seq_action.triggered.connect(self.__slot_import_sequence)
    self._panel.export_file_action.get_action().triggered.connect(
      self.__slot_export_file
    )
    self._panel.add_sequence_action.get_action().triggered.connect(
      self.__slot_add_sequence
    )
    self._panel.delete_object_action.get_action().triggered.connect(
      self.__slot_delete_object
    )

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
    self._current_snapshot = snapshot
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

  def __slot_add_sequence(self):
    from src.pyssa.controller import add_sequence_view_controller

    logger.log(
      log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
      "'Import sequence' button on the was clicked.",
    )
    self._external_controller = add_sequence_view_controller.AddSequenceViewController(
      the_app_state=self._app_state,
      a_parent=self._panel.window()
    )
    self._external_controller.restore_default_view()
    self._external_controller.get_view().show()

  def __slot_export_file(self) -> None:
    """Export selected objects to the filesystem.

    Sequences are exported as FASTA files; proteins are exported as PDB files.
    Protein pairs are silently skipped because they cannot be meaningfully
    exported as a single file type.

    The user is prompted to choose an output directory via a file dialog.
    All export operations run synchronously on the main thread.
    """
    logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Export file' button was clicked.")

    if not self._app_state.has_open_project():
      logger.warning("Export requested but no project is open.")
      return

    snapshot = self._current_snapshot
    if snapshot is None:
      logger.warning("Export requested but no selection snapshot is available.")
      return

    has_exportable = bool(snapshot.raw_sequences or snapshot.raw_standalone_proteins)
    if not has_exportable:
      logger.info("Export requested but no exportable objects (sequences or standalone proteins) are selected.")
      return

    # Ask the user for an output directory.
    output_dir = QtWidgets.QFileDialog.getExistingDirectory(
      self._panel.window(),
      "Select Export Directory",
      "",
    )
    if not output_dir:
      logger.info("Export cancelled by user (no directory selected).")
      return

    output_path = pathlib.Path(output_dir)
    project = self._app_state.require_project()
    exported_count = 0

    try:
      # Export sequences as FASTA.
      for seq_name in snapshot.raw_sequences:
        seq_record = project.search_sequence(seq_name)
        if seq_record is None:
          logger.warning("Sequence '%s' not found in project during export.", seq_name)
          continue
        fasta_filepath = output_path / f"{seq_name}.fasta"
        with open(fasta_filepath, "w") as fh:
          SeqIO.write(seq_record, fh, "fasta")
        logger.info("Exported sequence '%s' to '%s'.", seq_name, fasta_filepath)
        exported_count += 1

      # Export standalone proteins as PDB.
      for protein_obj in snapshot.raw_standalone_proteins:
        protein_name = protein_obj.get_molecule_object()
        pdb_filepath = output_path / f"{protein_name}.pdb"
        pdb_data = protein_obj.get_pdb_data()
        if not pdb_data:
          logger.warning("Protein '%s' has no PDB data, skipping export.", protein_name)
          continue
        bio_data.build_pdb_file(pdb_data, str(pdb_filepath))
        logger.info("Exported protein '%s' to '%s'.", protein_name, pdb_filepath)
        exported_count += 1

    except Exception as export_error:
      logger.error("Export failed with an unexpected error: %s", export_error, exc_info=True)
      msg_box = custom_message_box.CustomMessageBoxOk(
        f"An error occurred during export:\n{export_error}",
        "Export Failed",
        custom_message_box.CustomMessageBoxIcons.WARNING.value,
      )
      msg_box.exec()
      return

    if exported_count > 0:
      msg_box = custom_message_box.CustomMessageBoxOk(
        f"Successfully exported {exported_count} file(s) to:\n{output_dir}",
        "Export Successful",
        custom_message_box.CustomMessageBoxIcons.INFORMATION.value,
      )
      msg_box.exec()

  def __slot_delete_object(self) -> None:
    """Delete the currently selected objects from the model, project, and database.

    The deletion covers sequences, standalone proteins, and protein pairs.
    Before any deletion, a pre-flight check verifies that no selected
    standalone protein is still referenced by an existing protein pair.
    If such a conflict is found, the user is informed and the deletion is
    aborted without making any changes.

    The user must confirm the deletion via a custom message box before the
    operation proceeds.
    """
    logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Delete object' button was clicked.")

    if not self._app_state.has_open_project():
      logger.warning("Delete requested but no project is open.")
      return

    snapshot = self._current_snapshot
    if snapshot is None:
      logger.warning("Delete requested but no selection snapshot is available.")
      return

    has_deletable = bool(
      snapshot.raw_sequences or
      snapshot.raw_standalone_proteins or
      snapshot.raw_protein_pairs
    )
    if not has_deletable:
      logger.info("Delete requested but no deletable objects are selected.")
      return

    project = self._app_state.require_project()

    # Pre-flight check: prevent deletion of a protein that is still part of a pair.
    conflicting_proteins: list[str] = []
    for protein_obj in snapshot.raw_standalone_proteins:
      protein_name = protein_obj.get_molecule_object()
      if project.check_if_protein_is_in_any_protein_pair(protein_name):
        # Collect the name(s) of the pair(s) that reference this protein.
        pair_names = [
          pair.name
          for pair in project.protein_pairs
          if (
            pair.protein_1.get_molecule_object() == protein_name
            or pair.protein_2.get_molecule_object() == protein_name
          )
        ]
        for pair_name in pair_names:
          conflicting_proteins.append(f"'{protein_name}' (used in pair '{pair_name}')")

    if conflicting_proteins:
      conflict_list = "\n".join(f"  \u2022 {entry}" for entry in conflicting_proteins)
      msg_box = custom_message_box.CustomMessageBoxOk(
        f"The following protein(s) cannot be deleted because they are still part of a protein pair.\n"
        f"Please delete the protein pair first, then retry:\n\n{conflict_list}",
        "Cannot Delete Protein",
        custom_message_box.CustomMessageBoxIcons.WARNING.value,
      )
      msg_box.exec()
      return

    # Ask the user to confirm the deletion.
    confirm_box = custom_message_box.CustomMessageBoxDelete(
      "Are you sure you want to permanently delete the selected object(s)?\n"
      "This action cannot be undone.",
      "Confirm Deletion",
      custom_message_box.CustomMessageBoxIcons.DANGEROUS.value,
    )
    confirm_box.exec()
    if not confirm_box.response:
      logger.info("Deletion cancelled by user.")
      return

    hot_db = self._app_state.hot_db
    if hot_db is None:
      logger.error("Delete requested but no hot database is available.")
      return

    try:
      # Delete protein pairs first to avoid leaving orphaned pair records.
      for pair in list(snapshot.raw_protein_pairs):
        project.delete_specific_protein_pair(pair.name)
        hot_db.delete_protein_pair_full(pair.get_id())
        self._model.remove_protein_pair(pair)
        logger.info("Deleted protein pair '%s'.", pair.name)

      # Delete standalone proteins (conflicts were cleared by the pre-flight check above).
      for protein_obj in list(snapshot.raw_standalone_proteins):
        protein_name = protein_obj.get_molecule_object()
        logger.info("Attempting to delete protein '%s' with ID=%s", protein_name, protein_obj.get_id())

        # First delete from project (removes from in-memory list)
        project.delete_specific_protein(protein_name)
        logger.info("Deleted protein '%s' from project", protein_name)

        # Then delete from database
        hot_db.delete_protein_full(protein_obj.get_id())
        logger.info("Deleted protein '%s' (ID=%s) from database", protein_name, protein_obj.get_id())

        # Finally remove from UI model
        self._model.remove_protein(protein_obj)
        logger.info("Deleted protein '%s' from model/UI", protein_name)

      # Delete sequences.
      for seq_name in list(snapshot.raw_sequences):
        project.delete_specific_sequence(seq_name)
        hot_db.delete_sequence(seq_name)
        self._model.remove_sequence(seq_name)
        logger.info("Deleted sequence '%s'.", seq_name)

    except Exception as delete_error:
      logger.error("Deletion failed with an unexpected error: %s", delete_error, exc_info=True)
      msg_box = custom_message_box.CustomMessageBoxOk(
        f"An error occurred during deletion:\n{delete_error}",
        "Deletion Failed",
        custom_message_box.CustomMessageBoxIcons.WARNING.value,
      )
      msg_box.exec()
      return

    # Clear the cached snapshot to avoid stale object references after deletion.
    self._current_snapshot = None
    self.selectionSnapshotUpdated.emit(SelectionSnapshot())

    msg_box = custom_message_box.CustomMessageBoxOk(
      "Selected object(s) were successfully deleted.",
      "Deletion Successful",
      custom_message_box.CustomMessageBoxIcons.INFORMATION.value,
    )
    msg_box.exec()
