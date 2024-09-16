import logging
import pathlib
from typing import TYPE_CHECKING
from PyQt6 import QtWidgets
from PyQt6 import QtCore
from rdkit import Chem

from pyssa.gui.base_classes import base_controller
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging
from pyssa.gui.dialog import dialog_preview_import_ligand
from pyssa.gui.control import preview_import_ligand_controller
from pyssa.model.data_classes import ligand_preview

logger = default_logging.setup_logger(__file__)

__docformat__ = "google"

if TYPE_CHECKING:
  from pyssa.gui.dialog import dialog_import_ligand


class ImportLigandController(base_controller.BaseController):
  """Controller for the import ligand dialog."""

  component_task = QtCore.pyqtSignal(tuple)

  def __init__(self, a_dialog):
    """Constructor.
    
    Args:
      a_dialog: a dialog instance to be managed by the controller

    Raises:
      exception.NoneValueError: If `a_dialog` is None.
    """
    # <editor-fold desc="Checks">
    if a_dialog is None:
      default_logging.append_to_log_file(logger, "a_dialog is None.", logging.ERROR)
      raise exception.NoneValueError("a_dialog is None.")
    # </editor-fold>
    super().__init__()
    self._dialog: "dialog_import_ligand.DialogImportLigand" = a_dialog
    """The dialog that the controller should work with."""
    self.connect_all_signals()

  def connect_all_signals(self):
    """Connects all signals with their appropriate slots."""
    self._dialog.ui.btn_choose_ligand_file.clicked.connect(self.__slot_open_ligand_from_filesystem)

  def get_dialog(self) -> QtWidgets.QDialog:
    """Gets the dialog of the controller."""
    return self._dialog
  
  def restore_ui(self) -> None:
    """Restores the UI to default values."""
    self._dialog.setup_ui()

  def __slot_open_ligand_from_filesystem(self) -> None:
    """Loads a protein from the filesystem into the textbox."""
    default_logging.append_to_log_file(
      logger, "'Open ligand from filesystem' button was clicked."
    )
    try:
      file_name = QtWidgets.QFileDialog.getOpenFileName(
        self._dialog,
        "Open existing ligand",
        QtCore.QDir.homePath(),
        "SDF Files (*.sdf)",
      )
      if file_name == ("", ""):
        self._dialog.ui.lbl_status.setText("No file has been selected.")
      else:
        # display path in text box
        self._dialog.ui.txt_import_ligand.setText(str(file_name[0]))
        self._dialog.ui.btn_import_ligand.setEnabled(True)
        tmp_supplier = Chem.SDMolSupplier(r"C:\Users\student\Downloads\ATP_ideal.sdf")
        tmp_molecule: Chem.Mol = tmp_supplier[0]
        tmp_dialog = dialog_preview_import_ligand.DialogPreviewImportLigand()
        tmp_controller = preview_import_ligand_controller.PreviewImportLigandController(
          tmp_dialog,
          [
            ligand_preview.LigandPreview(tmp_molecule, pathlib.Path(r"C:\Users\student\Downloads\ATP_ideal.sdf")),
            ligand_preview.LigandPreview(tmp_molecule, pathlib.Path(r"C:\Users\student\Downloads\ATP_ideal.sdf")),
            ligand_preview.LigandPreview(tmp_molecule, pathlib.Path(r"C:\Users\student\Downloads\ATP_ideal.sdf"))
          ]
        )
        tmp_controller.get_dialog().exec()
    except FileNotFoundError:
      self._dialog.ui.lbl_status.setText("Opening the ligand structure failed!")
