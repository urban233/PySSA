import logging
from typing import TYPE_CHECKING
from PyQt6 import QtWidgets
from PyQt6 import QtCore
from pyssa.gui.base_classes import base_controller
from pyssa.model.data_classes import ligand_preview
from pyssa.model.util import exception
from pyssa.gui.custom_widgets import molecule_naming
from pyssa.model.pyssa_logging import default_logging


logger = default_logging.setup_logger(__file__)

__docformat__ = "google"

if TYPE_CHECKING:
  from pyssa.gui.dialog import dialog_preview_import_ligand


class PreviewImportLigandController(base_controller.BaseController):
  """Controller for the import ligand dialog."""

  component_task = QtCore.pyqtSignal(tuple)

  def __init__(self, a_dialog, a_list_of_ligand_previews: list["ligand_preview.LigandPreview"]):
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
    self._dialog: "dialog_preview_import_ligand.DialogPreviewImportLigand" = a_dialog
    """The dialog that the controller should work with."""
    self.connect_all_signals()
    self.ligand_previews: list["ligand_preview.LigandPreview"] = a_list_of_ligand_previews
    self._prepare_scroll_area_with_molecules()

  def _prepare_scroll_area_with_molecules(self) -> None:
    # Create a container widget that will hold the custom widgets
    self.container = QtWidgets.QFrame()
    self.container_layout = QtWidgets.QVBoxLayout(self.container)
    # Add multiple instances of CustomWidget to the container layout
    for tmp_molecule in self.ligand_previews:
      tmp_molecule.name = tmp_molecule.get_ligand_name()
      tmp_molecule.create_molecule_image()
      custom_widget = molecule_naming.MoleculeNaming(an_image_filepath=tmp_molecule.image_filepath, a_name=tmp_molecule.name, parent=self._dialog)
      self.container_layout.addWidget(custom_widget)
    # Set the container as the widget for the scroll area
    self._dialog.ui.scroll_area_molecules.setWidget(self.container)
    self._dialog.ui.scroll_area_molecules.setWidgetResizable(True)

  def connect_all_signals(self):
    """Connects all signals with their appropriate slots."""
    pass

  def get_dialog(self) -> QtWidgets.QDialog:
    """Gets the dialog of the controller."""
    return self._dialog
  
  def restore_ui(self) -> None:
    """Restores the UI to default values."""
    self._dialog.setup_ui()
