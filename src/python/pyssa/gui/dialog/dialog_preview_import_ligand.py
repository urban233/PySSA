from PyQt6 import QtCore
from PyQt6 import QtWidgets
from PyQt6.QtCore import Qt

from pyssa.gui.base_classes import base_dialog
from pyssa.model.util.gui_style import icons
from pyssa.model.preference import model_definitions
from pyssa.gui.dialog.forms.auto import auto_dialog_preview_import_ligand


class DialogPreviewImportLigand(base_dialog.BaseDialog):
  """Preview import ligand dialog"""

  dialogClosed = QtCore.pyqtSignal()
  """A signal indicating that the dialog is closed."""

  def __init__(self) -> None:
    """Constructor."""
    super().__init__()
    self.ui = auto_dialog_preview_import_ligand.Ui_Dialog()
    self.ui.setupUi(self)
    self.setup_ui()
    self.resize(450, 600)
    self.ui.btn_cancel.clicked.connect(self.close)
    self.setWindowModality(Qt.WindowModality.WindowModal)

  def setup_ui(self) -> None:
    """Sets up the initial ui."""
    self.ui.btn_import.setEnabled(False)
    #styles.color_bottom_frame_button(self.ui.btn_add_protein)
    icons.set_icon(self.ui.btn_help, model_definitions.IconsEnum.HELP)
    self.ui.btn_cancel.clicked.connect(self.close)
    self.setWindowTitle("Preview Import Ligand Structures")
    self.setWindowFlags(
      self.windowFlags() & ~QtCore.Qt.WindowType.WindowContextHelpButtonHint
    )

  def closeEvent(self, event) -> None:
    """Closes the dialog (with the closeEvent) and emits the 'dialogClosed' signal."""
    event.accept()
    self.dialogClosed.emit()

  def _close_dialog(self) -> None:
    """Closes the dialog and emits the 'dialogClosed' signal."""
    self.close()
    self.dialogClosed.emit()
