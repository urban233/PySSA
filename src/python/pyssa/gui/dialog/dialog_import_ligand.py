from PyQt6 import QtCore
from PyQt6 import QtWidgets
from PyQt6.QtCore import Qt

from pyssa.gui.base_classes import base_dialog
from pyssa.gui.dialog.forms.auto import auto_dialog_import_ligand


class DialogImportLigand(base_dialog.BaseDialog):
  """Import ligand dialog"""

  dialogClosed = QtCore.pyqtSignal()
  """A signal indicating that the dialog is closed."""

  def __init__(self) -> None:
    """Constructor."""
    super().__init__()
    self.ui = auto_dialog_import_ligand.Ui_Dialog()
    self.ui.setupUi(self)
    self.setup_ui()
    self.resize(450, 600)
    self.ui.btn_cancel.clicked.connect(self.close)
    self.setWindowModality(Qt.WindowModality.WindowModal)

  def setup_ui(self) -> None:
    """Sets up the initial ui."""
    self.ui.btn_import_ligand.setEnabled(False)
    self.ui.lbl_status.setText("")
    self.ui.lbl_status.setStyleSheet("""color: #ba1a1a; font-size: 11px;""")
    #styles.color_bottom_frame_button(self.ui.btn_add_protein)
    self.ui.btn_choose_ligand_file.setToolTip("Click to choose a .sdf file")
    #self.ui.btn_help.setIcon(QtGui.QIcon(":/icons/help_w200.png"))
    self.ui.btn_help.setIconSize(
      self.ui.btn_help.icon().actualSize(QtCore.QSize(30, 30))
    )
    self.ui.btn_help.setText("")
    self.ui.btn_cancel.clicked.connect(self.close)
    self.setWindowTitle("Import Ligand Structure")
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
