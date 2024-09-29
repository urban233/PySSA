from PyQt6 import QtCore
from PyQt6 import QtWidgets
from PyQt6.QtCore import Qt

from pyssa.gui.base_classes import base_dialog
from pyssa.model.util.gui_style import styles_utils
from pyssa.model.util.gui_style import icons
from pyssa.model.preference import model_definitions
from pyssa.gui.dialog.forms.auto import auto_dialog_import_sequence


class DialogFastaFileImportPreview(base_dialog.BaseDialog):
  """Fasta file import preview dialog"""

  dialogClosed = QtCore.pyqtSignal()
  """A signal indicating that the dialog is closed."""

  def __init__(self) -> None:
    """Constructor."""
    super().__init__()
    self.ui = auto_dialog_import_sequence.Ui_Dialog()
    self.ui.setupUi(self)
    self.setup_ui()
    self.setWindowModality(Qt.WindowModality.WindowModal)

  def setup_ui(self) -> None:
    """Sets up the initial ui."""
    self.ui.lbl_status.setText("")
    self.ui.lbl_status.setStyleSheet("""color: #ba1a1a; font-size: 11px;""")
    self.ui.btn_choose_fasta_file.setToolTip("Click to add a .fasta file")
    self.ui.btn_help.setText("")
    self.ui.btn_cancel.clicked.connect(self.close)
    self.setWindowTitle("FASTA File Import Preview")
    self.setWindowFlags(
      self.windowFlags() & ~QtCore.Qt.WindowType.WindowContextHelpButtonHint
    )
    icons.set_icon(self.ui.btn_help, model_definitions.IconsEnum.HELP)
    #styles_utils.set_stylesheet(self)
    self.resize(900, 600)
    #styles_utils.color_bottom_frame_button(self.ui.btn_import_sequence)

  def closeEvent(self, event) -> None:
    """Closes the dialog (with the closeEvent) and emits the 'dialogClosed' signal."""
    event.accept()
    self.dialogClosed.emit()

  def _close_dialog(self) -> None:
    """Closes the dialog and emits the 'dialogClosed' signal."""
    self.close()
    self.dialogClosed.emit()
