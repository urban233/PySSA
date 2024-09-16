from PyQt6 import QtCore
from PyQt6 import QtWidgets
from PyQt6 import QtGui
from PyQt6.QtCore import Qt

from pyssa.gui.base_classes import base_dialog
from pyssa.gui.dialog.forms.auto import auto_dialog_import_protein


class DialogImportProtein(base_dialog.BaseDialog):
  """Import protein dialog"""

  dialogClosed = QtCore.pyqtSignal()
  """A signal indicating that the dialog is closed."""

  def __init__(self) -> None:
    """Constructor."""
    super().__init__()
    self.ui = auto_dialog_import_protein.Ui_Dialog()
    self.ui.setupUi(self)
    self.setup_ui()
    self.resize(450, 600)
    self.ui.btn_cancel.clicked.connect(self.close)
    self.setWindowModality(Qt.WindowModality.WindowModal)

  def setup_ui(self) -> None:
    """Sets up the initial ui."""
    self.ui.btn_add_protein.setEnabled(False)
    self.ui.lbl_status.setText("")
    self.ui.lbl_status.setStyleSheet("""color: #ba1a1a; font-size: 11px;""")
    #styles.color_bottom_frame_button(self.ui.btn_add_protein)
    self.ui.btn_choose_protein.setToolTip("Click to add a .pdb file")
    #self.ui.btn_help.setIcon(QtGui.QIcon(":/icons/help_w200.png"))
    self.ui.btn_help.setIconSize(
      self.ui.btn_help.icon().actualSize(QtCore.QSize(30, 30))
    )
    self.ui.btn_help.setText("")
    self.ui.btn_cancel.clicked.connect(self.close)
    self.setWindowTitle("Import Protein Structure")
    self.setWindowFlags(
      self.windowFlags() & ~QtCore.Qt.WindowType.WindowContextHelpButtonHint
    )
    self.setWindowFlag(QtCore.Qt.WindowType.WindowContextHelpButtonHint, True)
    self.setWindowTitle(" ")

  def show_help(self):
    # This method is called when the help button is clicked
    QtWidgets.QMessageBox.information(self, 'Help', 'This is the help message.')

  def event(self, event):
    if event.type() == QtCore.QEvent.Type.WhatsThisClicked:
      self.show_help()
      return True  # Event is handled
    return super().event(event)

  def closeEvent(self, event) -> None:
    """Closes the dialog (with the closeEvent) and emits the 'dialogClosed' signal."""
    event.accept()
    self.dialogClosed.emit()
