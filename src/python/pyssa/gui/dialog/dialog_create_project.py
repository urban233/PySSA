from PyQt6 import QtCore, QtGui
from PyQt6 import QtWidgets
from PyQt6.QtCore import Qt

from pyssa.gui.base_classes import base_dialog
from pyssa.model.preference import model_definitions
from pyssa.model.util.gui_style import icons, styles_utils
from pyssa.gui.dialog.forms.auto import auto_dialog_create_project


class DialogCreateProject(base_dialog.BaseDialog):
  """Create project dialog"""

  dialogClosed = QtCore.pyqtSignal()
  """A signal indicating that the dialog is closed."""

  def __init__(self) -> None:
    """Constructor."""
    super().__init__()
    self.ui = auto_dialog_create_project.Ui_Dialog()
    self.ui.setupUi(self)
    self.setup_ui()
    self.ui.btn_cancel.clicked.connect(self.close)
    styles_utils.color_bottom_frame_button(self.ui.btn_new_create_project)
    self.setWindowModality(Qt.WindowModality.WindowModal)

  def setup_ui(self) -> None:
    """Sets up the initial ui."""
    self.ui.txt_new_project_name.clear()
    self.ui.txt_new_project_name.setStyleSheet("""QLineEdit {color: #000000; border-color: #DCDBE3;}""")
    self.ui.lbl_new_status_project_name.setText("")
    icons.set_icon(self.ui.btn_help, model_definitions.IconsEnum.HELP)
    self.ui.btn_new_create_project.setEnabled(False)
    self.setWindowIcon(QtGui.QIcon(model_definitions.ModelDefinitions.PLUGIN_LOGO_FILEPATH))
    self.setWindowTitle("Create Project")
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
