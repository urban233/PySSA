from PyQt6 import QtCore, QtGui
from PyQt6 import QtWidgets
from PyQt6.QtCore import Qt

from pyssa.gui.base_classes import base_dialog
from pyssa.model.preference import model_definitions
from pyssa.model.util.gui_style import icons
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
    self.setWindowModality(Qt.WindowModality.WindowModal)

  def setup_ui(self) -> None:
    """Sets up the initial ui."""
    self.ui.txt_new_project_name.clear()
    icons.set_icon(self.ui.btn_help, model_definitions.IconsEnum.HELP)
    self.ui.btn_new_create_project.setEnabled(False)
    self.setWindowIcon(QtGui.QIcon(model_definitions.ModelDefinitions.PLUGIN_LOGO_FILEPATH))
    #styles_util.color_bottom_frame_button(self.ui.btn_new_create_project)
    #styles_utils.set_stylesheet(self)
    # self.setStyleSheet(
    #   """
    #   QDialog { background: white; }
    #   QFrame#frame_bottom {
    #       background-color: #f7f8fa;
    #       border-style: solid;
    #       border-width: 1px;
    #       border-radius: 6px;
    #       border-color: qlineargradient(x1:0, y1:0, x2:0, y2:1, stop:0 #f9f9f9, stop:1 #f0f0f0);;
    #       border-top-color: #ebecf0;
    #       border-top-left-radius: 0px;
    #       border-top-right-radius: 0px;
    #   }
    #   QLabel#lbl_description {
    #       color: #367AF6;
    #       font: bold;
    #       font-size: 16px;
    #       font-family: "Segoe UI Variable";
    #   }
    #   """
    # )
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
