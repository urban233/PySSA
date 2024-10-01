from PyQt6 import QtCore
from PyQt6 import QtWidgets
from PyQt6 import QtGui
from PyQt6.QtCore import Qt
from pyssa.gui.base_classes import base_dialog
from .forms.auto import auto_dialog_open_project
from pyssa.model.preference import model_definitions
from pyssa.model.util.gui_style import icons
from pyssa.model.util.gui_style import styles_utils


class DialogOpenProject(base_dialog.BaseDialog):
  """Open project dialog"""

  dialogClosed = QtCore.pyqtSignal()
  """A signal indicating that the dialog is closed."""

  def __init__(self) -> None:
    """Constructor."""
    super().__init__()
    self.ui = auto_dialog_open_project.Ui_Dialog()
    self.ui.setupUi(self)
    self.setup_ui()
    self.resize(450, 600)
    self.ui.btn_cancel.clicked.connect(self.close)
    styles_utils.color_bottom_frame_button(self.ui.btn_open_project)
    self.setWindowModality(Qt.WindowModality.WindowModal)

  def setup_ui(self) -> None:
    """Sets up the initial ui."""
    self.ui.lbl_open_status_search.setText("")
    self.ui.txt_open_search.setText("")
    self.ui.txt_open_selected_project.setText("")
    self.ui.projects_list_view.setEditTriggers(
      QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
    )
    self.ui.projects_list_view.clearSelection()
    self.ui.btn_open_project.setEnabled(False)
    icons.set_icon(self.ui.btn_help, model_definitions.IconsEnum.HELP)
    self.setWindowIcon(QtGui.QIcon(model_definitions.ModelDefinitions.PLUGIN_LOGO_FILEPATH))
    self.setWindowTitle("Open Project")
    self.setWindowFlags(
      self.windowFlags() & ~QtCore.Qt.WindowType.WindowContextHelpButtonHint
    )

  def closeEvent(self, event) -> None:
    """Closes the dialog (with the closeEvent) and emits the 'dialogClosed' signal."""
    event.accept()
    self.dialogClosed.emit()
