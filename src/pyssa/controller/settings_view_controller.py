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
"""Module for the settings view controller."""
from typing import TYPE_CHECKING
import logging

from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.util import exception
from src.pyssa.util import gui_utils
from src.pyssa.logging_pyssa import log_handlers, log_levels
from src.pyssa.util import constants

from src.pyssa.gui.ui.views import settings_view
if TYPE_CHECKING:
  from src.pyssa.gui import app_state

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class SettingsViewController(QtCore.QObject):
  """Class for the SettingsViewController."""

  def __init__(
      self, the_app_state: "app_state.AppState", a_parent=None
  ) -> None:
    """Constructor.

    Args:
        the_app_state (app_state.AppState): The AppState object.
        a_parent: Parent widget to pass to the view.

    Raises:
        exception.IllegalArgumentError: If `the_app_state` is None.
    """
    # <editor-fold desc="Checks">
    if the_app_state is None:
      logger.error("the_app_state is None.")
      raise exception.IllegalArgumentError("the_app_state is None.")

    # </editor-fold>

    super().__init__()
    self._app_state = the_app_state
    self._parent = a_parent
    self._settings_manager = the_app_state._settings_manager
    self._view = settings_view.SettingsView(a_parent)
    self._initialize_ui()
    self.restore_default_view()
    self._connect_all_ui_elements_to_slot_functions()
    if self._app_state.job_scheduler.has_running_jobs():
      self._view.ui.btn_workspace_dir.setEnabled(False)
    else:
      self._view.ui.btn_workspace_dir.setEnabled(True)

  def get_view(self):
    return self._view

  # <editor-fold desc="Util methods">
  def _open_help_for_dialog(self) -> None:
    """Opens the help dialog."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Help' button was clicked."
    )
    # self._interface_manager.help_manager.open_pyssa_settings_page()

  def restore_default_view(self) -> None:
    """Restores the UI."""
    self._view.ui.tabWidget.setCurrentIndex(0)
    self._view.ui.txt_workspace_dir.setText(
      str(self._settings_manager.settings.get_workspace_path())
    )
    self._view.ui.spb_cycles.setValue(
      int(self._settings_manager.settings.get_cycles())
    )
    self._view.ui.dspb_cutoff.setValue(
      float(self._settings_manager.settings.get_cutoff())
    )
    self._view.ui.box_bg_color.setCurrentIndex(
      self._view.ui.box_bg_color.findText(
        self._settings_manager.settings.image_background_color
      ),
    )
    if self._settings_manager.settings.image_renderer == "0":
      self._view.ui.box_renderer.setCurrentIndex(0)
    else:
      self._view.ui.box_renderer.setCurrentIndex(1)
    self._view.ui.box_ray_trace_mode.setCurrentIndex(
      self._settings_manager.settings.image_ray_trace_mode
    )
    self._view.ui.box_ray_texture.setCurrentIndex(
      self._settings_manager.settings.image_ray_texture
    )
    self._view.ui.cb_pyssa_expert_mode.setChecked(
      self._settings_manager.settings.pyssa_expert_mode
    )

  # </editor-fold>

  def _initialize_ui(self) -> None:
    """Initializes the user interface by setting values and configuring options."""
    # combo box BgColor
    item_list_bg_color = [
        "black",
        "white",
    ]
    gui_utils.fill_combo_box(self._view.ui.box_bg_color, item_list_bg_color)
    # combo box Renderer
    item_list_renderer = [
        "default renderer",
        "PyMOL internal renderer",
    ]
    gui_utils.fill_combo_box(self._view.ui.box_renderer, item_list_renderer)
    # combo box RayTraceMode
    item_list_ray_trace_mode = [
        "normal color",
        "normal color + black outline",
        "black outline only",
        "quantized color + black outline",
    ]
    gui_utils.fill_combo_box(
        self._view.ui.box_ray_trace_mode, item_list_ray_trace_mode
    )
    # combo box Ray Texture
    item_list_ray_texture = [
        "None",
        "Matte 1",
        "Matte 2",
        "Swirl 1",
        "Fiber",
    ]
    gui_utils.fill_combo_box(
        self._view.ui.box_ray_texture, item_list_ray_texture
    )
    # customize spin boxes
    self._view.ui.spb_cycles.setMinimum(0)
    # self._view.ui.spb_cycles.setMaximum(20) # fixme: is a maximum needed?
    self._view.ui.spb_cycles.setSingleStep(1)
    self._view.ui.dspb_cutoff.setMinimum(0.00)
    self._view.ui.dspb_cutoff.setMaximum(20.00)
    self._view.ui.dspb_cutoff.setSingleStep(0.1)
    self._view.ui.txt_workspace_dir.setEnabled(False)
    # Hide the PyMOL settings tab for now
    self._view.ui.tabWidget.setTabVisible(3, False)

  def _connect_all_ui_elements_to_slot_functions(self) -> None:
    """Connects all UI elements to their corresponding slot functions in the class."""
    self._view.ui.btn_workspace_dir.clicked.connect(self.choose_workspace_dir)
    self._view.ui.btn_ok.clicked.connect(self.ok_dialog)
    self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)

  def choose_workspace_dir(self) -> None:
    """Opens a QFileDialog to choose a workspace directory."""
    logger.log(
        log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
        "'Choose workspace' button was clicked.",
    )
    gui_utils.choose_directory(
        self._view,
        self._view.ui.txt_workspace_dir,
    )

  def ok_dialog(self) -> None:
    """Sets all settings from the gui elements into the settings object and closes the dialog window."""
    logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'OK' button was clicked.")
    self._settings_manager.settings.set_workspace_path(
        self._view.ui.txt_workspace_dir.text()
    )
    self._settings_manager.settings.set_cycles(self._view.ui.spb_cycles.value())
    self._settings_manager.settings.set_cutoff(
        self._view.ui.dspb_cutoff.value()
    )
    # self._settings_manager.settings.color_vision_mode = self._view.ui.cb_color_vision_mode.currentText()
    self._settings_manager.settings.image_background_color = (
        self._view.ui.box_bg_color.currentText()
    )
    if self._view.ui.box_renderer.currentText() == "default renderer":
      self._settings_manager.settings.image_renderer = "0"
    else:
      self._settings_manager.settings.image_renderer = "-1"
    self._settings_manager.settings.image_ray_trace_mode = (
        self._view.ui.box_ray_trace_mode.currentIndex()
    )
    self._settings_manager.settings.image_ray_texture = (
        self._view.ui.box_ray_texture.currentIndex()
    )
    if self._settings_manager.settings.pyssa_expert_mode != self._view.ui.cb_pyssa_expert_mode.isChecked():
      tmp_activated_expert_mode = True
    else:
      tmp_activated_expert_mode = False
    self._settings_manager.settings.pyssa_expert_mode = (
      self._view.ui.cb_pyssa_expert_mode.isChecked()
    )

    self._settings_manager.settings.serialize_settings()
    logging.info("Settings were successfully saved.")
    self._view.close()

    if tmp_activated_expert_mode:
      tmp_msg = """
      You have activated the PySSA expert mode. You must keep this in mind:\n
      #1) Features can and will break data integrity if not used carefully.\n
      #2) Think before you use any of the features.\n
      #3) With great power comes great responsibility.\n
      """
      tmp_dialog = custom_message_box.CustomMessageBoxOk(
        tmp_msg,
        "Activate Expert Mode",
        custom_message_box.CustomMessageBoxIcons.WARNING.value,
      )
      tmp_dialog.exec()

    # Rebuild the workspace model inside the app state and trigger a UI update
    if hasattr(self._app_state, "_build_workspace_model"):
      self._app_state._build_workspace_model()
    self._app_state._on_state_changed()
