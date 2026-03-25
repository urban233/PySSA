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
"""Module for the help manager."""
import logging
import pathlib
from typing import Optional

from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.ui.views import help_view
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import constants
from src.pyssa.controller.help_dialog_registry import help_registry, HelpDialog

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class HelpManager:
  """Class for the HelpManager.

  This class manages help dialogs using a centralized registry system.
  Use open_help_dialog(dialog_id) to open any registered help page.
  """

  # <editor-fold desc="Class attributes">
  project_root_path = pathlib.Path(__file__).parent.absolute().parent.parent.parent
  """The Path to the project root directory."""

  docs_html_path = constants.DOCS_PATH
  """The Path to the docs html directory."""

  docs_help_html_path = pathlib.Path(docs_html_path, "help")
  """The Path to the docs help html directory."""

  index_html_filename = "index.html"
  """Generic file name for the index pages."""

  html_suffix = ".html"
  """Generic HTML file suffix."""
  # </editor-fold>

  def __init__(self) -> None:
    """Constructor."""
    super().__init__()
    self._view = help_view.BasicBrowserView()
    self._registry = help_registry

  def open_help_dialog(self, dialog_id: str) -> None:
    """Opens a help dialog by its registry ID.

    Args:
        dialog_id: The unique identifier of the help dialog to open

    Raises:
        ValueError: If the dialog_id is not found in the registry
    """
    html_path = self._registry.get_html_path(dialog_id)
    if html_path is None:
      logger.error(f"Help dialog with ID '{dialog_id}' not found in registry")
      raise ValueError(f"Help dialog with ID '{dialog_id}' not found")

    self.change_url(str(html_path))
    logger.debug(f"Opened help dialog: {dialog_id}")

  def get_registry(self) -> "help_registry":
    """Returns the help dialog registry for advanced usage.

    Returns:
        The HelpDialogRegistry instance
    """
    return self._registry

  def change_url(self, an_url: str) -> None:
    """Changes the URL of the browser window."""
    if an_url.find("\\"):
      tmp_url = an_url.replace("\\", "/")
    else:
      tmp_url = an_url

    self._view.browser.setUrl(QtCore.QUrl(tmp_url))
    self._view.show()
    self._view.activateWindow()

  # def display_help_over_model_dialog(self, the_model_dialog, an_url) -> QWebEngineView:
  #   tmp_view = help_view.BasicBrowserView(the_model_dialog)
  #   if an_url.find("\\"):
  #     tmp_url = an_url.replace("\\", "/")
  #   else:
  #     tmp_url = an_url
  #   tmp_view.browser.setUrl(QUrl(tmp_url))
  #   tmp_view.setModal(False)
  #   return tmp_view

  # <editor-fold desc="Main view help pages">
  def open_general_help_page(self) -> None:
    """Opens the help dialog on the general help page."""
    self.open_help_dialog("general")

  def open_sequences_tab_page(self) -> None:
    """Opens the help dialog on the sequences tab page."""
    self.open_help_dialog("sequences_tab")

  def open_additional_sequence_information_page(self) -> None:
    """Opens the help dialog on the additional sequence information page."""
    self.open_help_dialog("sequence_additional_info")

  def open_sequence_import_page(self) -> None:
    """Opens the help dialog on the import sequence page."""
    self.open_help_dialog("sequence_import")

  def open_sequence_add_page(self) -> None:
    """Opens the help dialog on the add sequence page."""
    self.open_help_dialog("sequence_add")

  def open_sequence_save_page(self) -> None:
    """Opens the help dialog on the save sequence page."""
    self.open_help_dialog("sequence_save")

  def open_sequence_delete_page(self) -> None:
    """Opens the help dialog on the delete sequence page."""
    self.open_help_dialog("sequence_delete")

  def open_proteins_tab_page(self) -> None:
    """Opens the help dialog on the proteins tab page."""
    self.open_help_dialog("proteins_tab")

  def open_protein_import_page(self) -> None:
    """Opens the help dialog on the import protein page."""
    self.open_help_dialog("protein_import")

  def open_protein_save_page(self) -> None:
    """Opens the help dialog on the save protein page."""
    self.open_help_dialog("protein_save")

  def open_protein_delete_page(self) -> None:
    """Opens the help dialog on the delete protein page."""
    self.open_help_dialog("protein_delete")

  def open_protein_pymol_scene_configuration_page(self) -> None:
    """Opens the help dialog on the protein pymol scene config page."""
    self.open_help_dialog("protein_pymol_scene_config")

  def open_protein_load_session_page(self) -> None:
    """Opens the help dialog on the open protein pymol session page."""
    self.open_help_dialog("protein_load_session")

  def open_protein_add_scene_page(self) -> None:
    """Opens the help dialog on the protein add scene page."""
    self.open_help_dialog("protein_add_scene")

  def open_protein_update_scene_page(self) -> None:
    """Opens the help dialog on the protein update scene page."""
    self.open_help_dialog("protein_update_scene")

  def open_protein_delete_scene_page(self) -> None:
    """Opens the help dialog on the protein delete scene page."""
    self.open_help_dialog("protein_delete_scene")

  def open_protein_pairs_tab_page(self) -> None:
    """Opens the help dialog on the protein pairs tab page."""
    self.open_help_dialog("protein_pairs_tab")

  def open_protein_pair_delete_page(self) -> None:
    """Opens the help dialog on the delete protein pair page."""
    self.open_help_dialog("protein_pair_delete")

  def open_protein_pair_pymol_scene_configuration_page(self) -> None:
    """Opens the help dialog on the protein pair pymol scene config page."""
    self.open_help_dialog("protein_pair_pymol_scene_config")

  def open_protein_pair_load_session_page(self) -> None:
    """Opens the help dialog on the open protein pair pymol session page."""
    self.open_help_dialog("protein_pair_load_session")

  def open_protein_pair_add_scene_page(self) -> None:
    """Opens the help dialog on the protein pair add scene page."""
    self.open_help_dialog("protein_pair_add_scene")

  def open_protein_pair_update_scene_page(self) -> None:
    """Opens the help dialog on the protein pair update scene page."""
    self.open_help_dialog("protein_pair_update_scene")

  def open_protein_pair_delete_scene_page(self) -> None:
    """Opens the help dialog on the protein pair delete scene page."""
    self.open_help_dialog("protein_pair_delete_scene")

  # </editor-fold>

  # <editor-fold desc="Project help pages">
  def open_create_project_page(self) -> None:
    """Opens the help dialog on the create new project page."""
    self.open_help_dialog("project_create")

  def open_delete_project_page(self) -> None:
    """Opens the help dialog on the delete project page."""
    self.open_help_dialog("project_delete")

  def open_open_project_page(self) -> None:
    """Opens the help dialog on the open project page."""
    self.open_help_dialog("project_open")

  def open_use_project_page(self) -> None:
    """Opens the help dialog on the use project page."""
    self.open_help_dialog("project_use")
  # </editor-fold>

  def open_advanced_prediction_configuration_page(self) -> None:
    """Opens the help dialog on the advanced prediction configuration page."""
    self.open_help_dialog("advanced_prediction_config")

  def open_colabfold_page(self) -> None:
    """Opens the help dialog on the colabfold page."""
    self.open_help_dialog("colabfold")

  def open_distance_analysis_page(self) -> None:
    """Opens the help dialog on the distance analysis page."""
    self.open_help_dialog("distance_analysis")

  def open_results_summary_page(self) -> None:
    """Opens the help dialog on the results summary page."""
    self.open_help_dialog("results_summary")

  def open_distance_data_visualizer_page(self) -> None:
    """Opens the help dialog on the distance data visualizer page."""
    self.open_help_dialog("distance_data_visualizer")

  def open_pyssa_settings_page(self) -> None:
    """Opens the help dialog on the pyssa settings page."""
    self.open_help_dialog("pyssa_settings")
