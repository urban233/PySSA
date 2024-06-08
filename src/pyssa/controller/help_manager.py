#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/zielesny/PySSA>
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
import subprocess

import pygetwindow
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtCore import QUrl

from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.gui.ui.views import help_view
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import constants

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class HelpManager:
  """Class for the HelpManager."""

  # <editor-fold desc="Class attributes">
  project_root_path = pathlib.Path(__file__).parent.absolute().parent.parent.parent
  """The Path to the project root directory."""

  docs_html_path = pathlib.Path(project_root_path, "docs", "build", "html")
  """The Path to the docs html directory."""

  docs_help_html_path = pathlib.Path(docs_html_path, "help")
  """The Path to the docs help html directory."""

  # </editor-fold>

  def __init__(self) -> None:
    """Constructor."""
    super().__init__()
    self._view = help_view.BasicBrowserView()

  def change_url(self, an_url: str) -> None:
    """Changes the URL of the browser window."""
    if an_url.find("\\"):
      tmp_url = an_url.replace("\\", "/")
    else:
      tmp_url = an_url

    self._view.browser.setUrl(QUrl(tmp_url))
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
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "help.html")))

  def open_sequences_tab_page(self) -> None:
    """Opens the help dialog on the sequences tab page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "sequences", "sequences_tab.html")))

  def open_additional_sequence_information_page(self) -> None:
    """Opens the help dialog on the additional sequence information page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "sequences", "additional_sequence_information.html")))

  def open_sequence_import_page(self) -> None:
    """Opens the help dialog on the import sequence page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "sequences", "sequence_import.html")))

  def open_sequence_add_page(self) -> None:
    """Opens the help dialog on the add sequence page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "sequences", "sequence_add.html")))

  def open_sequence_save_page(self) -> None:
    """Opens the help dialog on the save sequence page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "sequences", "sequence_save.html")))

  def open_sequence_delete_page(self) -> None:
    """Opens the help dialog on the delete sequence page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "sequences", "sequence_delete.html")))

  def open_proteins_tab_page(self) -> None:
    """Opens the help dialog on the proteins tab page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "proteins", "proteins_tab.html")))

  def open_protein_import_page(self) -> None:
    """Opens the help dialog on the import protein page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "proteins", "protein_import.html")))

  def open_protein_save_page(self) -> None:
    """Opens the help dialog on the save protein page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "proteins", "protein_save.html")))

  def open_protein_delete_page(self) -> None:
    """Opens the help dialog on the delete protein page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "proteins", "protein_delete.html")))

  def open_protein_pymol_scene_configuration_page(self) -> None:
    """Opens the help dialog on the protein pymol scene config page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "proteins", "protein_pymol_scene_configuration.html")))

  def open_protein_load_session_page(self) -> None:
    """Opens the help dialog on the open protein pymol session page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "proteins", "protein_load_session.html")))

  def open_protein_add_scene_page(self) -> None:
    """Opens the help dialog on the protein add scene page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "proteins", "protein_add_scene.html")))

  def open_protein_update_scene_page(self) -> None:
    """Opens the help dialog on the protein update scene page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "proteins", "protein_update_scene.html")))

  def open_protein_delete_scene_page(self) -> None:
    """Opens the help dialog on the protein delete scene page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "proteins", "protein_delete_scene.html")))

  def open_protein_pairs_tab_page(self) -> None:
    """Opens the help dialog on the protein pairs tab page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "protein_pairs", "protein_pairs_tab.html")))

  def open_protein_pair_delete_page(self) -> None:
    """Opens the help dialog on the delete protein pair page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "protein_pairs", "protein_pair_delete.html")))

  def open_protein_pair_pymol_scene_configuration_page(self) -> None:
    """Opens the help dialog on the protein pair pymol scene config page."""
    self.change_url(
      str(pathlib.Path(self.docs_help_html_path, "protein_pairs", "protein_pair_pymol_scene_configuration.html")))

  def open_protein_pair_load_session_page(self) -> None:
    """Opens the help dialog on the open protein pair pymol session page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "protein_pairs", "protein_pair_load_session.html")))

  def open_protein_pair_add_scene_page(self) -> None:
    """Opens the help dialog on the protein pair add scene page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "protein_pairs", "protein_pair_add_scene.html")))

  def open_protein_pair_update_scene_page(self) -> None:
    """Opens the help dialog on the protein pair update scene page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "protein_pairs", "protein_pair_update_scene.html")))

  def open_protein_pair_delete_scene_page(self) -> None:
    """Opens the help dialog on the protein pair delete scene page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "protein_pairs", "protein_pair_delete_scene.html")))

  # </editor-fold>

  # <editor-fold desc="Project help pages">
  def open_create_project_page(self) -> None:
    """Opens the help dialog on the create new project page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "project", "new_project.html")))

  def open_delete_project_page(self) -> None:
    """Opens the help dialog on the delete project page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "project", "delete_project.html")))

  def open_open_project_page(self) -> None:
    """Opens the help dialog on the open project page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "project", "open_project.html")))

  def open_use_project_page(self) -> None:
    """Opens the help dialog on the use project page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "project", "use_project.html")))
  # </editor-fold>

  def open_advanced_prediction_configuration_page(self) -> None:
    """Opens the help dialog on the advanced prediction configuration page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "protein_structure_prediction", "advanced_prediction_configuration.html")))

  def open_colabfold_page(self) -> None:
    """Opens the help dialog on the colabfold page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "protein_structure_prediction", "colabfold.html")))

  def open_distance_analysis_page(self) -> None:
    """Opens the help dialog on the distance analysis page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "protein_structure_analysis", "distance_analysis.html")))

  def open_results_summary_page(self) -> None:
    """Opens the help dialog on the results summary page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "results", "summary.html")))

  def open_distance_data_visualizer_page(self) -> None:
    """Opens the help dialog on the distance data visualizer page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "results", "distance_data_visualizer.html")))

  def open_pyssa_settings_page(self) -> None:
    """Opens the help dialog on the pyssa settings page."""
    self.change_url(str(pathlib.Path(self.docs_help_html_path, "settings", "pyssa_settings.html")))
