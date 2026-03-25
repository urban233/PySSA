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
"""Module for the help dialog registry.

This module provides a centralized registry for managing help dialogs in PySSA.
It allows easy registration and retrieval of help pages organized by categories.
"""
import logging
import pathlib
from typing import Optional
from dataclasses import dataclass

from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import constants

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


@dataclass
class HelpDialog:
    """Data class representing a help dialog configuration.

    Attributes:
        id: Unique identifier for the help dialog
        title: Display title for the help dialog
        category: Category grouping (e.g., "Project", "Sequences", "Proteins")
        html_path: Relative path to the HTML file from docs/help directory
        description: Optional short description of what the help page covers
    """
    id: str
    title: str
    category: str
    html_path: str
    description: str = ""


class HelpDialogRegistry:
    """Registry for managing help dialogs in PySSA.

    This class provides a centralized way to register and access help dialogs.
    Help dialogs are organized by categories and can be easily added or retrieved.
    """

    _instance: Optional["HelpDialogRegistry"] = None

    def __new__(cls):
        """Singleton pattern to ensure only one registry exists."""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance._initialized = False
        return cls._instance

    def __init__(self):
        """Initialize the help dialog registry."""
        if self._initialized:
            return

        self._dialogs: dict[str, HelpDialog] = {}
        self.docs_html_path = constants.DOCS_PATH
        self.docs_help_html_path = pathlib.Path(self.docs_html_path, "help")
        self.html_suffix = ".html"

        self._register_default_dialogs()
        self._initialized = True

    def register(self, dialog: HelpDialog) -> None:
        """Register a new help dialog.

        Args:
            dialog: HelpDialog instance to register
        """
        if dialog.id in self._dialogs:
            logger.warning(f"Help dialog with ID '{dialog.id}' already exists. Overwriting.")

        self._dialogs[dialog.id] = dialog
        logger.debug(f"Registered help dialog: {dialog.id}")

    def get_html_path(self, dialog_id: str) -> Optional[pathlib.Path]:
        """Get the full HTML file path for a help dialog.

        Args:
            dialog_id: Unique identifier of the help dialog

        Returns:
            Full path to the HTML file if dialog exists, None otherwise
        """
        dialog = self.get_dialog(dialog_id)
        if dialog:
            return pathlib.Path(self.docs_help_html_path, dialog.html_path)
        return None

    def _register_default_dialogs(self) -> None:
        """Register all default help dialogs for PySSA.

        This method registers all the built-in help pages. New help pages
        should be added here following the same pattern.
        """
        # General
        self.register(HelpDialog(
            id="general",
            title="General Help",
            category="General",
            html_path="index.html",
            description="General PySSA help and overview"
        ))

        # Project
        self.register(HelpDialog(
            id="project_create",
            title="Create Project",
            category="Project",
            html_path="project/new_project.html",
            description="How to create a new project"
        ))

        self.register(HelpDialog(
            id="project_delete",
            title="Delete Project",
            category="Project",
            html_path="project/delete_project.html",
            description="How to delete a project"
        ))

        self.register(HelpDialog(
            id="project_open",
            title="Open Project",
            category="Project",
            html_path="project/open_project.html",
            description="How to open an existing project"
        ))

        self.register(HelpDialog(
            id="project_use",
            title="Use Project",
            category="Project",
            html_path="project/use_project.html",
            description="How to work with a project"
        ))

        # Sequences
        self.register(HelpDialog(
            id="sequences_tab",
            title="Sequences Tab",
            category="Sequences",
            html_path="sequences/sequences_tab.html",
            description="Overview of the sequences tab"
        ))

        self.register(HelpDialog(
            id="sequence_additional_info",
            title="Additional Sequence Information",
            category="Sequences",
            html_path="sequences/additional_sequence_information.html",
            description="Viewing additional sequence information"
        ))

        self.register(HelpDialog(
            id="sequence_import",
            title="Import Sequence",
            category="Sequences",
            html_path="sequences/sequence_import.html",
            description="How to import sequences"
        ))

        self.register(HelpDialog(
            id="sequence_add",
            title="Add Sequence",
            category="Sequences",
            html_path="sequences/sequence_add.html",
            description="How to add a new sequence"
        ))

        self.register(HelpDialog(
            id="sequence_save",
            title="Save Sequence",
            category="Sequences",
            html_path="sequences/sequence_save.html",
            description="How to save a sequence"
        ))

        self.register(HelpDialog(
            id="sequence_delete",
            title="Delete Sequence",
            category="Sequences",
            html_path="sequences/sequence_delete.html",
            description="How to delete a sequence"
        ))

        # Proteins
        self.register(HelpDialog(
            id="proteins_tab",
            title="Proteins Tab",
            category="Proteins",
            html_path="proteins/proteins_tab.html",
            description="Overview of the proteins tab"
        ))

        self.register(HelpDialog(
            id="protein_import",
            title="Import Protein",
            category="Proteins",
            html_path="proteins/protein_import.html",
            description="How to import proteins"
        ))

        self.register(HelpDialog(
            id="protein_save",
            title="Save Protein",
            category="Proteins",
            html_path="proteins/protein_save.html",
            description="How to save a protein"
        ))

        self.register(HelpDialog(
            id="protein_delete",
            title="Delete Protein",
            category="Proteins",
            html_path="proteins/protein_delete.html",
            description="How to delete a protein"
        ))

        self.register(HelpDialog(
            id="protein_pymol_scene_config",
            title="PyMOL Scene Configuration",
            category="Proteins",
            html_path="proteins/protein_pymol_scene_configuration.html",
            description="Configure PyMOL scenes for proteins"
        ))

        self.register(HelpDialog(
            id="protein_load_session",
            title="Load PyMOL Session",
            category="Proteins",
            html_path="proteins/protein_load_session.html",
            description="Load a protein PyMOL session"
        ))

        self.register(HelpDialog(
            id="protein_add_scene",
            title="Add Scene",
            category="Proteins",
            html_path="proteins/protein_add_scene.html",
            description="Add a new PyMOL scene for a protein"
        ))

        self.register(HelpDialog(
            id="protein_update_scene",
            title="Update Scene",
            category="Proteins",
            html_path="proteins/protein_update_scene.html",
            description="Update an existing PyMOL scene"
        ))

        self.register(HelpDialog(
            id="protein_delete_scene",
            title="Delete Scene",
            category="Proteins",
            html_path="proteins/protein_delete_scene.html",
            description="Delete a PyMOL scene"
        ))

        # Protein Pairs
        self.register(HelpDialog(
            id="protein_pairs_tab",
            title="Protein Pairs Tab",
            category="Protein Pairs",
            html_path="protein_pairs/protein_pairs_tab.html",
            description="Overview of the protein pairs tab"
        ))

        self.register(HelpDialog(
            id="protein_pair_delete",
            title="Delete Protein Pair",
            category="Protein Pairs",
            html_path="protein_pairs/protein_pair_delete.html",
            description="How to delete a protein pair"
        ))

        self.register(HelpDialog(
            id="protein_pair_pymol_scene_config",
            title="Protein Pair PyMOL Scene Configuration",
            category="Protein Pairs",
            html_path="protein_pairs/protein_pair_pymol_scene_configuration.html",
            description="Configure PyMOL scenes for protein pairs"
        ))

        self.register(HelpDialog(
            id="protein_pair_load_session",
            title="Load Protein Pair Session",
            category="Protein Pairs",
            html_path="protein_pairs/protein_pair_load_session.html",
            description="Load a protein pair PyMOL session"
        ))

        self.register(HelpDialog(
            id="protein_pair_add_scene",
            title="Add Protein Pair Scene",
            category="Protein Pairs",
            html_path="protein_pairs/protein_pair_add_scene.html",
            description="Add a scene for a protein pair"
        ))

        self.register(HelpDialog(
            id="protein_pair_update_scene",
            title="Update Protein Pair Scene",
            category="Protein Pairs",
            html_path="protein_pairs/protein_pair_update_scene.html",
            description="Update a protein pair scene"
        ))

        self.register(HelpDialog(
            id="protein_pair_delete_scene",
            title="Delete Protein Pair Scene",
            category="Protein Pairs",
            html_path="protein_pairs/protein_pair_delete_scene.html",
            description="Delete a protein pair scene"
        ))

        # Prediction
        self.register(HelpDialog(
            id="advanced_prediction_config",
            title="Advanced Prediction Configuration",
            category="Protein Structure Prediction",
            html_path="protein_structure_prediction/advanced_prediction_configuration.html",
            description="Advanced configuration options for predictions"
        ))

        self.register(HelpDialog(
            id="colabfold",
            title="ColabFold",
            category="Protein Structure Prediction",
            html_path="protein_structure_prediction/colabfold.html",
            description="Using ColabFold for structure prediction"
        ))

        # Analysis
        self.register(HelpDialog(
            id="distance_analysis",
            title="Distance Analysis",
            category="Protein Structure Analysis",
            html_path="protein_structure_analysis/distance_analysis.html",
            description="Performing distance analysis on structures"
        ))

        # Results
        self.register(HelpDialog(
            id="results_summary",
            title="Results Summary",
            category="Results",
            html_path="results/summary.html",
            description="View and interpret results"
        ))

        self.register(HelpDialog(
            id="distance_data_visualizer",
            title="Distance Data Visualizer",
            category="Results",
            html_path="results/distance_data_visualizer.html",
            description="Visualize distance analysis data"
        ))

        # Settings
        self.register(HelpDialog(
            id="pyssa_settings",
            title="PySSA Settings",
            category="Settings",
            html_path="settings/pyssa_settings.html",
            description="Configure PySSA settings"
        ))


# Global registry instance
help_registry = HelpDialogRegistry()
