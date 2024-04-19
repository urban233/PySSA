#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2022
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
"""Module for the watcher class."""
import logging
import queue

from pyssa.internal.data_structures import job
from pyssa.internal.data_structures.data_classes import job_summary
from pyssa.logging_pyssa import log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class Watcher:
    """Holds information about the names that cannot be used in the current opened project."""

    def __init__(self):
        self.project = None
        self.protein_names_blacklist: list[str] = []
        self.protein_pair_names_blacklist: list[str] = []

    # <editor-fold desc="Private methods">
    def _retrieve_all_protein_names_from_project(self):
        """Retrieve all protein names from project."""
        tmp_all_protein_names = []
        for tmp_protein in self.project.proteins:
            tmp_all_protein_names.append(tmp_protein.get_molecule_object())
        return tmp_all_protein_names

    def _retrieve_all_protein_pair_names_from_project(self):
        """Retrieve all protein names from project."""
        tmp_all_protein_pair_names = []
        for tmp_protein_pair in self.project.protein_pairs:
            tmp_all_protein_pair_names.append(tmp_protein_pair.name)
        return tmp_all_protein_pair_names

    def _retrieve_all_prediction_names_of_the_current_project(
            self,
            a_prediction_queue: queue.Queue,
            the_current_prediction_job: "job_summary.PredictionJobSummary"
    ):
        # Retrieve all elements from the Queue
        elements = []
        while not a_prediction_queue.empty():
            tmp_job: "job.PredictionJob" = a_prediction_queue.get()
            logger.debug(f"Watcher project name: {self.project.get_project_name()} vs. Job project name: {tmp_job.frozen_project.get_project_name()}")
            if tmp_job.frozen_project.get_project_name() == self.project.get_project_name():
                for tmp_protein_info in tmp_job.prediction_protein_infos:
                    elements.append(tmp_protein_info.name)
        if the_current_prediction_job is not None:
            elements.extend(the_current_prediction_job.get_protein_names())
        return elements

    def _retrieve_all_distance_analysis_names_of_the_current_project(
            self,
            a_distance_analysis_queue: queue.Queue,
            the_current_distance_analysis_job
    ):
        # Retrieve all elements from the Queue
        elements = []
        while not a_distance_analysis_queue.empty():
            tmp_job: "job.DistanceAnalysisJob" = a_distance_analysis_queue.get()
            if tmp_job.frozen_project.get_project_name() == self.project.get_project_name():
                for tmp_analysis_name in tmp_job.list_with_analysis_names:
                    elements.append(tmp_analysis_name)
        if the_current_distance_analysis_job is not None:
            elements.extend(the_current_distance_analysis_job)
        return elements
    # </editor-fold>

    # <editor-fold desc="Public methods">
    def setup_blacklists(self,
                         a_project,
                         a_prediction_queue: queue.Queue,
                         a_distance_analysis_queue: queue.Queue,
                         the_current_prediction_job,
                         the_current_distance_analysis_job
                         ):
        # Reset lists
        self.protein_names_blacklist: list[str] = []
        self.protein_pair_names_blacklist: list[str] = []
        # Setup process
        self.project = a_project
        self.protein_names_blacklist: list[str] = self._retrieve_all_protein_names_from_project()
        self.protein_pair_names_blacklist: list[str] = self._retrieve_all_protein_pair_names_from_project()

        # Extend blacklist with job information
        self.protein_names_blacklist.extend(
            self._retrieve_all_prediction_names_of_the_current_project(
                a_prediction_queue, the_current_prediction_job
            )
        )
        self.protein_pair_names_blacklist.extend(
            self._retrieve_all_distance_analysis_names_of_the_current_project(
                a_distance_analysis_queue, the_current_distance_analysis_job
            )
        )
        logger.info(f"After setting up the blacklist, the protein blacklist is      {self.protein_names_blacklist}.")
        logger.info(f"After setting up the blacklist, the protein pair blacklist is {self.protein_pair_names_blacklist}.")
    # </editor-fold>

    # <editor-fold desc="Public methods">
    def add_proteins_from_new_job(self, prediction_protein_infos):
        logger.info(f"Adding the protein names from the new job {prediction_protein_infos} to the project's blacklist.")
        for tmp_protein_info in prediction_protein_infos:
            self.protein_names_blacklist.append(tmp_protein_info.name)
        logger.info(f"After adding the new names, the protein blacklist is {self.protein_names_blacklist}.")

    def add_protein(self, a_protein_name):
        logger.info(f"Adding the protein name {a_protein_name} to the project's blacklist.")
        self.protein_names_blacklist.append(a_protein_name)
        logger.info(f"After adding the name, the protein blacklist is {self.protein_names_blacklist}.")

    def remove_protein(self, a_protein_name):
        logger.info(f"Removing the protein name {a_protein_name} from the project's blacklist.")
        self.protein_names_blacklist.pop(self.protein_names_blacklist.index(a_protein_name))
        logger.info(f"After removing the name, the protein blacklist is {self.protein_names_blacklist}.")

    def add_protein_pairs_from_new_job(self, analysis_run_names):
        logger.info(f"Adding the protein pair names from the new job {analysis_run_names} to the project's blacklist.")
        for tmp_analysis_run_name in analysis_run_names:
            tmp_string = tmp_analysis_run_name.replace(";", "_")
            tmp_string_2 = tmp_string.replace(",", "_")
            self.protein_pair_names_blacklist.append(tmp_string_2)
        logger.info(f"After adding the new names, the protein pair blacklist is {self.protein_pair_names_blacklist}.")

    def add_protein_pair(self, a_protein_pair_name):
        logger.info(f"Adding the protein pair name {a_protein_pair_name} to the project's blacklist.")
        self.protein_pair_names_blacklist.append(a_protein_pair_name)
        logger.info(f"After adding the name, the protein pair blacklist is {self.protein_pair_names_blacklist}.")

    def remove_protein_pair(self, a_protein_pair_name):
        logger.info(f"Removing the protein pair name {a_protein_pair_name} from the project's blacklist.")
        self.protein_pair_names_blacklist.pop(self.protein_pair_names_blacklist.index(a_protein_pair_name))
        logger.info(f"After removing the name, the protein pair blacklist is {self.protein_pair_names_blacklist}.")

    def is_protein_name_on_blacklist(self, a_protein_name) -> bool:
        if a_protein_name in self.protein_names_blacklist:
            return True
        return False

    # </editor-fold>
