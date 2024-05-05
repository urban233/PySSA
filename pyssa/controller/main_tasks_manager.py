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
"""Module for the main tasks manager."""
# TODO: check if this class is necessary
from pyssa.internal.thread import tasks
from pyssa.presenter import main_presenter_async
from pyssa.util import constants


class MainTasksManager:
    """A manager class which manages the threads for the main tasks such as prediction and analysis."""

    prediction_task: "tasks.Task"
    distance_analysis_task: "tasks.Task"

    def __init__(self):
        self.prediction_task = None
        self.distance_analysis_task = None

    def start_prediction_task(self, the_prediction_task: "tasks.Task"):
        """Runs a structure prediction"""
        constants.PYSSA_LOGGER.info("Running only a prediction.")
        # No analysis after prediction
        self.prediction_task = the_prediction_task
        self.prediction_task.start()

    def check_if_prediction_task_is_finished(self) -> bool:
        if self.prediction_task is None:
            raise ValueError("Cannot check the state of the prediction task if the task object is None! Please check first if the task object is None!")
        return self.prediction_task.is_finished()

    def start_distance_analysis_task(self, the_distance_analysis_task: "tasks.Task"):
        """Runs a structure distance_analysis"""
        constants.PYSSA_LOGGER.info("Running only a distance_analysis.")
        # No analysis after distance_analysis
        self.distance_analysis_task = the_distance_analysis_task
        self.distance_analysis_task.start()

    def check_if_distance_analysis_task_is_finished(self) -> bool:
        if self.distance_analysis_task is None:
            raise ValueError("Cannot check the state of the distance_analysis task if the task object is None! Please check first if the task object is None!")
        return self.distance_analysis_task.is_finished()
