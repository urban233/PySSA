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
"""Module for the job manager class."""
import logging
import queue
import subprocess
import time
from typing import Union

from PyQt5 import QtCore

from pyssa.controller import interface_manager
from pyssa.gui.ui.custom_widgets import job_entry
from pyssa.internal.data_structures import project, structure_prediction, job
from pyssa.internal.data_structures.data_classes import job_summary
from pyssa.internal.thread import tasks
from pyssa.logging_pyssa import log_handlers
from pyssa.util import enums, exception, exit_codes

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class JobManager:

    def __init__(self):
        self._prediction_queue: "queue.Queue" = queue.Queue()
        self._prediction_queue_lock = QtCore.QMutex()
        self._is_prediction_queue_running = False
        self.current_prediction_job: "job_summary.PredictionJobSummary" = None
        self._distance_analysis_queue: "queue.Queue" = queue.Queue()
        self._is_distance_analysis_queue_running = False
        self.current_distance_analysis_job: "job_summary.DistanceAnalysisJobSummary" = None
        self._is_prediction_and_distance_analysis_queue_running = False
        self._prediction_and_distance_analysis_queue = queue.Queue()
        self._ray_tracing_queue: "queue.Queue" = queue.Queue()
        self._is_ray_tracing_queue_running = False

    # general approach
    def put_job_into_queue(self, a_job: Union["job.PredictionJob", "job.DistanceAnalysisJob", "job.RayTracingJob"]):
        if a_job.type == enums.JobType.PREDICTION:
            if self._is_prediction_queue_running:
                logger.debug("Prediction queue is already running.")
                self._prediction_queue.put(a_job)
            else:
                logger.debug("Prediction queue needs to be started.")
                self._prediction_queue.put(a_job)
                self._prediction_queue_thread = tasks.Task(
                    target=self._execute_prediction_job_queue,
                    args=(0, 0),
                    post_func=self._prediction_queue_finished,
                )
                self._prediction_queue_thread.start()
        elif a_job.type == enums.JobType.DISTANCE_ANALYSIS:
            if self._is_distance_analysis_queue_running:
                logger.debug("Distance analysis queue is already running.")
                self._distance_analysis_queue.put(a_job)
            else:
                logger.debug("Distance analysis queue needs to be started.")
                self._distance_analysis_queue.put(a_job)
                self._distance_analysis_queue_thread = tasks.Task(
                    target=self._execute_distance_analysis_job_queue,
                    args=(0, 0),
                    post_func=self._distance_analysis_queue_finished,
                )
                self._distance_analysis_queue_thread.start()
        elif a_job.type == enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS:
            # if self._is_prediction_queue_running:
            #     logger.debug("Prediction queue is already running.")
            #     self._prediction_queue.put(a_job)
            # else:
            #     logger.debug("Prediction queue needs to be started.")
            #     self._prediction_queue.put(a_job)
            #     self._prediction_queue_thread = tasks.Task(
            #         target=self._execute_prediction_job_queue,
            #         args=(0, 0),
            #         post_func=self._prediction_queue_finished,
            #     )
            #     self._prediction_queue_thread.start()

            if self._is_prediction_and_distance_analysis_queue_running:
                logger.debug("Distance analysis queue is already running.")
                self._prediction_and_distance_analysis_queue.put(a_job)
            else:
                logger.debug("Distance analysis queue needs to be started.")
                self._prediction_and_distance_analysis_queue.put(a_job)
                self._prediction_and_distance_analysis_queue_thread = tasks.Task(
                    target=self._execute_prediction_and_distance_analysis_job_queue,
                    args=(0, 0),
                    post_func=self._prediction_and_distance_analysis_queue_finished,
                )
                self._prediction_and_distance_analysis_queue_thread.start()
        elif a_job.type == enums.JobType.RAY_TRACING:
            if self._is_ray_tracing_queue_running:
                logger.debug("Ray-tracing queue is already running.")
                self._ray_tracing_queue.put(a_job)
            else:
                logger.debug("Ray-tracing queue needs to be started.")
                self._ray_tracing_queue.put(a_job)
                self._ray_tracing_queue_thread = tasks.Task(
                    target=self._execute_ray_tracing_job_queue,
                    args=(0, 0),
                    post_func=self._ray_tracing_queue_finished,
                )
                self._ray_tracing_queue_thread.start()
        else:
            pass

    def pop_job_from_queue(self, a_job: Union["job.PredictionJob", "job.DistanceAnalysisJob", "job.RayTracingJob"]):
        if a_job.type == enums.JobType.PREDICTION:
            self._prediction_queue_lock.lock()
            tmp_queue_elements = []
            while not self._prediction_queue.empty():
                tmp_queue_elements.append(self._prediction_queue.get())
            if len(tmp_queue_elements) == 0:
                self._prediction_queue_lock.unlock()
                return
            tmp_queue_elements.pop(tmp_queue_elements.index(a_job))
            for element in tmp_queue_elements:
                self._prediction_queue.put(element)
            self._prediction_queue_lock.unlock()

    def is_prediction_job_currently_running(self, a_job: "job.PredictionJob") -> bool:
        if a_job.type == enums.JobType.PREDICTION:
            self._prediction_queue_lock.lock()
            tmp_queue_elements = []
            while not self._prediction_queue.empty():
                tmp_queue_elements.append(self._prediction_queue.get())
            for element in tmp_queue_elements:
                self._prediction_queue.put(element)
            self._prediction_queue_lock.unlock()
            if a_job in tmp_queue_elements:
                return False
            return True

    def stop_prediction_queue(self):
        # Stop the active queue
        logger.info("Stopping prediction queue.")
        self._prediction_queue_thread = None
        self._is_prediction_queue_running = False
        # Restart queue if it is not empty
        if not self._prediction_queue.empty():
            logger.info("Restating prediction queue.")
            self._prediction_queue_thread = tasks.Task(
                target=self._execute_prediction_job_queue,
                args=(0, 0),
                post_func=self._prediction_queue_finished,
            )
            self._prediction_queue_thread.start()
        else:
            logger.info("Prediction queue is empty.")

    def get_queue(self, a_type: enums.JobType):
        if a_type == enums.JobType.PREDICTION:
            return self._prediction_queue
        elif a_type == enums.JobType.DISTANCE_ANALYSIS:
            return self._distance_analysis_queue
        elif a_type == enums.JobType.RAY_TRACING:
            return self._ray_tracing_queue
        else:
            return None

    # <editor-fold desc="Prediction job">
    def create_prediction_job(
            self,
            a_project: "project.Project",
            the_prediction_protein_infos,
            the_prediction_configuration,
            the_project_lock: QtCore.QMutex,
            the_interface_manager: "interface_manager.InterfaceManager",
    ):
        tmp_prediction_job = job.PredictionJob(
            a_project,
            the_prediction_protein_infos,
            the_prediction_configuration,
            the_project_lock
        )
        tmp_protein_names = []
        for tmp_protein_info in the_prediction_protein_infos:
            tmp_protein_names.append(tmp_protein_info.name)
        tmp_prediction_job.update_status_bar_signal.connect(the_interface_manager.status_bar_manager.update_job_entry)
        tmp_prediction_job.job_entry_widget = job_entry.JobEntry(
            "Running ColabFold prediction",
            a_project.get_project_name(),
            tmp_protein_names,
            []
        )
        tmp_prediction_job.job_entry_widget.btn_cancel_job.clicked.connect(tmp_prediction_job.cancel_job)
        tmp_prediction_job.cancel_job_signal.connect(the_interface_manager.cancel_job)
        return tmp_prediction_job, tmp_prediction_job.job_entry_widget

    def put_prediction_job_into_queue(self, a_prediction_job: "job.PredictionJob"):
        if self._is_prediction_queue_running:
            logger.debug("Prediction queue is already running.")
            self._prediction_queue.put(a_prediction_job)
        else:
            logger.debug("Prediction queue needs to be started.")
            self._prediction_queue.put(a_prediction_job)
            self._prediction_queue_thread = tasks.Task(
                target=self._execute_prediction_job_queue,
                args=(0, 0),
                post_func=self._prediction_queue_finished,
            )
            self._prediction_queue_thread.start()

    def _execute_prediction_job_queue(self, placeholder_1, placeholder_2):
        logger.debug("Starting prediction queue ...")
        while True:
            self._is_prediction_queue_running = True
            self._prediction_queue_lock.lock()
            tmp_prediction_job: "job.PredictionJob" = self._prediction_queue.get()
            self._prediction_queue_lock.unlock()
            if tmp_prediction_job is None:
                self._is_prediction_queue_running = False
                break
            self.current_prediction_job = job_summary.PredictionJobSummary(tmp_prediction_job.prediction_protein_infos)
            tmp_prediction_job.run_job()
            if self._prediction_queue.empty():
                logger.info("The prediction queue is empty and will now end execution.")
                break
            self._prediction_queue.task_done()
        self._is_prediction_queue_running = False
        self.current_prediction_job = None
        logger.info("Prediction queue is finished.")
        return "Finished.", 0

    def _prediction_queue_finished(self):
        logger.info("Prediction queue is empty and thread is no longer running.")

    def stop_prediction_queue_execution(self):
        """Stops the queue after gracefully after all items until the None are executed."""
        logger.debug("Stopping prediction queue ...")
        self._prediction_queue.put(None)
    # </editor-fold>

    # <editor-fold desc="Distance analysis job">
    def create_distance_analysis_job(
            self,
            a_project: "project.Project",
            the_project_lock: QtCore.QMutex,
            the_interface_manager: "interface_manager.InterfaceManager",
            a_list_with_analysis_names,
            a_cutoff: float,
            cycles: int
    ):
        tmp_distance_analysis_job = job.DistanceAnalysisJob(
            a_project,
            the_project_lock,
            a_list_with_analysis_names,
            a_cutoff,
            cycles
        )
        tmp_distance_analysis_job.update_status_bar_signal.connect(
            the_interface_manager.status_bar_manager.update_job_entry)
        tmp_distance_analysis_job.job_entry_widget = job_entry.JobEntry(
            "Running distance analysis", a_project.get_project_name()
        )
        return tmp_distance_analysis_job, tmp_distance_analysis_job.job_entry_widget

    def put_distance_analysis_job_into_queue(self, a_distance_analysis_job: "job.PredictionJob"):
        if self._is_distance_analysis_queue_running:
            logger.debug("Distance analysis queue is already running.")
            self._distance_analysis_queue.put(a_distance_analysis_job)
        else:
            logger.debug("Distance analysis queue needs to be started.")
            self._distance_analysis_queue.put(a_distance_analysis_job)
            self._distance_analysis_queue_thread = tasks.Task(
                target=self._execute_distance_analysis_job_queue,
                args=(0, 0),
                post_func=self._distance_analysis_queue_finished,
            )
            self._distance_analysis_queue_thread.start()

    def _execute_distance_analysis_job_queue(self, placeholder_1, placeholder_2):
        logger.debug("Starting distance analysis queue ...")
        while True:
            self._is_distance_analysis_queue_running = True
            tmp_distance_analysis_job: "job.DistanceAnalysisJob" = self._distance_analysis_queue.get()
            if tmp_distance_analysis_job is None:
                self._is_distance_analysis_queue_running = False
                break
            self.current_distance_analysis_job = job_summary.DistanceAnalysisJobSummary(
                tmp_distance_analysis_job.list_with_analysis_names
            )
            tmp_distance_analysis_job.run_job()
            if self._distance_analysis_queue.empty():
                logger.info("The distance analysis queue is empty and will now end execution.")
                break
            self._distance_analysis_queue.task_done()
        self._is_distance_analysis_queue_running = False
        self.current_distance_analysis_job = None
        logger.info("Distance analysis queue is finished.")
        return "Finished.", 0

    def _distance_analysis_queue_finished(self):
        logger.info("Distance analysis queue is empty and thread is no longer running.")

    def stop_distance_analysis_queue_execution(self):
        """Stops the queue after gracefully after all items until the None are executed."""
        logger.debug("Stopping distance analysis queue ...")
        self._distance_analysis_queue.put(None)
    # </editor-fold>

    # <editor-fold desc="Distance analysis job">
    def create_prediction_and_distance_analysis_job(
            self,
            a_prediction_job: "job.PredictionJob",
            a_distance_analysis_job: "job.DistanceAnalysisJob",
            the_interface_manager
    ):
        tmp_prediction_and_distance_analysis_job = job.PredictionAndDistanceAnalysisJob(
            a_prediction_job,
            a_distance_analysis_job
        )
        tmp_prediction_and_distance_analysis_job.update_status_bar_signal.connect(
            the_interface_manager.status_bar_manager.update_job_entry)
        tmp_prediction_and_distance_analysis_job.job_entry_widget = job_entry.JobEntry(
            "Running ColabFold prediction + distance analysis", a_prediction_job.frozen_project.get_project_name()
        )
        return tmp_prediction_and_distance_analysis_job, tmp_prediction_and_distance_analysis_job.job_entry_widget

    def _execute_prediction_and_distance_analysis_job_queue(self, placeholder_1, placeholder_2):
        logger.debug("Starting prediction and distance analysis queue ...")

        while self._is_prediction_queue_running:
            time.sleep(60)

        while True:
            self._is_prediction_and_distance_analysis_queue_running = True
            tmp_prediction_and_distance_analysis_job: "job.DistanceAnalysisJob" = self._prediction_and_distance_analysis_queue.get()
            if tmp_prediction_and_distance_analysis_job is None:
                self._is_prediction_and_distance_analysis_queue_running = False
                break
            # self.current_prediction_and_distance_analysis_job = job_summary.DistanceAnalysisJobSummary(
            #     tmp_prediction_and_distance_analysis_job.list_with_analysis_names
            # )
            tmp_prediction_and_distance_analysis_job.run_job()
            if self._prediction_and_distance_analysis_queue.empty():
                logger.info("The prediction and distance analysis queue is empty and will now end execution.")
                break
            self._prediction_and_distance_analysis_queue.task_done()
        self._is_prediction_and_distance_analysis_queue_running = False
        self.current_prediction_and_distance_analysis_job = None
        logger.info("Prediction and distance analysis queue is finished.")
        return "Finished.", 0

    def _prediction_and_distance_analysis_queue_finished(self):
        logger.info("Prediction and distance analysis queue is empty and thread is no longer running.")

    def stop_prediction_and_distance_analysis_queue_execution(self):
        """Stops the queue after gracefully after all items until the None are executed."""
        logger.debug("Stopping prediction and distance analysis queue ...")
        self._prediction_and_distance_analysis_queue.put(None)
    # </editor-fold>

    # <editor-fold desc="Ray-tracing job">
    def create_ray_tracing_job(
            self,
            the_destination_image_filepath,
            the_cached_session_filepath,
            image_ray_trace_mode,
            image_ray_texture,
            image_renderer,
            the_interface_manager,
            a_project_name: str
    ):
        tmp_ray_tracing_job = job.RayTracingJob(
            the_destination_image_filepath,
            the_cached_session_filepath,
            image_ray_trace_mode,
            image_ray_texture,
            image_renderer
        )
        tmp_ray_tracing_job.update_status_bar_signal.connect(the_interface_manager.status_bar_manager.update_job_entry)
        tmp_ray_tracing_job.job_entry_widget = job_entry.JobEntry("Creating ray-traced image", a_project_name)
        return tmp_ray_tracing_job, tmp_ray_tracing_job.job_entry_widget

    def put_ray_tracing_job_into_queue(self, a_ray_tracing_job: "job.RayTracingJob"):
        if self._is_ray_tracing_queue_running:
            logger.debug("Ray-tracing queue is already running.")
            self._ray_tracing_queue.put(a_ray_tracing_job)
        else:
            logger.debug("Ray-tracing queue needs to be started.")
            self._ray_tracing_queue.put(a_ray_tracing_job)
            self._ray_tracing_queue_thread = tasks.Task(
                target=self._execute_ray_tracing_job_queue,
                args=(0, 0),
                post_func=self._ray_tracing_queue_finished,
            )
            self._ray_tracing_queue_thread.start()

    def _execute_ray_tracing_job_queue(self, placeholder_1, placeholder_2):
        logger.debug("Starting ray-tracing queue ...")
        while True:
            self._is_ray_tracing_queue_running = True
            tmp_ray_tracing_job: "job.RayTracingJob" = self._ray_tracing_queue.get()
            if tmp_ray_tracing_job is None:
                self._is_ray_tracing_queue_running = False
                break
            # Do something with the task
            tmp_ray_tracing_job.run_job()
            if self._ray_tracing_queue.empty():
                logger.info("The ray-tracing queue is empty and will now end execution.")
                break
            self._ray_tracing_queue.task_done()
        logger.info("Ray-tracing queue is finished.")
        self._is_ray_tracing_queue_running = False
        return "Finished.", 0

    def _ray_tracing_queue_finished(self):
        logger.info("Ray-tracing queue is empty and thread is no longer running.")
    # </editor-fold>
