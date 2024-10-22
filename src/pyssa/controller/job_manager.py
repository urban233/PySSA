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
"""Module for the job manager class."""
import logging
import multiprocessing
import queue
import subprocess
import time
from typing import Union, Optional

from PyQt5 import QtCore

from src.pyssa.gui.ui.custom_widgets import job_entry
from src.pyssa.internal.data_structures import job
from src.pyssa.internal.data_structures.data_classes import job_summary
from src.pyssa.internal.thread import tasks
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import enums, constants, exception
import zmq

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class JobManager:
  """Manages all major jobs that can be run with PySSA."""

  def __init__(self) -> None:
    """Constructor."""
    context = zmq.Context()
    self._main_socket = context.socket(zmq.REQ)
    self._main_socket.connect("tcp://127.0.0.1:8070")

    self._prediction_queue: "queue.Queue" = queue.Queue()
    self._prediction_queue_lock = QtCore.QMutex()
    self._is_prediction_queue_running = False
    self.current_prediction_job: Optional[
        "job_summary.PredictionJobSummary"
    ] = None
    self._prediction_socket = context.socket(zmq.REQ)
    self._prediction_socket.connect("tcp://127.0.0.1:8071")

    self._distance_analysis_queue: "queue.Queue" = queue.Queue()
    self._is_distance_analysis_queue_running = False
    self.current_distance_analysis_job: Optional[
        "job_summary.DistanceAnalysisJobSummary"
    ] = None
    self._distance_analysis_socket = context.socket(zmq.REQ)
    self._distance_analysis_socket.connect(
        "tcp://127.0.0.1:8072"
    )  # Connecting to the server on port 5555

    self._is_prediction_and_distance_analysis_queue_running = False
    self._prediction_and_distance_analysis_queue = queue.Queue()

    self._ray_tracing_queue: "queue.Queue" = queue.Queue()
    self._is_ray_tracing_queue_running = False
    self._ray_tracing_socket = context.socket(zmq.REQ)
    self._ray_tracing_socket.connect(
        "tcp://127.0.0.1:8074"
    )  # Connecting to the server on port 5555

    self._general_purpose_socket = context.socket(zmq.REQ)
    self._general_purpose_socket.connect(
        "tcp://127.0.0.1:8075"
    )  # Connecting to the server on port 5555

    # poller = zmq.Poller()
    # poller.register(self._ray_tracing_socket, zmq.POLLOUT)  # Register the socket for outgoing events
    # # Check if the socket is ready for sending data
    # socks = dict(poller.poll(timeout=1000))  # Timeout is in milliseconds
    # if self._ray_tracing_socket in socks and socks[self._ray_tracing_socket] == zmq.POLLOUT:
    #     print("Socket is ready for sending data")
    # else:
    #     print("Socket is not ready for sending data")

  def stop_auxiliary_pymol(self) -> None:
    """Sends an "Abort" message to the main_socket, receives a response, and then sends a JSON message with job_type "Abort" to the main_socket again."""
    self._main_socket.send_string("Abort")
    response = self._main_socket.recv_string()
    logger.debug(f"Received response: {response}")
    message = {
        "job_type": "Abort",
    }
    self._main_socket.send_json(message)
    response = self._main_socket.recv_string()
    logger.debug(f"Received response: {response}")

  def start_auxiliary_pymol(self) -> None:
    """Starts the auxiliary Pymol instance.

    This method is used to start the auxiliary Pymol instance for advanced rendering.
    If the DEBUGGING constant is set to True, debugging mode will be activated and the method returns.
    """
    if constants.DEBUGGING:
      logger.debug("Debugging activated.")
      return
    process = subprocess.Popen(
        [
            constants.PYTHON_FILEPATH,
            f"{constants.PROGRAM_SRC_PATH}\\auxiliary_pymol\\main.py",
        ],
        creationflags=subprocess.CREATE_NO_WINDOW,
    )
    if process.poll() is None:
      logger.debug("main.py of auxiliary pymol started correctly.")
    else:
      logger.debug("main.py of auxiliary pymol failed to start.")

  def get_general_purpose_socket_pair(self) -> tuple:
    """Gets a tuple containing the main socket and the general purpose socket.

    Returns:
        The tuple contains the main socket and the general purpose socket.
    """
    return self._main_socket, self._general_purpose_socket

  # general approach
  def put_job_into_queue(
      self,
      a_job: Union[
          "job.PredictionJob", "job.DistanceAnalysisJob", "job.RayTracingJob"
      ],
  ) -> None:
    """Puts a job into the corresponding queue based on its type.

    Args:
        a_job (Union["job.PredictionJob", "job.DistanceAnalysisJob", "job.RayTracingJob"]): A job object.

    Raises:
        exception.IllegalArgumentError: If `a_job` is None.
    """
    # <editor-fold desc="Checks">
    if a_job is None:
      logger.error("a_job is None.")
      raise exception.IllegalArgumentError("a_job is None.")

    # </editor-fold>

    if a_job.type == enums.JobType.PREDICTION:
      if self._is_prediction_queue_running:
        logger.debug("Prediction queue is already running.")
        self._prediction_queue.put(a_job)
      else:
        logger.debug("Prediction queue needs to be started.")
        self._prediction_queue.put(a_job)
        self._prediction_queue_thread = tasks.LegacyTask(
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
        self._distance_analysis_queue_thread = tasks.LegacyTask(
            target=self._execute_distance_analysis_job_queue,
            args=(0, 0),
            post_func=self._distance_analysis_queue_finished,
        )
        self._distance_analysis_queue_thread.start()
    elif a_job.type == enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS:
      if self._is_prediction_and_distance_analysis_queue_running:
        logger.debug("Distance analysis queue is already running.")
        self._prediction_and_distance_analysis_queue.put(a_job)
      else:
        logger.debug("Distance analysis queue needs to be started.")
        self._prediction_and_distance_analysis_queue.put(a_job)
        self._prediction_and_distance_analysis_queue_thread = tasks.LegacyTask(
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
        self._ray_tracing_queue_thread = tasks.LegacyTask(
            target=self._execute_ray_tracing_job_queue,
            args=(0, 0),
            post_func=self._ray_tracing_queue_finished,
        )
        self._ray_tracing_queue_thread.start()
    else:
      pass

  def pop_job_from_queue(
      self,
      a_job: Union[
          "job.PredictionJob", "job.DistanceAnalysisJob", "job.RayTracingJob"
      ],
  ) -> None:
    """Removes a job from the queue.

    Args:
        a_job (Union["job.PredictionJob", "job.DistanceAnalysisJob", "job.RayTracingJob"]): The job to be removed from the queue.

    Raises:
        exception.IllegalArgumentError: If `a_job` is None.
    """
    # <editor-fold desc="Checks">
    if a_job is None:
      logger.error("a_job is None.")
      raise exception.IllegalArgumentError("a_job is None.")

    # </editor-fold>

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

  def is_prediction_job_currently_running(
      self, a_job: "job.PredictionJob"
  ) -> bool:
    """Checks if a prediction job is currently running.

    Args:
        a_job (job.PredictionJob): An instance of job.PredictionJob.

    Returns:
        A boolean value indicating whether the given prediction job is currently running.

    Raises:
        exception.IllegalArgumentError: If `a_job` is None.
        ValueError: If `a_job.type` is not 'structure prediction'.
    """
    # <editor-fold desc="Checks">
    if a_job is None:
      logger.error("a_job is None.")
      raise exception.IllegalArgumentError("a_job is None.")
    if a_job.type != enums.JobType.PREDICTION:
      logger.error("a_job.type is not 'structure prediction'.")
      raise ValueError("a_job.type is not 'structure prediction'.")

    # </editor-fold>

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

  def stop_prediction_queue(self) -> None:
    """Stops the prediction queue."""
    # Stop the active queue
    logger.info("Stopping prediction queue.")
    self._is_prediction_queue_running = False
    # Restart queue if it is not empty
    if not self._prediction_queue.empty():
      logger.info("Restating prediction queue.")
      self._prediction_queue_thread = tasks.LegacyTask(
          target=self._execute_prediction_job_queue,
          args=(0, 0),
          post_func=self._prediction_queue_finished,
      )
      self._prediction_queue_thread.start()
    else:
      logger.info("Prediction queue is empty.")

  def get_queue(self, a_type: "enums.JobType") -> Optional[queue.Queue]:
    """Gets the queue based on the given job type.

    Args:
        a_type (enums.JobType): The type of job for which the queue is required.

    Returns:
        The corresponding queue based on the given job type. If no matching queue is found, returns None.

    Raises:
        exception.IllegalArgumentError: If `a_type` is None.
    """
    # <editor-fold desc="Checks">
    if a_type is None:
      logger.error("a_type is None.")
      raise exception.IllegalArgumentError("a_type is None.")

    # </editor-fold>

    if a_type == enums.JobType.PREDICTION:
      return self._prediction_queue
    elif a_type == enums.JobType.DISTANCE_ANALYSIS:
      return self._distance_analysis_queue
    elif a_type == enums.JobType.RAY_TRACING:
      return self._ray_tracing_queue
    else:
      return None

  def there_are_jobs_running(self) -> bool:
    """Check if there are any jobs currently running.

    Returns:
        True if there are jobs running, False otherwise.
    """
    if self._is_prediction_queue_running:
      return True
    elif self._is_distance_analysis_queue_running:
      return True
    elif self._is_prediction_and_distance_analysis_queue_running:
      return True
    elif self._is_ray_tracing_queue_running:
      return True
    else:
      return False

  # <editor-fold desc="Prediction job">
  def create_prediction_job(
      self,
      a_project: "project.Project",
      the_prediction_protein_infos,
      the_prediction_configuration,
      the_project_lock: QtCore.QMutex,
      the_interface_manager: "interface_manager.InterfaceManager",
  ) -> tuple[job.PredictionJob, job_entry.JobEntryWidget]:
    """Creates a new prediction job.

    Args:
        a_project (project.Project): The project in which the prediction job will be created.
        the_prediction_protein_infos (Needs to be checked): A list of `ProteinInfo` objects containing information about the proteins for which the prediction is to be made.
        the_prediction_configuration (Needs to be checked): The configuration object for the prediction job.
        the_project_lock (QtCore.QMutex): A `QtCore.QMutex` object for ensuring thread safety when accessing the project.
        the_interface_manager (interface_manager.InterfaceManager): An `InterfaceManager` object for managing the interface.

    Returns:
        A tuple containing the created prediction job object and its corresponding job entry widget.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_project is None:
      logger.error("a_project is None.")
      raise exception.IllegalArgumentError("a_project is None.")
    if the_prediction_protein_infos is None:
      logger.error("the_prediction_protein_infos is None.")
      raise exception.IllegalArgumentError(
          "the_prediction_protein_infos is None."
      )
    if the_prediction_configuration is None:
      logger.error("the_prediction_configuration is None.")
      raise exception.IllegalArgumentError(
          "the_prediction_configuration is None."
      )
    if the_project_lock is None:
      logger.error("the_project_lock is None.")
      raise exception.IllegalArgumentError("the_project_lock is None.")
    if the_interface_manager is None:
      logger.error("the_interface_manager is None.")
      raise exception.IllegalArgumentError("the_interface_manager is None.")

    # </editor-fold>

    tmp_prediction_job = job.PredictionJob(
        self._main_socket,
        self._prediction_socket,
        self._general_purpose_socket,
        a_project,
        the_prediction_protein_infos,
        the_prediction_configuration,
        the_project_lock,
    )
    tmp_protein_names = []
    for tmp_protein_info in the_prediction_protein_infos:
      tmp_protein_names.append(tmp_protein_info.name)
    tmp_prediction_job.update_job_entry_signal.connect(
        the_interface_manager.update_job_entry
    )

    tmp_prediction_job.job_entry_widget = job_entry.JobEntryWidget(
        "Running ColabFold prediction",
        job_summary.JobBaseInformation(
            tmp_prediction_job.type,
            a_project.get_project_name(),
            tmp_protein_names,
            [],
            enums.JobProgress.WAITING,
        ),
    )
    tmp_prediction_job.job_entry_widget.ui.btn_cancel_job.clicked.connect(
        tmp_prediction_job.cancel_job
    )
    tmp_prediction_job.cancel_job_signal.connect(
        the_interface_manager.cancel_job
    )
    return tmp_prediction_job, tmp_prediction_job.job_entry_widget

  def put_prediction_job_into_queue(
      self, a_prediction_job: "job.PredictionJob"
  ) -> None:
    """Puts a prediction job into the appropriate queue.

    Args:
        a_prediction_job (job.PredictionJob): A prediction job to be added to the queue.

    Raises:
        exception.IllegalArgumentError: If `a_prediction_job` is None.
    """
    # <editor-fold desc="Checks">
    if a_prediction_job is None:
      logger.error("a_prediction_job is None.")
      raise exception.IllegalArgumentError("a_prediction_job is None.")

    # </editor-fold>

    if self._is_prediction_queue_running:
      logger.debug("Prediction queue is already running.")
      self._prediction_queue.put(a_prediction_job)
    else:
      logger.debug("Prediction queue needs to be started.")
      self._prediction_queue.put(a_prediction_job)
      self._prediction_queue_thread = tasks.LegacyTask(
          target=self._execute_prediction_job_queue,
          args=(0, 0),
          post_func=self._prediction_queue_finished,
      )
      self._prediction_queue_thread.start()

  def _execute_prediction_job_queue(
      self, placeholder_1: int, placeholder_2: int
  ) -> tuple[str, int]:
    """Executes the prediction job queue.

    Args:
        placeholder_1: Placeholder is needed for the LegacyTask.
        placeholder_2: Placeholder is needed for the LegacyTask.

    Returns:
        A tuple containing a status message ("Finished.") and a status code (0).
    """
    logger.debug("Starting prediction queue ...")

    while self._is_prediction_and_distance_analysis_queue_running:
      # Waiting for a prediction and distance analysis combi job to finish.
      time.sleep(60)

    self._is_prediction_queue_running = True
    while self._is_prediction_queue_running is True:
      self._prediction_queue_lock.lock()
      tmp_prediction_job: "job.PredictionJob" = self._prediction_queue.get()
      self._prediction_queue_lock.unlock()
      if tmp_prediction_job is None:
        self._is_prediction_queue_running = False
        break
      self.current_prediction_job = job_summary.PredictionJobSummary(
          tmp_prediction_job.prediction_protein_infos
      )
      tmp_prediction_job.run_job()
      if self._prediction_queue.empty():
        logger.info("The prediction queue is empty and will now end execution.")
        break
      self._prediction_queue.task_done()
    self._is_prediction_queue_running = False
    self.current_prediction_job = None
    logger.info("Prediction queue is finished.")
    return "Finished.", 0

  def _prediction_queue_finished(self) -> None:
    """Logs when the thread is no longer running."""
    logger.info("Prediction queue is empty and thread is no longer running.")

  def stop_prediction_queue_execution(self) -> None:
    """Stops the execution of the prediction queue.

    This method stops the execution of the prediction queue by adding a None value to the prediction queue.
    Once the None value is retrieved by the prediction queue consumer, the prediction queue execution will be stopped.
    """
    logger.debug("Stopping prediction queue ...")
    self._prediction_queue.put(None)

  # </editor-fold>

  # <editor-fold desc="Distance analysis job">
  def create_distance_analysis_job(
      self,
      a_project: "project.Project",
      the_project_lock: QtCore.QMutex,
      the_interface_manager: "interface_manager.InterfaceManager",
      a_list_with_analysis_names: list,
      a_cutoff: float,
      cycles: int,
  ) -> tuple[job.DistanceAnalysisJob, job_entry.JobEntryWidget]:
    """Creates a new distance analysis job.

    Args:
        a_project (project.Project): The project object to perform distance analysis on.
        the_project_lock (QtCore.QMutex): The lock object used for thread synchronization.
        the_interface_manager (interface_manager.InterfaceManager): The interface manager object that handles the user interface.
        a_list_with_analysis_names (list): A list of analysis names to be used for distance analysis.
        a_cutoff (float): The cutoff value for distance analysis.
        cycles (int): The number of cycles to perform for distance analysis.

    Returns:
        A tuple containing the distance analysis job object and the job entry widget object.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_project is None:
      logger.error("a_project is None.")
      raise exception.IllegalArgumentError("a_project is None.")
    if the_project_lock is None:
      logger.error("the_project_lock is None.")
      raise exception.IllegalArgumentError("the_project_lock is None.")
    if the_interface_manager is None:
      logger.error("the_interface_manager is None.")
      raise exception.IllegalArgumentError("the_interface_manager is None.")
    if a_list_with_analysis_names is None:
      logger.error("a_list_with_analysis_names is None.")
      raise exception.IllegalArgumentError(
          "a_list_with_analysis_names is None."
      )
    if a_cutoff is None:
      logger.error("a_cutoff is None.")
      raise exception.IllegalArgumentError("a_cutoff is None.")
    if cycles is None:
      logger.error("cycles is None.")
      raise exception.IllegalArgumentError("cycles is None.")

    # </editor-fold>

    tmp_distance_analysis_job = job.DistanceAnalysisJob(
        self._main_socket,
        self._distance_analysis_socket,
        a_project,
        the_project_lock,
        a_list_with_analysis_names,
        a_cutoff,
        cycles,
    )
    tmp_protein_pair_names = []
    for tmp_analysis_name in a_list_with_analysis_names:
      tmp_protein_1_name, tmp_protein_2_name = self.get_protein_names(
          tmp_analysis_name
      )
      print(tmp_protein_1_name)
      print(tmp_protein_2_name)
      if tmp_protein_1_name == tmp_protein_2_name:
        print(
            f"This analysis run name contains two identical protein names: {tmp_analysis_name}"
        )
        tmp_protein_pair_name = tmp_analysis_name.replace(";", "_1_", 1)
        tmp_protein_pair_name = tmp_protein_pair_name.replace(";", "_2_", 1)
        tmp_protein_pair_name = tmp_protein_pair_name.replace(",", "_")
      else:
        tmp_protein_pair_name = tmp_analysis_name.replace(";", "_")
        tmp_protein_pair_name = tmp_protein_pair_name.replace(",", "_")
      tmp_protein_pair_names.append(tmp_protein_pair_name)
    tmp_distance_analysis_job.update_job_entry_signal.connect(
        the_interface_manager.update_job_entry
    )
    tmp_distance_analysis_job.job_entry_widget = job_entry.JobEntryWidget(
        "Running distance analysis",
        job_summary.JobBaseInformation(
            tmp_distance_analysis_job.type,
            a_project.get_project_name(),
            [],
            tmp_protein_pair_names,
            enums.JobProgress.WAITING,
        ),
    )
    return tmp_distance_analysis_job, tmp_distance_analysis_job.job_entry_widget

  @staticmethod
  def get_protein_names(string: str) -> tuple[str, str]:
    """Gets the two protein names of the protein pair analysis name.

    Args:
        string (str): A string containing protein names separated by semicolons, with a specific format.

    Returns:
        A tuple of two protein names extracted from the input string. The first protein name is obtained by extracting the substring before the first semicolon.
        The second protein name is obtained by extracting the substring after the first occurrence of "_vs_" and before the second semicolon.

    Raises:
        exception.IllegalArgumentError: If `string` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if string is None or string == "":
      logger.error("string is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "string is either None or an empty string."
      )

    # </editor-fold>

    index_of_first_semicolon = string.find(";")
    sub1 = string[:index_of_first_semicolon]
    tmp_string = string.replace(";", "", 1)
    sub2 = tmp_string[tmp_string.find("_vs_") + 4 :]
    index_of_second_semicolon = sub2.find(";")
    sub3 = sub2[:index_of_second_semicolon]
    return sub1, sub3

  def put_distance_analysis_job_into_queue(
      self, a_distance_analysis_job: "job.PredictionJob"
  ) -> None:
    """Puts distance analysis job into queue.

    Args:
        a_distance_analysis_job: An instance of job.PredictionJob that represents the distance analysis job to be put into the queue.

    Raises:
        exception.IllegalArgumentError: If `a_distance_analysis_job` is None.
    """
    # <editor-fold desc="Checks">
    if a_distance_analysis_job is None:
      logger.error("a_distance_analysis_job is None.")
      raise exception.IllegalArgumentError("a_distance_analysis_job is None.")

    # </editor-fold>

    if self._is_distance_analysis_queue_running:
      logger.debug("Distance analysis queue is already running.")
      self._distance_analysis_queue.put(a_distance_analysis_job)
    else:
      logger.debug("Distance analysis queue needs to be started.")
      self._distance_analysis_queue.put(a_distance_analysis_job)
      self._distance_analysis_queue_thread = tasks.LegacyTask(
          target=self._execute_distance_analysis_job_queue,
          args=(0, 0),
          post_func=self._distance_analysis_queue_finished,
      )
      self._distance_analysis_queue_thread.start()

  def _execute_distance_analysis_job_queue(
      self, placeholder_1: int, placeholder_2: int
  ) -> tuple[str, int]:
    """Executes the distance analysis job queue.

    Args:
        placeholder_1: Placeholder is needed for the LegacyTask.
        placeholder_2: Placeholder is needed for the LegacyTask.

    Returns:
        A tuple containing a status message ("Finished.") and a status code (0).
    """
    logger.debug("Starting distance analysis queue ...")

    self._is_distance_analysis_queue_running = True
    while self._is_distance_analysis_queue_running is True:
      tmp_distance_analysis_job: "job.DistanceAnalysisJob" = (
        self._distance_analysis_queue.get()
      )
      if tmp_distance_analysis_job is None:
        self._is_distance_analysis_queue_running = False
        break
      self.current_distance_analysis_job = (
        job_summary.DistanceAnalysisJobSummary(
          tmp_distance_analysis_job.list_with_analysis_names,
        )
      )
      tmp_distance_analysis_job.job_entry_widget.job_base_information.job_progress = (
        enums.JobProgress.RUNNING
      )
      tmp_distance_analysis_job.run_job()
      if self._distance_analysis_queue.empty():
        logger.info(
          "The distance analysis queue is empty and will now end execution."
        )
        break
      self._distance_analysis_queue.task_done()
    self._is_distance_analysis_queue_running = False
    self.current_distance_analysis_job = None
    logger.info("Distance analysis queue is finished.")
    return "Finished.", 0

  def _distance_analysis_queue_finished(self) -> None:
    """Logs when the distance analysis queue is finished."""
    logger.info(
        "Distance analysis queue is empty and thread is no longer running."
    )

  def stop_distance_analysis_queue_execution(self) -> None:
    """Stops the queue after gracefully after all items until the None are executed."""
    logger.debug("Stopping distance analysis queue ...")
    self._distance_analysis_queue.put(None)

  # </editor-fold>

  # <editor-fold desc="Prediction and distance analysis job">
  def create_prediction_and_distance_analysis_job(
      self,
      a_prediction_job: "job.PredictionJob",
      a_distance_analysis_job: "job.DistanceAnalysisJob",
      the_interface_manager: "interface_manager.InterfaceManager",
  ) -> tuple[job.PredictionAndDistanceAnalysisJob, job_entry.JobEntryWidget]:
    """Creates a new prediction and distance analysis job.

    Args:
        a_prediction_job (job.PredictionJob): A PredictionJob object representing the prediction job to be included in the PredictionAndDistanceAnalysisJob.
        a_distance_analysis_job (job.DistanceAnalysisJob): A DistanceAnalysisJob object representing the distance analysis job to be included in the PredictionAndDistanceAnalysisJob.
        the_interface_manager (interface_manager.InterfaceManager): An InterfaceManager object responsible for managing the user interface.

    Returns:
        A tuple containing the created PredictionAndDistanceAnalysisJob object and its associated JobEntryWidget object.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_prediction_job is None:
      logger.error("a_prediction_job is None.")
      raise exception.IllegalArgumentError("a_prediction_job is None.")
    if a_distance_analysis_job is None:
      logger.error("a_distance_analysis_job is None.")
      raise exception.IllegalArgumentError("a_distance_analysis_job is None.")
    if the_interface_manager is None:
      logger.error("the_interface_manager is None.")
      raise exception.IllegalArgumentError("the_interface_manager is None.")

    # </editor-fold>

    tmp_prediction_and_distance_analysis_job = (
        job.PredictionAndDistanceAnalysisJob(
            a_prediction_job,
            a_distance_analysis_job,
        )
    )
    tmp_prediction_and_distance_analysis_job.update_job_entry_signal.connect(
        the_interface_manager.update_job_entry
    )
    tmp_job_entry_widget = job_entry.JobEntryWidget(
        "Running ColabFold prediction + distance analysis",
        job_summary.JobBaseInformation(
            tmp_prediction_and_distance_analysis_job.type,
            a_prediction_job.frozen_project.get_project_name(),
            a_prediction_job.job_entry_widget.job_base_information.protein_names,
            a_distance_analysis_job.job_entry_widget.job_base_information.protein_pair_names,
            enums.JobProgress.WAITING,
        ),
    )
    tmp_prediction_and_distance_analysis_job.job_entry_widget = (
        tmp_job_entry_widget
    )
    tmp_prediction_and_distance_analysis_job.prediction_job.job_entry_widget.job_base_information.job_type = (
        enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS
    )
    tmp_prediction_and_distance_analysis_job.distance_analysis_job.job_entry_widget.job_base_information.job_type = (
        enums.JobType.PREDICTION_AND_DISTANCE_ANALYSIS
    )
    return (
        tmp_prediction_and_distance_analysis_job,
        tmp_prediction_and_distance_analysis_job.job_entry_widget,
    )

  def _execute_prediction_and_distance_analysis_job_queue(
      self, placeholder_1: int, placeholder_2: int
  ):
    """Executes the prediction and distance analysis job queue.

    Args:
        placeholder_1: Placeholder is needed for the LegacyTask.
        placeholder_2: Placeholder is needed for the LegacyTask.

    Returns:
        A tuple containing a status message ("Finished.") and a status code (0).
    """
    logger.debug("Starting prediction and distance analysis queue ...")

    while self._is_prediction_queue_running:
      time.sleep(60)

    self._is_prediction_and_distance_analysis_queue_running = True
    while self._is_prediction_and_distance_analysis_queue_running is True:
      tmp_prediction_and_distance_analysis_job: (
        "job.PredictionAndDistanceAnalysisJob"
      ) = self._prediction_and_distance_analysis_queue.get()
      if tmp_prediction_and_distance_analysis_job is None:
        self._is_prediction_and_distance_analysis_queue_running = False
        break
      # self.current_prediction_and_distance_analysis_job = job_summary.DistanceAnalysisJobSummary(
      #     tmp_prediction_and_distance_analysis_job.list_with_analysis_names
      # )
      tmp_prediction_and_distance_analysis_job.run_job()
      if self._prediction_and_distance_analysis_queue.empty():
        logger.info(
          "The prediction and distance analysis queue is empty and will now end execution."
        )
        break
      self._prediction_and_distance_analysis_queue.task_done()
    self._is_prediction_and_distance_analysis_queue_running = False
    self.current_prediction_and_distance_analysis_job = None
    logger.info("Prediction and distance analysis queue is finished.")
    return "Finished.", 0

  def _prediction_and_distance_analysis_queue_finished(self) -> None:
    """Logs when the queue is finished."""
    logger.info(
        "Prediction and distance analysis queue is empty and thread is no longer running."
    )

  def stop_prediction_and_distance_analysis_queue_execution(self) -> None:
    """Stops the queue after gracefully after all items until the None are executed."""
    logger.debug("Stopping prediction and distance analysis queue ...")
    self._prediction_and_distance_analysis_queue.put(None)

  # </editor-fold>

  # <editor-fold desc="Ray-tracing job">
  def create_ray_tracing_job(
      self,
      the_destination_image_filepath: str,
      the_cached_session_filepath: str,
      image_ray_trace_mode: int,
      image_ray_texture: int,
      image_renderer: str,
      the_interface_manager: "interface_manager.InterfaceManager",
      a_project_name: str,
  ) -> tuple[job.RayTracingJob, job_entry.JobEntryWidget]:
    """Creates a new ray-tracing job.

    Args:
        the_destination_image_filepath (str): The file path for the destination image.
        the_cached_session_filepath (str): The file path for the cached session.
        image_ray_trace_mode (int): The ray trace mode for the image.
        image_ray_texture (int): The ray texture for the image.
        image_renderer (str): The renderer for the image.
        the_interface_manager (interface_manager.InterfaceManager): An instance of the InterfaceManager class.
        a_project_name (str): The name of the project.

    Returns:
        A tuple containing an instance of the RayTracingJob class and its corresponding JobEntryWidget instance.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None or if `the_destination_image_filepath`, `the_cached_session_filepath`, `image_renderer` or `a_project_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if (
        the_destination_image_filepath is None
        or the_destination_image_filepath == ""
    ):
      logger.error(
          "the_destination_image_filepath is either None or an empty string."
      )
      raise exception.IllegalArgumentError(
          "the_destination_image_filepath is either None or an empty string."
      )
    if the_cached_session_filepath is None or the_cached_session_filepath == "":
      logger.error(
          "the_cached_session_filepath is either None or an empty string."
      )
      raise exception.IllegalArgumentError(
          "the_cached_session_filepath is either None or an empty string."
      )
    if image_ray_trace_mode is None:
      logger.error("image_ray_trace_mode is None.")
      raise exception.IllegalArgumentError("image_ray_trace_mode is None.")
    if image_ray_texture is None:
      logger.error("image_ray_texture is None.")
      raise exception.IllegalArgumentError("image_ray_texture is None.")
    if image_renderer is None or image_renderer == "":
      logger.error("image_renderer is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_value_to_check is either None or an empty string."
      )
    if the_interface_manager is None:
      logger.error("the_interface_manager is None.")
      raise exception.IllegalArgumentError("the_interface_manager is None.")
    if a_project_name is None or a_project_name == "":
      logger.error("a_project_name is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_project_name is either None or an empty string."
      )

    # </editor-fold>

    tmp_ray_tracing_job = job.RayTracingJob(
        self._main_socket,
        self._ray_tracing_socket,
        the_destination_image_filepath,
        the_cached_session_filepath,
        image_ray_trace_mode,
        image_ray_texture,
        image_renderer,
    )
    tmp_ray_tracing_job.update_job_entry_signal.connect(
        the_interface_manager.update_job_entry
    )
    tmp_ray_tracing_job.job_entry_widget = job_entry.JobEntryWidget(
        "Creating ray-traced image",
        job_summary.JobBaseInformation(
            enums.JobType.RAY_TRACING,
            a_project_name,
            [],
            [],
            enums.JobProgress.WAITING,
        ),
    )
    tmp_ray_tracing_job.job_entry_widget.job_base_information.add_image_filepath(
        the_destination_image_filepath
    )
    return tmp_ray_tracing_job, tmp_ray_tracing_job.job_entry_widget

  def put_ray_tracing_job_into_queue(
      self, a_ray_tracing_job: "job.RayTracingJob"
  ) -> None:
    """Puts a given ray tracing job into the ray tracing queue and starts the queue if not already running.

    Args:
        a_ray_tracing_job (job.RayTracingJob): The ray tracing job to be added to the queue.

    Raises:
        exception.IllegalArgumentError: If `a_ray_tracing_job` is None.
    """
    # <editor-fold desc="Checks">
    if a_ray_tracing_job is None:
      logger.error("a_ray_tracing_job is None.")
      raise exception.IllegalArgumentError("a_ray_tracing_job is None.")

    # </editor-fold>

    if self._is_ray_tracing_queue_running:
      logger.debug("Ray-tracing queue is already running.")
      self._ray_tracing_queue.put(a_ray_tracing_job)
    else:
      logger.debug("Ray-tracing queue needs to be started.")
      self._ray_tracing_queue.put(a_ray_tracing_job)
      self._ray_tracing_queue_thread = tasks.LegacyTask(
          target=self._execute_ray_tracing_job_queue,
          args=(0, 0),
          post_func=self._ray_tracing_queue_finished,
      )
      self._ray_tracing_queue_thread.start()

  def _execute_ray_tracing_job_queue(
      self, placeholder_1: int, placeholder_2: int
  ) -> tuple[str, int]:
    """Executes the ray-tracing job queue.

    Args:
        placeholder_1: Placeholder is needed for the LegacyTask.
        placeholder_2: Placeholder is needed for the LegacyTask.

    Returns:
        A tuple containing a status message ("Finished.") and a status code (0).
    """
    logger.debug("Starting ray-tracing queue ...")

    self._is_ray_tracing_queue_running = True
    while self._is_ray_tracing_queue_running is True:
      tmp_ray_tracing_job: "job.RayTracingJob" = self._ray_tracing_queue.get()
      if tmp_ray_tracing_job is None:
        self._is_ray_tracing_queue_running = False
        break
      # Do something with the task
      tmp_ray_tracing_job.run_job()
      if self._ray_tracing_queue.empty():
        logger.info(
          "The ray-tracing queue is empty and will now end execution."
        )
        break
      self._ray_tracing_queue.task_done()
    logger.info("Ray-tracing queue is finished.")
    self._is_ray_tracing_queue_running = False
    return "Finished.", 0

  def _ray_tracing_queue_finished(self) -> None:
    """Logs when the queue is finished."""
    logger.info("Ray-tracing queue is empty and thread is no longer running.")

  # </editor-fold>
