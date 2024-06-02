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
"""Module contains the database thread class."""
import logging
import queue
from PyQt5 import QtCore
from src.pyssa.controller import database_manager
from src.pyssa.internal.thread import tasks
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import constants, enums, exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class DatabaseThread(QtCore.QObject):
  """Custom class that can handle the database from a different thread."""

  _database_filepath: str
  """The filepath of the database file."""

  def __init__(self, the_database_filepath: str) -> None:
    """Constructor.

    Args:
        the_database_filepath: The filepath of the database used by this database thread.

    Raises:
        exception.IllegalArgumentError: If `the_database_filepath` is None.
    """
    # <editor-fold desc="Checks">
    if the_database_filepath is None:
      logger.error("the_database_filepath is None")
      raise exception.IllegalArgumentError("the_database_filepath is None.")

    # </editor-fold>

    super(DatabaseThread, self).__init__()
    self._database_filepath = the_database_filepath
    self._queue: "queue.Queue" = queue.Queue()
    self._is_queue_running = False
    self._stop_thread = False
    self._setup_operations_mapping()
    self._thread = None

  def _setup_operations_mapping(self) -> None:
    """Sets up the operations mapping for the class."""
    self._operations_mapping = {
        enums.SQLQueryType.INSERT_NEW_PROTEIN: self.__wrapper_insert_new_protein,
        enums.SQLQueryType.DELETE_EXISTING_PROTEIN: self.__wrapper_delete_existing_protein,
        enums.SQLQueryType.DELETE_SPECIFIC_CHAIN: self.__wrapper_delete_specific_chain,
        enums.SQLQueryType.INSERT_NEW_PROTEIN_PAIR: self.__wrapper_insert_new_protein_pair,
        enums.SQLQueryType.DELETE_EXISTING_PROTEIN_PAIR: self.__wrapper_delete_existing_protein_pair,
        enums.SQLQueryType.UPDATE_PYMOL_SESSION_PROTEIN: self.__wrapper_update_pymol_session_of_protein,
        enums.SQLQueryType.UPDATE_PYMOL_SESSION_PROTEIN_PAIR: self.__wrapper_update_pymol_session_of_protein_pair,
        enums.SQLQueryType.INSERT_NEW_SEQUENCE: self.__wrapper_insert_new_sequence,
        enums.SQLQueryType.DELETE_EXISTING_SEQUENCE: self.__wrapper_delete_existing_sequence,
        enums.SQLQueryType.UPDATE_SEQUENCE_NAME: self.__wrapper_update_sequence_name,
    }

  def set_database_filepath(self, a_filepath: str) -> None:
    """Sets the filepath of the database.

    Args:
        a_filepath (str): The filepath of the database.

    Raises:
        exception.IllegalArgumentError: If `a_filepath` is either None or an empty string.
    """
    # <editor-fold desc="Checks">
    if a_filepath is None or a_filepath == "":
      logger.error("a_filepath is either None or an empty string.")
      raise exception.IllegalArgumentError(
          "a_filepath is either None or an empty string."
      )

    # </editor-fold>

    self._database_filepath = a_filepath

  def put_database_operation_into_queue(
      self, a_database_operation: "database_operation.DatabaseOperation"
  ) -> None:
    """Puts a database operation object into the queue for execution.

    Args:
        a_database_operation (database_operation.DatabaseOperation): The database operation object to be put into the queue.

    Raises:
        exception.IllegalArgumentError: If `a_database_operation` is None.
    """
    # <editor-fold desc="Checks">
    if a_database_operation is None:
      logger.error("a_database_operation is None.")
      raise exception.IllegalArgumentError("a_database_operation is None.")

    # </editor-fold>

    if self._is_queue_running:
      self._queue.put(a_database_operation)
    else:
      self._queue.put(a_database_operation)
      self._thread = tasks.LegacyTask(
          target=self._execute_queue,
          args=(0, 0),
          post_func=self._queue_finished,
      )
      self._thread.start()

  def queue_is_running(self) -> bool:
    """Check if the queue is currently running.

    Returns:
        True if the queue is running, False otherwise.
    """
    return self._is_queue_running

  def _execute_queue(
      self, placeholder_1: int, placeholder_2: int
  ) -> tuple[str, int]:
    """Executes the operations in the queue.

    Args:
        placeholder_1 (int): The first placeholder.
        placeholder_2 (int): The second placeholder.

    Returns:
        tuple[str, int]: A tuple with the message "Finished." and the integer 0.
    """
    try:
      with database_manager.DatabaseManager(
          self._database_filepath, "database_thread"
      ) as tmp_database_manager:
        while self._stop_thread is False:
          self._is_queue_running = True
          tmp_database_operation: "database_operation.DatabaseOperation" = (
              self._queue.get()
          )
          if (
              tmp_database_operation.sql_query_type
              is enums.SQLQueryType.CLOSE_PROJECT
          ):
            logger.info("Received request to close the project.")
            self._stop_thread = True
            break
          logger.info(
              f"Running {tmp_database_operation.sql_query_type} database operation."
          )
          self._process_work(tmp_database_manager, tmp_database_operation)
          logger.info(
              f"Finished {tmp_database_operation.sql_query_type} database operation."
          )
          self._queue.task_done()
          if self._queue.empty():
            logger.info("The queue is empty and will now end execution.")
      self._is_queue_running = False
      self._stop_thread = False
    except Exception as e:
      logger.error(e)
      return "", 1
    else:
      return "Finished.", 0

  def _queue_finished(self) -> None:
    """Check if the queue is empty and the thread is no longer running."""
    logger.info("Queue is empty and thread is no longer running.")

  def _process_work(
      self,
      the_db_manager: "database_manager.DatabaseManager",
      a_database_operation: "database_operation.DatabaseOperation",
  ) -> None:
    """Processes the operations in the queue.

    Args:
        the_db_manager (database_manager.DatabaseManager): The database manager that will be used to perform the database operation.
        a_database_operation (database_operation.DatabaseOperation): The database operation that needs to be processed.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_db_manager is None:
      logger.error("the_db_manager is None.")
      raise exception.IllegalArgumentError("the_db_manager is None.")
    if a_database_operation is None:
      logger.error("a_database_operation is None.")
      raise exception.IllegalArgumentError("a_database_operation is None.")

    # </editor-fold>

    # Get the function based on the operation, defaulting to a generic function
    operation_function = self._operations_mapping.get(
        a_database_operation.sql_query_type, self._default_operation
    )
    # Call the selected function with the provided data
    operation_function(the_db_manager, a_database_operation.buffered_data)

  @staticmethod
  def _default_operation(placeholder_1: int, placeholder_2: int) -> None:
    """Performs a default operation and logs a warning message.

    Args:
        placeholder_1 (int): The first placeholder value.
        placeholder_2 (int): The second placeholder value.
    """
    constants.PYSSA_LOGGER.warning("This operation does not exists!")

  @staticmethod
  def __wrapper_insert_new_protein(
      the_db_manager: "database_manager.DatabaseManager",
      the_buffered_data: tuple,
  ) -> None:
    """Wrapper function for inserting a new protein into the database.

    Args:
        the_db_manager (database_manager.DatabaseManager): An instance of the class "database_manager.DatabaseManager" that provides access to the database and methods to manipulate data.
        the_buffered_data (tuple): A tuple containing the data to be inserted into the database.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_db_manager is None:
      logger.error("the_db_manager is None.")
      raise exception.IllegalArgumentError("the_db_manager is None.")
    if the_buffered_data is None:
      logger.error("the_buffered_data is None.")
      raise exception.IllegalArgumentError("the_buffered_data is None.")
    # </editor-fold>

    _, tmp_protein = the_buffered_data
    the_db_manager.insert_new_protein(tmp_protein)

  @staticmethod
  def __wrapper_delete_existing_protein(
      the_db_manager: "database_manager.DatabaseManager",
      the_buffered_data: tuple,
  ) -> None:
    """Wrapper function for deleting an existing protein.

    Args:
        the_db_manager (database_manager.DatabaseManager): An instance of the class "database_manager.DatabaseManager" that provides access to the database and methods to manipulate data.
        the_buffered_data (tuple): A tuple containing the data to be inserted into the database.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_db_manager is None:
      logger.error("the_db_manager is None.")
      raise exception.IllegalArgumentError("the_db_manager is None.")
    if the_buffered_data is None:
      logger.error("the_buffered_data is None.")
      raise exception.IllegalArgumentError("the_buffered_data is None.")
    # </editor-fold>

    _, tmp_protein_id = the_buffered_data
    the_db_manager.delete_existing_protein(tmp_protein_id)

  @staticmethod
  def __wrapper_delete_specific_chain(
      the_db_manager: "database_manager.DatabaseManager",
      the_buffered_data: tuple,
  ) -> None:
    """Wrapper function for deleting a specific chain.

    Args:
        the_db_manager (database_manager.DatabaseManager): An instance of the class "database_manager.DatabaseManager" that provides access to the database and methods to manipulate data.
        the_buffered_data (tuple): A tuple containing the data to be inserted into the database.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_db_manager is None:
      logger.error("the_db_manager is None.")
      raise exception.IllegalArgumentError("the_db_manager is None.")
    if the_buffered_data is None:
      logger.error("the_buffered_data is None.")
      raise exception.IllegalArgumentError("the_buffered_data is None.")
    # </editor-fold>

    tmp_protein_id, tmp_chain_id = the_buffered_data
    the_db_manager.delete_specific_chain(tmp_protein_id, tmp_chain_id)

  @staticmethod
  def __wrapper_insert_new_protein_pair(
      the_db_manager: "database_manager.DatabaseManager",
      the_buffered_data: tuple,
  ) -> None:
    """Wrapper function for inserting a new protein pair into the database.

    Args:
        the_db_manager (database_manager.DatabaseManager): An instance of the class "database_manager.DatabaseManager" that provides access to the database and methods to manipulate data.
        the_buffered_data (tuple): A tuple containing the data to be inserted into the database.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_db_manager is None:
      logger.error("the_db_manager is None.")
      raise exception.IllegalArgumentError("the_db_manager is None.")
    if the_buffered_data is None:
      logger.error("the_buffered_data is None.")
      raise exception.IllegalArgumentError("the_buffered_data is None.")
    # </editor-fold>

    _, tmp_protein_pair = the_buffered_data
    the_db_manager.insert_new_protein_pair(the_buffered_data[1])

  @staticmethod
  def __wrapper_delete_existing_protein_pair(
      the_db_manager: "database_manager.DatabaseManager",
      the_buffered_data: tuple,
  ) -> None:
    """Wrapper function for deleting an existing protein pair.

    Args:
        the_db_manager (database_manager.DatabaseManager): An instance of the class "database_manager.DatabaseManager" that provides access to the database and methods to manipulate data.
        the_buffered_data (tuple): A tuple containing the data to be inserted into the database.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_db_manager is None:
      logger.error("the_db_manager is None.")
      raise exception.IllegalArgumentError("the_db_manager is None.")
    if the_buffered_data is None:
      logger.error("the_buffered_data is None.")
      raise exception.IllegalArgumentError("the_buffered_data is None.")
    # </editor-fold>

    _, tmp_protein_pair_id = the_buffered_data
    the_db_manager.delete_existing_protein_pair(tmp_protein_pair_id)

  @staticmethod
  def __wrapper_update_pymol_session_of_protein(
      the_db_manager: "database_manager.DatabaseManager",
      the_buffered_data: tuple,
  ) -> None:
    """Wrapper function for updating a protein pymol session.

    Args:
        the_db_manager (database_manager.DatabaseManager): An instance of the class "database_manager.DatabaseManager" that provides access to the database and methods to manipulate data.
        the_buffered_data (tuple): A tuple containing the data to be inserted into the database.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_db_manager is None:
      logger.error("the_db_manager is None.")
      raise exception.IllegalArgumentError("the_db_manager is None.")
    if the_buffered_data is None:
      logger.error("the_buffered_data is None.")
      raise exception.IllegalArgumentError("the_buffered_data is None.")
    # </editor-fold>

    _, tmp_protein = the_buffered_data
    the_db_manager.update_pymol_session_of_protein(
        tmp_protein.get_id(), tmp_protein.pymol_session
    )

  @staticmethod
  def __wrapper_update_pymol_session_of_protein_pair(
      the_db_manager: "database_manager.DatabaseManager",
      the_buffered_data: tuple,
  ) -> None:
    """Wrapper function for updating a protein pair pymol session.

    Args:
        the_db_manager (database_manager.DatabaseManager): An instance of the class "database_manager.DatabaseManager" that provides access to the database and methods to manipulate data.
        the_buffered_data (tuple): A tuple containing the data to be inserted into the database.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_db_manager is None:
      logger.error("the_db_manager is None.")
      raise exception.IllegalArgumentError("the_db_manager is None.")
    if the_buffered_data is None:
      logger.error("the_buffered_data is None.")
      raise exception.IllegalArgumentError("the_buffered_data is None.")
    # </editor-fold>

    _, tmp_new_pymol_session, tmp_protein_pair = the_buffered_data
    the_db_manager.update_pymol_session_of_protein_pair(
        tmp_protein_pair.get_id(), tmp_protein_pair.pymol_session
    )

  @staticmethod
  def __wrapper_insert_new_sequence(
      the_db_manager: "database_manager.DatabaseManager",
      the_buffered_data: tuple,
  ) -> None:
    """Wrapper function for inserting a new sequence into the database.

    Args:
        the_db_manager (database_manager.DatabaseManager): An instance of the class "database_manager.DatabaseManager" that provides access to the database and methods to manipulate data.
        the_buffered_data (tuple): A tuple containing the data to be inserted into the database.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_db_manager is None:
      logger.error("the_db_manager is None.")
      raise exception.IllegalArgumentError("the_db_manager is None.")
    if the_buffered_data is None:
      logger.error("the_buffered_data is None.")
      raise exception.IllegalArgumentError("the_buffered_data is None.")
    # </editor-fold>

    _, tmp_seq_record_obj = the_buffered_data
    the_db_manager.insert_new_sequence(tmp_seq_record_obj)

  @staticmethod
  def __wrapper_delete_existing_sequence(
      the_db_manager: "database_manager.DatabaseManager",
      the_buffered_data: tuple,
  ) -> None:
    """Wrapper function for deleting an existing sequence.

    Args:
        the_db_manager (database_manager.DatabaseManager): An instance of the class "database_manager.DatabaseManager" that provides access to the database and methods to manipulate data.
        the_buffered_data (tuple): A tuple containing the data to be inserted into the database.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_db_manager is None:
      logger.error("the_db_manager is None.")
      raise exception.IllegalArgumentError("the_db_manager is None.")
    if the_buffered_data is None:
      logger.error("the_buffered_data is None.")
      raise exception.IllegalArgumentError("the_buffered_data is None.")
    # </editor-fold>

    _, tmp_seq_record = the_buffered_data
    the_db_manager.delete_existing_sequence(tmp_seq_record.name)

  @staticmethod
  def __wrapper_update_sequence_name(
      the_db_manager: "database_manager.DatabaseManager",
      the_buffered_data: tuple,
  ) -> None:
    """Wrapper function for updating a sequence name.

    Args:
        the_db_manager (database_manager.DatabaseManager): An instance of the class "database_manager.DatabaseManager" that provides access to the database and methods to manipulate data.
        the_buffered_data (tuple): A tuple containing the data to be inserted into the database.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if the_db_manager is None:
      logger.error("the_db_manager is None.")
      raise exception.IllegalArgumentError("the_db_manager is None.")
    if the_buffered_data is None:
      logger.error("the_buffered_data is None.")
      raise exception.IllegalArgumentError("the_buffered_data is None.")
    # </editor-fold>

    _, tmp_new_seq_name, tmp_old_seq_name, tmp_seq = the_buffered_data
    the_db_manager.update_sequence_name(
        tmp_new_seq_name, tmp_old_seq_name, tmp_seq
    )

  def stop(self) -> None:
    """Stops the thread.

    This method sets a stop event that can be used to gracefully stop the thread.
    """
    self._stop_thread = True
    self._stop_event.set()


# class DatabaseThread(threading.Thread):
#
#     _database_filepath: str
#
#     def __init__(self, the_database_filepath: str):
#         super(DatabaseThread, self).__init__()
#         self._database_filepath = the_database_filepath
#         self._queue: "queue.Queue" = queue.Queue()
#         self._setup_operations_mapping()
#         self._stop_event = threading.Event()
#
#     def _setup_operations_mapping(self):
#         self._operations_mapping = {
#             enums.SQLQueryType.INSERT_NEW_PROTEIN: self.__wrapper_insert_new_protein,
#             enums.SQLQueryType.DELETE_EXISTING_PROTEIN: self.__wrapper_delete_existing_protein,
#             enums.SQLQueryType.INSERT_NEW_PROTEIN_PAIR: self.__wrapper_insert_new_protein_pair,
#             enums.SQLQueryType.DELETE_EXISTING_PROTEIN_PAIR: self.__wrapper_delete_existing_protein_pair,
#             enums.SQLQueryType.UPDATE_PYMOL_SESSION_PROTEIN: self.__wrapper_update_pymol_session_of_protein,
#             enums.SQLQueryType.UPDATE_PYMOL_SESSION_PROTEIN_PAIR: self.__wrapper_update_pymol_session_of_protein_pair,
#             enums.SQLQueryType.INSERT_NEW_SEQUENCE: self.__wrapper_insert_new_sequence,
#             enums.SQLQueryType.DELETE_EXISTING_SEQUENCE: self.__wrapper_delete_existing_sequence,
#             enums.SQLQueryType.UPDATE_SEQUENCE_NAME: self.__wrapper_update_sequence_name,
#         }
#
#     def set_database_filepath(self, a_filepath: str) -> None:
#         """Sets the database filepath into the object."""
#         self._database_filepath = a_filepath
#
#     def put_database_operation_into_queue(self, a_database_operation: database_operation.DatabaseOperation):
#         """Puts a DatabaseOperation object into the queue."""
#         self._queue.put(a_database_operation)
#
#     def run(self):
#         #tmp_database_manager = database_manager.DatabaseManager(self._database_filepath)
#         #tmp_database_manager.open_project_database()
#         # with database_manager.DatabaseManager(self._database_filepath) as db_manager:
#         #     while not self._stop_event.is_set():
#         #         try:
#         #             tmp_database_operation: database_operation.DatabaseOperation = self._queue.get(timeout=1)
#         #             if tmp_database_operation.sql_query_type is enums.SQLQueryType.CLOSE_PROJECT:
#         #                 logger.info("Received request to close the project.")
#         #                 break
#         #             logger.info(f"Running {tmp_database_operation.sql_query_type} database operation.")
#         #             self._process_work(db_manager, tmp_database_operation)
#         #             logger.info(f"Finished {tmp_database_operation.sql_query_type} database operation.")
#         #         except queue.Empty:
#         #             pass  # Continue checking for new tasks
#         #     self._stop_event.set()
#         #     self._queue.join()
#         with database_manager.DatabaseManager(self._database_filepath, "database_thread") as tmp_database_manager:
#             while not self._stop_event.is_set():
#                 try:
#                     tmp_database_operation: database_operation.DatabaseOperation = self._queue.get(timeout=1)
#                     if tmp_database_operation.sql_query_type is enums.SQLQueryType.CLOSE_PROJECT:
#                         logger.info("Received request to close the project.")
#                         break
#                     logger.info(f"Running {tmp_database_operation.sql_query_type} database operation.")
#                     self._process_work(tmp_database_manager, tmp_database_operation)
#                     logger.info(f"Finished {tmp_database_operation.sql_query_type} database operation.")
#                 except queue.Empty:
#                     pass  # Continue checking for new tasks
#             self._stop_event.set()
#         self._queue.join()
#
#     def join(self, timeout=None):
#         """ Stop the thread. """
#         threading.Thread.join(self, timeout)
#
#     def _process_work(self, the_db_manager, a_database_operation: database_operation.DatabaseOperation):
#         # Get the function based on the operation, defaulting to a generic function
#         operation_function = self._operations_mapping.get(a_database_operation.sql_query_type, self._default_operation)
#         # Call the selected function with the provided data
#         operation_function(the_db_manager, a_database_operation.buffered_data)
#
#     @staticmethod
#     def _default_operation(placeholder_1, placeholder_2):
#         constants.PYSSA_LOGGER.warning("This operation does not exists!")
#
#     @staticmethod
#     def __wrapper_insert_new_protein(the_db_manager, the_buffered_data: tuple):
#         _, tmp_protein = the_buffered_data
#         the_db_manager.insert_new_protein(tmp_protein)
#
#     @staticmethod
#     def __wrapper_delete_existing_protein(the_db_manager, the_buffered_data: tuple):
#         _, tmp_protein_id = the_buffered_data
#         the_db_manager.delete_existing_protein(tmp_protein_id)
#
#     @staticmethod
#     def __wrapper_insert_new_protein_pair(the_db_manager, the_buffered_data: tuple):
#         _, tmp_protein_pair = the_buffered_data
#         the_db_manager.insert_new_protein_pair(the_buffered_data[1])
#
#     @staticmethod
#     def __wrapper_delete_existing_protein_pair(the_db_manager, the_buffered_data: tuple):
#         _, tmp_protein_pair_id = the_buffered_data
#         the_db_manager.delete_existing_protein_pair(tmp_protein_pair_id)
#
#     @staticmethod
#     def __wrapper_update_pymol_session_of_protein(the_db_manager, the_buffered_data: tuple):
#         _, tmp_protein = the_buffered_data
#         #tmp_protein.save_pymol_session_as_base64_string()
#         the_db_manager.update_pymol_session_of_protein(tmp_protein.get_id(), tmp_protein.pymol_session)
#
#     @staticmethod
#     def __wrapper_update_pymol_session_of_protein_pair(the_db_manager, the_buffered_data: tuple):
#         _, tmp_new_pymol_session, tmp_protein_pair = the_buffered_data
#         #tmp_protein_pair.save_session_of_protein_pair()
#         the_db_manager.update_pymol_session_of_protein_pair(tmp_protein_pair.get_id(), tmp_protein_pair.pymol_session)
#
#     @staticmethod
#     def __wrapper_insert_new_sequence(the_db_manager: "database_manager.DatabaseManager", the_buffered_data: tuple):
#         _, tmp_seq_record_obj = the_buffered_data
#         the_db_manager.insert_new_sequence(tmp_seq_record_obj)
#
#     @staticmethod
#     def __wrapper_delete_existing_sequence(the_db_manager: "database_manager.DatabaseManager", the_buffered_data: tuple):
#         _, tmp_seq_record = the_buffered_data
#         the_db_manager.delete_existing_sequence(tmp_seq_record.name)
#
#     @staticmethod
#     def __wrapper_update_sequence_name(the_db_manager, the_buffered_data: tuple):
#         _, tmp_new_seq_name, tmp_old_seq_name, tmp_seq = the_buffered_data
#         the_db_manager.update_sequence_name(tmp_new_seq_name, tmp_old_seq_name, tmp_seq)
#
#     def stop(self):
#         # Stop the thread
#         self._stop_event.set()
