import logging
import threading
import queue
from PyQt5 import QtCore
from pyssa.controller import database_manager
from pyssa.internal.thread import tasks
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants, enums
from pyssa.internal.data_structures.data_classes import database_operation

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class DatabaseThread(QtCore.QObject):

    _database_filepath: str

    def __init__(self, the_database_filepath: str):
        super(DatabaseThread, self).__init__()
        self._database_filepath = the_database_filepath
        self._queue: "queue.Queue" = queue.Queue()
        self._is_queue_running = False
        self._setup_operations_mapping()
        self._thread = None

    def _setup_operations_mapping(self):
        self._operations_mapping = {
            enums.SQLQueryType.INSERT_NEW_PROTEIN: self.__wrapper_insert_new_protein,
            enums.SQLQueryType.DELETE_EXISTING_PROTEIN: self.__wrapper_delete_existing_protein,
            enums.SQLQueryType.INSERT_NEW_PROTEIN_PAIR: self.__wrapper_insert_new_protein_pair,
            enums.SQLQueryType.DELETE_EXISTING_PROTEIN_PAIR: self.__wrapper_delete_existing_protein_pair,
            enums.SQLQueryType.UPDATE_PYMOL_SESSION_PROTEIN: self.__wrapper_update_pymol_session_of_protein,
            enums.SQLQueryType.UPDATE_PYMOL_SESSION_PROTEIN_PAIR: self.__wrapper_update_pymol_session_of_protein_pair,
            enums.SQLQueryType.INSERT_NEW_SEQUENCE: self.__wrapper_insert_new_sequence,
            enums.SQLQueryType.DELETE_EXISTING_SEQUENCE: self.__wrapper_delete_existing_sequence,
            enums.SQLQueryType.UPDATE_SEQUENCE_NAME: self.__wrapper_update_sequence_name,
        }

    def set_database_filepath(self, a_filepath: str) -> None:
        """Sets the database filepath into the object."""
        self._database_filepath = a_filepath

    def put_database_operation_into_queue(self, a_database_operation: "database_operation.DatabaseOperation"):
        """Puts a DatabaseOperation object into the queue."""
        if self._is_queue_running:
            self._queue.put(a_database_operation)
        else:
            self._queue.put(a_database_operation)
            self._thread = tasks.Task(
                target=self._execute_queue,
                args=(0, 0),
                post_func=self._queue_finished
            )
            self._thread.start()

    def queue_is_running(self):
        return self._is_queue_running

    def _execute_queue(self, placeholder_1, placeholder_2):
        with database_manager.DatabaseManager(self._database_filepath, "database_thread") as tmp_database_manager:
            while True:
                self._is_queue_running = True
                tmp_database_operation: "database_operation.DatabaseOperation" = self._queue.get()
                # if tmp_database_operation is None:
                #     break
                if tmp_database_operation.sql_query_type is enums.SQLQueryType.CLOSE_PROJECT:
                    logger.info("Received request to close the project.")
                    break
                logger.info(f"Running {tmp_database_operation.sql_query_type} database operation.")
                self._process_work(tmp_database_manager, tmp_database_operation)
                logger.info(f"Finished {tmp_database_operation.sql_query_type} database operation.")
                self._queue.task_done()
                if self._queue.empty():
                    logger.info("The queue is empty and will now end execution.")
        self._is_queue_running = False
        return "Finished.", 0

    def _queue_finished(self):
        logger.info("Queue is empty and thread is no longer running.")

    def _process_work(self, the_db_manager, a_database_operation: database_operation.DatabaseOperation):
        # Get the function based on the operation, defaulting to a generic function
        operation_function = self._operations_mapping.get(a_database_operation.sql_query_type, self._default_operation)
        # Call the selected function with the provided data
        operation_function(the_db_manager, a_database_operation.buffered_data)

    @staticmethod
    def _default_operation(placeholder_1, placeholder_2):
        constants.PYSSA_LOGGER.warning("This operation does not exists!")

    @staticmethod
    def __wrapper_insert_new_protein(the_db_manager, the_buffered_data: tuple):
        _, tmp_protein = the_buffered_data
        the_db_manager.insert_new_protein(tmp_protein)

    @staticmethod
    def __wrapper_delete_existing_protein(the_db_manager, the_buffered_data: tuple):
        _, tmp_protein_id = the_buffered_data
        the_db_manager.delete_existing_protein(tmp_protein_id)

    @staticmethod
    def __wrapper_insert_new_protein_pair(the_db_manager, the_buffered_data: tuple):
        _, tmp_protein_pair = the_buffered_data
        the_db_manager.insert_new_protein_pair(the_buffered_data[1])

    @staticmethod
    def __wrapper_delete_existing_protein_pair(the_db_manager, the_buffered_data: tuple):
        _, tmp_protein_pair_id = the_buffered_data
        the_db_manager.delete_existing_protein_pair(tmp_protein_pair_id)

    @staticmethod
    def __wrapper_update_pymol_session_of_protein(the_db_manager, the_buffered_data: tuple):
        _, tmp_protein = the_buffered_data
        the_db_manager.update_pymol_session_of_protein(tmp_protein.get_id(), tmp_protein.pymol_session)

    @staticmethod
    def __wrapper_update_pymol_session_of_protein_pair(the_db_manager, the_buffered_data: tuple):
        _, tmp_new_pymol_session, tmp_protein_pair = the_buffered_data
        the_db_manager.update_pymol_session_of_protein_pair(tmp_protein_pair.get_id(), tmp_protein_pair.pymol_session)

    @staticmethod
    def __wrapper_insert_new_sequence(the_db_manager: "database_manager.DatabaseManager", the_buffered_data: tuple):
        _, tmp_seq_record_obj = the_buffered_data
        the_db_manager.insert_new_sequence(tmp_seq_record_obj)

    @staticmethod
    def __wrapper_delete_existing_sequence(the_db_manager: "database_manager.DatabaseManager", the_buffered_data: tuple):
        _, tmp_seq_record = the_buffered_data
        the_db_manager.delete_existing_sequence(tmp_seq_record.name)

    @staticmethod
    def __wrapper_update_sequence_name(the_db_manager, the_buffered_data: tuple):
        _, tmp_new_seq_name, tmp_old_seq_name, tmp_seq = the_buffered_data
        the_db_manager.update_sequence_name(tmp_new_seq_name, tmp_old_seq_name, tmp_seq)

    def stop(self):
        # Stop the thread
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
