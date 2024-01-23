import os.path
import threading
import queue
import sqlite3
from pyssa.controller import database_manager
from pyssa.util import constants, enums
from pyssa.internal.data_structures.data_classes import database_operation

class DatabaseThread(threading.Thread):

    _database_filepath: str

    def __init__(self, the_database_filepath: str):
        super(DatabaseThread, self).__init__()
        self._database_filepath = the_database_filepath
        self._queue: "queue.Queue" = queue.Queue()
        self._setup_operations_mapping()
        self._stop_event = threading.Event()

    def _setup_operations_mapping(self):
        self._operations_mapping = {
            enums.SQLQueryType.INSERT_NEW_PROTEIN: self.__wrapper_insert_new_protein,
            enums.SQLQueryType.DELETE_EXISTING_PROTEIN: self.__wrapper_delete_existing_protein
        }

    def set_database_filepath(self, a_filepath: str) -> None:
        """Sets the database filepath into the object."""
        self._database_filepath = a_filepath

    def put_database_operation_into_queue(self, a_database_operation: database_operation.DatabaseOperation):
        """Puts a DatabaseOperation object into the queue."""
        self._queue.put(a_database_operation)

    def run(self):
        with database_manager.DatabaseManager(self._database_filepath) as db_manager:
            while not self._stop_event.is_set():
                    try:
                        tmp_database_operation = self._queue.get(timeout=1)
                        self._process_work(db_manager, tmp_database_operation)
                    except queue.Empty:
                        pass  # Continue checking for new tasks

    def _process_work(self, the_db_manager, a_database_operation: database_operation.DatabaseOperation):
        # Get the function based on the operation, defaulting to a generic function
        operation_function = self._operations_mapping.get(a_database_operation.sql_query_type, self._default_operation)
        # Call the selected function with the provided data
        operation_function(the_db_manager, a_database_operation.buffered_data)

    @staticmethod
    def _default_operation():
        constants.PYSSA_LOGGER.warning("This operation does not exists!")

    @staticmethod
    def __wrapper_insert_new_protein(the_db_manager, the_buffered_data: tuple):
        _, tmp_protein = the_buffered_data
        the_db_manager.insert_new_protein(tmp_protein)

    @staticmethod
    def __wrapper_delete_existing_protein(the_db_manager, the_buffered_data: tuple):
        _, tmp_protein_id = the_buffered_data
        the_db_manager.delete_existing_protein(tmp_protein_id)

    def stop(self):
        # Stop the thread
        self._stop_event.set()
