import logging

from PyQt6 import QtCore
from PyQt6 import QtSql

from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)


class UnifiedDatabaseAccess:
  """Class for accessing the database in a unified way."""

  def __init__(self, a_database_path: str, a_connection_name: str) -> None:
    """Constructor.

    Args:
        a_database_path (str): The path to the SQLite database file.
        a_connection_name (str): The name of the database connection.

    Raises:
        exception.NoneValueError: If any of the arguments are None.
        exception.NoneValueError: If any of the arguments are an empty string.

    Notes:
      The connection name must be unique for the thread where it is used!!!
    """
    # <editor-fold desc="Checks">
    if a_database_path is None:
      default_logging.append_to_log_file(logger, "a_database_path is None.", logging.ERROR)
      raise exception.NoneValueError("a_database_path is None.")
    if a_database_path == "":
      default_logging.append_to_log_file(logger, "a_database_path is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_database_path is an empty string.")
    if a_connection_name is None:
      default_logging.append_to_log_file(logger, "a_connection_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_connection_name is None.")
    if a_connection_name == "":
      default_logging.append_to_log_file(logger, "a_connection_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_connection_name is an empty string.")
    # </editor-fold>
    self.mutex = QtCore.QMutex()
    self.db = QtSql.QSqlDatabase.addDatabase("QSQLITE", a_connection_name)
    self.db.setDatabaseName(a_database_path)

  def connect(self) -> bool:
    """Connects to the database.

    This method acquires a lock on the mutex, attempts to open the database connection,
    and returns a boolean value indicating whether the connection was successful.

    Returns:
        True if the connection was successful, False otherwise.
    """
    self.mutex.lock()
    try:
      if not self.db.open():
        # Handle connection error
        logger.error("Error connecting to database!")
        return False
      return True
    finally:
      self.mutex.unlock()

  def disconnect(self) -> None:
    """Disconnects from the database.

    This method closes the database connection and releases the associated mutex lock.
    """
    self.mutex.lock()
    try:
      self.db.close()
    finally:
      self.mutex.unlock()

  def execute_query(
          self,
          a_qsql_query: QtSql.QSqlQuery,
          params: tuple
  ) -> QtSql.QSqlQuery:
    """Executes a QSQL query with given parameters.

    Notes:
        Be aware this function has side effects!
        The given QSqlQuery will carry the results of the query after this function finished!

    Args:
        a_qsql_query (QtSql.QSqlQuery): The QSQL query to execute.
        params (tuple): The parameters to bind to the query.

    Returns:
        QtSql.QSqlQuery: The executed QSQL query.

    Raises:
        exception.NoneValueError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_qsql_query is None:
      default_logging.append_to_log_file(logger, "a_qsql_query is None.", logging.ERROR)
      raise exception.NoneValueError("a_qsql_query is None.")
    if params is None:
      default_logging.append_to_log_file(logger, "params is None.", logging.ERROR)
      raise exception.NoneValueError("params is None.")
    # </editor-fold>
    self.mutex.lock()
    try:
      i = 0
      for tmp_param in params:
        if isinstance(tmp_param, bytes):
          a_qsql_query.bindValue(i, QtCore.QByteArray(tmp_param), type=QtSql.QSql.ParamTypeFlag.In | QtSql.QSql.ParamTypeFlag.Binary)
        else:
          a_qsql_query.bindValue(i, tmp_param)
        i += 1
      if not a_qsql_query.exec():
        # Handle query execution error
        default_logging.append_to_log_file(
          logger,f"Error executing query: {a_qsql_query.lastError().text()}"
        )
    except Exception as e:
      logger.error(e)
    finally:
      self.mutex.unlock()
    return a_qsql_query
