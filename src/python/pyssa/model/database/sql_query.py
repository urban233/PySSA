import logging
from typing import TYPE_CHECKING, Optional, Union
from PyQt6 import QtSql
from pyssa.model.pyssa_logging import default_logging
from pyssa.model.util import exception

if TYPE_CHECKING:
  from pyssa.model.database import common_query_statements

logger = default_logging.setup_logger(__file__)


class SQLQuery:
  """Contains all SQL queries needed for the PyDD application in form of static functions."""

  def __init__(
          self,
          a_database: Optional[QtSql.QSqlDatabase] = None,
          a_raw_query_string: Optional[str] = None
  ) -> None:
    """Constructor.

    Args:
      a_database: A QSqlDatabase instance

    Notes:
      The query object must be created in the thread where it will be used!!!
    """
    self.db: QtSql.QSqlDatabase = a_database
    self.raw_query_string: Optional[str] = a_raw_query_string
    self.query: Optional[QtSql.QSqlQuery] = None
    self.query_data: Optional[tuple] = None

  def prepare_qsql_query(self, a_query: str) -> bool:
    """Prepares the QSQL query to be run on the database.

    Args:
        a_query: A string representation of the query

    Returns:
        True if the preparation was successful, Otherwise: false

    Raises:
        exception.NoneValueError: If `a_query` or `self.db` is None.
        exception.IllegalArgumentError: If `a_query` is an empty string.
        ConnectionError: If the database is not open.
    """
    # <editor-fold desc="Checks">
    if a_query is None:
      default_logging.append_to_log_file(logger, "a_query is None.", logging.ERROR)
      raise exception.NoneValueError("a_query is None.")
    if a_query == "":
      default_logging.append_to_log_file(logger, "a_query is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_query is an empty string.")
    if self.db is None:
      default_logging.append_to_log_file(logger, "self.db is None.", logging.ERROR)
      raise exception.NoneValueError("self.db is None.")
    # </editor-fold>
    if not self.db.isOpen():
      raise ConnectionError(
        "The database must be open to perform this operation!"
      )
    tmp_qsql_query = QtSql.QSqlQuery(self.db)
    if tmp_qsql_query.prepare(a_query):
      self.query = tmp_qsql_query
      return True
    return False

  def get_all_results(self) -> list:
    """Gets all results of the QSqlQuery object."""
    results = []
    while self.query.next():
      record = self.query.record()
      row = {}
      for i in range(record.count()):
        row[record.fieldName(i)] = self.query.value(i)
      results.append(row)
    return results

  @staticmethod
  def from_query(
          a_query_statement: str,
          the_query_data: tuple = ()
  ) -> "SQLQuery":
    """Alternative constructor.

    Args:
        a_query_statement: The SQL query statement.
        the_query_data: A tuple containing the data necessary for the query

    Returns:
        A SQLQuery object prepared with the SQL statement if the preparation was successful, Otherwise: None

    Raises:
        exception.NoneValueError: If any of the arguments are None.
        exception.IllegalArgumentError: If `a_query_statement` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_query_statement is None:
      default_logging.append_to_log_file(logger, "a_query_statement is None.", logging.ERROR)
      raise exception.NoneValueError("a_query_statement is None.")
    if a_query_statement == "":
      default_logging.append_to_log_file(logger, "a_query_statement is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_query_statement is an empty string.")
    if the_query_data is None:
      default_logging.append_to_log_file(logger, "the_query_data is None.", logging.ERROR)
      raise exception.NoneValueError("the_query_data is None.")
    # </editor-fold>
    tmp_sql_query: SQLQuery = SQLQuery(a_raw_query_string=a_query_statement)
    tmp_sql_query.query_data = the_query_data
    return tmp_sql_query

  @staticmethod
  def from_query_with_database(
          a_database: QtSql.QSqlDatabase,
          a_query_statement: str,
          the_query_data: tuple = ()
  ) -> "SQLQuery":
    """Alternative constructor.

    Args:
        a_database: The QSqlDatabase object representing the database connection.
        a_query_statement: The SQL query statement.
        the_query_data: A tuple containing the data necessary for the query

    Returns:
        A SQLQuery object prepared with the SQL statement if the preparation was successful, Otherwise: None

    Raises:
        exception.NoneValueError: If any of the arguments are None.
        exception.IllegalArgumentError: If `a_query_statement` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_database is None:
      default_logging.append_to_log_file(logger, "a_database is None.", logging.ERROR)
      raise exception.NoneValueError("a_database is None.")
    if a_query_statement is None:
      default_logging.append_to_log_file(logger, "a_query_statement is None.", logging.ERROR)
      raise exception.NoneValueError("a_query_statement is None.")
    if a_query_statement == "":
      default_logging.append_to_log_file(logger, "a_query_statement is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_query_statement is an empty string.")
    if the_query_data is None:
      default_logging.append_to_log_file(logger, "the_query_data is None.", logging.ERROR)
      raise exception.NoneValueError("the_query_data is None.")
    # </editor-fold>
    tmp_sql_query: SQLQuery = SQLQuery(a_database=a_database)
    tmp_sql_query.query_data = the_query_data
    if tmp_sql_query.prepare_qsql_query(a_query_statement):
      default_logging.append_to_log_file(logger, "Preparing the QSQL query was successful.", logging.DEBUG)
    else:
      default_logging.append_to_log_file(logger, "QSQL query could not be prepared.", logging.WARNING)
    return tmp_sql_query


class SQLQueryResult:
  """Stores the results of an SQL query."""
  def __init__(self, results: list) -> None:
    """Constructor.

    Args:
      results: A list containing the query results

    Raises:
      exception.NoneValueError: If `results` is None.
    """
    # <editor-fold desc="Checks">
    if results is None:
      default_logging.append_to_log_file(logger, "results is None.", logging.ERROR)
      raise exception.NoneValueError("results is None.")
    # </editor-fold>
    self.results = results
