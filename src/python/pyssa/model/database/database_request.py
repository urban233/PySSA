import logging

from pyssa.model.database import sql_query
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)


class DatabaseRequest:
  """Defines the request used by the database manager."""

  def __init__(self, a_list_of_queries: list["sql_query.SQLQuery"]) -> None:
    """Constructor.

    Args:
      a_list_of_queries: A list of SQL queries.

    Raises:
      exception.NoneValueError: If `a_list_of_queries` is None.
      exception.IllegalArgumentError: If `a_list_of_queries` is a empty list.
    """
    # <editor-fold desc="Checks">
    if a_list_of_queries is None:
      default_logging.append_to_log_file(logger, "a_list_of_queries is None.", logging.ERROR)
      raise exception.NoneValueError("a_list_of_queries is None.")
    if a_list_of_queries == "":
      default_logging.append_to_log_file(logger, "a_list_of_queries is a empty list.", logging.ERROR)
      raise exception.IllegalArgumentError("a_list_of_queries is a empty list.")
    # </editor-fold>
    # <editor-fold desc="Instance attributes">
    self.queries: list["sql_query.SQLQuery"] = a_list_of_queries
    """List of queries that will be executed in order"""
    self.await_results: bool = False
    """Flag to indicate if the request should be awaited or not"""
    # </editor-fold>
