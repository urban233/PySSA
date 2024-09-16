import logging
import time

from pyssa.model.database import database_handler
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)


class DatabasesManager:
  """Manages all database handlers."""

  def __init__(self) -> None:
    """Constructor."""
    self.handlers: dict[str, "database_handler.DatabaseHandler"] = {}

  def add_handler(self, a_key, a_database_handler: "database_handler.DatabaseHandler") -> None:
    """Adds a database handler to the dict."""
    self.handlers[a_key] = a_database_handler
    default_logging.append_to_log_file(
      logger,
      f"Register new database handler: {a_key}. After registration these are the handlers: {self.handlers.keys()}",
      logging.INFO
    )

  def remove_handler(self, a_key) -> None:
    """Removes a database handler from the dict."""
    del self.handlers[a_key]

  def remove_unused_handlers(self) -> None:
    """Remove all database handlers from the dict that are no longer been used."""
    for a_key in self.handlers.keys():
      if self.handlers[a_key].is_current_project == False and self.handlers[a_key].is_working == False:
        default_logging.append_to_log_file(logger, f"Remove unused database handler {a_key}.", logging.INFO)
        del self.handlers[a_key]

  def get_handler(self, a_key) -> "database_handler.DatabaseHandler":
    """Gets a database handler from the dict."""
    return self.handlers[a_key]

  def stop_all_handlers(self) -> None:
    """Stops all handlers."""
    for a_key in self.handlers.keys():
      self.handlers[a_key].is_current_project = False
      while self.handlers[a_key].is_working is True:
        default_logging.append_to_log_file(logger, f"Database handler {a_key} is still working ...", logging.INFO)
        time.sleep(1)
