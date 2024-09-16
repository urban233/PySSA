import logging
import pathlib
from typing import Optional
from pyssa.model.preference import model_definitions

__docformat__ = "google"


def setup_logger(
        a_name: str,
        level: int = logging.DEBUG,
        log_path: pathlib.Path = model_definitions.ModelDefinitions.DEFAULT_LOG_PATH,
        add_console_handler: bool = True
) -> Optional[logging.Logger]:
  """Sets up a logger with a FileHandler directing output to the specified file.

  Args:
      a_name (str): Name for the logger.
      level (int, optional): Logging level (Default: logging.INFO).
      log_path (str): The path for the log file. (Default: The default PySSA log path.)
      add_console_handler (bool): A flag for whether to add a console handler. (Default: True).

  Returns:
      A logger object or None if any of the arguments are None or `a_name` is an empty string.
  """
  # <editor-fold desc="Checks">
  if a_name is None or a_name == "":
    return None
  if log_path is None or log_path == "":
    return None
  if add_console_handler is None:
    return None
  # </editor-fold>

  # <editor-fold desc="Definitions of important paths">
  if log_path == model_definitions.ModelDefinitions.DEFAULT_LOG_PATH:
    if not model_definitions.ModelDefinitions.DEFAULT_LOG_PATH.exists():
      model_definitions.ModelDefinitions.DEFAULT_LOG_PATH.mkdir(parents=True)
  else:
    if not log_path.exists():
      log_path.mkdir(parents=True)

  # </editor-fold>

  tmp_logger = logging.getLogger(a_name)
  tmp_logger.setLevel(level)

  # File handler for separate log files
  tmp_handler = logging.FileHandler(model_definitions.ModelDefinitions.LOG_FILEPATH)
  tmp_log_formatter = logging.Formatter(
    "%(asctime)s: %(name)s %(levelname)s - %(message)s"
  )
  tmp_handler.setFormatter(tmp_log_formatter)
  tmp_logger.addHandler(tmp_handler)

  if add_console_handler:
    # Console handler for logging to console
    tmp_console_handler = logging.StreamHandler()
    tmp_console_formatter = logging.Formatter(
      "%(levelname)s: %(message)s"
    )  # Simpler console format
    tmp_console_handler.setFormatter(tmp_console_formatter)
    tmp_logger.addHandler(tmp_console_handler)
  return tmp_logger


def append_to_log_file(a_logger: logging.Logger, a_message: str, a_level: int = logging.DEBUG) -> bool:
  """Appends a message to the log file.

  Args:
    a_logger: A logger instance
    a_message: A message to append
    a_level: A log level to use for the message (default: logging.DEBUG)

  Notes:
    If a console handler is added during the logger setup, then the message
    is also appended to the console log.
  """
  # <editor-fold desc="Checks">
  if a_logger is None:
    return False
  if a_message is None:
    return False
  if a_level is None:
    return False
  # TODO: Add check for the available log levels
  # </editor-fold>
  try:
    a_logger.log(a_level, a_message)
    return True
  except Exception as e:
    print(e)  # Logging is here not possible if an error occurs during the log process
    return False
