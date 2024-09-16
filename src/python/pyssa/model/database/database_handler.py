import logging
import pathlib
import pickle
import queue
import time
from typing import Optional
from tea.concurrent import task
from tea.concurrent import action
from pyssa.model.central_objects import project
from pyssa.model.central_objects import sequence
from pyssa.model.central_objects import protein
from pyssa.model.central_objects import ligand
from pyssa.model.central_objects import protein_ligand_complex
from pyssa.model.central_objects import protein_protein_complex
from pyssa.model.central_objects.dummy import dummy_protein
from pyssa.model.central_objects.dummy import dummy_protein_ligand_complex
from pyssa.model.central_objects.dummy import dummy_protein_protein_complex
from pyssa.model.database import unified_database_access
from pyssa.model.database import common_query_statements
from pyssa.model.database import sql_query
from pyssa.model.database import database_request
from pyssa.model.pyssa_logging import default_logging
from pyssa.model.util import exception

logger = default_logging.setup_logger(__file__)


class DatabaseHandler:
  """Class for handling all database operations."""

  def __init__(self, a_filepath: pathlib.Path) -> None:
    """Constructor.

    Args:
      a_filepath: A filepath to the .db file

    Raises:
      exception.NoneValueError: If `a_filepath` is None.
      exception.IllegalArgumentError: If `a_filepath` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_filepath is None:
      default_logging.append_to_log_file(logger, "a_filepath is None.", logging.ERROR)
      raise exception.NoneValueError("a_filepath is None.")
    if a_filepath == "":
      default_logging.append_to_log_file(logger, "a_filepath is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_filepath is an empty string.")
    # </editor-fold>
    # self._database_access = unified_database_access.UnifiedDatabaseAccess(
    #   str(a_filepath), a_filepath.name.replace(".db", "")
    # )
    # <editor-fold desc="Instance attributes">
    self._database_path = a_filepath
    """The filepath to the database file."""
    self._push_queue = queue.Queue()
    """A queue that is used to store database requests that can run without feedback."""
    self._await_request_queue = queue.Queue()
    """A queue that gets database requests that have priority and must be awaited."""
    self._await_request_results_queue = queue.Queue()
    """A queue that has the results of the database requests that have priority and must be awaited."""
    self._stop_work: bool = False
    """A flag to stop the working method from outside."""
    self.is_working: bool = False
    """A flag indicating whether the handler is working."""
    self.is_current_project: bool = False
    """A flag indicating if the handler, handles the currently opened project."""
    # </editor-fold>

  # <editor-fold desc="Database work setup">
  # def __async_do_work(self) -> None:
  #   """Does the work which are in the queues.
  #
  #   Notes:
  #     This method must run in a separate thread!
  #   """
  #   self._database_access_async = unified_database_access.UnifiedDatabaseAccess(
  #     str(self._database_path), "AsyncConnection"
  #   )
  #   self._stop_work = False
  #   while self._stop_work is False:
  #     self.is_working = True
  #     if self._read_queue.empty() is False:
  #       # <editor-fold desc="Database operation: Reads & direct writes">
  #       while self._read_queue.empty() is False:
  #         try:
  #           # Reading operation
  #           tmp_sql_query: "sql_query.SQLQuery" = self._run_async_single_query(self._read_queue)
  #           self._read_results_queue.put(
  #             sql_query.SQLQueryResult(tmp_sql_query.get_all_results())
  #           )
  #           # Writing operation
  #           if self._direct_write_queue.empty() is False:
  #             self._run_async_single_query(self._direct_write_queue)
  #         except queue.Empty:
  #           default_logging.append_to_log_file(logger,
  #                                              "The read operation queue is empty, which shouldn't be possible.",
  #                                              logging.WARNING)
  #         # Writing operation
  #         while self._direct_write_queue.empty() is False:
  #           self._run_async_single_query(self._direct_write_queue)
  #       # </editor-fold>
  #     elif self._direct_write_queue.empty() is False:
  #       # <editor-fold desc="Database operation: Only direct writes">
  #       while self._direct_write_queue.empty() is False:
  #         self._run_async_single_query(self._direct_write_queue)
  #       # </editor-fold>
  #     else:
  #       time.sleep(1)
  #     if self._write_queue.empty() is False:
  #       # <editor-fold desc="Database operation: Writes">
  #       try:
  #         self._run_async_single_query(self._write_queue)
  #       except queue.Empty:
  #         default_logging.append_to_log_file(logger, "The write operation queue is empty, which shouldn't be possible.", logging.WARNING)
  #       # </editor-fold>
  #     else:
  #       time.sleep(1)
  #     default_logging.append_to_log_file(logger, "Database: Working ...")
  #   self.is_working = False

  def __async_do_work(self) -> None:
    """Does the work which are in the queues.

    Notes:
      This method must run in a separate thread!
    """
    self._database_access_async = unified_database_access.UnifiedDatabaseAccess(
      str(self._database_path), f"AsyncConnectionOf{self._database_path.name.replace('.db', '')}"
    )
    self._stop_work = False
    while self._stop_work is False:
      self.is_working = True
      if self._database_access_async.connect():
        if self._await_request_queue.empty() is False:
          # <editor-fold desc="Database operation: Await">
          while self._await_request_queue.empty() is False:
            try:
              tmp_database_request: "database_request.DatabaseRequest" = self._await_request_queue.get()
              for tmp_request_query in tmp_database_request.queries:
                tmp_sql_query: "sql_query.SQLQuery" = self._run_async_single_query(tmp_request_query)
                tmp_results = sql_query.SQLQueryResult(tmp_sql_query.get_all_results())
                self._await_request_results_queue.put(tmp_results)
            except queue.Empty:
              default_logging.append_to_log_file(logger,
                                                 "The await operation queue is empty, which should be impossible here.",
                                                 logging.WARNING)
          # </editor-fold>
        else:
          time.sleep(1)
        if self._push_queue.empty() is False:
          # <editor-fold desc="Database operation: Async">
          try:
            tmp_database_request: "database_request.DatabaseRequest" = self._push_queue.get()
            for tmp_request_query in tmp_database_request.queries:
              self._run_async_single_query(tmp_request_query)
          except queue.Empty:
            default_logging.append_to_log_file(logger, "The async operation queue is empty, which should be impossible here.", logging.WARNING)
          # </editor-fold>
        else:
          time.sleep(1)
        default_logging.append_to_log_file(logger, f"Database handler {self._database_path.name.replace('.db', '')}: Working ...")
      if self.is_current_project is False and self._push_queue.empty() is True:
        default_logging.append_to_log_file(
          logger,
          f"Database handler {self._database_path.name.replace('.db', '')} does not belong to the currently opened project and has no insert operations left. "
          "Therefore the handler will shutdown.",
          logging.INFO
        )
        self._database_access_async.disconnect()
        self._stop_work = True
    self.is_working = False

  def _run_async_single_query(self, a_sql_query: "sql_query.SQLQuery") -> "sql_query.SQLQuery":
    """Runs a single query asynchronously from a given queue.

    Args:
      a_sql_query: The query to run

    Returns:
      A SQLQuery object with the result of the query

    Raises:
      exception.NoneValueError: If `a_queue` is None.
    """
    # <editor-fold desc="Checks">
    if a_sql_query is None:
      default_logging.append_to_log_file(logger, "a_sql_query is None.", logging.ERROR)
      raise exception.NoneValueError("a_sql_query is None.")
    # </editor-fold>
    a_sql_query.db = self._database_access_async.db
    a_sql_query.prepare_qsql_query(a_sql_query.raw_query_string)
    self._database_access_async.execute_query(
      a_sql_query.query, params=a_sql_query.query_data
    )
    return a_sql_query

  def create_do_work_task(self) -> task.Task:
    """Creates a task that runs the __async_do_work method."""
    return task.Task().from_action(
      action.Action(
        self.__async_do_work
      )
    )

  def stop_work(self) -> None:
    """Stops the do_work method."""
    self._stop_work = True

  def put_to_push_query(self, a_database_request: "database_request.DatabaseRequest") -> None:
    """Puts data to the write queue."""
    self._push_queue.put(a_database_request)

  def put_to_await_request_query(self, a_database_request: "database_request.DatabaseRequest") -> None:
    """Puts data to the read queue."""
    self._await_request_queue.put(a_database_request)

  def get_await_request_query_result(self) -> "sql_query.SQLQueryResult":
    """Gets results from the read queue."""
    return self._await_request_results_queue.get()

  def build_new_database_and_change_to_it(self, a_database_filepath: pathlib.Path) -> task.Task:
    """Builds a new database and changes teh current active database to it.

    Args:
      a_database_filepath: The filepath to the database to change to

    Raises:
      exception.NoneValueError: If `a_database_filepath` is None.
    """
    # <editor-fold desc="Checks">
    if a_database_filepath is None:
      default_logging.append_to_log_file(logger, "a_database_filepath is None.", logging.ERROR)
      raise exception.NoneValueError("a_database_filepath is None.")
    # </editor-fold>
    # Stop work of old database
    self.stop_work()
    # Set filepath to new database
    while self.is_working is True:
      time.sleep(1)
    self.set_database_filepath(a_database_filepath)
    self.build_new_database(a_database_filepath.name.replace(".db", ""))
    # Create work loop task for new database
    return self.create_do_work_task()

  # </editor-fold>

  def set_database_filepath(self, a_filepath: pathlib.Path) -> None:
    """Sets the filepath of the database file.

    Args:
      a_filepath: A filepath to the .db file

    Raises:
      exception.NoneValueError: If `a_filepath` is None.
      exception.IllegalArgumentError: If `a_filepath` is an empty string.

    Notes:
      This will create a new database access object which overrides the old one!
    """
    # <editor-fold desc="Checks">
    if a_filepath is None:
      default_logging.append_to_log_file(logger, "a_filepath is None.", logging.ERROR)
      raise exception.NoneValueError("a_filepath is None.")
    if a_filepath == "":
      default_logging.append_to_log_file(logger, "a_filepath is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_filepath is an empty string.")
    # </editor-fold>
    # self._database_access = unified_database_access.UnifiedDatabaseAccess(
    #   str(a_filepath), a_filepath.name.replace(".db", "")
    # )
    self._database_path = a_filepath

  # def build_new_database(self) -> None:
  #   """Builds a new database for a new project."""
  #   if self._database_access.connect():
  #     tmp_queries = [
  #       sql_query.SQLQuery.from_query_with_database(
  #         self._database_access.db,
  #         common_query_statements.CommonQueryStatements.CREATE_TABLE_PROJECT
  #       ),
  #       sql_query.SQLQuery.from_query_with_database(
  #         self._database_access.db,
  #         common_query_statements.CommonQueryStatements.CREATE_TABLE_SEQUENCE
  #       ),
  #       sql_query.SQLQuery.from_query_with_database(
  #         self._database_access.db,
  #         common_query_statements.CommonQueryStatements.CREATE_TABLE_PROTEIN
  #       ),
  #       sql_query.SQLQuery.from_query_with_database(
  #         self._database_access.db,
  #         common_query_statements.CommonQueryStatements.CREATE_TABLE_LIGAND
  #       ),
  #       sql_query.SQLQuery.from_query_with_database(
  #         self._database_access.db,
  #         common_query_statements.CommonQueryStatements.CREATE_TABLE_PROTEIN_LIGAND_COMPLEX
  #       ),
  #       sql_query.SQLQuery.from_query_with_database(
  #         self._database_access.db,
  #         common_query_statements.CommonQueryStatements.CREATE_TABLE_PROTEIN_PROTEIN_COMPLEX
  #       ),
  #       sql_query.SQLQuery.from_query_with_database(
  #         self._database_access.db,
  #         common_query_statements.CommonQueryStatements.CREATE_TABLE_DUMMY_PROTEIN
  #       ),
  #       sql_query.SQLQuery.from_query_with_database(
  #         self._database_access.db,
  #         common_query_statements.CommonQueryStatements.CREATE_TABLE_DUMMY_PROTEIN_LIGAND_COMPLEX
  #       ),
  #       sql_query.SQLQuery.from_query_with_database(
  #         self._database_access.db,
  #         common_query_statements.CommonQueryStatements.CREATE_TABLE_DUMMY_PROTEIN_PROTEIN_COMPLEX
  #       ),
  #       sql_query.SQLQuery.from_query_with_database(
  #         self._database_access.db,
  #         common_query_statements.CommonQueryStatements.CREATE_TABLE_JOB
  #       )
  #     ]
  #     for tmp_query in tmp_queries:
  #       self._database_access.execute_query(tmp_query.query, ())
  #
  #     # --- Locking the database
  #     self._database_access.disconnect()

  # <editor-fold desc="Database operations">
  def build_new_database(self, a_project_name: str) -> None:
    """Builds a new database for a new project.

    Args:
      a_project_name: The name of the project

    Raises:
      exception.NoneValueError: If `a_project_name` is None.
    """
    # <editor-fold desc="Checks">
    if a_project_name is None:
      default_logging.append_to_log_file(logger, "a_project_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_project_name is None.")
    # </editor-fold>
    self.put_to_await_request_query(
      database_request.DatabaseRequest(
        [
          sql_query.SQLQuery.from_query(
            common_query_statements.CommonQueryStatements.CREATE_TABLE_PROJECT
          ),
          sql_query.SQLQuery.from_query(
            common_query_statements.CommonQueryStatements.CREATE_TABLE_SEQUENCE
          ),
          sql_query.SQLQuery.from_query(
            common_query_statements.CommonQueryStatements.CREATE_TABLE_PROTEIN
          ),
          sql_query.SQLQuery.from_query(
            common_query_statements.CommonQueryStatements.CREATE_TABLE_LIGAND
          ),
          sql_query.SQLQuery.from_query(
            common_query_statements.CommonQueryStatements.CREATE_TABLE_PROTEIN_LIGAND_COMPLEX
          ),
          sql_query.SQLQuery.from_query(
            common_query_statements.CommonQueryStatements.CREATE_TABLE_PROTEIN_PROTEIN_COMPLEX
          ),
          sql_query.SQLQuery.from_query(
            common_query_statements.CommonQueryStatements.CREATE_TABLE_DUMMY_PROTEIN
          ),
          sql_query.SQLQuery.from_query(
            common_query_statements.CommonQueryStatements.CREATE_TABLE_DUMMY_PROTEIN_LIGAND_COMPLEX
          ),
          sql_query.SQLQuery.from_query(
            common_query_statements.CommonQueryStatements.CREATE_TABLE_DUMMY_PROTEIN_PROTEIN_COMPLEX
          ),
          sql_query.SQLQuery.from_query(
            common_query_statements.CommonQueryStatements.CREATE_TABLE_JOB
          ),
          sql_query.SQLQuery.from_query(
            a_query_statement=common_query_statements.CommonQueryStatements.INSERT_PROJECT,
            the_query_data=(a_project_name,)
          )
        ]
      )
    )

  # <editor-fold desc="Project table operations">
  def get_complete_project(self, a_name: str) -> "project.Project":
    """Creates a new project object with all information from the database.

    Args:
      a_name: The name of the project

    Returns:
      A project object

    Raises:
      exception.NoneValueError: If `a_name` is None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    # </editor-fold>
    default_logging.append_to_log_file(logger, "Getting complete project from database ...")
    tmp_project = project.Project(a_name)
    default_logging.append_to_log_file(logger, "Getting all sequences from database ...")
    tmp_project.sequences = self.awaits_get_all_sequences()
    default_logging.append_to_log_file(logger, "Getting all proteins from database ...")
    tmp_project.proteins = self.awaits_get_all_proteins()
    default_logging.append_to_log_file(logger, "Getting all ligands from database ...")
    tmp_project.ligands = self.get_all_ligands()
    default_logging.append_to_log_file(logger, "Getting all protein-ligand complexes from database ...")
    tmp_project.protein_ligand_complexes = self.get_all_protein_ligand_complexes(tmp_project.proteins,
                                                                                 tmp_project.ligands)
    default_logging.append_to_log_file(logger, "Getting all protein-protein complexes from database ...")
    tmp_project.protein_protein_complexes = self.get_all_protein_protein_complexes(tmp_project.proteins)
    # tmp_project.dummy_proteins = self.get_all_dummy_proteins()
    # tmp_project.dummy_protein_ligand_complexes = self.get_all_dummy_protein_ligand_complexs()
    # tmp_project.dummy_protein_protein_complexes = self.get_all_dummy_protein_protein_complexs()
    default_logging.append_to_log_file(
      logger,
      "Retrieved all necessary data from the database and now return the project."
    )
    return tmp_project
  # </editor-fold>

  # <editor-fold desc="Sequence table operations">
  # <editor-fold desc="Writing">
  def insert_new_sequence(
          self, a_sequence: "sequence.Sequence", direct_write: bool = False
  ) -> None:
    """Inserts a new sequence in the database.

    Args:
      a_sequence (sequence.Sequence): The sequence object to be inserted into the database.
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_sequence is None:
      logger.error("a_sequence is None.")
      raise exception.NoneValueError("a_sequence is None.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.INSERT_SEQUENCE,
      the_query_data=(
        a_sequence.name,
        a_sequence.get_chain_sequence_map_as_binary(),
        a_sequence.type,
        0
      )
    )
    if direct_write:
      self.put_to_await_request_query(
        database_request.DatabaseRequest([tmp_sql_query])
      )
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  def update_chain_sequence_map(self, a_name: str, a_chain_sequence_map: bytes, direct_write: bool = False) -> None:
    """Updates an existing Chain Sequence map object of an existing sequence.

    Args:
      a_name: The name of the sequence to update the SeqRecord of
      a_chain_sequence_map: The chain sequence map as binary object
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if a_chain_sequence_map is None:
      default_logging.append_to_log_file(logger, "a_chain_sequence_map is None.", logging.ERROR)
      raise exception.NoneValueError("a_chain_sequence_map is None.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.UPDATE_CHAIN_SEQUENCE_MAP,
      the_query_data=(
        a_chain_sequence_map,
        a_name
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  def delete_sequence_by_name(self, a_name: str, direct_write: bool = False) -> None:
    """Deletes a sequence by its name.

    Args:
      a_name: The name of the sequence to delete
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.DELETE_SEQUENCE,
      the_query_data=(
        a_name,
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  # </editor-fold>

  # <editor-fold desc="Reading">
  def awaits_get_sequence_by_name(self, a_name: str) -> Optional["sequence.Sequence"]:
    """Gets a sequence from the database by its name.

    Args:
      a_name: The name of the sequence to get

    Returns:
      A sequence object or None if the sequence could not be found.

    Raises:
      exception.NoneValueError: If `a_name` is None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.GET_SEQUENCE,
      the_query_data=(a_name,)
    )
    self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    tmp_query_results: sql_query.SQLQueryResult = self.get_await_request_query_result()
    tmp_sequence = sequence.Sequence.from_chain_sequence_map(
      tmp_query_results.results[0]["name"], pickle.loads(tmp_query_results.results[0]["chain_sequence_map"])
    )
    return tmp_sequence

  def awaits_get_all_sequences(self) -> list["sequence.Sequence"]:
    """Gets all sequences from the database.

    Returns:
      A list of all sequences or an empty list if no sequences are in the table.
    """
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.GET_ALL_SEQUENCES,
      the_query_data=(0,)
    )
    self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    tmp_query_results: sql_query.SQLQueryResult = self.get_await_request_query_result()
    tmp_sequences: list = []
    for tmp_result in tmp_query_results.results:
      tmp_sequences.append(
        sequence.Sequence.from_chain_sequence_map(
          tmp_result["name"], pickle.loads(tmp_result["chain_sequence_map"])
        )
      )
    return tmp_sequences

  # </editor-fold>
  # </editor-fold>

  # <editor-fold desc="Protein table operations">
  # <editor-fold desc="Writing">
  def insert_new_protein(self,
                         a_protein: "protein.Protein",
                         direct_write: bool = False) -> None:
    """Inserts a new protein in the database.

    Args:
      a_protein: The protein object to be inserted into the database.
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Returns:
      The ID of the inserted protein or -1 if operation failed.

    Raises:
      exception.NoneValueError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_protein is None:
      logger.error("a_protein is None.")
      raise exception.IllegalArgumentError("a_protein is None.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    try:
      tmp_sql_query = sql_query.SQLQuery.from_query(
        a_query_statement=common_query_statements.CommonQueryStatements.INSERT_PROTEIN,
        the_query_data=(
          a_protein.name,
          a_protein.get_structure_as_binary(),
          a_protein.get_chain_sequence_map_as_binary(),
          a_protein.get_session_as_binary(),
          0
        )
      )
      if direct_write:
        self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
      else:
        self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))
      default_logging.append_to_log_file(logger, f"The protein {a_protein} will be added to the database {self._database_path.name}.",
                                         logging.INFO)
    except Exception as e:
      default_logging.append_to_log_file(
        logger,
        f"An error occurred while adding the protein {a_protein} to the database {self._database_path.name}. The error message is {e.__str__()}.",
        logging.ERROR
      )

  def update_protein_structure(self, a_name: str, a_structure: bytes, direct_write: bool = False) -> None:
    """Updates the protein structure.

    Args:
      a_name: The name of the protein to update the structure for
      a_structure: The structure as binary data
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if a_structure is None:
      default_logging.append_to_log_file(logger, "a_structure is None.", logging.ERROR)
      raise exception.NoneValueError("a_structure is None.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.UPDATE_PROTEIN_STRUCTURE,
      the_query_data=(
        a_structure,
        a_name
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  def update_protein_chain_sequence_map(
          self,
          a_name: str,
          a_chain_sequence_map: bytes,
          direct_write: bool = False
  ) -> None:
    """Updates the chain sequence map of an existing protein.

    Args:
      a_name: The name of the protein to update the chain sequence map for
      a_chain_sequence_map: The chain sequence map as binary data
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if a_chain_sequence_map is None:
      default_logging.append_to_log_file(logger, "a_chain_sequence_map is None.", logging.ERROR)
      raise exception.NoneValueError("a_chain_sequence_map is None.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.UPDATE_PROTEIN_CHAIN_SEQ_MAP,
      the_query_data=(
        a_chain_sequence_map,
        a_name
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  def update_protein_session(self, a_name: str, a_session: bytes, direct_write: bool = False) -> None:
    """Updates the protein session of an existing protein.

    Args:
      a_name: The name of the protein to update the session for
      a_session: The session as binary data
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if a_session is None:
      default_logging.append_to_log_file(logger, "a_session is None.", logging.ERROR)
      raise exception.NoneValueError("a_session is None.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.UPDATE_PROTEIN_SESSION,
      the_query_data=(
        a_session,
        a_name
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  def delete_protein_by_name(self, a_name: str, direct_write: bool = False) -> None:
    """Deletes a protein by its name.

    Args:
      a_name: The name of the protein to delete
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.DELETE_PROTEIN,
      the_query_data=(
        a_name,
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))
  # </editor-fold>

  # <editor-fold desc="Reading">
  def awaits_get_protein_by_name(self, a_name: str) -> Optional["protein.Protein"]:
    """Gets a protein from the database by its name.

    Args:
      a_name: The name of the protein to get

    Returns:
      A protein or None if the protein could not be found.

    Raises:
      exception.NoneValueError: If `a_name` is None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.GET_PROTEIN,
      the_query_data=(a_name,)
    )
    self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    tmp_query_results: sql_query.SQLQueryResult = self.get_await_request_query_result()
    tmp_protein = protein.Protein(
      a_name, pickle.loads(tmp_query_results.results[0]["structure"])
    )
    tmp_protein.id = tmp_query_results.results[0]["id"]
    tmp_protein.chain_sequence_map = pickle.loads(tmp_query_results.results[0]["chain_sequence_map"])
    tmp_protein.session = tmp_query_results.results[0]["session"]
    return tmp_protein

  def awaits_get_all_proteins(self) -> list["protein.Protein"]:
    """Gets all proteins from the database.

    Returns:
      A list of all proteins from the database or an empty list if no protein could be found.
    """
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.GET_ALL_PROTEINS,
      the_query_data=(0,)
    )
    self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    tmp_query_results: sql_query.SQLQueryResult = self.get_await_request_query_result()
    tmp_proteins: list = []
    for tmp_result in tmp_query_results.results:
      tmp_protein = protein.Protein(
        tmp_result["name"], pickle.loads(tmp_result["structure"])
      )
      tmp_protein.id = tmp_result["id"]
      #tmp_protein.chain_sequence_map = pickle.loads(tmp_result["chain_sequence_map"])
      tmp_protein.session = tmp_result["session"]
      tmp_proteins.append(tmp_protein)
    return tmp_proteins

  # </editor-fold>
  # </editor-fold>

  # <editor-fold desc="Ligand table operations">
  # <editor-fold desc="Writing">
  def insert_new_ligand(self, a_ligand: "ligand.Ligand", direct_write: bool = False) -> None:
    """Inserts a new ligand in the database.

    Args:
      a_ligand: The ligand object to be inserted into the database.
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_ligand is None:
      logger.error("a_ligand is None.")
      raise exception.NoneValueError("a_ligand is None.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.INSERT_LIGAND,
      the_query_data=(
        a_ligand.name,
        a_ligand.get_structure_as_binary(),
        a_ligand.session,
        0
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  def update_ligand_structure(self, a_name: str, a_structure: bytes, direct_write: bool = False) -> None:
    """Updates the ligand structure of an existing ligand.

    Args:
      a_name: The name of the ligand to update the structure for
      a_structure: The structure as binary data
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if a_structure is None:
      default_logging.append_to_log_file(logger, "a_structure is None.", logging.ERROR)
      raise exception.NoneValueError("a_structure is None.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.UPDATE_LIGAND_STRUCTURE,
      the_query_data=(
        a_structure,
        a_name
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  def update_ligand_session(self, a_name: str, a_session: bytes, direct_write: bool = False) -> None:
    """Updates the ligand session of an existing ligand.

    Args:
      a_name: The name of the ligand to update the structure for
      a_session: The session as binary data
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if a_session is None:
      default_logging.append_to_log_file(logger, "a_session is None.", logging.ERROR)
      raise exception.NoneValueError("a_session is None.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.UPDATE_LIGAND_SESSION,
      the_query_data=(
        a_session,
        a_name
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  def delete_ligand_by_name(self, a_name: str, direct_write: bool = False) -> None:
    """Deletes a ligand by its name.

    Args:
      a_name: The name of the ligand to delete
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.DELETE_LIGAND,
      the_query_data=(
        a_name,
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  # </editor-fold>

  # <editor-fold desc="Reading">
  def get_ligand_by_name(self, a_name: str) -> Optional["ligand.Ligand"]:
    """Gets a ligand from the database by its name.

    Args:
      a_name: The name of the ligand to get

    Returns:
      A ligand object or None if the ligand could not be found

    Raises:
      exception.NoneValueError: If `a_name` is None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.GET_LIGAND,
      the_query_data=(a_name,)
    )
    self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    tmp_query_results: sql_query.SQLQueryResult = self.get_await_request_query_result()
    tmp_ligand = ligand.Ligand(pickle.loads(tmp_query_results.results[0]["structure"]))
    tmp_ligand.id = tmp_query_results.results[0]["id"]
    tmp_ligand.name = tmp_query_results.results[0]["name"]
    tmp_ligand.session = tmp_query_results.results[0]["session"]
    return tmp_ligand

  def get_all_ligands(self) -> list["ligand.Ligand"]:
    """Gets all ligands from the database.

    Returns:
      A list of all ligands or an empty list if no ligands could be found.
    """
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.GET_LIGANDS,
      the_query_data=(0,)
    )
    self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    tmp_query_results: sql_query.SQLQueryResult = self.get_await_request_query_result()
    tmp_ligands = []
    for tmp_result in tmp_query_results.results:
      tmp_ligand = ligand.Ligand(pickle.loads(tmp_result["structure"]))
      tmp_ligand.id = tmp_result["id"]
      tmp_ligand.name = tmp_result["name"]
      tmp_ligand.session = tmp_result["session"]
      tmp_ligands.append(tmp_ligand)
    return tmp_ligands
  # </editor-fold>
  # </editor-fold>

  # <editor-fold desc="Protein-ligand complex table operations">
  # <editor-fold desc="Writing">
  def insert_new_protein_ligand_complex(
          self,
          a_protein_ligand_complex: "protein_ligand_complex.ProteinLigandComplex",
          direct_write: bool = False
  ) -> None:
    """Inserts a new protein-ligand complex in the database.

    Args:
      a_protein_ligand_complex: The protein-ligand complex object to be inserted into the database.
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_protein_ligand_complex is None:
      default_logging.append_to_log_file(logger, "a_protein_ligand_complex is None.", logging.ERROR)
      raise exception.NoneValueError("a_protein_ligand_complex is None.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.INSERT_PROTEIN_LIGAND_COMPLEX,
      the_query_data=(
        a_protein_ligand_complex.name,
        a_protein_ligand_complex.protein.id,
        a_protein_ligand_complex.ligand.id,
        a_protein_ligand_complex.get_results_as_binary(),
        a_protein_ligand_complex.session,
        0
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  def update_protein_ligand_complex_session(self, a_name: str, a_session: bytes, direct_write: bool = False) -> None:
    """Updates the protein-ligand complex session of an existing protein-ligand complex.

    Args:
      a_name: The name of the protein-ligand complex to update the session for
      a_session: The session as binary data
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if a_session is None:
      default_logging.append_to_log_file(logger, "a_session is None.", logging.ERROR)
      raise exception.NoneValueError("a_session is None.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.UPDATE_PROTEIN_LIGAND_COMPLEX_SESSION,
      the_query_data=(
        a_session,
        a_name
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  def delete_protein_ligand_complex_by_name(self, a_name: str, direct_write: bool = False) -> None:
    """Deletes a protein-ligand complex by its name.

    Args:
      a_name: The name of the protein-ligand complex to delete
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.DELETE_PROTEIN_LIGAND_COMPLEX,
      the_query_data=(
        a_name,
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  # </editor-fold>

  # <editor-fold desc="Reading">
  def get_protein_ligand_complex_by_name(
          self,
          a_name: str,
          a_protein_list: list["protein.Protein"],
          a_ligand_list: list["ligand.Ligand"]
  ) -> Optional["protein_ligand_complex.ProteinLigandComplex"]:
    """Gets a ligand from the database by its name.

    Args:
      a_name: The name of the ligand to get
      a_protein_list: A list of proteins to use for the complex
      a_ligand_list: A list of ligands to use for the complex

    Returns:
      A protein-ligand complex object or None if the protein-ligand complex could not be found

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string or if one of the lists are empty.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if a_protein_list is None:
      default_logging.append_to_log_file(logger, "a_protein_list is None.", logging.ERROR)
      raise exception.NoneValueError("a_protein_list is None.")
    if a_protein_list == []:
      default_logging.append_to_log_file(logger, "a_protein_list is a empty list.", logging.ERROR)
      raise exception.IllegalArgumentError("a_protein_list is a empty list")
    if a_ligand_list is None:
      default_logging.append_to_log_file(logger, "a_ligand_list is None.", logging.ERROR)
      raise exception.NoneValueError("a_ligand_list is None.")
    if a_ligand_list == []:
      default_logging.append_to_log_file(logger, "a_ligand_list is a empty list.", logging.ERROR)
      raise exception.IllegalArgumentError("a_ligand_list is a empty list")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.GET_PROTEIN_LIGAND_COMPLEX,
      the_query_data=(a_name,)
    )
    self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    tmp_query_results: sql_query.SQLQueryResult = self.get_await_request_query_result()

    tmp_complex_protein: Optional["protein.Protein"] = None
    for tmp_protein in a_protein_list:
      if tmp_protein.id == tmp_query_results.results[0]["protein_id"]:
        tmp_complex_protein = tmp_protein
        break

    tmp_complex_ligand: Optional["ligand.Ligand"] = None
    for tmp_ligand in a_ligand_list:
      if tmp_ligand.id == tmp_query_results.results[0]["ligand_id"]:
        tmp_complex_ligand = tmp_ligand
        break

    tmp_protein_ligand_complex = protein_ligand_complex.ProteinLigandComplex(
      tmp_query_results.results[0]["name"],
      tmp_complex_protein,
      tmp_complex_ligand
    )
    tmp_protein_ligand_complex.results = pickle.loads(tmp_query_results.results[0]["results"])
    tmp_protein_ligand_complex.session = tmp_query_results.results[0]["session"]
    return tmp_protein_ligand_complex

  def get_all_protein_ligand_complexes(
          self,
          a_protein_list: list["protein.Protein"],
          a_ligand_list: list["ligand.Ligand"]
  ) -> list["protein_ligand_complex.ProteinLigandComplex"]:
    """Gets all protein_ligand_complexes from the database.

    Args:
      a_protein_list: A list of proteins to use for the complex
      a_ligand_list: A list of ligands to use for the complex

    Returns:
      A list of all protein_ligand_complexes or an empty list if no protein_ligand_complexes could be found.

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If one of the lists are empty.
    """
    # <editor-fold desc="Checks">
    if a_protein_list is None:
      default_logging.append_to_log_file(logger, "a_protein_list is None.", logging.ERROR)
      raise exception.NoneValueError("a_protein_list is None.")
    if a_protein_list == []:
      default_logging.append_to_log_file(logger, "a_protein_list is a empty list.", logging.WARNING)
      return []
    if a_ligand_list is None:
      default_logging.append_to_log_file(logger, "a_ligand_list is None.", logging.ERROR)
      raise exception.NoneValueError("a_ligand_list is None.")
    if a_ligand_list == []:
      default_logging.append_to_log_file(logger, "a_ligand_list is a empty list.", logging.WARNING)
      return []
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.GET_PROTEIN_LIGAND_COMPLEXES,
      the_query_data=(0,)
    )
    self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    tmp_query_results: sql_query.SQLQueryResult = self.get_await_request_query_result()

    tmp_protein_ligand_complexes = []
    for tmp_result in tmp_query_results.results:
      tmp_complex_protein: Optional["protein.Protein"] = None
      for tmp_protein in a_protein_list:
        if tmp_protein.id == tmp_result["protein_id"]:
          tmp_complex_protein = tmp_protein
          break

      tmp_complex_ligand: Optional["ligand.Ligand"] = None
      for tmp_ligand in a_ligand_list:
        if tmp_ligand.id == tmp_result["ligand_id"]:
          tmp_complex_ligand = tmp_ligand
          break

      tmp_protein_ligand_complex = protein_ligand_complex.ProteinLigandComplex(
        tmp_result["name"],
        tmp_complex_protein,
        tmp_complex_ligand
      )
      tmp_protein_ligand_complex.results = pickle.loads(tmp_result["results"])
      tmp_protein_ligand_complex.session = tmp_result["session"]
      tmp_protein_ligand_complexes.append(tmp_protein_ligand_complex)
    return tmp_protein_ligand_complexes

  # </editor-fold>
  # </editor-fold>

  # <editor-fold desc="Protein-protein complex table operations">
  # <editor-fold desc="Writing">
  def insert_new_protein_protein_complex(
          self,
          a_protein_protein_complex: "protein_protein_complex.ProteinProteinComplex",
          direct_write: bool = False
  ) -> None:
    """Inserts a new protein-protein complex in the database.

    Args:
      a_protein_protein_complex: The protein-protein complex object to be inserted into the database.
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_protein_protein_complex is None:
      default_logging.append_to_log_file(logger, "a_protein_protein_complex is None.", logging.ERROR)
      raise exception.NoneValueError("a_protein_protein_complex is None.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.INSERT_PROTEIN_PROTEIN_COMPLEX,
      the_query_data=(
        a_protein_protein_complex.name,
        a_protein_protein_complex.protein_1.id,
        a_protein_protein_complex.protein_2.id,
        a_protein_protein_complex.get_results_as_binary(),
        a_protein_protein_complex.session,
        0
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  def update_protein_protein_complex_session(self, a_name: str, a_session: bytes, direct_write: bool = False) -> None:
    """Updates the protein-protein complex session of an existing protein-protein complex.

    Args:
      a_name: The name of the protein-protein complex to update the session for
      a_session: The session as binary data
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if a_session is None:
      default_logging.append_to_log_file(logger, "a_session is None.", logging.ERROR)
      raise exception.NoneValueError("a_session is None.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.UPDATE_PROTEIN_PROTEIN_COMPLEX_SESSION,
      the_query_data=(
        a_session,
        a_name
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  def delete_protein_protein_complex_by_name(self, a_name: str, direct_write: bool = False) -> None:
    """Deletes a protein-protein complex by its name.

    Args:
      a_name: The name of the protein-protein complex to delete
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.DELETE_PROTEIN_PROTEIN_COMPLEX,
      the_query_data=(
        a_name,
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  # </editor-fold>

  # <editor-fold desc="Reading">
  def get_protein_protein_complex_by_name(
          self,
          a_name: str,
          a_protein_list: list["protein.Protein"]
  ) -> Optional["protein_protein_complex.ProteinProteinComplex"]:
    """Gets a protein from the database by its name.

    Args:
      a_name: The name of the protein to get
      a_protein_list: A list of proteins to use for the complex

    Returns:
      A protein-protein complex object or None if the protein-protein complex could not be found

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string or `a_protein_list` is an empty list.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if a_protein_list is None:
      default_logging.append_to_log_file(logger, "a_protein_list is None.", logging.ERROR)
      raise exception.NoneValueError("a_protein_list is None.")
    if a_protein_list == []:
      default_logging.append_to_log_file(logger, "a_protein_list is a empty list.", logging.ERROR)
      raise exception.IllegalArgumentError("a_protein_list is a empty list")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.GET_PROTEIN_PROTEIN_COMPLEX,
      the_query_data=(a_name,)
    )
    self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    tmp_query_results: sql_query.SQLQueryResult = self.get_await_request_query_result()

    tmp_complex_protein_1: Optional["protein.Protein"] = None
    for tmp_protein in a_protein_list:
      if tmp_protein.id == tmp_query_results.results[0]["protein_1_id"]:
        tmp_complex_protein_1 = tmp_protein
        break

    tmp_complex_protein_2: Optional["protein.Protein"] = None
    for tmp_protein in a_protein_list:
      if tmp_protein.id == tmp_query_results.results[0]["protein_2_id"]:
        tmp_complex_protein_2 = tmp_protein
        break

    tmp_protein_protein_complex = protein_protein_complex.ProteinProteinComplex(
      tmp_query_results.results[0]["name"],
      tmp_complex_protein_1,
      tmp_complex_protein_2
    )
    tmp_protein_protein_complex.results = pickle.loads(tmp_query_results.results[0]["results"])
    tmp_protein_protein_complex.session = tmp_query_results.results[0]["session"]
    return tmp_protein_protein_complex

  def get_all_protein_protein_complexes(self, a_protein_list: list["protein.Protein"]) -> list["protein_protein_complex.ProteinProteinComplex"]:
    """Gets all protein_protein_complexes from the database.

    Args:
      a_protein_list: A list of proteins to use for the complex

    Returns:
      A list of all protein_protein_complexes or an empty list if no protein_protein_complexes could be found.

    Raises:
      exception.NoneValueError: If `a_protein_list` is None.
      exception.IllegalArgumentError: If `a_protein_list` is an empty list.
    """
    # <editor-fold desc="Checks">
    if a_protein_list is None:
      default_logging.append_to_log_file(logger, "a_protein_list is None.", logging.ERROR)
      raise exception.NoneValueError("a_protein_list is None.")
    if a_protein_list == []:
      default_logging.append_to_log_file(logger, "a_protein_list is a empty list.", logging.ERROR)
      return []
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.GET_PROTEIN_PROTEIN_COMPLEXES,
      the_query_data=(0,)
    )
    self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    tmp_query_results: sql_query.SQLQueryResult = self.get_await_request_query_result()

    tmp_protein_protein_complexes = []
    for tmp_result in tmp_query_results.results:
      tmp_complex_protein_1: Optional["protein.Protein"] = None
      for tmp_protein in a_protein_list:
        if tmp_protein.id == tmp_result["protein_1_id"]:
          tmp_complex_protein_1 = tmp_protein
          break

      tmp_complex_protein_2: Optional["protein.Protein"] = None
      for tmp_protein in a_protein_list:
        if tmp_protein.id == tmp_result["protein_2_id"]:
          tmp_complex_protein_2 = tmp_protein
          break

      tmp_protein_protein_complex = protein_protein_complex.ProteinProteinComplex(
        tmp_result["name"],
        tmp_complex_protein_1,
        tmp_complex_protein_2
      )
      tmp_protein_protein_complex.results = pickle.loads(tmp_result["results"])
      tmp_protein_protein_complex.session = tmp_result["session"]
      tmp_protein_protein_complexes.append(tmp_protein_protein_complex)
    return tmp_protein_protein_complexes
  # </editor-fold>
  # </editor-fold>

  # <editor-fold desc="Dummy protein table operations">
  # <editor-fold desc="Writing">
  def insert_new_dummy_protein(self, a_dummy_protein: "dummy_protein.DummyProtein", direct_write: bool = False) -> None:
    """Inserts a new dummy_protein in the database.

    Args:
      a_dummy_protein (dummy_protein.DummyProtein): The dummy_protein object to be inserted into the database.
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_dummy_protein is None:
      logger.error("a_dummy_protein is None.")
      raise exception.NoneValueError("a_dummy_protein is None.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.INSERT_DUMMY_PROTEIN,
      the_query_data=(
        a_dummy_protein.name,
        a_dummy_protein.chain_sequence_map,
        0
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  def delete_dummy_protein_by_name(self, a_name: str, direct_write: bool = False) -> None:
    """Deletes a dummy_protein by its name.

    Args:
      a_name: The name of the dummy_protein to delete
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.DELETE_DUMMY_PROTEIN,
      the_query_data=(
        a_name,
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))
  # </editor-fold>

  # <editor-fold desc="Reading">
  def get_dummy_protein_by_name(self, a_name: str) -> Optional["dummy_protein.DummyProtein"]:
    """Gets a dummy_protein from the database by its name.

    Args:
      a_name: The name of the dummy_protein to get

    Returns:
      A dummy_protein object or None if the dummy_protein could not be found.

    Raises:
      exception.NoneValueError: If `a_name` is None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.GET_DUMMY_PROTEIN,
      the_query_data=(a_name,)
    )
    self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    tmp_query_results: sql_query.SQLQueryResult = self.get_await_request_query_result()
    return dummy_protein.DummyProtein(
      tmp_query_results.results[0]["name"],
      pickle.loads(tmp_query_results.results[0]["chain_sequence_map"])
    )

  def get_all_dummy_proteins(self) -> list["dummy_protein.DummyProtein"]:
    """Gets all dummy_proteins from the database.

    Returns:
      A list of all dummy_proteins or an empty list if no dummy_proteins are in the table.
    """
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.GET_DUMMY_PROTEINS,
      the_query_data=(0,)
    )
    self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    tmp_query_results: sql_query.SQLQueryResult = self.get_await_request_query_result()
    tmp_dummy_proteins = []
    for tmp_result in tmp_query_results.results:
      tmp_dummy_proteins.append(
        dummy_protein.DummyProtein(
          tmp_result["name"],
          pickle.loads(tmp_result["chain_sequence_map"])
        )
      )
    return tmp_dummy_proteins

  # </editor-fold>
  # </editor-fold>

  # <editor-fold desc="Dummy protein-ligand complex table operations">
  # <editor-fold desc="Writing">
  def insert_new_dummy_protein_ligand_complex(
          self,
          a_dummy_protein_ligand_complex: "dummy_protein_ligand_complex.DummyProteinLigandComplex",
          direct_write: bool = False
  ) -> None:
    """Inserts a new dummy_protein_ligand_complex in the database.

    Args:
      a_dummy_protein_ligand_complex: The dummy_protein_ligand_complex object to be inserted into the database.
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_dummy_protein_ligand_complex is None:
      logger.error("a_dummy_protein_ligand_complex is None.")
      raise exception.NoneValueError("a_dummy_protein_ligand_complex is None.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.INSERT_DUMMY_PROTEIN_LIGAND_COMPLEX,
      the_query_data=(
        a_dummy_protein_ligand_complex.name,
        0
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  def delete_dummy_protein_ligand_complex_by_name(self, a_name: str, direct_write: bool = False) -> None:
    """Deletes a dummy_protein_ligand_complex by its name.

    Args:
      a_name: The name of the dummy_protein_ligand_complex to delete
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.DELETE_DUMMY_PROTEIN_LIGAND_COMPLEX,
      the_query_data=(
        a_name,
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))
  # </editor-fold>

  # <editor-fold desc="Reading">
  def get_dummy_protein_ligand_complex_by_name(self, a_name: str) -> Optional["dummy_protein_ligand_complex.DummyProteinLigandComplex"]:
    """Gets a dummy_protein_ligand_complex from the database by its name.

    Args:
      a_name: The name of the dummy_protein_ligand_complex to get

    Returns:
      A dummy_protein_ligand_complex object or None if the dummy_protein_ligand_complex could not be found.

    Raises:
      exception.NoneValueError: If `a_name` is None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    # </editor-fold>
    raise NotImplementedError()

  def get_all_dummy_protein_ligand_complexs(self) -> list["dummy_protein_ligand_complex.DummyProteinLigandComplex"]:
    """Gets all dummy_protein_ligand_complexs from the database.

    Returns:
      A list of all dummy_protein_ligand_complexs or an empty list if no dummy_protein_ligand_complexs are in the table.
    """
    raise NotImplementedError()
  # </editor-fold>
  # </editor-fold>

  # <editor-fold desc="Dummy protein-ligand complex table operations">
  # <editor-fold desc="Writing">
  def insert_new_dummy_protein_protein_complex(
          self,
          a_dummy_protein_protein_complex: "dummy_protein_protein_complex.DummyProteinProteinComplex",
          direct_write: bool = False
  ) -> None:
    """Inserts a new dummy_protein_protein_complex in the database.

    Args:
      a_dummy_protein_protein_complex: The dummy_protein_protein_complex object to be inserted into the database.
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_dummy_protein_protein_complex is None:
      logger.error("a_dummy_protein_protein_complex is None.")
      raise exception.NoneValueError("a_dummy_protein_protein_complex is None.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.INSERT_DUMMY_PROTEIN_PROTEIN_COMPLEX,
      the_query_data=(
        a_dummy_protein_protein_complex.name,
        0
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  def delete_dummy_protein_protein_complex_by_name(self, a_name: str, direct_write: bool = False) -> None:
    """Deletes a dummy_protein_protein_complex by its name.

    Args:
      a_name: The name of the dummy_protein_protein_complex to delete
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.DELETE_DUMMY_PROTEIN_PROTEIN_COMPLEX,
      the_query_data=(
        a_name,
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))
  # </editor-fold>

  # <editor-fold desc="Reading">
  def get_dummy_protein_protein_complex_by_name(self, a_name: str) -> Optional[
    "dummy_protein_protein_complex.DummyProteinProteinComplex"]:
    """Gets a dummy_protein_protein_complex from the database by its name.

    Args:
      a_name: The name of the dummy_protein_protein_complex to get

    Returns:
      A dummy_protein_protein_complex object or None if the dummy_protein_protein_complex could not be found.

    Raises:
      exception.NoneValueError: If `a_name` is None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    # </editor-fold>
    NotImplementedError()

  def get_all_dummy_protein_protein_complexs(self) -> list["dummy_protein_protein_complex.DummyProteinProteinComplex"]:
    """Gets all dummy_protein_protein_complexs from the database.

    Returns:
      A list of all dummy_protein_protein_complexs or an empty list if no dummy_protein_protein_complexs are in the table.
    """
    raise NotImplementedError()
  # </editor-fold>
  # </editor-fold>

  # <editor-fold desc="Job table operations">
  # <editor-fold desc="Writing">
  def insert_new_job(self, a_job: "job.Job", direct_write: bool = False) -> None:
    """Inserts a new job in the database.

    Args:
      a_job: The job object to be inserted into the database.
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_job is None:
      logger.error("a_job is None.")
      raise exception.NoneValueError("a_job is None.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.INSERT_JOB,
      the_query_data=(
        a_job.name,
        a_job.type,
        a_job.status,
        a_job.dummy_id,
        a_job.dummy_type,
        0
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  def update_job_status(self, a_name: str, a_status: bytes, direct_write: bool = False) -> None:
    """Updates the job status of an existing job.

    Args:
      a_name: The name of the job to update the status for
      a_status: The status as binary data
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if a_status is None:
      default_logging.append_to_log_file(logger, "a_status is None.", logging.ERROR)
      raise exception.NoneValueError("a_status is None.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.UPDATE_JOB_STATUS,
      the_query_data=(
        a_status,
        a_name
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  def delete_job_by_name(self, a_name: str, direct_write: bool = False) -> None:
    """Deletes a job by its name.

    Args:
      a_name: The name of the job to delete
      direct_write: Flag which indicates whether the write operation should be done directly. (Default: False)

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    if direct_write is None:
      default_logging.append_to_log_file(logger, "direct_write is None.", logging.ERROR)
      raise exception.NoneValueError("direct_write is None.")
    # </editor-fold>
    tmp_sql_query = sql_query.SQLQuery.from_query(
      a_query_statement=common_query_statements.CommonQueryStatements.DELETE_JOB,
      the_query_data=(
        a_name,
      )
    )
    if direct_write:
      self.put_to_await_request_query(database_request.DatabaseRequest([tmp_sql_query]))
    else:
      self.put_to_push_query(database_request.DatabaseRequest([tmp_sql_query]))

  # </editor-fold>

  # <editor-fold desc="Reading">
  def get_job_by_name(self, a_name: str) -> Optional["job.Job"]:
    """Gets a job from the database by its name.

    Args:
      a_name: The name of the job to get

    Returns:
      A job object or None if the job could not be found.

    Raises:
      exception.NoneValueError: If `a_name` is None.
      exception.IllegalArgumentError: If `a_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_name is None:
      default_logging.append_to_log_file(logger, "a_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_name is None.")
    if a_name == "":
      default_logging.append_to_log_file(logger, "a_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_name is an empty string.")
    # </editor-fold>
    NotImplementedError()

  def get_all_jobs(self) -> list["job.Job"]:
    """Gets all jobs from the database.

    Returns:
      A list of all jobs or an empty list if no jobs are in the table.
    """
    raise NotImplementedError()
  # </editor-fold>
  # </editor-fold>
  # </editor-fold>
