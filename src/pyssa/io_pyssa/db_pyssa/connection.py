# connection.py
"""
Thread-safe SQLite connection management for PySSA.

Each ProjectDatabase owns one ProjectConnectionPool. A pool maintains one
QSqlDatabase connection *per thread* (keyed by thread-id), which is the
correct model for Qt's SQL driver: connections must not be shared across
threads, but within a thread a single connection is reused.
"""
from __future__ import annotations

import logging
import threading
from contextlib import contextmanager
from typing import Generator

from src.pyssa.gui.qt import QtSql

logger = logging.getLogger(__name__)


class ProjectConnectionPool:
    """Maintains one QSqlDatabase connection per calling thread.

    Args:
        db_path:  Absolute path to the SQLite file.
        pool_id:  A short, stable identifier for the project (used as the
                  QSqlDatabase connection-name prefix so names are unique
                  across multiple open projects).
    """

    def __init__(self, db_path: str, pool_id: str) -> None:
        self._db_path = db_path
        self._pool_id = pool_id
        self._lock = threading.Lock()
        self._connections: dict[int, QtSql.QSqlDatabase] = {}

    # ------------------------------------------------------------------
    # Public helpers
    # ------------------------------------------------------------------

    @contextmanager
    def connection(self, max_retries: int = 3, retry_delay: float = 0.1) -> Generator[QtSql.QSqlDatabase, None, None]:
        """Context manager that yields an open, thread-local QSqlDatabase.

        Retries on transient errors (SQLITE_BUSY, connection failures).

        Args:
            max_retries: Maximum number of retry attempts.
            retry_delay: Delay in seconds between retries.

        Example::

            with pool.connection() as db:
                query = QtSql.QSqlQuery(db)
                query.prepare("SELECT …")
                query.exec()

        Raises:
            ConnectionError: If connection cannot be opened after retries.
        """
        import time

        db = self._get_or_create()
        last_error = None

        for attempt in range(max_retries + 1):
            if not db.isOpen():
                if db.open():
                    logger.debug("Connection opened for thread %d", threading.get_ident())
                    break
                last_error = db.lastError().text()
                if attempt < max_retries:
                    logger.warning(
                        "[%s] Connection attempt %d/%d failed: %s. Retrying in %.1fs...",
                        self._pool_id, attempt + 1, max_retries + 1, last_error, retry_delay
                    )
                    time.sleep(retry_delay)
                    # Exponential backoff
                    retry_delay *= 2
                else:
                    raise ConnectionError(
                        f"[{self._pool_id}] Cannot open database after {max_retries + 1} attempts: {last_error}"
                    )
            else:
                break

        try:
            yield db
        finally:
            # Keep the connection alive for reuse – closing is deferred to
            # close_all(), called when the project is unloaded.
            pass

    def close_all(self) -> None:
        """Close every thread-local connection and remove it from Qt's registry."""
        with self._lock:
            for tid, db in list(self._connections.items()):
                conn_name = db.connectionName()
                db.close()
                QtSql.QSqlDatabase.removeDatabase(conn_name)
                logger.debug("Closed connection %s (thread %d)", conn_name, tid)
            self._connections.clear()

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    def _get_or_create(self) -> QtSql.QSqlDatabase:
        tid = threading.get_ident()
        with self._lock:
            if tid not in self._connections:
                conn_name = f"{self._pool_id}_t{tid}"
                db = QtSql.QSqlDatabase.addDatabase("QSQLITE", conn_name)
                db.setDatabaseName(self._db_path)
                self._connections[tid] = db
                logger.debug("Created connection %s", conn_name)
            return self._connections[tid]
