# queries.py
"""
Thin, stateless helpers for building and running QSqlQuery objects.

These functions do *not* manage connections – callers obtain a connection
from ProjectConnectionPool.connection() and pass the db handle in.
"""
from __future__ import annotations

import logging
from typing import Any

from src.pyssa.gui.qt import QtSql

logger = logging.getLogger(__name__)


class QueryError(RuntimeError):
    """Raised when a SQL statement cannot be prepared or executed."""


def prepare(db: QtSql.QSqlDatabase, sql: str) -> QtSql.QSqlQuery:
    """Prepare *sql* against *db* and return the query object.

    Raises:
        QueryError: If the statement cannot be prepared.
    """
    if not db.isOpen():
        raise QueryError("Database is not open.")
    q = QtSql.QSqlQuery(db)
    if not q.prepare(sql):
        raise QueryError(f"Prepare failed: {q.lastError().text()}\nSQL: {sql}")
    return q


def execute(q: QtSql.QSqlQuery, *params: Any, max_retries: int = 3) -> QtSql.QSqlQuery:
    """Bind positional *params* and execute *q*.

    Retries on SQLITE_BUSY errors (database locked).

    Args:
        q: The prepared query.
        *params: Positional parameters to bind.
        max_retries: Maximum number of retry attempts for SQLITE_BUSY.

    Returns:
        The executed query so callers can iterate results immediately.

    Raises:
        QueryError: If execution fails after retries.
    """
    import time

    for i, value in enumerate(params):
        q.bindValue(i, value)

    last_error = None
    for attempt in range(max_retries + 1):
        if q.exec():
            return q

        last_error = q.lastError().text()
        error_number = q.lastError().nativeErrorCode()

        # SQLite error code 5 is SQLITE_BUSY (database is locked)
        if error_number == "5" and attempt < max_retries:
            retry_delay = 0.1 * (2 ** attempt)  # Exponential backoff
            logger.warning(
                "Query execution failed with SQLITE_BUSY. Retry %d/%d in %.2fs",
                attempt + 1, max_retries, retry_delay
            )
            time.sleep(retry_delay)
        else:
            break

    raise QueryError(f"Exec failed after {max_retries + 1} attempts: {last_error}")


def run(db: QtSql.QSqlDatabase, sql: str, *params: Any) -> QtSql.QSqlQuery:
    """Prepare and execute *sql* in one call. Convenience wrapper."""
    return execute(prepare(db, sql), *params)


def scalar(db: QtSql.QSqlDatabase, sql: str, *params: Any) -> Any:
    """Execute *sql* and return the single value in column 0 of the first row."""
    q = run(db, sql, *params)
    if q.next():
        return q.value(0)
    return None


def rows(db: QtSql.QSqlDatabase, sql: str, *params: Any) -> list[tuple]:
    """Execute *sql* and return all result rows as a list of tuples."""
    q = run(db, sql, *params)
    result = []
    while q.next():
        result.append(tuple(q.value(i) for i in range(q.record().count())))
    return result


def last_insert_id(db: QtSql.QSqlDatabase) -> int:
    """Return the rowid of the most recently inserted row."""
    value = scalar(db, "SELECT last_insert_rowid()")
    return int(value) if value is not None else -1
