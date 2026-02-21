# __init__.py
"""
pyssa_db – persistence layer for PySSA.

Public surface
--------------
ProjectDatabase     Read/write API for one project database.
                    Owns a write_queue (ProjectWriteQueue) for
                    fire-and-forget async writes via ThreadRuntime.

WriteOperation      A unit of deferred write work.
OperationType       Enum of all supported deferred operations.

DatabaseThread is intentionally absent from this package — it has been
replaced by ProjectWriteQueue, which uses the application's shared
ThreadRuntime rather than a dedicated thread per project.
"""

from .project_database import ProjectDatabase
from .write_queue import WriteOperation, OperationType
from .cold_project_handle import ColdProjectHandle

__all__ = [
    "ProjectDatabase",
    "WriteOperation",
    "OperationType",
    "ColdProjectHandle"
]
