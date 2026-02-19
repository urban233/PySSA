"""Defines the command protocol between host and PyMOL worker processes."""
from __future__ import annotations

import enum
from dataclasses import dataclass, field
from typing import Any


class CommandType(str, enum.Enum):
    """Strongly-typed core command identifiers sent to the worker process.

    Attributes:
        SHUTDOWN: Instructs the worker to stop its event loop and exit cleanly.
        LOAD_SESSION: Instructs the worker to load a ``.pse`` session file.
        DO: Instructs the worker to execute a PyMOL command string.
    """

    SHUTDOWN = "shutdown"
    LOAD_SESSION = "load_session"
    DO = "do"


@dataclass(frozen=True, slots=True)
class WorkerCommand:
    """Immutable value object representing a single instruction sent to a PyMOL worker.

    This class is the sole data contract between the host process and the
    worker subprocess. It is intentionally frozen to prevent accidental
    mutation after dispatch and uses ``__slots__`` for minimal memory overhead.

    Prefer the factory class-methods (shutdown(), load_session(), do()) over
    constructing instances directly.

    Attributes:
        command_type: The high-level operation the worker should perform.
        session_path: Absolute path to the ``.pse`` file the worker must load
            before executing the command. None is only valid for SHUTDOWN.
        pymol_command: The raw PyMOL command string (e.g. ``"zoom"``). Only
            meaningful for DO commands.
        args: Positional arguments forwarded to the PyMOL command.
        sync: When True the host blocks and waits for an acknowledgement reply.

    Example:
        Using the factory methods::

            cmd = WorkerCommand.do("select", args=("sele", "resi 50"), sync=True)
            cmd = WorkerCommand.shutdown()
    """

    command_type: CommandType
    session_path: str | None = None
    pymol_command: str | None = None
    args: tuple[Any, ...] = field(default_factory=tuple)
    sync: bool = False

    # ------------------------------------------------------------------
    # Factories
    # ------------------------------------------------------------------

    @classmethod
    def shutdown(cls) -> WorkerCommand:
        """Creates a shutdown command that instructs the worker to exit.

        Returns:
            A WorkerCommand with command_type SHUTDOWN.
        """
        return cls(command_type=CommandType.SHUTDOWN, sync=False)

    @classmethod
    def load_session(cls, session_path: str, *, sync: bool = False) -> WorkerCommand:
        """Creates a command that instructs the worker to load a session file.

        Args:
            session_path: Absolute path to the ``.pse`` file to load.
            sync: When True, the host blocks until the worker confirms the
                session has been loaded.

        Returns:
            A WorkerCommand with command_type LOAD_SESSION.
        """
        return cls(
            command_type=CommandType.LOAD_SESSION,
            session_path=session_path,
            sync=sync,
        )

    @classmethod
    def do(
            cls,
            session_path: str,
            pymol_command: str,
            args: tuple[Any, ...] = (),
            *,
            sync: bool = False,
    ) -> WorkerCommand:
        """Creates a command that instructs the worker to execute a PyMOL command.

        Args:
            session_path: Absolute path to the ``.pse`` file. The worker will
                load it if it differs from the currently loaded session.
            pymol_command: A valid PyMOL command string, e.g. ``"zoom"`` or
                ``"select"``.
            args: Positional arguments appended to the command string.
            sync: When True, the host blocks until the worker sends back an
                acknowledgement.

        Returns:
            A WorkerCommand with command_type DO.
        """
        return cls(
            command_type=CommandType.DO,
            session_path=session_path,
            pymol_command=pymol_command,
            args=args,
            sync=sync,
        )

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def build_pymol_string(self) -> str:
        """Constructs the full PyMOL command string including its arguments.

        The resulting string is ready to be passed directly to ``cmd.do()``.

        Returns:
            The fully assembled PyMOL command string, e.g. ``"select sele, resi 50"``.

        Raises:
            ValueError: If ``pymol_command`` is None.
        """
        if self.pymol_command is None:
            raise ValueError("pymol_command must be set to build a PyMOL string.")
        if not self.args:
            return self.pymol_command
        args_str = ", ".join(str(a) for a in self.args)
        return f"{self.pymol_command} {args_str}"
