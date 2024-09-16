import sys
import threading
import time
import zmq
from typing import TYPE_CHECKING
from PyQt6.QtCore import pyqtSignal
from PyQt6 import QtCore
from chimerax.core import commands
from chimerax.core import session
from tea.concurrent import task_result, action, task_scheduler, task_manager

if TYPE_CHECKING:
  from . import adapter_frame


class AdapterFrameController(QtCore.QObject):
  """Controller for managing the AdapterFrame (gui tool for ChimeraX)"""

  command_data = pyqtSignal(tuple)

  def __init__(self, the_tool_instance,
               a_view: "adapter_frame.AdapterFrame", a_session) -> None:
    """Constructor."""
    super().__init__()
    self.tool_instance = the_tool_instance
    self._view = a_view
    self.session = a_session
    self.setup_connection()
    self.is_checking = True
    self.connect_signals()
    self.task_scheduler = task_scheduler.TaskScheduler()
    self.task_manager = task_manager.TaskManager()

  def setup_connection(self) -> None:
    """Sets up the connection to PyDD."""
    context = zmq.Context()
    self._recv_socket = context.socket(zmq.PULL)
    self._recv_socket.bind("tcp://127.0.0.1:9071")
    self._sender_socket = context.socket(zmq.PUSH)
    self._sender_socket.bind("tcp://127.0.0.1:9072")

  def connect_signals(self):
    """Connects all signals to their slots."""
    self._view.ui.btn_run.clicked.connect(self.run_command)
    self._view.ui.btn_close.clicked.connect(self.close)
    self.command_data.connect(self.run_chimerax_command)

  def run_command(self) -> None:
    """Runs a command triggered through the button in the bundle.

    ONLY use for debug purposes!!!

    """
    data = self._recv_socket.recv_json()
    commands.run(self.session, f"{data['command']} {data['arguments']}")

  def close(self):
    self.tool_instance.display(False)

  def start_checking(self) -> None:
    """Starts the checking routine in a separate thread."""
    tmp_task_result = task_result.TaskResult.from_action(
      action.Action(self.handle_requests),
      self.__await_handle_requests
    )
    self.task_manager.append_task_result(tmp_task_result)
    self.task_scheduler.schedule(tmp_task_result)

  def handle_requests(self) -> None:
    """Handles all requests that are sent from PyDD.

    Notes:
      Be aware that this method is only used to get the data from PyDD
      (command, arguments) and NOT for any functionality.
      EVERY api call from ChimeraX MUST be in the main thread!!! -> check
      run_chimerax_command() as an example
    """
    self.is_checking = True
    while self.is_checking is True:
      print("Checking requests ...")
      data = self._recv_socket.recv_json()
      if data["command"] == "Stop":
        print("Stopping checking requests...")
        self.is_checking = False
      else:
        self.command_data.emit((data["command"], data["arguments"]))
    print("Exiting while loop ...")

  def __await_handle_requests(self) -> None:
    """Awaits the handle requests method."""
    print("Thread closed for checking requests.")

  def run_chimerax_command(self, data: tuple):
    """Runs a chimerax command.

    Note:
      data => tuple of command (0) and arguments (1)
    """
    if threading.main_thread():
      commands.run(self.session, f"{data[0]} {data[1]}")
    else:
      raise PermissionError("You cannot use API commands from a separate thread!")
