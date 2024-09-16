import threading
import zmq
from PyQt5.QtCore import pyqtSignal
from PyQt5 import QtCore
from pymol import cmd


class Adapter(QtCore.QObject):
  command_data = pyqtSignal(tuple)

  def __init__(self):
    super().__init__()
    self.setup_connection()
    self.connect_signals()
    self.start_checking()

  def start_checking(self):
    self.thread = threading.Thread(target=self.handle_requests)
    self.thread.daemon = True
    self.thread.start()

  def setup_connection(self) -> None:
    """Sets up the connection to PyDD."""
    context = zmq.Context()
    self._recv_socket = context.socket(zmq.PULL)
    self._recv_socket.bind("tcp://127.0.0.1:9071")
    self._sender_socket = context.socket(zmq.PUSH)
    self._sender_socket.bind("tcp://127.0.0.1:9072")

  def connect_signals(self):
    """Connects all signals to their slots."""
    self.command_data.connect(self.run_pymol_command)

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

  def run_pymol_command(self, data: tuple) -> None:
    print(f"Coloring in {data[1]}")
    cmd.color(data[1])
    print(f"Finished coloring in {data[1]}")
    cmd.scene(key="new", action="append")
    print(f"Finished saving scene.")
