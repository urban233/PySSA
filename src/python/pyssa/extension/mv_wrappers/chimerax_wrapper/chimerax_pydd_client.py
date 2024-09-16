import time
from typing import Optional

import zmq
import subprocess
from src.python.pydd.mv_wrappers.abstract_wrapper import imv_client


class ChimeraxPyddClient(imv_client.IMvClient):
  """Class that implements the mv_client interface for communication with ChimeraX."""

  def __init__(self) -> None:
    """Constructor."""
    self.process: Optional[subprocess.Popen[bytes]] = None

  def start_viewer(self) -> None:
    """Starts ChimeraX in a separate process."""
    self.process = subprocess.Popen(r"C:\ProgramData\IBCI\PyDD\bin\ChimeraX\bin\ChimeraX.exe")
    if self.process.pid is None:
      time.sleep(2)  # Only for debug purposes

  def connect_to_viewer(self) -> None:
    """Connects to ChimeraX."""
    context = zmq.Context()
    self._sender_socket = context.socket(zmq.PUSH)
    self._sender_socket.connect("tcp://127.0.0.1:9071")
    self._recv_socket = context.socket(zmq.PULL)
    self._recv_socket.connect("tcp://127.0.0.1:9072")

  def stop_viewer(self) -> None:
    """Stops the checking requests routine in the ChimeraX bundle."""
    # TODO: This is only a working solution, the command should be wrapped in an enum for better handling
    self._sender_socket.send_json({
      "command": "Stop",
      "arguments": (),
    })
