import time

import zmq

from PyQt6 import QtWidgets

from sandbox.pymol_adapter import command_exec
from sandbox.gui.forms.auto import auto_simple_test


class SimpleFrame(QtWidgets.QMainWindow):
  """A simple frame for testing purposes and prototyping."""
  def __init__(self, a_process):
    """Constructor."""
    super().__init__()
    # build ui object
    self.ui = auto_simple_test.Ui_MainWindow()
    self.ui.setupUi(self)
    self.chimerax_process = a_process
    self.connect_signals()
    time.sleep(5)
    self.connect_to_chimerax()

  # <editor-fold desc="Util methods">
  def closeEvent(self, event) -> None:  # noqa: ANN001
    """Overrides the closeEvent of the QMainWindow class.

    Args:
        event: The event object representing the close event.
    """
    self.chimerax_process.terminate()

  def connect_signals(self):
    self.ui.btn_run.clicked.connect(self.run_action)
  # </editor-fold>

  def connect_to_chimerax(self):
    context = zmq.Context()
    self._sender_socket = context.socket(zmq.PUSH)
    self._sender_socket.connect("tcp://127.0.0.1:9071")
    self._recv_socket = context.socket(zmq.PULL)
    self._recv_socket.connect("tcp://127.0.0.1:9072")

  def run_action(self):
    print("Running single command")
    switch = False
    for i in range(10000):
      time.sleep(0.05)
      if switch:
        command_exec.PyMOLCommandExec().run_single_command(
          self._sender_socket, "color", "red"
        )
        switch = False
      else:
        command_exec.PyMOLCommandExec().run_single_command(
          self._sender_socket, "color", "blue"
        )
        switch = True
      print(f"The {i}th command finished.")
    print("Finished")
