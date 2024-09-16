import subprocess
from abc import ABC, abstractmethod


class IMvClient(ABC):
  """Interface that unifies the communication between PyDD and any molecular visualization software."""

  process: subprocess.Popen[bytes]
  """The name of the active PyMOL session."""

  @abstractmethod
  def start_viewer(self):
    """Starts the molecular viewer application."""
    raise NotImplementedError()

  @abstractmethod
  def connect_to_viewer(self):
    """Connects PyDD the molecular viewer application."""
    raise NotImplementedError()

  @abstractmethod
  def stop_viewer(self):
    """Stops the molecular viewer application."""
    raise NotImplementedError()
