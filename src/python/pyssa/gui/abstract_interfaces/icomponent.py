from abc import ABC, abstractmethod
from typing import Callable, Optional


class IComponent(ABC):
  """Interface for all components."""

  @abstractmethod
  def run_tool(self, a_callable: Optional[Callable] = None):
    """Starts and runs the tool.

    Args:
      a_callable: The function to execute after the tool is closed. (Default: None)
    """
    raise NotImplementedError()
