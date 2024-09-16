from abc import ABC, abstractmethod
from typing import Union, Optional, Callable
from PyQt6 import QtCore
from src.python.pydd.mv_wrappers.abstract_wrapper import imv_client


class ISessionManager(ABC):
  """Interface that unifies the session management in different molecular visualization software."""

  # <editor-fold desc="Class attributes">
  user_mv_client: "imv_client.IMvClient"
  """A client that can be used to connect to the molecular visualization software."""
  
  user_mv_client_lock: QtCore.QMutex
  """A mutex lock that allows us to lock the `mv_client`."""
  
  session_name: str
  """The name of the active PyMOL session."""

  session_object_type: str
  """The object type (protein or protein_pair) of the active PyMOL session."""

  session_objects: list
  """A list of all objects in the active PyMOL session."""

  current_scene_name: str
  """The name of the current scene."""

  all_scenes: list[str]
  """A list of all scenes in the active PyMOL session."""
  
  # </editor-fold>

  @abstractmethod
  def _check_session_integrity(self, a_protein_name: str) -> bool:
    raise NotImplementedError()
    
  @abstractmethod
  def cmd(self, a_command_name: "pymol_enums.CommandEnum", the_args: tuple = ()) -> dict:
    """Runs a command.

    Args:
      a_command_name (pymol_enums.CommandEnum): The name of the command to be executed.
      the_args (tuple): (Optional) The arguments to be passed to the command. Defaults to ().

    Returns:
      The result of the command execution as a dictionary.
    
    Raises:
      exception.IllegalArgumentError: If any of the arguments are None.
      exception.PyMOLCommandFailedError: If the command failed.
    """
    raise NotImplementedError()
  
  @abstractmethod
  def async_cmd(self,
                the_task_manager: "task_manager.TaskManager",
                the_task_scheduler: "task_scheduler.TaskScheduler",
                a_command_name: "pymol_enums.CommandEnum",
                the_args: tuple = (),
                an_await_function: Optional[Callable] = None) -> None:
    """Runs a single command in an asynchronous manner.

    The results will be given to the await function as argument.
    Therefore, the await function needs an argument to use the results.

    Args:
      the_task_manager (task_manager.TaskManager): The task manager instance of the interface manager.
      the_task_scheduler (task_scheduler.TaskScheduler): The task scheduler instance of the interface manager.
      a_command_name (pymol_enums.CommandEnum): The name of the command to be executed.
      the_args (tuple): (Optional) The arguments to be passed to the command. Defaults to ().
      an_await_function (Callable): An await function to call after the async function finished. Defaults to None.
    
    Raises:
      exception.IllegalArgumentError: If any of the arguments is None except `an_await_function`.
    """
    raise NotImplementedError()
  
  @abstractmethod
  def async_cmds(self,
                 the_task_manager: "task_manager.TaskManager",
                 the_task_scheduler: "task_scheduler.TaskScheduler",
                 the_command_names: tuple["pymol_enums.CommandEnum"],
                 the_args: tuple[tuple],
                 an_await_function: Optional[Callable] = None) -> None:
    """Runs multiple commands sequentially in an asynchronous manner.

    The results will be given to the await function as argument.
    Therefore, the await function needs an argument to use the results.

    Args:
      the_task_manager (task_manager.TaskManager): The task manager instance of the interface manager.
      the_task_scheduler (task_scheduler.TaskScheduler): The task scheduler instance of the interface manager.
      the_command_names (tuple[pymol_enums.CommandEnum]): The name of the command to be executed.
      the_args tuple[tuple]: The arguments to be passed to the command.
      an_await_function (Callable): An await function to call after the async function finished. Defaults to None.
    
    Raises:
      exception.IllegalArgumentError: If any of the arguments is None except `an_await_function` or `the_command_names` or `the_args` are an empty list.
    """
    raise NotImplementedError()
  
  def reinitialize_session(self) -> None:
    """Reinitialize the pymol session and class attributes."""
    # reset class attributes
    self.session_name = ""
    self.session_object_type = ""
    self.session_objects: list = []
    # reset actual pymol session
    self.user_mv_client.reinitialize_session()


class SessionManagerAdapter:
  """Class that functions as a general session manager which is independent of the used molecular visualization software."""
  
  def __init__(self,
               a_session_manager: Union[
                   "pymol_session_manager.PymolSessionManager",
                   "chimerax_session_manager.ChimeraXSessionManager"
               ]) -> None:
    self.session_manager = a_session_manager
