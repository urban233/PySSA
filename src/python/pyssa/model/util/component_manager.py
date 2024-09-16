import logging
from typing import TYPE_CHECKING
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)

__docformat__ = "google"

if TYPE_CHECKING:
  from src.python.pydd.components.internal.abstract_interfaces import icomponent


class ComponentManager:
  """Class for managing all internal and external components.

  Note:
    Adding new components: To add new components it is necessary to go to the
      main frame controller class and add the new component to the dict in the
      component_manager.ComponentManager constructor. In addition, a new dict
      key must be defined in the model_definitions.ComponentsEnum!

  """

  def __init__(self, the_components: dict[int, "icomponent.IComponent"]) -> None:
    """Constructor.

    Args:
      the_components: A dictionary of all available components

    Raises:
      exception.NoneValueError: If `the_components` is None.

    """
    # <editor-fold desc="Checks">
    if the_components is None:
      default_logging.append_to_log_file(logger, "the_components is None.", logging.ERROR)
      raise exception.NoneValueError("the_components is None.")
    # </editor-fold>
    # <editor-fold desc="Instance attributes">
    self.components: dict[int, "icomponent.IComponent"] = the_components
    """Contains all components used in PyDD"""
    # </editor-fold>
