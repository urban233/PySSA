import pathlib
from dataclasses import dataclass
import datetime
from typing import Callable, Optional

from pyssa.gui.custom_widgets import custom_button


@dataclass
class WorkspaceProject:
  """Dataclass representing a project as a combination of its filepath and modification date
  and not as 'Project' object.
  """
  filepath: pathlib.Path
  date_modified: datetime.datetime

  def __init__(self, a_filepath, a_date_modified, a_callable):
    self.filepath = a_filepath
    self.date_modified = a_date_modified
    self.home_button = self.get_as_button(a_callable)
    self.open_button = self.get_as_button(a_callable)
    self.delete_button = self.get_as_button_with_checkbox()

  def get_as_button(self, a_callable: Callable):
    """Get the dataclass in the form of a Button"""
    return custom_button.ProjectOverviewButton(
      self.get_name(),
      self.date_modified.strftime("%d.%m.%Y %I:%M %p"),
      a_callable
    )

  def get_as_button_with_checkbox(self, a_callable: Optional[Callable] = None):
    """Get the dataclass in the form of a Button with a checkbox"""
    return custom_button.ProjectOverviewButtonWithCheckbox(
      self.get_name(),
      self.date_modified.strftime("%d.%m.%Y %I:%M %p"),
      a_callable
    )

  def get_name(self) -> str:
    """Gets the name of the project based on its filepath."""
    return self.filepath.name.replace(".db", "")
