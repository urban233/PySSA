import pathlib

from src.pyssa.gui.qt import QtGui


class Workspace:

  def __init__(self, a_path: pathlib.Path):
    self.path = a_path
    self._model = QtGui.QStandardItemModel()

  def get_model(self):
    return self._model

  def construct_project_db_path(self, a_project_name: str):
    return pathlib.Path(self.path / f"{a_project_name}.db")

  def delete_project(self, a_project_name: str):
    self._model.removeRow(self._model.findItems(a_project_name)[0].row())
    self.construct_project_db_path(a_project_name).unlink()

  def get_projects_as_string_list(self) -> list[str]:
    """Converts the internal QStandardItemModel into a list of project names.

    Returns:
      list[str]: A list containing the names of all projects in the workspace.
    """
    return [
        self._model.item(row).text()
        for row in range(self._model.rowCount())
        if self._model.item(row) is not None
    ]
