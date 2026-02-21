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
