import os
from typing import Union

from PyQt6 import QtWidgets, QtGui
from PyQt6 import QtCore

from pyssa.gui.custom_widgets import custom_button
from pyssa.model.qmodel import project_model
from pyssa.model.data_classes import workspace_project
from pyssa.model.preference import model_definitions


class LineEditSearch(QtWidgets.QLineEdit):

  def __init__(self, parent):
    super().__init__(parent)
    self.setPlaceholderText("Search")
    self.setGeometry(50, 20, 200, 30)
    self.setStyleSheet("""
      QLineEdit {
        font-size: 13px;
        padding: 3px;
      }
    """)
    self.model: Union["project_model.ProjectModel", None] = None
    self.textChanged.connect(self.on_text_changed)
    self.search_dialog = _SearchDialog(parent)
    self.setFixedSize(275, 30)

  def add_model_to_search(self, a_model: Union["project_model.ProjectModel"]):
    self.model: Union["project_model.ProjectModel"] = a_model

  def connect_search_dialog_with_line_edit(self, a_callable):
    self.search_dialog.item_clicked.connect(a_callable)

  def on_text_changed(self, text):
    if text:
      # Perform the search with a delay to avoid processing on every keystroke
      QtCore.QTimer.singleShot(300, lambda: self.search(text))
    else:
      self.search_dialog.hide()

  def search(self, text):
    if not text:
      self.search_dialog.hide()
      return

    results = self.model.findItems(text, QtCore.Qt.MatchFlag.MatchContains)

    if results:
      self.search_dialog.update_search_results(results[0:6])
      self.search_dialog.show_below_widget(self)
    else:
      self.search_dialog.hide()

    # Keep the focus on the QLineEdit
    self.setFocus()


class _SearchDialog(QtWidgets.QDialog):
  item_clicked = QtCore.pyqtSignal(tuple)

  def __init__(self, parent=None):
    super().__init__(parent)
    # <editor-fold desc="GUI">
    self.list_widget = QtWidgets.QListWidget()
    layout = QtWidgets.QVBoxLayout()
    layout.addWidget(self.list_widget)
    self.setLayout(layout)
    self.setStyleSheet("""
      QDialog {
        background: white;
        border: 1px solid #636363;
        border-radius: 4px;
      }
      QListWidget {
        border: none;
      }
      QListWidget::item {
        background: white;
        padding: 5px;
        font-size: 20px;
      }
      QListWidget::item:hover {
        background: #e3ecfa;
      }
      QListView::item:selected {
        color: black;
        background: #a5b9d1;
        border: 2px solid #4a78b0;
      }
    """)
    # Set the dialog to behave like a tooltip, so it doesn't steal focus
    self.setWindowFlags(QtCore.Qt.WindowType.Tool | QtCore.Qt.WindowType.FramelessWindowHint)
    self.setFixedWidth(275)
    self.setSizePolicy(QtWidgets.QSizePolicy.Policy.Minimum, QtWidgets.QSizePolicy.Policy.Minimum)
    self.list_widget.setSizePolicy(QtWidgets.QSizePolicy.Policy.Minimum, QtWidgets.QSizePolicy.Policy.Minimum)
    # </editor-fold>
    self.list_widget.clicked.connect(self.emit_item_clicked)

  def update_search_results(self, results: list[QtGui.QStandardItem]):
    self.list_widget.clear()
    for tmp_result in results:
      tmp_workspace_project: "workspace_project.WorkspaceProject" = tmp_result.data(
        model_definitions.RolesEnum.OBJECT_ROLE)
      item = QtWidgets.QListWidgetItem(tmp_workspace_project.get_name())
      item.setData(model_definitions.RolesEnum.OBJECT_ROLE, tmp_workspace_project)
      self.list_widget.addItem(item)
    self.setFixedHeight(self._calculate_height())

  def _calculate_height(self):
    return 50 + (self.list_widget.count() - 1) * 30

  def show_below_widget(self, widget):
    pos = widget.mapToGlobal(widget.rect().bottomLeft())
    self.move(pos)
    self.show()

  def mousePressEvent(self, event):
    # Override the mousePressEvent to prevent closing when clicking inside the dialog
    event.accept()

  def emit_item_clicked(self) -> None:
    tmp_workspace_project: "workspace_project.WorkspaceProject" = self.list_widget.currentItem().data(model_definitions.RolesEnum.OBJECT_ROLE)
    self.item_clicked.emit((tmp_workspace_project.get_name(), tmp_workspace_project.date_modified))
    self.close()
