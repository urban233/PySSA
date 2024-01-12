from PyQt5.QtCore import Qt, QAbstractTableModel, QModelIndex
from pyssa.util import enums


class ApplicationModel:
    application_state: dict
    app_enums = enums.ApplicationModelEnum

    def __init__(self, a_project):
        self.application_state: dict = {
            self.app_enums.PROJECT: a_project
        }


class ProteinTableModel(QAbstractTableModel):
    def __init__(self, proteins, parent=None):
        super().__init__(parent)
        self.proteins = proteins

    def rowCount(self, parent=QModelIndex()):
        return len(self.proteins)

    def columnCount(self, parent=QModelIndex()):
        # Assuming two columns: attribute name and attribute value
        return 2

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid() or role != Qt.DisplayRole:
            return None

        protein = self.proteins[index.row()]

        if index.column() == 0:  # Attribute name column
            attribute_names = list(protein.__annotations__.keys())
            return attribute_names[index.column()]
        elif index.column() == 1:  # Attribute value column
            attribute_name = list(protein.__annotations__.keys())[index.row()]
            return getattr(protein, attribute_name)

        return None

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            if section == 0:
                return "Attribute Name"
            elif section == 1:
                return "Attribute Value"

        return None
