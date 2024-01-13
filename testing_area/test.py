#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2022
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/urban233/PySSA>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
from PyQt5.QtCore import Qt, QModelIndex
from PyQt5.QtGui import QStandardItem
from PyQt5.QtWidgets import QApplication, QMainWindow, QTreeView, QVBoxLayout, QWidget, QPushButton
from PyQt5 import QtCore


class TreeModel(QtCore.QAbstractItemModel):
    def __init__(self, parent=None):
        super(TreeModel, self).__init__(parent)
        self.root_item = QStandardItem("Root Item")

    def rowCount(self, parent_index):
        if not parent_index.isValid():
            return self.root_item.rowCount()
        else:
            parent_item = parent_index.internalPointer()
            return parent_item.rowCount()

    def columnCount(self, parent_index):
        return self.root_item.columnCount()

    def data(self, index, role):
        if not index.isValid():
            return None

        item = index.internalPointer()

        if role == Qt.DisplayRole:
            return item.data(index.column())

        return None

    def index(self, row, column, parent_index):
        if not self.hasIndex(row, column, parent_index):
            return QModelIndex()

        if not parent_index.isValid():
            parent_item = self.root_item
        else:
            parent_item = parent_index.internalPointer()

        child_item = parent_item.child(row)
        if child_item:
            return self.createIndex(row, column, child_item)
        else:
            return QModelIndex()

    def parent(self, index):
        if not index.isValid():
            return QModelIndex()

        child_item = index.internalPointer()
        parent_item = child_item.parent()

        if parent_item == self.root_item:
            return QModelIndex()

        return self.createIndex(parent_item.row(), 0, parent_item)

    def flags(self, index):
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable

    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.root_item.data(section)

        return None

    def addItem(self, parent_index, text):
        parent_item = parent_index.internalPointer() if parent_index.isValid() else self.root_item
        new_item = QStandardItem(text)
        parent_item.appendRow(new_item)

        # Emit dataChanged signal to notify views that the data has changed
        self.dataChanged.emit(self.index(0, 0, parent_index), self.index(parent_item.rowCount() - 1, self.columnCount(parent_index) - 1, parent_index))

    def getItemIndex(self, item):
        if item:
            return self.createIndex(item.row(), 0, item)

        return QModelIndex()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main_window = QMainWindow()
    main_window.setWindowTitle('Custom Tree Model Example')
    main_window.setGeometry(100, 100, 800, 600)

    central_widget = QWidget(main_window)
    main_window.setCentralWidget(central_widget)

    layout = QVBoxLayout(central_widget)
    tree_view = QTreeView(central_widget)
    tree_model = TreeModel()
    tree_view.setModel(tree_model)

    # Button to add new item
    add_button = QPushButton('Add Item', central_widget)
    add_button.clicked.connect(lambda: tree_model.addItem(tree_view.currentIndex(), 'New Item'))

    # Button to retrieve data
    retrieve_button = QPushButton('Retrieve Data', central_widget)
    retrieve_button.clicked.connect(lambda: print(tree_model.getItemData(tree_model.getItemIndex(tree_model.root_item))))

    layout.addWidget(tree_view)
    layout.addWidget(add_button)
    layout.addWidget(retrieve_button)

    main_window.show()
    sys.exit(app.exec_())
