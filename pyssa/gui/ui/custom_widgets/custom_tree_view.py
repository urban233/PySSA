from PyQt5.QtWidgets import QApplication, QMainWindow, QTreeView, QMenu, QAction, QVBoxLayout, QWidget, QAbstractItemView, QFileSystemModel
import sys


class CustomTreeView(QTreeView):
    def __init__(self, parent=None):
        super().__init__(parent)

        # Enable the context menu for the tree view
        self.setContextMenuPolicy(3)  # 3 corresponds to Qt.CustomContextMenu
        self.customContextMenuRequested.connect(self.showContextMenu)

    def contextMenuEvent(self, event):
        # Override the contextMenuEvent method to prevent the default menu
        pass

    def showContextMenu(self, point):
        index = self.indexAt(point)
        if index.isValid():
            menu = QMenu(self)

            # Example actions
            action1 = QAction("Option 1", self)
            action1.triggered.connect(lambda: self.onContextMenuSelection(index, "Option 1"))

            action2 = QAction("Option 2", self)
            action2.triggered.connect(lambda: self.onContextMenuSelection(index, "Option 2"))

            menu.addAction(action1)
            menu.addAction(action2)

            menu.exec_(self.mapToGlobal(point))

    def onContextMenuSelection(self, index, option):
        print(f"Selected {option} for item at row {index.row()}, column {index.column()}")
