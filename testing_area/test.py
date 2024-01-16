import sys
from PyQt5.QtWidgets import QApplication, QTableView, QItemDelegate, QComboBox
from PyQt5.QtGui import QStandardItemModel
from PyQt5.QtCore import Qt

class ComboBoxDelegate(QItemDelegate):
    def createEditor(self, parent, option, index):
        combo_box = QComboBox(parent)
        # Add items to the combo box if needed
        combo_box.addItems(["Option 1", "Option 2", "Option 3"])
        return combo_box

    def setEditorData(self, editor, index):
        current_text = index.model().data(index, role=Qt.DisplayRole)
        combo_box = editor
        combo_box.setCurrentText(current_text)

    def setModelData(self, editor, model, index):
        combo_box = editor
        model.setData(index, combo_box.currentText(), role=Qt.EditRole)

class MainWindow(QTableView):
    def __init__(self):
        super().__init__()

        # Create a simple model
        model = QStandardItemModel(4, 3)
        self.setModel(model)

        # Set the combo box delegate for a specific column (e.g., column 1)
        combo_delegate = ComboBoxDelegate()
        self.setItemDelegateForColumn(1, combo_delegate)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
