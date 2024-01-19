import sys
from PyQt5.QtCore import Qt, QEvent
from PyQt5.QtGui import QStandardItemModel, QStandardItem
from PyQt5.QtWidgets import QApplication, QTableView, QStyledItemDelegate, QComboBox, QWidget, QVBoxLayout, QStyleOptionComboBox, QStyle

class ComboBoxDelegate(QStyledItemDelegate):
    def __init__(self, parent=None):
        super(ComboBoxDelegate, self).__init__(parent)
        self.current_index = None

    def createEditor(self, parent, option, index):
        combo_box = QComboBox(parent)
        combo_box.addItems(["Option 1", "Option 2", "Option 3"])
        combo_box.activated.connect(self.commitAndCloseEditor)
        return combo_box

    def updateEditorGeometry(self, editor, option, index):
        editor.setGeometry(option.rect)

    def setEditorData(self, editor, index):
        value = index.data(Qt.DisplayRole)
        editor.setCurrentText(value)

    def setModelData(self, editor, model, index):
        model.setData(index, editor.currentText(), role=Qt.EditRole)

    def paint(self, painter, option, index):
        value = index.data(Qt.DisplayRole)
        style_option = QStyleOptionComboBox()
        style_option.rect = option.rect
        style_option.currentText = value if value is not None else ""

        # Draw the ComboBox
        QApplication.style().drawComplexControl(QStyle.CC_ComboBox, style_option, painter)

        # Draw the selected item text
        painter.drawText(option.rect, Qt.AlignLeft | Qt.AlignVCenter, str(value))

    def editorEvent(self, event, model, option, index):
        if event.type() == QEvent.MouseButtonPress and event.button() == Qt.LeftButton:
            self.current_index = index
            view = self.parent()  # Get the view from the parent
            if view is not None:
                view.openPersistentEditor(index)
            return True
        return super().editorEvent(event, model, option, index)

    def commitAndCloseEditor(self):
        editor = self.sender()
        if isinstance(editor, QComboBox):
            self.commitData.emit(editor)
            self.closeEditor.emit(editor)

class TableViewWithComboBox(QWidget):
    def __init__(self):
        super(TableViewWithComboBox, self).__init__()

        self.initUI()

    def initUI(self):
        self.table_view = QTableView(self)
        self.model = QStandardItemModel(self)
        self.model.setHorizontalHeaderLabels(["Column 1", "Column 2"])
        self.table_view.setModel(self.model)

        # Add some sample data
        for row in range(5):
            item1 = QStandardItem(f"Item {row + 1}")
            item2 = QStandardItem()
            self.model.setItem(row, 0, item1)
            self.model.setItem(row, 1, item2)

        # Set the custom delegate for the second column
        self.table_view.setItemDelegateForColumn(1, ComboBoxDelegate(self))

        layout = QVBoxLayout(self)
        layout.addWidget(self.table_view)
        self.setLayout(layout)

        self.setGeometry(100, 100, 500, 300)
        self.setWindowTitle('TableView with ComboBox')
        self.show()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = TableViewWithComboBox()
    sys.exit(app.exec_())
