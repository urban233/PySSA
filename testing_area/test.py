from PyQt5.QtWidgets import QApplication, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget, QLineEdit
from PyQt5.QtCore import Qt

class CustomLineEdit(QLineEdit):
    def keyPressEvent(self, event):
        # Get the key code
        key = event.key()

        # Check if the key is allowed (e.g., disallow 'ä' in any case)
        if key == Qt.Key_Adiaeresis:
            # Ignore the key event
            return

        # Call the base class implementation to handle other keys
        super(CustomLineEdit, self).keyPressEvent(event)

class MainWindow(QWidget):
    def __init__(self):
        super(MainWindow, self).__init__()

        # Create a table widget with 3 rows and 3 columns
        self.table_widget = QTableWidget(3, 3)

        # Use the custom line edit in the table widget
        for row in range(self.table_widget.rowCount()):
            for col in range(self.table_widget.columnCount()):
                line_edit = CustomLineEdit()
                self.table_widget.setCellWidget(row, col, line_edit)

        # Add the table widget to a layout
        layout = QVBoxLayout()
        layout.addWidget(self.table_widget)

        # Set the layout
        self.setLayout(layout)

if __name__ == '__main__':
    app = QApplication([])

    # Create the main window
    main_window = MainWindow()
    main_window.setWindowTitle('Table Widget Example')
    main_window.show()

    app.exec_()

