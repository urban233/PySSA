# import sys
# from PyQt5.QtWidgets import QApplication, QDialog, QVBoxLayout, QSplitter, QTableWidget, QTableWidgetItem, QTextEdit, QFrame
#
#
# class ResizableDialog(QDialog):
#     def __init__(self):
#         super().__init__()
#
#         self.setWindowTitle("Resizable Dialog")
#         self.resize(800, 600)  # Set initial size
#
#         main_layout = QVBoxLayout()
#         self.setLayout(main_layout)
#
#         # First splitter to divide the dialog into two sections
#         first_splitter = QSplitter()
#         main_layout.addWidget(first_splitter)
#
#         # Left side (plot area)
#         plot_area = QFrame()
#         first_splitter.addWidget(plot_area)
#
#         # Right side (table)
#         table_widget = QTableWidget(10, 3)
#         first_splitter.addWidget(table_widget)
#
#         # Second splitter within the plot area to split it horizontally
#         second_splitter = QSplitter()
#         plot_area_layout = QVBoxLayout()
#         plot_area.setLayout(plot_area_layout)
#         plot_area_layout.addWidget(second_splitter)
#
#         # Left part of the plot area
#         plot_left = QTextEdit()
#         second_splitter.addWidget(plot_left)
#
#         # Right part of the plot area
#         plot_right = QTextEdit()
#         second_splitter.addWidget(plot_right)
#         second_splitter.setOrientation(0)  # Set orientation to horizontal
#
#
# if __name__ == "__main__":
#     app = QApplication(sys.argv)
#     dialog = ResizableDialog()
#     dialog.show()
#     sys.exit(app.exec_())

import sys
from PyQt5.QtWidgets import QMainWindow, QApplication, QVBoxLayout, QWidget
from PyQt5.QtCore import QTimer, Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
import threading

class PlotThread(threading.Thread):
    def __init__(self, parent, width, height):
        super().__init__()
        self.parent = parent
        self.width = width
        self.height = height

    def run(self):
        aspect_ratio = self.parent.ax.get_aspect()
        aspect_ratio = float(aspect_ratio)  # Convert to float
        if self.width / self.height > aspect_ratio:
            new_width = self.height * aspect_ratio
            new_height = self.height
        else:
            new_width = self.width
            new_height = self.width / aspect_ratio
        self.parent.figure.set_size_inches(new_width / self.parent.figure.dpi,
                                           new_height / self.parent.figure.dpi)
        self.parent.canvas.draw()

class MyMainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.layout = QVBoxLayout(self.central_widget)

        self.figure, self.ax = plt.subplots()
        self.canvas = FigureCanvas(self.figure)
        self.layout.addWidget(self.canvas)

        self.aspect_ratio = self.ax.get_aspect()
        self.aspect_ratio = float(self.aspect_ratio)  # Convert to float

        self.resize_timer = QTimer(self)
        self.resize_timer.setInterval(500)  # Delay in milliseconds
        self.resize_timer.setSingleShot(True)
        self.resize_timer.timeout.connect(self.delayed_resize)

    def delayed_resize(self):
        width = self.centralWidget().width()
        height = self.centralWidget().height()
        plot_thread = PlotThread(self, width, height)
        plot_thread.start()

    def resizeEvent(self, event):
        self.resize_timer.start()
        super().resizeEvent(event)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    mainWindow = MyMainWindow()
    mainWindow.setGeometry(100, 100, 800, 600)
    mainWindow.show()
    sys.exit(app.exec_())

