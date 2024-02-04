import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import numpy as np

# Sample data
data = np.random.randn(100)


class MyMainWindow(QMainWindow):
    def __init__(self):
        super(MyMainWindow, self).__init__()

        self.central_widget = QWidget(self)
        self.setCentralWidget(self.central_widget)

        # Create a Matplotlib figure and canvas
        self.figure = Figure()
        self.ax1 = self.figure.add_subplot(211)
        self.ax2 = self.figure.add_subplot(212)

        self.canvas = FigureCanvas(self.figure)
        self.ax1.plot(np.arange(len(data)), data, color='blue', label='Line Plot')
        self.ax1.set_title('Top Subplot (Line Plot)')
        self.ax1.legend()

        self.ax2.hist(data, bins=20, color='green', alpha=0.7, label='Histogram')
        self.ax2.set_title('Bottom Subplot (Histogram)')
        self.ax2.legend()

        # Add some space between subplots
        self.figure.subplots_adjust(hspace=0.4)

        # Layout setup
        layout = QVBoxLayout(self.centralWidget())
        layout.addWidget(self.canvas)


def main():
    app = QApplication(sys.argv)
    main_window = MyMainWindow()
    main_window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
