import sys
import numpy as np
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QMainWindow, QApplication, QVBoxLayout, QWidget
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

class MyMainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Horizontal Bar Chart in PyQt5")

        # Create a QWidget to hold the Matplotlib plot
        widget = QWidget()
        self.setCentralWidget(widget)
        layout = QVBoxLayout()
        widget.setLayout(layout)

        # Create a Matplotlib figure and canvas
        self.figure = plt.Figure()
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas)

        # Generate sample data
        data = np.random.uniform(0, 10, 1000)  # Generate 1000 random distances

        # Define the bins and calculate the frequencies
        bins = np.arange(0, 5, 1)
        frequencies, _ = np.histogram(data, bins=bins)

        # Create the horizontal bar chart with dedicated space between bars
        ax = self.figure.add_subplot(111)
        bar_height = 0.2  # Adjust as needed
        space_height = 0.1  # Adjust as needed
        bars = []
        for i, freq in enumerate(frequencies):
            y = bins[i] + space_height * i  # Adjusted y-position
            bar = ax.barh(y, freq, height=bar_height, align='center')
            bars.append(bar)

        # Add labels for each bar with their respective frequency
        for i, bar in enumerate(bars):
            for rect in bar:
                width = rect.get_width()
                ax.text(width, rect.get_y() + rect.get_height()/2, f'{int(width)}',
                        va='center', ha='left')

        # Customize the plot
        ax.set_xlabel('Frequency')
        ax.set_ylabel('Distance Intervals')

        # Set y-tick labels to indicate intervals
        ax.set_yticks(bins[:-1] + space_height * np.arange(len(frequencies)))
        ax.set_yticklabels([f'{bins[i]} - {bins[i+1]}' for i in range(len(bins)-1)])

        # Adjust layout
        self.figure.tight_layout()

        # Show the plot
        self.canvas.draw()

def main():
    app = QApplication(sys.argv)
    window = MyMainWindow()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()

