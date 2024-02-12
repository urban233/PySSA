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

import matplotlib.pyplot as plt
import numpy as np


def decreased_space_bar_chart(data):
    num_bars = len(data)
    bar_width = 2 / (num_bars + 1)  # Adjusting bar width to decrease space between bars

    fig, ax = plt.subplots()

    # Calculate x positions for bars
    x_positions = [i * (1 + bar_width) for i in range(num_bars)]

    # Plot bars
    ax.bar(x_positions, data, width=bar_width)

    ax.set_xticks([i + 0.1 * bar_width for i in x_positions])  # Adjusting x ticks position

    plt.show()


if __name__ == "__main__":
    # Example data
    data = np.random.randint(1, 10, size=5)

    # Generate bar chart with decreased space between bars
    decreased_space_bar_chart(data)
