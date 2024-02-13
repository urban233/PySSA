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

import pymol
from pymol import cmd


def get_residue_colors( sele ):
    pymol.stored.colors = []
    cmd.iterate(sele, "stored.colors.append((chain, resi, name, color))")
    res_colors = {}
    for chain, resi, name, color in pymol.stored.colors:
        if name == 'CA':  # c-alpha atom
            res_colors[(chain, resi, name)] = color
    return res_colors


if __name__ == "__main__":
    PYMOL_COLORS = [
        "red",
        "green",
        "limegreen",
        "blue",
        "skyblue",
        "yellow",
        "limon",
        "magenta",
        "hotpink",
        "violet",
        "cyan",
        "greencyan",
        "orange",
        "lightorange",
        "white",
    ]
    cmd.fetch("3bmp")
    for tmp_color in PYMOL_COLORS:
        cmd.color(tmp_color)
        print(get_residue_colors("3bmp"))
