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
            print(color)
            break
            # res_colors[(chain, resi, name)] = color
    return res_colors

def get_residue_repr( sele ):
    pymol.stored.reps = []
    cmd.iterate(sele, "stored.reps.append((chain, resi, name, reps))")
    res_colors = {}
    for chain, resi, name, reps in pymol.stored.reps:
        if name == 'CA':  # c-alpha atom
            print(reps)
            res_colors[(chain, resi, name)] = reps
    return res_colors

def get_chain_repr_index(a_selection_string: str, chain_letter: str) -> int:
    pymol.stored.reps = []
    cmd.iterate(a_selection_string, "stored.reps.append((chain, resi, name, reps))")
    for chain, resi, name, reps_index in pymol.stored.reps:
        if chain == chain_letter:  # c-alpha atom
            return reps_index

def get_all_combis_as_dict() -> list:
    from itertools import product

    elements = ['sticks', 'spheres', 'surface', 'cartoon', 'ribbon', 'lines', 'mesh', 'dots']

    # Generate all combinations of 0s and 1s for the given number of elements
    binary_combinations = list(product([0, 1], repeat=len(elements)))

    results: list = []
    # Print all combinations
    for combination in binary_combinations:
        results.append({elements[i]: bit for i, bit in enumerate(combination)})
    return results

def make_dict(keys: list, values: list) -> dict:
    my_dict = {k: v for k, v in zip(keys, values)}
    return my_dict


if __name__ == "__main__":
    from pymol import preset

    PYMOL_COLORS = [
        "red", "tv_red", "salmon", "raspberry",
        "green", "tv_green", "palegreen", "forest",
        "blue", "tv_blue", "lightblue", "skyblue",
        "yellow", "tv_yellow", "paleyellow", "sand",
        "magenta", "purple", "pink", "hotpink",
        "cyan", "aquamarine", "palecyan", "teal",
        "orange", "tv_orange", "lightorange", "olive",
        "white", "grey70", "grey30", "black"
    ]
    cmd.fetch("3bmp")
    for tmp_color in PYMOL_COLORS:
        print(tmp_color)
        cmd.color(tmp_color, "3bmp")
        get_residue_colors("3bmp")
    # PYMOL_REPS = {
    #     1: 'sticks',
    #     2: 'spheres',
    #     3: 'ball and stick',  # Nothing is known on how to hide this repr
    #     4: 'surface',
    #     32: 'cartoon',
    #     64: 'ribbon',
    #     128: 'lines',
    #     256: 'mesh',
    #     512: 'dots',
    # }
    # cmd.fetch("3bmp")
    # #cmd.hide("cartoon", "3bmp")
    # cmd.show("lines", "3bmp")
    # # preset.ball_and_stick("3bmp")
    # print(get_chain_repr_index("3bmp", "A"))
    #
    # combis = get_all_combis_as_dict()
    # indices = range(512)
    # my_dict = {k: v for k, v in zip(indices, combis)}
    # print(my_dict)
    #
    # # _---------------------------------------------------
    # from itertools import product
    #
    # # List of elements
    # elements = ['sticks', 'spheres', 'surface', 'cartoon', 'ribbon', 'lines', 'mesh', 'dots']
    #
    # # Generate all combinations of 0s and 1s for the given number of elements
    # binary_combinations = list(product([0, 1], repeat=len(elements)))
    #
    # # Define the filename for the Python script
    # filename = "pymol_commands.py"
    #
    # # Open the file for writing
    # with open(filename, "w") as file:
    #     file.write("from test import get_chain_repr_index\n")
    #     file.write("from test import get_all_combis_as_dict\n")
    #     file.write("from test import make_dict\n")
    #
    #     file.write("combis_dict = get_all_combis_as_dict()\n")
    #     file.write("indices_list = []\n")
    #     # Loop through each combination
    #     for index, combination in enumerate(binary_combinations):
    #         # Generate PyMOL commands for each element based on the combination
    #         cmd_list = []
    #         for i, bit in enumerate(combination):
    #             if bit:
    #                 cmd_list.append(f'cmd.show("{elements[i]}")')
    #             else:
    #                 cmd_list.append(f'cmd.hide("{elements[i]}")')
    #
    #         # Write the code for this combination to the file
    #         file.write(f"# Combination {index + 1}\n")
    #         file.write("from pymol import cmd\n")
    #         file.write("cmd.fetch('3bmp')\n")
    #         file.write("\n")
    #         for cmd in cmd_list:
    #             file.write(f"{cmd}\n")
    #         file.write("\n")
    #         file.write("indices_list.append(get_chain_repr_index('3bmp', 'A'))\n")
    #         file.write("cmd.reinitialize()\n\n")
    #
    #     file.write("print(make_dict(indices_list, combis_dict))")
