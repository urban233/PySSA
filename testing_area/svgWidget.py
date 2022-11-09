from PyQt5 import QtWidgets
from PyQt5 import QtSvg
import sys

app = QtWidgets.QApplication(sys.argv)

viewer = QtSvg.QSvgWidget()
viewer.load("/home/matt/Documents/test_pymol/Results/plots/distance_plot/distance_plot_selected_prediction_1.svg")
viewer.show()

if __name__ == '__main__':
    app.exec()
