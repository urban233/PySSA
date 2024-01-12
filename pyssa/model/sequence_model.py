from PyQt5 import QtGui


class SequenceModel(QtGui.QStandardItemModel):
    def __init__(self):
        super(SequenceModel, self).__init__()

    def add_sequence(self, a_sequence: str) -> None:
        i = self.rowCount()
        j = 0
        for tmp_amino_acid in a_sequence:
            self.setItem(i, j, QtGui.QStandardItem(tmp_amino_acid))
            print(self.item(i,j).text())