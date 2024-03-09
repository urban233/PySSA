from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
from PyQt5 import QtCore


class FastaFileImportPreviewTableDelegate(QtWidgets.QItemDelegate):
    def __init__(self, parent=None):
        super(FastaFileImportPreviewTableDelegate, self).__init__(parent)

    def createEditor(self, parent, option, index):
        editor = super(FastaFileImportPreviewTableDelegate, self).createEditor(parent, option, index)
        return editor

    def setModelData(self, editor, model, index):
        tmp_item_text = editor.text()
        tmp_column = index.column()

        if tmp_column == 0:
            allowed_chars = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_-"
            tmp_item_new_text = ''.join(char for char in tmp_item_text if char in allowed_chars)
        elif tmp_column == 1:
            allowed_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
            tmp_item_new_text = ''.join(char for char in tmp_item_text if char in allowed_chars)
            tmp_item_new_text = ''.join(char for char in tmp_item_new_text if char not in self._get_all_chain_letters(model, index))
        # elif tmp_column == 2:
        #     allowed_chars = {"C", "D", "S", "Q", "K", "I", "P", "T", "F", "N", "G", "H", "L", "R", "W", "A", "V", "E", "Y", "M"}
        #     tmp_item_new_text = ''.join(char for char in tmp_item_text if char in allowed_chars)
        else:
            model.setData(index, "")
            return
        if tmp_item_new_text == "":
            model.setData(index, "Invalid input!")
        else:
            model.setData(index, tmp_item_new_text)

    def _get_all_chain_letters(self, model: QtCore.QAbstractTableModel, index):
        # TODO: this needs to be fixed! (Problem: all chain letters are considered for one sequence, which is wrong!)
        tmp_all_chains = []
        for tmp_row_no in range(model.rowCount()):
            tmp_all_chains.append(model.index(tmp_row_no, 1).data(Qt.DisplayRole))
        return tmp_all_chains
