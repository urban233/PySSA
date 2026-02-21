import sys
from PyQt6.QtGui import QStandardItemModel, QStandardItem
from PyQt6.QtCore import Qt
model = QStandardItemModel()
sec = QStandardItem('Proteins')
model.appendRow(sec)
prot = QStandardItem('1abc')
prot.setData('protein', Qt.ItemDataRole.UserRole + 1)
sec.appendRow(prot)
header = QStandardItem('Chains')
prot.appendRow(header)
chain = QStandardItem('A')
chain.setData('chain', Qt.ItemDataRole.UserRole + 1)
header.appendRow(chain)
res = QStandardItem('42 - GLY')
res.setData('residue', Qt.ItemDataRole.UserRole + 1)
chain.appendRow(res)
atom = QStandardItem('CA')
atom.setData('atom', Qt.ItemDataRole.UserRole + 1)
res.appendRow(atom)
def test(idx):
  if idx.data(Qt.ItemDataRole.UserRole + 1) == 'atom':
    r = idx.parent()
    c = r.parent()
    p = c.parent().parent()
    return f'({p.data(Qt.ItemDataRole.DisplayRole)} and chain {c.data(Qt.ItemDataRole.DisplayRole)})'
print('Atom Result:', test(atom.index()))
