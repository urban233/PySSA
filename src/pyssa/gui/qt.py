#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/urban233/PySSA>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

try:
    from PyQt6 import QtCore, QtGui, QtWidgets, QtSql
    from PyQt6.QtCore import pyqtSignal
    from PyQt6.QtCore import Qt
    from PyQt6.QtWebEngineWidgets import QWebEngineView

    IS_PYQT6 = True
except ImportError:
    from PyQt5 import QtCore, QtGui, QtWidgets, QtSql
    from PyQt5.QtCore import pyqtSignal
    from PyQt5.QtCore import Qt

    IS_PYQT6 = False

__docformat__ = "google"

# Handle enum differences (PyQt6 uses enum classes, PyQt5 uses constants)
if not IS_PYQT6:
    QtCore.Qt.AlignmentFlag = QtCore.Qt
    QtCore.Qt.WindowType = QtCore.Qt
    QtGui.QAction = QtWidgets.QAction
    QtGui.QActionGroup = QtWidgets.QActionGroup

# Alias `exec` method for QDialog, QApplication, etc.
if not IS_PYQT6:
    QtWidgets.QDialog.exec = QtWidgets.QDialog.exec_
    QtWidgets.QApplication.exec = QtWidgets.QApplication.exec_
# else:
#     QtWidgets.QDialog.exec = QtWidgets.QDialog.exec
#     QtWidgets.QApplication.exec = QtWidgets.QApplication.exec

# Handle moved classes (e.g., QStringListModel moved from QtGui to QtCore in PyQt6)
if IS_PYQT6:
    QtGui.QStringListModel = QtCore.QStringListModel
    QtWidgets.QAbstractItemView.NoEditTriggers = (
        QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
    )
    Qt.WindowModal = Qt.WindowModality.WindowModal
    Qt.DisplayRole = Qt.ItemDataRole.DisplayRole
    Qt.WaitCursor = Qt.CursorShape.WaitCursor

    QtWidgets.QActionGroup = QtGui.QActionGroup
    QtWidgets.QAction = QtGui.QAction
    QtWidgets.QShortcut = QtGui.QShortcut
    QtCore.QSortFilterProxyModel.setFilterRegExp = QtCore.QSortFilterProxyModel.setFilterRegularExpression
    QtGui.QFont.Monospace = QtGui.QFont.StyleHint.Monospace

    def copy_attributes(target_class, source_class):
        for attr in dir(source_class):
            if not attr.startswith('_'):
                setattr(target_class, attr, getattr(source_class, attr))

    copy_attributes(QtCore.QEvent, QtCore.QEvent.Type)
    copy_attributes(QtCore.Qt, QtCore.Qt.AlignmentFlag)
    copy_attributes(QtCore.Qt, QtCore.Qt.CaseSensitivity)
    copy_attributes(QtCore.Qt, QtCore.Qt.CheckState)
    copy_attributes(QtCore.Qt, QtCore.Qt.ContextMenuPolicy)
    copy_attributes(QtCore.Qt, QtCore.Qt.DockWidgetArea)
    copy_attributes(QtCore.Qt, QtCore.Qt.FocusPolicy)
    copy_attributes(QtCore.Qt, QtCore.Qt.GestureType)
    copy_attributes(QtCore.Qt, QtCore.Qt.ItemFlag)
    copy_attributes(QtCore.Qt, QtCore.Qt.Key)
    copy_attributes(QtCore.Qt, QtCore.Qt.KeyboardModifier)
    copy_attributes(QtCore.Qt, QtCore.Qt.MouseButton)
    copy_attributes(QtCore.Qt, QtCore.Qt.Orientation)
    copy_attributes(QtCore.Qt, QtCore.Qt.WindowType)
    copy_attributes(QtGui.QFont, QtGui.QFont.StyleHint)
    copy_attributes(QtWidgets.QAbstractItemView, QtWidgets.QAbstractItemView.ScrollHint)
    copy_attributes(QtWidgets.QAbstractItemView, QtWidgets.QAbstractItemView.SelectionBehavior)
    copy_attributes(QtWidgets.QAbstractItemView, QtWidgets.QAbstractItemView.SelectionMode)
    copy_attributes(QtWidgets.QBoxLayout, QtWidgets.QBoxLayout.Direction)
    copy_attributes(QtWidgets.QMainWindow, QtWidgets.QMainWindow.DockOption)
    # copy_attributes(QtWidgets.QOpenGLWidget, QtOpenGLWidgets.QOpenGLWidget.UpdateBehavior)
    copy_attributes(QtWidgets.QSizePolicy, QtWidgets.QSizePolicy.Policy)
    copy_attributes(QtWidgets.QTreeWidgetItem, QtWidgets.QTreeWidgetItem.ChildIndicatorPolicy)

    QtCore.Qt.MidButton = QtCore.Qt.MiddleButton
    QtCore.Qt.WA_LayoutUsesWidgetRect = QtCore.Qt.WidgetAttribute.WA_LayoutUsesWidgetRect
