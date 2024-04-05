#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2022
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
"""Module for functions which change the ui in some way."""
from PyQt5 import QtCore
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
from pyssa.util import input_validator


def select_matching_string_in_q_list_view(a_string_to_match: str,
                                          a_q_list_view_to_select: QtWidgets.QListView,
                                          a_q_line_edit: QtWidgets.QLineEdit):
    """Selects a QListView item from a QListView with a string to match and sets it into a q_line_edit.

    Notes:
        This function must run in the MAIN THREAD!
    """
    # Clears any active selections
    # if a_q_list_view_to_select.currentIndex().isValid():
    #     a_q_list_view_to_select.selectionModel().clearCurrentIndex()

    if a_string_to_match == "":
        a_q_line_edit.clear()
        a_q_list_view_to_select.selectionModel().clearCurrentIndex()
        a_q_list_view_to_select.selectionModel().clearSelection()
        return

    tmp_match_items = input_validator.find_match_in_model(a_q_list_view_to_select.model(), a_string_to_match)
    if len(tmp_match_items) == 0:
        return
    # Sets the selection
    a_q_list_view_to_select.selectionModel().setCurrentIndex(
        tmp_match_items[0].index(),
        QtCore.QItemSelectionModel.Clear | QtCore.QItemSelectionModel.Select
    )
    a_q_line_edit.setText(tmp_match_items[0].data(Qt.DisplayRole))
