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
"""Module for the input validator class."""
from pymol import Qt
from PyQt5 import QtGui
from PyQt5 import QtWidgets
from pyssa.gui.ui.styles import styles


class InputValidator:
    """Class for validating any input from the user."""

    def __init__(self) -> None:
        """Empty constructor."""
        pass

    @staticmethod
    def validate_project_name(
        list_of_projects: QtWidgets.QListWidget,
        txt_for_project_name: QtWidgets.QTextEdit,
        lbl_for_status_project_name: QtWidgets.QLabel,
        btn_for_next_step: QtWidgets.QPushButton,
        cb_for_add_reference: QtWidgets.QCheckBox = None,
    ) -> None:
        """This function validates the input of the project name in real-time.

        Args:
            list_of_projects:
                list widget which holds all projects from the workspace
            txt_for_project_name:
                line edit widget which is used to enter the project name
            lbl_for_status_project_name:
                label which is used to give feedback if the input is legal or not
            btn_for_next_step:
                push button which is used to execute either the next step or to create a project
            cb_for_add_reference (optional):
                checkbox widget which is used to add a reference
        """
        if list_of_projects.currentItem() is not None:
            list_of_projects.currentItem().setSelected(False)
        # set color for lineEdit
        txt_for_project_name.setStyleSheet("color: #f44336")
        if len(txt_for_project_name.text()) == 0:
            lbl_for_status_project_name.setText("")
            if cb_for_add_reference is not None:
                cb_for_add_reference.setCheckable(False)
                cb_for_add_reference.setStyleSheet("color: #E1E1E1;")
            btn_for_next_step.setEnabled(False)
            styles.color_button_not_ready(btn_for_next_step)
        elif len(txt_for_project_name.text()) > 20:
            lbl_for_status_project_name.setText("Project name is too long (max. 20 characters).")
            btn_for_next_step.setEnabled(False)
            styles.color_button_not_ready(btn_for_next_step)
            return
        else:
            regex = Qt.QtCore.QRegularExpression()
            regex.setPattern("(([a-z])|([A-Z])|([0-9])|(-)|(_)){0,20}")
            validator = QtGui.QRegularExpressionValidator(regex)
            for i in range(len(txt_for_project_name.text())):
                result = validator.validate(txt_for_project_name.text(), i)
                if result[0] > 0:
                    txt_for_project_name.setStyleSheet("color: #000000")
                    lbl_for_status_project_name.setText("")
                    if cb_for_add_reference is not None:
                        cb_for_add_reference.setCheckable(True)
                        cb_for_add_reference.setStyleSheet("color: black;")
                    btn_for_next_step.setEnabled(True)
                    styles.color_button_ready(btn_for_next_step)
                else:
                    txt_for_project_name.setStyleSheet("color: #f44336")
                    lbl_for_status_project_name.setText("Invalid character.")
                    if cb_for_add_reference is not None:
                        cb_for_add_reference.setCheckable(False)
                        cb_for_add_reference.setStyleSheet("color: #E1E1E1;")
                    btn_for_next_step.setEnabled(False)
                    styles.color_button_not_ready(btn_for_next_step)
                    return
            item = list_of_projects.findItems(
                txt_for_project_name.text(),
                Qt.QtCore.Qt.MatchContains | Qt.QtCore.Qt.MatchExactly,
            )
            if len(item) != 0:
                list_of_projects.setCurrentItem(item[0])
                txt_for_project_name.setStyleSheet("color: #f44336")
                lbl_for_status_project_name.setText("Project name already exists.")
                if cb_for_add_reference is not None:
                    cb_for_add_reference.setCheckable(False)
                    cb_for_add_reference.setStyleSheet("color: #E1E1E1;")
                btn_for_next_step.setEnabled(False)
                styles.color_button_not_ready(btn_for_next_step)

    @staticmethod
    def validate_search_input(
        list_for_projects: QtWidgets.QListWidget,
        txt_for_search: QtWidgets.QTextEdit,
        lbl_for_status_search: QtWidgets.QLabel,
        txt_for_selected_project: QtWidgets.QTextEdit = None,
        status_message: str = None,
    ) -> None:
        """This function validates the input of the project name in real-time.

        Args:
            list_for_projects:
                list widget where all projects from the workspace are stored
            txt_for_search:
                line edit widget which ii used for entering the search term
            lbl_for_status_search:
                label which gives feedback
            txt_for_selected_project:
                line edit widget which is used to display the selected project
            status_message:
                the message which should get displayed.
        """
        if list_for_projects.currentItem() is not None:
            list_for_projects.currentItem().setSelected(False)
        # set color for lineEdit
        txt_for_search.setStyleSheet("background-color: white")
        if len(txt_for_search.text()) == 0:
            lbl_for_status_search.setText("")
            if txt_for_selected_project is not None:
                txt_for_selected_project.setText("")
        else:
            item = list_for_projects.findItems(
                txt_for_search.text(),
                Qt.QtCore.Qt.MatchContains | Qt.QtCore.Qt.MatchExactly,
            )
            if len(item) != 0:
                list_for_projects.setCurrentItem(item[0])
                lbl_for_status_search.setText("")
                if txt_for_selected_project is not None:
                    txt_for_selected_project.setText(list_for_projects.currentItem().text())
            else:
                txt_for_search.setStyleSheet("color: #f44336")
                if status_message is None:
                    lbl_for_status_search.setText("Project name does not exists.")
                else:
                    lbl_for_status_search.setText(status_message)
                if txt_for_selected_project is not None:
                    txt_for_selected_project.setText("")

    @staticmethod
    def validate_protein_name(
        txt_for_protein_name: QtWidgets.QTextEdit,
        lbl_for_status_protein_name: QtWidgets.QLabel,
        btn_next: QtWidgets.QPushButton = None,
    ) -> None:
        """This function validates the input of the protein name in real-time.

        Args:
            txt_for_protein_name:
                line edit widget which ii used for entering the protein name
            lbl_for_status_protein_name:
                label which gives feedback if the protein name is legal
            btn_next:
                push button which is used for the next step
        """
        # set color for lineEdit
        txt_for_protein_name.setStyleSheet("color: #f44336")
        if len(txt_for_protein_name.text()) == 0:
            lbl_for_status_protein_name.setText("")
            if btn_next is not None:
                btn_next.setEnabled(False)
                styles.color_button_not_ready(btn_next)
        else:
            regex = Qt.QtCore.QRegularExpression()
            regex.setPattern("(([a-z])|([A-Z])|([0-9])|(-)|(_)){0,20}")
            validator = QtGui.QRegularExpressionValidator(regex)
            for i in range(len(txt_for_protein_name.text())):
                result = validator.validate(txt_for_protein_name.text(), i)
                if result[0] > 0:
                    txt_for_protein_name.setStyleSheet("color: #000000")
                    lbl_for_status_protein_name.setText("")
                    if btn_next is not None:
                        btn_next.setEnabled(True)
                        styles.color_button_ready(btn_next)
                else:
                    txt_for_protein_name.setStyleSheet("color: #f44336")
                    lbl_for_status_protein_name.setText("Invalid character.")
                    if btn_next is not None:
                        btn_next.setEnabled(False)
                        styles.color_button_not_ready(btn_next)
                    return

    @staticmethod
    def validate_protein_sequence(
        txt_protein_sequence: QtWidgets.QTextEdit,
        lbl_status_protein_sequence: QtWidgets.QLabel,
        btn_next: QtWidgets.QPushButton,
    ) -> None:
        """This function validates the input of the protein sequence in real-time.

        Args:
            txt_protein_sequence:
                line edit widget which is used to enter the protein sequence
            lbl_status_protein_sequence:
                label which gives feedback if the protein sequence is legal or not
            btn_next:
                button which is used to get to the next step
        """
        # set color for lineEdit
        txt_protein_sequence.setStyleSheet("color: #f44336")
        if len(txt_protein_sequence.toPlainText()) == 0:
            lbl_status_protein_sequence.setText("")
            btn_next.setEnabled(False)
            styles.color_button_not_ready(btn_next)
        else:
            regex = Qt.QtCore.QRegularExpression()
            regex.setPattern("(([A])|([C-I])|([K-N])|([P-T])|([V-W])|([Y]))+")
            validator = QtGui.QRegularExpressionValidator(regex)
            for i in range(len(txt_protein_sequence.toPlainText())):
                result = validator.validate(txt_protein_sequence.toPlainText(), i)
                if result[0] > 0:
                    txt_protein_sequence.setStyleSheet("color: #000000")
                    lbl_status_protein_sequence.setText("")
                    btn_next.setEnabled(True)
                    styles.color_button_ready(btn_next)
                else:
                    txt_protein_sequence.setStyleSheet("color: #f44336")
                    lbl_status_protein_sequence.setText("Invalid character.")
                    btn_next.setEnabled(False)
                    styles.color_button_not_ready(btn_next)
                    return
