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
"""Module for the input validator class."""
import logging
import re
from PyQt5.QtCore import Qt
from PyQt5 import QtGui
from PyQt5 import QtWidgets
from PyQt5 import QtCore

from src.pyssa.internal.thread import thread_util
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = 'google'


def validate_input_for_project_name(
    the_current_entered_text: str, the_current_projects: set
) -> tuple[bool, str, str]:
  """Validates the input for a project name.

  Args:
      the_current_entered_text (str): The text entered by the user.
      the_current_projects (set): The set of existing project names.

  Returns:
      A tuple with three elements:
          a boolean indicating if the input is valid
          a string that contains the stylesheet for the line edit
          a string that contains a message

  Raises:
      exception.IllegalArgumentError: If either `the_current_entered_text` or `the_current_projects` is None.
  """
  # <editor-fold desc="Checks">
  if the_current_entered_text is None:
    logger.error('the_current_entered_text is None.')
    raise exception.IllegalArgumentError('the_current_entered_text is None.')
  if the_current_projects is None:
    logger.error('the_current_projects is None.')
    raise exception.IllegalArgumentError('the_current_projects is None.')

  # </editor-fold>

  allowed_chars = {
      '0',
      '1',
      '2',
      '3',
      '4',
      '5',
      '6',
      '7',
      '8',
      '9',
      'a',
      'b',
      'c',
      'd',
      'e',
      'f',
      'g',
      'h',
      'i',
      'j',
      'k',
      'l',
      'm',
      'n',
      'o',
      'p',
      'q',
      'r',
      's',
      't',
      'u',
      'v',
      'w',
      'x',
      'y',
      'z',
      'A',
      'B',
      'C',
      'D',
      'E',
      'F',
      'G',
      'H',
      'I',
      'J',
      'K',
      'L',
      'M',
      'N',
      'O',
      'P',
      'Q',
      'R',
      'S',
      'T',
      'U',
      'V',
      'W',
      'X',
      'Y',
      'Z',
      '-',
      '_',
  }
  for char in the_current_entered_text:
    if char not in allowed_chars:
      return (
          False,
          """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}""",
          'Invalid character.',
      )
  if the_current_entered_text in the_current_projects:
    return (
        False,
        """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}""",
        'Project name already exists!',
    )
  if the_current_entered_text == '':
    return (
        False,
        """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}""",
        'Please enter a project name.',
    )
  else:
    return True, """QLineEdit {color: #000000; border-color: #DCDBE3;}""", ''


def find_match_in_model(
    a_model, a_text_to_search_for
) -> list[QtGui.QStandardItem]:
  """Find match in model.

  Args:
      a_model: The model to search in.
      a_text_to_search_for: The text to search for.

  Returns:
      A list of exactly or partially matched items in the model.
  """
  tmp_exactly_matched_items = a_model.findItems(
      a_text_to_search_for,
      QtCore.Qt.MatchExactly,
  )
  tmp_partial_matched_items = a_model.findItems(
      a_text_to_search_for,
      QtCore.Qt.MatchContains,
  )
  if len(tmp_exactly_matched_items) == 1:
    return tmp_exactly_matched_items
  return tmp_partial_matched_items


class InputValidator:
  """Class for validating any input from the user."""

  line_edit: QtWidgets.QLineEdit
  """The line edit widget to validate."""

  def __init__(self, a_line_edit: QtWidgets.QLineEdit) -> None:
    """Constructor.

    Args:
        a_line_edit (QtWidgets.QLineEdit): The line edit widget to validate.

    Raises:
        exception.IllegalArgumentError: If a_line_edit is None.
        exception.NotMainThreadError: If function is called not from the main thread.
    """
    # <editor-fold desc="Checks">
    if a_line_edit is None:
      logger.error('a_line_edit is None.')
      raise exception.IllegalArgumentError('a_line_edit is None.')
    if not thread_util.is_main_thread():
      raise exception.NotMainThreadError()

    # </editor-fold>

    self.line_edit = a_line_edit

  def validate_input_for_sequence_name(
      self, the_current_entered_text: str, the_current_sequences: set
  ) -> tuple[bool, str]:
    """Validates the input for a protein name.

    Args:
        the_current_entered_text (str): The text entered for the sequence name.
        the_current_sequences (set): A set of existing sequence names.

    Returns:
        A tuple containing a boolean indicating whether the input is valid and a string representing an error message if applicable.

    Raises:
        exception.IllegalArgumentError: If either `the_current_entered_text` or `the_current_sequences` is None.
        exception.NotMainThreadError: If function is called not from the main thread.
    """
    # <editor-fold desc="Checks">
    if the_current_entered_text is None:
      logger.error('the_current_entered_text is None.')
      raise exception.IllegalArgumentError('the_current_entered_text is None.')
    if the_current_sequences is None:
      logger.error('the_current_sequences is None.')
      raise exception.IllegalArgumentError('the_current_sequences is None.')
    if not thread_util.is_main_thread():
      raise exception.NotMainThreadError()

    # </editor-fold>

    allowed_chars = {
        '0',
        '1',
        '2',
        '3',
        '4',
        '5',
        '6',
        '7',
        '8',
        '9',
        'a',
        'b',
        'c',
        'd',
        'e',
        'f',
        'g',
        'h',
        'i',
        'j',
        'k',
        'l',
        'm',
        'n',
        'o',
        'p',
        'q',
        'r',
        's',
        't',
        'u',
        'v',
        'w',
        'x',
        'y',
        'z',
        'A',
        'B',
        'C',
        'D',
        'E',
        'F',
        'G',
        'H',
        'I',
        'J',
        'K',
        'L',
        'M',
        'N',
        'O',
        'P',
        'Q',
        'R',
        'S',
        'T',
        'U',
        'V',
        'W',
        'X',
        'Y',
        'Z',
        '-',
        '_',
    }
    for char in the_current_entered_text:
      if char not in allowed_chars:
        self.line_edit.setStyleSheet(
            """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}""",
        )
        return False, 'Invalid character!'
    if the_current_entered_text in the_current_sequences:
      self.line_edit.setStyleSheet(
          """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}""",
      )
      return False, 'Sequence name already exists!'
    if the_current_entered_text == '':
      self.line_edit.setStyleSheet(
          """QLineEdit {color: #ba1a1a; border-color: #ba1a1a;}""",
      )
      return False, 'Please enter a sequence name!'
    self.line_edit.setStyleSheet(
        """QLineEdit {color: #000000; border-color: #DCDBE3;}""",
    )
    return True, ''

  def validate_input_for_protein_sequence(
      self, the_current_entered_text: str
  ) -> tuple[bool, str]:
    """Validates the input for a protein sequence.

    Args:
        the_current_entered_text (str): The text entered by the user.

    Returns:
        A tuple containing a boolean value indicating whether the input is valid and a string message.
        If the input is valid, the boolean value will be True and the message will be an empty string.
        If the input is invalid, the boolean value will be False and the message will describe the error.

    Raises:
        exception.IllegalArgumentError: If the_current_entered_text is None.
        exception.NotMainThreadError: If function is called not from the main thread.
    """
    # <editor-fold desc="Checks">
    if the_current_entered_text is None:
      logger.error('the_current_entered_text is None.')
      raise exception.IllegalArgumentError('the_current_entered_text is None.')
    if not thread_util.is_main_thread():
      raise exception.NotMainThreadError()

    # </editor-fold>

    allowed_chars = {
        'C',
        'D',
        'S',
        'Q',
        'K',
        'I',
        'P',
        'T',
        'F',
        'N',
        'G',
        'H',
        'L',
        'R',
        'W',
        'A',
        'V',
        'E',
        'Y',
        'M',
        ',',
    }
    for char in the_current_entered_text:
      if char not in allowed_chars:
        self.line_edit.setStyleSheet(
            """QTextEdit {color: #ba1a1a; border-color: #ba1a1a;}""",
        )
        return False, 'Invalid character!'
    if the_current_entered_text == '':
      self.line_edit.setStyleSheet(
          """QTextEdit {color: #ba1a1a; border-color: #ba1a1a;}""",
      )
      return False, 'Please enter a protein sequence!'
    elif (
        the_current_entered_text[-1] == ','
        and the_current_entered_text[-2] == ','
    ):
      self.line_edit.setStyleSheet(
          """QTextEdit {color: #ba1a1a; border-color: #ba1a1a;}""",
      )
      return False, 'You already entered a protein sequence delimter!'
    elif self._detect_multiple_commas(the_current_entered_text) > 0:
      self.line_edit.setStyleSheet(
          """QTextEdit {color: #ba1a1a; border-color: #ba1a1a;}""",
      )
      return False, 'You already entered a protein sequence delimter!'
    else:
      self.line_edit.setStyleSheet(
          """QTextEdit {color: #000000; border-color: #DCDBE3;}""",
      )
      return True, ''

  def _detect_multiple_commas(self, text: str) -> int:
    """Detects the number of commas in the text.

    Args:
        text (str): The text to search for multiple commas.

    Returns:
        The number of occurrences of multiple commas in the text.

    Raises:
        exception.IllegalArgumentError: If text is None.
    """
    # <editor-fold desc="Checks">
    if text is None:
      logger.error('text is None.')
      raise exception.IllegalArgumentError('text is None.')

    # </editor-fold>

    pattern = r',{2,}'
    matches = re.findall(pattern, text)
    return len(matches)

  @staticmethod
  def validate_search_input(
      list_for_projects: QtCore.QStringListModel,
      txt_for_search: QtWidgets.QLineEdit,
      lbl_for_status_search: QtWidgets.QLabel,
      txt_for_selected_project: QtWidgets.QLineEdit = None,
      status_message: str = None,
  ) -> None:
    """Validates the input of the project name in real-time.

    Args:
        list_for_projects (QtCore.QStringListModel): A list widget where all projects from the workspace are stored.
        txt_for_search (QtWidgets.QLineEdit): A line edit widget which is used for entering the search term.
        lbl_for_status_search (QtWidgets.QLabel): A label which gives feedback.
        txt_for_selected_project (QtWidgets.QLineEdit): A line edit widget which is used to display the selected project.
        status_message (str): The message which should get displayed. Default is None.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
        exception.NotMainThreadError: If function is called not from the main thread.
    """
    # <editor-fold desc="Checks">
    if list_for_projects is None:
      logger.error('list_for_projects is None.')
      raise exception.IllegalArgumentError('list_for_projects is None.')
    if txt_for_search is None:
      logger.error('txt_for_search is None.')
      raise exception.IllegalArgumentError('txt_for_search is None.')
    if lbl_for_status_search is None:
      logger.error('lbl_for_status_search is None.')
      raise exception.IllegalArgumentError('lbl_for_status_search is None.')
    if txt_for_selected_project is None:
      logger.error('txt_for_selected_project is None.')
      raise exception.IllegalArgumentError('txt_for_selected_project is None.')
    if status_message is None:
      logger.error('status_message is None.')
      raise exception.IllegalArgumentError('status_message is None.')
    if not thread_util.is_main_thread():
      raise exception.NotMainThreadError()

    # </editor-fold>

    list_model = list_for_projects
    txt_for_search = txt_for_search
    lbl_for_status_search = lbl_for_status_search
    txt_for_selected_project = txt_for_selected_project

    # Reset styles for lineEdit
    txt_for_search.setStyleSheet('background-color: white; color: black')

    if len(txt_for_search.text()) == 0:
      lbl_for_status_search.setText('')
      if txt_for_selected_project is not None:
        txt_for_selected_project.setText('')
    else:
      string_list = list_model.stringList()
      matching_items = [
          item for item in string_list if txt_for_search.text() in item
      ]

      if matching_items:
        index = string_list.index(matching_items[0])
        list_model.setData(
            list_model.index(index), Qt.Checked, Qt.CheckStateRole
        )
        lbl_for_status_search.setText('')
        if txt_for_selected_project is not None:
          txt_for_selected_project.setText(matching_items[0])
          # fixme: configure color settings for selected project in list_model
          # list_model.setData(txt_for_search.setStyleSheet("background-color: blue"))

      else:
        txt_for_search.setStyleSheet('background-color: white; color: #f44336')
        if status_message is None:
          lbl_for_status_search.setText('Project name does not exist.')
          lbl_for_status_search.setStyleSheet('color: #f44336')
        else:
          lbl_for_status_search.setText(status_message)
        if txt_for_selected_project is not None:
          txt_for_selected_project.setText('')

  @staticmethod
  def validate_protein_name(
      txt_for_protein_name: QtWidgets.QTextEdit,
      lbl_for_status_protein_name: QtWidgets.QLabel,
      btn_next: QtWidgets.QPushButton = None,
  ) -> None:
    """This function validates the input of the protein name in real-time.

    Args:
        txt_for_protein_name (QtWidgets.QTextEdit): A line edit widget which ii used for entering the protein name.
        lbl_for_status_protein_name (QtWidgets.QLabel): A label which gives feedback if the protein name is legal.
        btn_next (QtWidgets.QPushButton): A push button which is used for the next step.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
        exception.NotMainThreadError: If function is called not from the main thread.
    """
    # <editor-fold desc="Checks">
    if txt_for_protein_name is None:
      logger.error('txt_for_protein_name is None.')
      raise exception.IllegalArgumentError('txt_for_protein_name is None.')
    if lbl_for_status_protein_name is None:
      logger.error('lbl_for_status_protein_name is None.')
      raise exception.IllegalArgumentError(
          'lbl_for_status_protein_name is None.'
      )
    if btn_next is None:
      logger.error('btn_next is None.')
      raise exception.IllegalArgumentError('btn_next is None.')
    if not thread_util.is_main_thread():
      raise exception.NotMainThreadError()

    # </editor-fold>

    # set color for lineEdit
    txt_for_protein_name.setStyleSheet('color: #f44336')
    if len(txt_for_protein_name.text()) == 0:
      lbl_for_status_protein_name.setText('')
      if btn_next is not None:
        btn_next.setEnabled(False)
    else:
      regex = Qt.QtCore.QRegularExpression()
      regex.setPattern('(([a-z])|([A-Z])|([0-9])|(-)|(_)){0,20}')
      validator = QtGui.QRegularExpressionValidator(regex)
      for i in range(len(txt_for_protein_name.text())):
        result = validator.validate(txt_for_protein_name.text(), i)
        if result[0] > 0:
          txt_for_protein_name.setStyleSheet('color: #000000')
          lbl_for_status_protein_name.setText('')
          if btn_next is not None:
            btn_next.setEnabled(True)
        else:
          txt_for_protein_name.setStyleSheet('color: #f44336')
          lbl_for_status_protein_name.setText('Invalid character.')
          if btn_next is not None:
            btn_next.setEnabled(False)
          return
