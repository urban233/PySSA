#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/zielesny/PySSA>
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
"""Module for functions which reduce code duplicates in the main module."""
import logging
import os
import pathlib

from PyQt5 import QtWidgets

from src.pyssa.gui.ui.custom_dialogs import custom_message_box
from src.pyssa.internal.thread import thread_util
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import exception
from src.pyssa.util.void import rvoid

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


def fill_combo_box(combo_box: QtWidgets.QComboBox, item_list: list) -> None:
  """This function fills a pyqt combobox.

  Args:
      combo_box (QtWidgets.QComboBox):
           pyqt combo box which should be filled
      item_list (list):
          list of items which should be placed in the combo box

  Raises:
      exception.IllegalArgumentError: If combo_box is None.
      exception.IllegalArgumentError: If item_list is either None or an empty list.
      exception.NotMainThreadError: If function is called not from the main thread.
  """
  # <editor-fold desc="Checks">
  if combo_box is None:
    logger.error("combo_box is None.")
    raise exception.IllegalArgumentError("combo_box is None.")
  if item_list is None or len(item_list) == 0:
    logger.error("item_list is either None or an empty list.")
    raise exception.IllegalArgumentError(
        "item_list is either None or an empty list."
    )
  if not thread_util.is_main_thread():
    raise exception.NotMainThreadError()

  # </editor-fold>

  combo_box.clear()
  for item in item_list:
    combo_box.addItem(item)


def choose_directory(
    self, txt_box_dir: QtWidgets.QLineEdit
) -> None:  # noqa: ANN001
  """Open s a QFileDialog to choose a filepath.

  Args:
      self: The main_window object.
      txt_box_dir (QtWidgets.QLineEdit): The QLineEdit widget containing the current directory path.

  Raises:
      exception.IllegalArgumentError: If either `self` or `txt_box_dir` is None.
      exception.NotMainThreadError: If function is called not from the main thread.
  """
  # <editor-fold desc="Checks">
  if self is None:
    logger.error("self is None.")
    raise exception.IllegalArgumentError("self is None.")
  if txt_box_dir is None:
    logger.error("txt_box_dir is None.")
    raise exception.IllegalArgumentError("txt_box_dir is None.")
  if not thread_util.is_main_thread():
    raise exception.NotMainThreadError()

  # </editor-fold>

  current_file_path = pathlib.Path(txt_box_dir.text())
  new_file_path = pathlib.Path(
      QtWidgets.QFileDialog.getExistingDirectory(
          self,
          "Open Directory",
          str(current_file_path),
          options=QtWidgets.QFileDialog.ShowDirsOnly,
      ),
  )
  rvoid(os.access(new_file_path, os.W_OK))
  if new_file_path != pathlib.Path("") and os.access(new_file_path, os.W_OK):
    txt_box_dir.setText(str(new_file_path))
  elif new_file_path != pathlib.Path("") and not os.access(
      new_file_path, os.W_OK
  ):
    tmp_dialog = custom_message_box.CustomMessageBoxOk(
        "You do not have write permissions for this directroy. Please choose another one.",
        "Permission Error",
        custom_message_box.CustomMessageBoxIcons.WARNING.value,
    )
    tmp_dialog.exec_()
    txt_box_dir.setText(str(current_file_path))
  else:
    txt_box_dir.setText(str(current_file_path))


def hide_gui_elements(gui_elements: list) -> None:
  """Hides gui elements.

  Args:
      gui_elements (list): A list of pyqt gui elements.

  Raises:
      exception.IllegalArgumentError: If gui_elements is None.
      exception.NotMainThreadError: If function is called not from the main thread.
  """
  # <editor-fold desc="Checks">
  if gui_elements is None:
    logger.error("gui_elements is None.")
    raise exception.IllegalArgumentError("gui_elements is None.")
  if not thread_util.is_main_thread():
    raise exception.NotMainThreadError()

  # </editor-fold>

  for gui_element in gui_elements:
    if gui_element is not None:
      gui_element.hide()


def show_gui_elements(gui_elements: list) -> None:
  """Shows gui elements.

  Args:
      gui_elements (list): A list of pyqt gui elements.

  Raises:
      exception.IllegalArgumentError: If gui_elements is None.
      exception.NotMainThreadError: If function is called not from the main thread.
  """
  # <editor-fold desc="Checks">
  if gui_elements is None:
    logger.error("gui_elements is None.")
    raise exception.IllegalArgumentError("gui_elements is None.")
  if not thread_util.is_main_thread():
    raise exception.NotMainThreadError()

  # </editor-fold>

  for gui_element in gui_elements:
    if gui_element is not None:
      gui_element.show()


def manage_gui_visibility(
    gui_elements_to_show: list, gui_elements_to_hide: list
) -> None:
  """Manages a combination of "show_gui_elements" and "hide_gui_elements" manage the visibility of gui elements.

  Args:
      gui_elements_to_show (list): A list which contains the gui elements which should be displayed.
      gui_elements_to_hide (list): A list which contains the gui elements which should be hidden.

  Raises:
      exception.IllegalArgumentError: If either `gui_elements_to_show` or `gui_elements_to_hide` is None.
      exception.NotMainThreadError: If function is called not from the main thread.
  """
  # <editor-fold desc="Checks">
  if gui_elements_to_show is None:
    logger.error("gui_elements_to_show is None.")
    raise exception.IllegalArgumentError("gui_elements_to_show is None.")
  if gui_elements_to_hide is None:
    logger.error("gui_elements_to_hide is None.")
    raise exception.IllegalArgumentError("gui_elements_to_hide is None.")
  if not thread_util.is_main_thread():
    raise exception.NotMainThreadError()

  # </editor-fold>

  show_gui_elements(gui_elements_to_show)
  hide_gui_elements(gui_elements_to_hide)


def disable_text_box(
    text_box: QtWidgets.QLineEdit, text_box_label: QtWidgets.QLabel
) -> None:
  """Disables a text box and grays out the corresponding label.

  Args:
      text_box (QtWidgets.QLineEdit): A PyQt line edit.
      text_box_label (QtWidgets.QLabel): A PyQt label which describes the text box.

  Raises:
      exception.IllegalArgumentError: If either `text_box` or `text_box_label` is None.
      exception.NotMainThreadError: If function is called not from the main thread.
  """
  # <editor-fold desc="Checks">
  if text_box is None:
    logger.error("text_box is None.")
    raise exception.IllegalArgumentError("text_box is None.")
  if text_box_label is None:
    logger.error("text_box_label is None.")
    raise exception.IllegalArgumentError("text_box_label is None.")
  if not thread_util.is_main_thread():
    raise exception.NotMainThreadError()

  # </editor-fold>

  text_box_label.setStyleSheet("color: #E1E1E1")
  text_box.setStyleSheet("background-color: white")
  text_box.setEnabled(False)


def enable_text_box(
    text_box: QtWidgets.QLineEdit, text_box_label: QtWidgets.QLabel
) -> None:
  """This function enables a text box and colors the corresponding label black.

  Args:
      text_box (QtWidgets.QLineEdit): A PyQt line edit.
      text_box_label (QtWidgets.QLabel): A PyQt label which describes the text box.

  Raises:
      exception.IllegalArgumentError: If either `text_box` or `text_box_label` is None.
      exception.NotMainThreadError: If function is called not from the main thread.
  """
  # <editor-fold desc="Checks">
  if text_box is None:
    logger.error("text_box is None.")
    raise exception.IllegalArgumentError("text_box is None.")
  if text_box_label is None:
    logger.error("text_box_label is None.")
    raise exception.IllegalArgumentError("text_box_label is None.")
  if not thread_util.is_main_thread():
    raise exception.NotMainThreadError()

  # </editor-fold>

  text_box_label.setStyleSheet("color: black")
  text_box.setStyleSheet("background-color: white")
  text_box.setEnabled(True)
