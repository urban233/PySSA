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
"""Module for the status bar manager."""
import logging

from PyQt5 import QtCore
from PyQt5 import QtWidgets
from pyssa.gui.ui.custom_widgets import custom_label
from pyssa.logging_pyssa import log_handlers
from pyssa.util import constants
from pyssa.util import exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class StatusBarManager:
  """A class to manage the statusbar style and messages."""

  def __init__(self, the_main_view: QtWidgets.QMainWindow) -> None:
    """Constructor.

    Args:
        the_main_view (QMainWindow): The main view of the application.

    Raises:
        exception.IllegalArgumentError: If `the_main_view` is None.
    """
    # <editor-fold desc="Checks">
    if the_main_view is None:
      logger.error("the_main_view is None.")
      raise exception.IllegalArgumentError("the_main_view is None.")

    # </editor-fold>

    self._view = the_main_view
    self._update_signal = None

    self._progress_bar = QtWidgets.QProgressBar()
    self._permanent_message = custom_label.PermanentMessageLabel()

    self._menu_task = QtWidgets.QMenu()
    self._is_menu_open = False
    self._abort_action = QtWidgets.QAction("Abort Job")
    self._menu_task.addAction(self._abort_action)

    self._view.status_bar.addPermanentWidget(self._progress_bar)
    self._view.status_bar.addPermanentWidget(self._permanent_message)
    self._view.status_bar.addPermanentWidget(self._view.btn_open_job_overview)
    self._view.status_bar.addPermanentWidget(
        self._view.btn_open_job_notification
    )
    self._progress_bar.hide()
    self._view.btn_open_job_overview.show()
    self._view.btn_open_job_notification.show()
    self.temp_message_timer = QtCore.QTimer()

    # self._connect_ui_elements()

  # <editor-fold desc="Util methods">

  # <editor-fold desc="Methods for styling the status bar">
  def _style_status_bar_for_normal_message(self) -> None:
    """Sets custom style sheet for a normal message."""
    self._view.status_bar.setStyleSheet(
        """
            QStatusBar {
                background-color: #F2F2F2;
                border-style: solid;
                border-width: 2px;
                border-radius: 4px;
                border-color: #DCDBE3;
            }
        """
    )

  def _style_status_bar_for_long_running_task_message(self) -> None:
    """Sets custom style sheet for a long-running message."""
    self._view.status_bar.setStyleSheet(
        """
            QStatusBar {
                background-color: #ff9000;
                border-style: solid;
                border-width: 2px;
                border-radius: 4px;
                border-color: #5b5b5b;
            }
        """
    )

  def _style_status_bar_for_error_message(self) -> None:
    """Sets custom style sheet for an error message."""
    self._view.status_bar.setStyleSheet(
        """
            QStatusBar {
                background-color: #ff9000;
                border-style: solid;
                border-width: 2px;
                border-radius: 4px;
                border-color: #5b5b5b;
            }
        """
    )

  # </editor-fold>

  def _setup_status_bar_message_timer(
      self, running_task: bool = False, the_long_running_task_message: str = ""
  ) -> None:
    """Connects the timer to reset the status bar to the long-running task message.

    Args:
        running_task (bool): Flag for indicating if a long-running task is currently running.
        the_long_running_task_message (str): Message to be displayed when a long-running task is running.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if running_task is None:
      logger.error("running_task is None.")
      raise exception.IllegalArgumentError("running_task is None.")
    if the_long_running_task_message is None:
      logger.error("the_long_running_task_message is None.")
      raise exception.IllegalArgumentError(
          "the_long_running_task_message is None."
      )

    # </editor-fold>

    if self.temp_message_timer:
      self.temp_message_timer.stop()  # Stop previous timer if exists
    self.temp_message_timer.setSingleShot(True)
    if running_task:
      self.temp_message_timer.timeout.connect(
          lambda a_long_running_task_message=the_long_running_task_message: self._switch_to_long_running_task_message(
              a_long_running_task_message
          )
      )
    else:
      self.temp_message_timer.timeout.connect(self._restore_status_bar)
    self.temp_message_timer.start(
        5000
    )  # Display temporary message for 5 seconds

  def _switch_to_long_running_task_message(
      self, a_long_running_task_message: str
  ) -> None:
    """Shows a long-running task message as a permanent message."""
    self.show_permanent_message(a_long_running_task_message)

  def _restore_status_bar(self) -> None:
    """Restores the statusbar."""
    self._style_status_bar_for_normal_message()
    self._view.status_bar.showMessage("")

  # </editor-fold>

  def _connect_ui_elements(self) -> None:
    """Connects all UI elements to their corresponding slot functions in the class."""
    # self._permanent_message.textChanged.connect(self._manage_status_bar_ui)
    raise NotImplementedError()

  # <editor-fold desc="Public methods">
  def show_permanent_message(self, a_message: str) -> None:
    """Shows a permanent message in the statusbar.

    Args:
        a_message: A string representing the message that will be displayed as a permanent message.

    Raises:
        exception.IllegalArgumentError: If `a_message` is None.
    """
    # <editor-fold desc="Checks">
    if a_message is None:
      logger.error("a_message is None.")
      raise exception.IllegalArgumentError("a_message is None.")

    # </editor-fold>

    self._permanent_message.setText(a_message)

  def show_error_message(
      self, a_message: str, overwrite_permanent_message: bool = True
  ) -> None:
    """Shows an error message in the statusbar.

    Args:
        a_message (str): The error message to be displayed.
        overwrite_permanent_message (bool, optional): Flag indicating whether to overwrite the permanent message. Defaults to True.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_message is None:
      logger.error("a_message is None.")
      raise exception.IllegalArgumentError("a_message is None.")
    if overwrite_permanent_message is None:
      logger.error("overwrite_permanent_message is None.")
      raise exception.IllegalArgumentError(
          "overwrite_permanent_message is None."
      )

    # </editor-fold>

    self._style_status_bar_for_error_message()
    self._view.status_bar.showMessage("")
    if overwrite_permanent_message is True:
      self._permanent_message.show()
      self._permanent_message.setText(a_message)
    else:
      self._view.status_bar.showMessage(a_message, 999999)

  def show_temporary_message(
      self,
      a_temporary_message: str,
      a_with_timeout_flag: bool = True,
      a_timeout: int = constants.STATUS_MESSAGE_TIMEOUT,
  ) -> None:
    """Shows a temporary message in the statusbar.

    Args:
        a_temporary_message (str): The message to be displayed temporarily in the status bar.
        a_with_timeout_flag (bool): Optional parameter that specifies whether the message should be displayed for a limited time. Defaults to True.
        a_timeout (int): Optional parameter that specifies the amount of time (in milliseconds) the message should be displayed if a_with_timeout_flag is set to True. Defaults to the value of constants.STATUS_MESSAGE_TIMEOUT.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if a_temporary_message is None:
      logger.error("a_temporary_message is None.")
      raise exception.IllegalArgumentError("a_temporary_message is None.")
    if a_with_timeout_flag is None:
      logger.error("a_with_timeout_flag is None.")
      raise exception.IllegalArgumentError("a_with_timeout_flag is None.")
    if a_timeout is None:
      logger.error("a_timeout is None.")
      raise exception.IllegalArgumentError("a_timeout is None.")

    # </editor-fold>

    self._style_status_bar_for_normal_message()
    self._permanent_message.setText("")
    if a_with_timeout_flag:
      self._view.status_bar.showMessage(a_temporary_message, a_timeout)
    else:
      self._view.status_bar.showMessage(a_temporary_message, 999999)

  def update_progress_bar(self, a_message_value_tuple: tuple) -> None:
    """Updates the progress bar with the given message and value.

    Args:
        a_message_value_tuple (tuple): A tuple containing the message and value to be displayed on the progress bar. The message should be a string, and the value should be an integer between 0 and 100 (inclusive).

    Raises:
        exception.IllegalArgumentError: If `a_message_value_tuple` is None.
        ValueError: If the value is less than 0 or greater than 100.
    """
    # <editor-fold desc="Checks">
    if a_message_value_tuple is None:
      logger.error("a_message_value_tuple is None.")
      raise exception.IllegalArgumentError("a_message_value_tuple is None.")
    tmp_message, tmp_value = a_message_value_tuple
    if tmp_value < 0 or tmp_value > 100:
      raise ValueError("Value for progress bar must be between 0 and 100!")

    # </editor-fold>

    self._progress_bar.show()
    self._progress_bar.setFormat(f"{tmp_value}%")
    self._progress_bar.setValue(tmp_value)
    self._permanent_message.show()
    self._permanent_message.setText(tmp_message)

  def hide_progress_bar(self) -> None:
    """Hides the progress bar and reset the permanent message."""
    self._progress_bar.hide()
    self._permanent_message.hide()
    self._permanent_message.setText("")

  # </editor-fold>
