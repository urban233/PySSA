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
"""Module for the status bar manager."""
import logging

from src.pyssa.gui import main_window
from src.pyssa.gui.qt import QtCore
from src.pyssa.gui.qt import QtWidgets
from src.pyssa.gui.ui.custom_widgets import custom_label
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import constants
from src.pyssa.util import exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class StatusBarManager:
  """A class to manage the statusbar style and messages."""

  def __init__(self, the_main_view: main_window.MainWindow) -> None:
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

    self._progress_bar = GradientProgressBar()
    # self._progress_bar.setRange(0, 100)
    # self._progress_bar.setValue(0)
    self._permanent_message = custom_label.PermanentMessageLabel()

    self._view.status_bar.addPermanentWidget(self._progress_bar)
    self._view.status_bar.addPermanentWidget(self._permanent_message)
    self._progress_bar.setMaximumWidth(100)
    self._progress_bar.hide()
    self._permanent_message.hide()
    self.temp_message_timer = QtCore.QTimer()
    self._restore_status_bar()

  # <editor-fold desc="Util methods">

  def _setup_progress_bar_animation(self):
    self.timer = QtCore.QTimer()
    self.timer.timeout.connect(self.update_progress)
    self.timer.start(100)  # 100 ms per step (slower than default)

  def update_progress(self):
    value = (self._progress_bar.value() + 2) % 101  # loop from 0→100
    self._progress_bar.setValue(value)

  # <editor-fold desc="Methods for styling the status bar">
  def _style_progress_bar(self):
    self._progress_bar.setStyleSheet(
      """
      QProgressBar {
          border-style: solid;
          border-width: 2px;
          border-radius: 4px;
          border-color: #DCDBE3;
          background-color: #efefef;
          max-height: 17px;
          max-width: 80px;
          text-align: center;
          color: black;
      }
      QProgressBar::chunk {
          background-color: #367af6;
          width: 50px;
      }
      """
    )

  def _style_status_bar_for_normal_message(self) -> None:
    """Sets custom style sheet for a normal message."""
    self._view.status_bar.setStyleSheet(
        """
            QStatusBar {
                background-color: #eeeff0;
                min-height: 1.3em;
                max-height: 1.3em;
            }
        """
    )

  def _style_status_bar_for_long_running_task_message(self) -> None:
    """Sets custom style sheet for a long-running message."""
    self._view.status_bar.setStyleSheet(
        """
            QStatusBar {
                background-color: #ff9000;
            }
        """
    )

  def _style_status_bar_for_error_message(self) -> None:
    """Sets custom style sheet for an error message."""
    self._view.status_bar.setStyleSheet(
        """
            QStatusBar {
                background-color: #ff9000;
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
  def show_permanent_message(
          self,
          a_message: str,
          a_with_progress_bar: bool = False
  ) -> None:
    """Shows a permanent message in the statusbar.

    Args:
        a_message: A string representing the message that will be displayed as a permanent message.
        a_with_progress_bar: A boolean flag for whether to show a progress bar along with the message.

    Raises:
        exception.IllegalArgumentError: If `a_message` is None.
    """
    # <editor-fold desc="Checks">
    if a_message is None:
      logger.error("a_message is None.")
      raise exception.IllegalArgumentError("a_message is None.")

    # </editor-fold>

    if a_message == "":
      self._permanent_message.hide()
    else:
      self._permanent_message.show()
      self._permanent_message.setText(a_message)

    if a_with_progress_bar:
      self._progress_bar.show()
    else:
      self._progress_bar.hide()

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
    self._progress_bar.hide()

  def show_temporary_message(
      self,
      a_temporary_message: str,
      a_with_timeout_flag: bool = True,
      a_timeout: int = constants.STATUS_MESSAGE_TIMEOUT
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

  def hide_progress_bar(self) -> None:
    """Hides the progress bar and reset the permanent message."""
    self._progress_bar.hide()
    self._permanent_message.hide()
    self._permanent_message.setText("")

  # </editor-fold>


class GradientProgressBar(QtWidgets.QProgressBar):
  def __init__(self):
    super().__init__()
    self.setRange(0, 100)
    self.setValue(0)
    self.setTextVisible(False)
    self.setFixedWidth(180)
    self.setFixedHeight(10)
    self.setStyleSheet(self.get_stylesheet(0))

    self.offset = 0
    self.timer = QtCore.QTimer()
    self.timer.timeout.connect(self.update_gradient)
    self.timer.start(50)  # slower movement: increase interval to slow down

  def update_gradient(self):
    # Move the gradient slowly
    self.offset = (self.offset + 2) % 100
    self.setStyleSheet(self.get_stylesheet(self.offset))
    # Increment value for continuous bar fill illusion
    self.setValue((self.value() + 1) % 101)

  def get_stylesheet(self, offset):
    # Use offset to shift gradient position
    return f"""
        QProgressBar {{
            border: 1px solid #ccc;
            border-radius: 6px;
            background-color: #f0f0f0;
        }}
        QProgressBar::chunk {{
            border-radius: 6px;
            background: qlineargradient(
                x1:0, y1:0, x2:1, y2:0,
                stop:0.0 rgba(90, 173, 226, 255),
                stop:{0.3 + offset/200:.2f} rgba(150, 200, 255, 255),
                stop:{0.7 + offset/200:.2f} rgba(90, 173, 226, 255),
                stop:1.0 rgba(90, 173, 226, 255)
            );
        }}
        """
