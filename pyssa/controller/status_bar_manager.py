from PyQt5 import QtCore

from pyssa.controller import main_tasks_manager
from pyssa.util import enums


class StatusBarManager:
    """A class to manage the statusbar style and messages."""
    def __init__(self, the_main_view):
        self._view = the_main_view
        self.temp_message_timer = QtCore.QTimer()

    # <editor-fold desc="Util methods">

    # <editor-fold desc="Methods for styling the status bar">
    def _style_status_bar_for_normal_message(self):
        self._view.status_bar.setStyleSheet("""
            QStatusBar {
                background-color: white;
                border-style: solid;
                border-width: 2px;
                border-radius: 4px;
                border-color: #DCDBE3;
            }
        """)

    def _style_status_bar_for_long_running_task_message(self):
        self._view.status_bar.setStyleSheet("""
            QStatusBar {
                background-color: #ff9000;
                border-style: solid;
                border-width: 2px;
                border-radius: 4px;
                border-color: #5b5b5b;
            }
        """)

    def _style_status_bar_for_error_message(self):
        self._view.status_bar.setStyleSheet("""
            QStatusBar {
                background-color: #ba1a1a;
                color: white;
                border-style: solid;
                border-width: 2px;
                border-radius: 4px;
                border-color: #5b5b5b;
            }
        """)

    # </editor-fold>

    def _setup_status_bar_message_timer(self, running_task=False, the_long_running_task_message=""):
        """Connects the timer to reset the status bar to the long-runnnig task message."""
        if self.temp_message_timer:
            self.temp_message_timer.stop()  # Stop previous timer if exists
        self.temp_message_timer.setSingleShot(True)
        if running_task:
            self.temp_message_timer.timeout.connect(lambda a_long_running_task_message=the_long_running_task_message: self._switch_to_long_running_task_message(a_long_running_task_message))
        else:
            self.temp_message_timer.timeout.connect(self._restore_status_bar)
        self.temp_message_timer.start(5000)  # Display temporary message for 5 seconds

    def _switch_to_long_running_task_message(self, a_long_running_task_message):
        self.show_long_running_task_message(a_long_running_task_message)

    def _restore_status_bar(self):
        self._style_status_bar_for_normal_message()
        self._view.status_bar.showMessage("")

    # </editor-fold>

    def show_long_running_task_message(self, a_message):
        self._style_status_bar_for_long_running_task_message()
        self._view.status_bar.showMessage(a_message)

    def show_error_message(self, a_message):
        self._style_status_bar_for_error_message()
        self._view.status_bar.showMessage(a_message)

    def show_temporary_message(self, a_temporary_message, the_main_tasks_manager: "main_tasks_manager.MainTasksManager"):
        self._style_status_bar_for_normal_message()
        self._view.status_bar.showMessage(a_temporary_message)

        if the_main_tasks_manager.prediction_task is not None:
            # A prediction task might be running
            if the_main_tasks_manager.check_if_prediction_task_is_finished() is False:
                # Revert back to long-running tasks message
                self._setup_status_bar_message_timer(
                    running_task=True,
                    the_long_running_task_message=enums.StatusMessages.PREDICTION_IS_RUNNING.value
                )
        elif the_main_tasks_manager.distance_analysis_task is not None:
            # A distance analysis task might be running
            pass
        else:
            self._setup_status_bar_message_timer(running_task=False)

    def _update_progress_bar(self, a_message_value_tuple: tuple):
        tmp_message, tmp_value = a_message_value_tuple
        if tmp_value < 0 or tmp_value > 100:
            raise ValueError("Value for progress bar must be between 0 and 100!")
        self._view.progress_bar.show()
        self._view.progress_bar.setFormat(tmp_message)
        self._view.progress_bar.setValue(tmp_value)
