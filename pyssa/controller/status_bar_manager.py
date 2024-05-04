from PyQt5 import QtCore
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from pyssa.controller import main_tasks_manager
from pyssa.gui.ui.custom_widgets import custom_label, job_entry
from pyssa.util import enums, constants


class StatusBarManager:
    """A class to manage the statusbar style and messages."""
    def __init__(self, the_main_view, the_main_task_manager: "main_tasks_manager.MainTasksManager"):
        self._view = the_main_view
        self._main_task_manager = the_main_task_manager
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
        self._view.status_bar.addPermanentWidget(self._view.btn_open_job_notification)
        self._progress_bar.hide()
        self._view.btn_open_job_overview.show()
        self._view.btn_open_job_notification.show()
        self.temp_message_timer = QtCore.QTimer()

        self._connect_ui_elements()

    # <editor-fold desc="Util methods">

    # <editor-fold desc="Methods for styling the status bar">
    def _style_status_bar_for_normal_message(self):
        self._view.status_bar.setStyleSheet("""
            QStatusBar {
                background-color: #F2F2F2;
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
                background-color: #ff9000;
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
            self.temp_message_timer.timeout.connect(lambda
                                                        a_long_running_task_message=the_long_running_task_message: self._switch_to_long_running_task_message(
                a_long_running_task_message))
        else:
            self.temp_message_timer.timeout.connect(self._restore_status_bar)
        self.temp_message_timer.start(5000)  # Display temporary message for 5 seconds

    def _switch_to_long_running_task_message(self, a_long_running_task_message):
        self.show_permanent_message(a_long_running_task_message)

    def _restore_status_bar(self):
        self._style_status_bar_for_normal_message()
        self._view.status_bar.showMessage("")

    # </editor-fold>

    def _connect_ui_elements(self):
        self._permanent_message.textChanged.connect(self._manage_status_bar_ui)

    def _manage_status_bar_ui(self, the_current_text):
        if self._main_task_manager.prediction_task is None:
            #self._btn_task.setEnabled(False)
            return
        if not self._main_task_manager.prediction_task.is_finished():
            #self._btn_task.setEnabled(True)
            pass
        else:
            #self._btn_task.setEnabled(False)
            pass

    # <editor-fold desc="Public methods">
    def show_permanent_message(self, a_message):
        #self._style_status_bar_for_long_running_task_message()
        #self._view.status_bar.showMessage(a_message)
        self._permanent_message.setText(a_message)

    def show_error_message(self, a_message):
        self._style_status_bar_for_error_message()
        if self._permanent_message.text() == "":
            self._permanent_message.show()
            self._permanent_message.setText(a_message)
        else:
            self._view.status_bar.showMessage(a_message, 999999)

    def show_temporary_message(self, a_temporary_message, a_with_timeout_flag: bool = True, a_timeout: int = constants.STATUS_MESSAGE_TIMEOUT):
        self._style_status_bar_for_normal_message()
        self._permanent_message.setText("")
        if a_with_timeout_flag:
            self._view.status_bar.showMessage(a_temporary_message, a_timeout)
        else:
            self._view.status_bar.showMessage(a_temporary_message, 999999)

        # if the_main_tasks_manager.prediction_task is not None:
        #     # A prediction task might be running
        #     if the_main_tasks_manager.check_if_prediction_task_is_finished() is False:
        #         # Revert back to long-running tasks message
        #         self._setup_status_bar_message_timer(
        #             running_task=True,
        #             the_long_running_task_message=enums.StatusMessages.PREDICTION_IS_RUNNING.value
        #         )
        # elif the_main_tasks_manager.distance_analysis_task is not None:
        #     # A distance analysis task might be running
        #     pass
        # else:
        #     self._setup_status_bar_message_timer(running_task=False)

    def update_progress_bar(self, a_message_value_tuple: tuple):
        tmp_message, tmp_value = a_message_value_tuple
        if tmp_value < 0 or tmp_value > 100:
            raise ValueError("Value for progress bar must be between 0 and 100!")

        self._progress_bar.show()
        self._progress_bar.setFormat(f"{tmp_value}%")
        self._progress_bar.setValue(tmp_value)
        self._permanent_message.show()
        self._permanent_message.setText(tmp_message)

    def hide_progress_bar(self):
        self._progress_bar.hide()
        self._permanent_message.hide()
        self._permanent_message.setText("")
    # </editor-fold>
