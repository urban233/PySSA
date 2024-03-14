from PyQt5 import QtCore
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from pyssa.controller import main_tasks_manager
from pyssa.util import enums


class StatusBarManager:
    """A class to manage the statusbar style and messages."""
    def __init__(self, the_main_view):
        self._view = the_main_view
        self._progress_bar = QtWidgets.QProgressBar()
        self._permanent_message = QtWidgets.QLabel()
        self._btn_task = QtWidgets.QPushButton()
        self._menu_task = QtWidgets.QMenu()
        self._abort_action = QtWidgets.QAction("Abort Task")

        add_sequence_icon = QtGui.QIcon(QtGui.QPixmap(":icons/keyboard_arrow_down_w400.svg"))
        self._btn_task.setIcon(add_sequence_icon)
        self._btn_task.setText("")
        self._btn_task.setIconSize(add_sequence_icon.actualSize(QtCore.QSize(30, 30)))
        self._btn_task.setStyleSheet("""
            QPushButton {
                background-color: white;
                border: none;
                border-width: 2px;
                border-radius: 10px;
                padding: 2px;
                min-width: 20px;
                max-width: 20px;
                min-height: 20px;
                max-height: 20px
            }
            QPushButton::hover {
                background: qlineargradient(x1: 0, y1: 0, x2: 0, y2: 1,
                                            stop: 0 #4B91F7, stop: 0.4 #367AF6,
                                            stop: 0.5 #367AF6, stop: 1.0 #4B91F7);
                background: white;
                color: white;
                color: #4B91F7;
                border: 2px solid #DCDBE3;
            }
        """)
        self._view.status_bar.addPermanentWidget(self._progress_bar)
        self._view.status_bar.addPermanentWidget(self._permanent_message)
        self._view.status_bar.addPermanentWidget(self._btn_task)
        self._progress_bar.hide()
        self._btn_task.hide()
        self.temp_message_timer = QtCore.QTimer()

        self._connect_ui_elements()

    def _connect_ui_elements(self):
        self._btn_task.clicked.connect(self.__slot_open_menu)

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
                border-color: #DCDBE3;
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
        #self._style_status_bar_for_long_running_task_message()
        #self._view.status_bar.showMessage(a_message)
        self._permanent_message.setText(a_message)

    def show_error_message(self, a_message):
        self._style_status_bar_for_error_message()
        self._permanent_message.setText(a_message)
        #self._view.status_bar.showMessage(a_message)

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

    def update_progress_bar(self, a_message_value_tuple: tuple):
        tmp_message, tmp_value = a_message_value_tuple
        if tmp_value < 0 or tmp_value > 100:
            raise ValueError("Value for progress bar must be between 0 and 100!")
        # self._view.progress_bar.show()
        # self._view.progress_bar.setFormat(f"{tmp_value}%")
        # self._view.progress_bar.setValue(tmp_value)

        self._progress_bar.show()
        self._progress_bar.setFormat(f"{tmp_value}%")
        self._progress_bar.setValue(tmp_value)
        self._permanent_message.setText(tmp_message)
        self._btn_task.show()

    def __slot_open_menu(self):
        self._menu_task.addAction(self._abort_action)
        self._menu_task.exec_(self._btn_task.mapToGlobal(QtCore.QPoint(0, -25)))
