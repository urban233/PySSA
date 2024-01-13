import logging
import os
import pathlib
import sys

from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt
from PyQt5 import QtGui

from pyssa.gui.ui.dialogs import dialog_startup
from pyssa.internal.data_structures import project, settings
from pyssa.internal.thread import tasks
from pyssa.io_pyssa import safeguard
from pyssa.logging_pyssa import log_handlers
from pyssa.presenter import main_presenter_async
from pyssa.util import constants, enums, exception, main_window_util
from pyssa.gui.ui.views import main_view
from pyssa.gui.ui.views import open_project_view
from pyssa.model import application_model
from pyssa.controller import interface_manager


logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class MainViewController:
    """Class for main presenter of the pyssa plugin."""

    """
    The main view of the pyssa plugin.
    """
    _view: "main_view.MainView"

    """
    
    """
    _app_model: "application_model.ApplicationModel"

    """
    The path to the active workspace.
    """
    _workspace_path: pathlib.Path

    """
    The application settings object.
    """
    _application_settings: "settings.Settings"

    """
    A watcher for the project state.
    """
    _interface_manager: "interface_manager.InterfaceManager"

    """
    The active task of the application.
    """
    _active_task: tasks.Task

    def __init__(self, a_view: "main_view.MainView") -> None:
        """Constructor.

        Args:
            a_view: the main view of the pyssa plugin.

        Raises:
            IllegalArgumentException: If an argument is illegal.
        """
        # <editor-fold desc="Checks">
        safeguard.Safeguard.check_if_value_is_not_none(a_view, logger)

        # </editor-fold>

        self._view: "main_view.MainView" = a_view
        self._app_model: "application_model.ApplicationModel" = application_model.ApplicationModel(project.Project())
        self._external_view = None
        self._protein_model = QtGui.QStandardItemModel()
        self._workspace_path = constants.DEFAULT_WORKSPACE_PATH
        self._workspace_status = f"Current workspace: {str(self._workspace_path)}"
        self._workspace_label = QtWidgets.QLabel(f"Current Workspace: {self._workspace_path}")

        # Configure settings
        # <editor-fold desc="Setup App Settings">
        self._application_settings = settings.Settings(constants.SETTINGS_DIR, constants.SETTINGS_FILENAME)
        if not os.path.exists(constants.SETTINGS_FULL_FILEPATH):
            constants.PYSSA_LOGGER.info("Settings file not found, open configuration dialog.")
            # Configuration dialog to setup setting file
            dialog = dialog_startup.DialogStartup()
            dialog.exec_()

            # checks if the cancel button was pressed
            if dialog_startup.global_var_terminate_app == 1:
                os.remove(constants.SETTINGS_FULL_FILEPATH)
                constants.PYSSA_LOGGER.info("Configuration dialog closed, and removed new settings file.")
                sys.exit()

            self._application_settings.app_launch = 1
            self._application_settings.workspace_path = pathlib.Path(dialog_startup.global_var_startup_workspace)

            constants.PYSSA_LOGGER.info("Demo projects are getting downloaded and extracted ...")
            import zipfile

            with zipfile.ZipFile(pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects.zip"), "r") as zip_ref:
                zip_ref.extractall(pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects"))
            constants.PYSSA_LOGGER.info(
                "Demo projects are downloaded and extracted.\n Import of demo projects started ...",
            )

            path_of_demo_projects = pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects")
            tmp_project = project.Project("", dialog_startup.global_var_startup_workspace)
            for tmp_filename in os.listdir(path_of_demo_projects):
                try:
                    tmp_project = tmp_project.deserialize_project(
                        pathlib.Path(f"{path_of_demo_projects}/{tmp_filename}"),
                        self._application_settings,
                    )
                except exception.IllegalArgumentError:
                    constants.PYSSA_LOGGER.warning(
                        "The workspace path does not exist on this system, " "but this is due to the demo projects.",
                    )
                tmp_project.set_workspace_path(dialog_startup.global_var_startup_workspace)
                new_filepath = pathlib.Path(f"{dialog_startup.global_var_startup_workspace}/{tmp_filename}")
                tmp_project.serialize_project(new_filepath)
            constants.PYSSA_LOGGER.info("Import process of demo projects finished.")
            try:
                os.remove(pathlib.Path(f"{constants.SETTINGS_DIR}/demo-projects.zip"))
            except FileNotFoundError:
                constants.PYSSA_LOGGER.warning("Zip archive of demo projects could not be found!")
            constants.PYSSA_LOGGER.info("Serialize settings ...")
            self._application_settings.serialize_settings()
            constants.PYSSA_LOGGER.info("Serialize settings finished.")

            QtWidgets.QApplication.restoreOverrideCursor()
        self._application_settings = main_window_util.setup_app_settings(self._application_settings)
        # </editor-fold>

        self._interface_manager: "interface_manager.InterfaceManager"
        self._setup_statusbar()
        self._connect_all_ui_elements_with_slot_functions()

    def _setup_statusbar(self) -> None:
        """Sets up the status bar and fills it with the current workspace."""
        self._view.setStatusBar(self._view.status_bar)
        self._interface_manager = interface_manager.InterfaceManager(project.Project(),
                                                                     self._application_settings)

    def _connect_all_ui_elements_with_slot_functions(self):
        self._view.ui.action_open_project.triggered.connect(self._open_project)
        self._view.ui.action_close_project.triggered.connect(self._close_project)
        self._view.ui.proteins_tree_view.clicked.connect(self._show_protein_information)

    def _close_project(self):
        """Closes the current project"""
        self._interface_manager.set_new_project(project.Project())
        self._interface_manager.refresh_main_view(self._view)

    def _open_project(self) -> None:
        self._external_view = open_project_view.OpenProjectView()
        self._external_view.return_value.connect(self._post_open_project)
        self._external_view.show()

    def _post_open_project(self, return_value: str):
        self._interface_manager.start_wait_spinner(self._view)
        tmp_filepath = pathlib.Path(f"{return_value}.xml")
        self._active_task = tasks.Task(
            target=main_presenter_async.open_project,
            args=(
                self._workspace_path,
                tmp_filepath,
                self._application_settings,
            ),
            post_func=self.__await_open_project,
        )
        self._active_task.start()
        self._interface_manager.update_status_bar("Opening existing project ...", self._view)

    def __await_open_project(self, a_result: tuple) -> None:
        self._interface_manager.set_new_project(a_result[1])
        self._interface_manager.update_status_bar(self._workspace_status, self._view)
        self._interface_manager.refresh_main_view(self._view)
        self._interface_manager.stop_wait_spinner(self._view)

    def _show_protein_information(self) -> None:
        tmp_type = self._view.ui.proteins_tree_view.model().data(
            self._view.ui.proteins_tree_view.currentIndex(), enums.ModelEnum.TYPE_ROLE
        )
        if tmp_type == "protein":
            tmp_first_chain = self._view.ui.proteins_tree_view.currentIndex().child(0, 0)
            self._interface_manager.show_chain_pymol_parameters(tmp_first_chain, self._view)
        elif tmp_type == "chain":
            self._interface_manager.show_chain_pymol_parameters(self._protein_model.itemFromIndex(self._view.ui.proteins_tree_view.currentIndex()), self._view)
        else:
            raise ValueError("Unknown type!")
        self._view.status_bar.showMessage(f"Active protein structure: {tmp_type}")

    def update_status(self, message: str) -> None:
        """Updates the status bar of the main view with a custom message."""
        self._view.status_bar.showMessage(message)
