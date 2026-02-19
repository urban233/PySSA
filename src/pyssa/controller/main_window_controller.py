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
"""Main window controller class for the PySSA frontend application.

Authors: Martin Urban, Hannah Kullik

Version: 2.0.0
"""
import logging
import pathlib
import shutil
from multiprocessing import Pipe, Process

from pyssa.controller import create_project_view_controller, interface_manager
from pyssa.controller import database_manager
from pyssa.controller import open_project_view_controller
from pyssa.controller import use_project_view_controller
from pyssa.controller import delete_project_view_controller
from pyssa.internal.data_structures import project
from pyssa.internal.thread import database_thread
from pyssa.internal.thread.async_pyssa import project_async
from pyssa.logging_pyssa import log_handlers, log_levels
from src.pyssa.controller import status_bar_manager
from src.pyssa.gui.qt import QtCore, QtWidgets, QtGui
from src.pyssa.gui import main_window
# from PySSA.model.qmodel import protein_model
from src.pyssa.internal.pymol import pml_worker, worker_command
from src.pyssa.util import enums
from tea.thread import task_result_factory, task_result, action

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class MainWindowController:
    """Controller class for the main window."""

    def __init__(self, a_main_window: "main_window.MainWindow") -> None:
        """Constructor.

        Args:
          a_main_window: Main window instance
        """
        # <editor-fold desc="Checks">
        # psa_comm_api.macros.REQUIRE_NOT_NONE(a_main_window, "a_main_window is None.")
        # </editor-fold>
        # <editor-fold desc="Instance attributes">
        # <editor-fold desc="Private">
        self._main_window = a_main_window
        self._interface_manager = interface_manager.InterfaceManager(a_main_window)
        # self._pyssa_api = api.CoreAPI.instance()
        # self._context = communication.ZmqContext()
        # self._push_endpoint = communication.PushEndpoint.connect(
        #   self._context, "tcp://127.0.0.1:8070"
        # )
        self._user_pymol = self._main_window.user_pymol
        # self._PySSA_worker = w_worker.WWorker()
        # self._PySSA_worker.start()
        self._pymol_worker_connection, child_conn = Pipe()
        self._pymol_worker_process = Process(
            target=pml_worker.start_pml_worker, args=(child_conn,)
        )
        self._pymol_worker_process.start()
        self.feedback_timer = QtCore.QTimer()
        self._status_bar_manager = status_bar_manager.StatusBarManager(
            self._main_window
        )
        # </editor-fold>
        # </editor-fold>
        # self._init_main_window()
        self._connect_all_signals_with_their_slots()
        # self.aux_pymol_client = aux_pymol_client.AuxPyMOLClient(self.context)
        # self.protein_model = protein_model.ProteinModel()
        # self.protein_model.add_protein(self._user_pymol.get_cmd_module().get_model())
        self._init_main_window()

    # <editor-fold desc="Private methods">
    def _init_main_window(self) -> None:
        """Sets up the main window."""
        self._set_models_to_views()

    def _set_models_to_views(self) -> None:
        """Sets all models to their corresponding views."""
        # self._main_window.side_panel_molecule_objects.tree_view.setModel(self.protein_model)
        pass

    def _connect_all_signals_with_their_slots(self) -> None:
        """Connects all relevant widget signals with their appropriate slots."""
        # self._main_window.dialogClosed.connect(self.__slot_close_application)

        self._main_window.action_new_project.triggered.connect(self.__slot_create_project)


        # # <editor-fold desc="Session ribbon slots">
        # # <editor-fold desc="Session slots">
        self._main_window.viewer_toolbar_actions.get("open_session").get_action().triggered.connect(
            self.__slot_open_session
        )
        # # </editor-fold>
        # # <editor-fold desc="Scene slots">
        self._main_window.viewer_toolbar_actions.get("create_scene").get_action().triggered.connect(
            self.__slot_save_scene
        )
        self._main_window.viewer_toolbar_actions.get("save_scene").get_action().triggered.connect(
            self.__slot_update_scene
        )
        self._main_window.viewer_toolbar_actions.get("delete_scene").get_action().triggered.connect(
            self.__slot_delete_scene
        )
        # # </editor-fold>
        # # <editor-fold desc="Representation slots">
        self._main_window.viewer_toolbar_actions.get("cartoon").get_action().triggered.connect(
            self.__slot_display_as_cartoon
        )
        self._main_window.cartoon_show_action.triggered.connect(self.__slot_show_as_cartoon)
        self._main_window.cartoon_hide_action.triggered.connect(self.__slot_hide_cartoon)

        self._main_window.viewer_toolbar_actions.get("sticks").get_action().triggered.connect(
            self.__slot_display_as_sticks
        )
        self._main_window.sticks_show_action.triggered.connect(self.__slot_show_as_sticks)
        self._main_window.sticks_hide_action.triggered.connect(self.__slot_hide_sticks)

        self._main_window.viewer_toolbar_actions.get("ribbon").get_action().triggered.connect(
            self.__slot_display_as_ribbon
        )
        self._main_window.ribbon_show_action.triggered.connect(self.__slot_show_as_ribbon)
        self._main_window.ribbon_hide_action.triggered.connect(self.__slot_hide_ribbon)

        self._main_window.viewer_toolbar_actions.get("lines").get_action().triggered.connect(
            self.__slot_display_as_lines
        )
        self._main_window.lines_show_action.triggered.connect(self.__slot_show_as_lines)
        self._main_window.lines_hide_action.triggered.connect(self.__slot_hide_lines)

        self._main_window.viewer_toolbar_actions.get("spheres").get_action().triggered.connect(
            self.__slot_display_as_spheres
        )
        self._main_window.spheres_show_action.triggered.connect(self.__slot_show_as_spheres)
        self._main_window.spheres_hide_action.triggered.connect(self.__slot_hide_spheres)

        self._main_window.viewer_toolbar_actions.get("dots").get_action().triggered.connect(
            self.__slot_display_as_dots
        )
        self._main_window.dots_show_action.triggered.connect(self.__slot_show_as_dots)
        self._main_window.dots_hide_action.triggered.connect(self.__slot_hide_dots)

        self._main_window.viewer_toolbar_actions.get("mesh").get_action().triggered.connect(
            self.__slot_display_as_mesh
        )
        self._main_window.mesh_show_action.triggered.connect(self.__slot_show_as_mesh)
        self._main_window.mesh_hide_action.triggered.connect(self.__slot_hide_mesh)
        self._main_window.viewer_toolbar_actions.get("surface").get_action().triggered.connect(
            self.__slot_display_as_surface
        )
        self._main_window.surface_show_action.triggered.connect(self.__slot_show_as_surface)
        self._main_window.surface_hide_action.triggered.connect(self.__slot_hide_surface)
        # # </editor-fold>
        # # <editor-fold desc="Color slots">
        self._main_window.viewer_toolbar_actions.get("color").get_action().triggered.connect(
            self.__slot_display_color_grid
        )
        # # <editor-fold desc="Color pads">
        self._main_window.color_grid.c_red.clicked.connect(lambda: self.__slot_apply_color("red"))
        self._main_window.color_grid.c_tv_red.clicked.connect(
            lambda: self.__slot_apply_color("tv_red")
        )
        self._main_window.color_grid.c_salomon.clicked.connect(
            lambda: self.__slot_apply_color("salmon")
        )
        self._main_window.color_grid.c_raspberry.clicked.connect(
            lambda: self.__slot_apply_color("raspberry")
        )

        self._main_window.color_grid.c_green.clicked.connect(
            lambda: self.__slot_apply_color("green")
        )
        self._main_window.color_grid.c_tv_green.clicked.connect(
            lambda: self.__slot_apply_color("tv_green")
        )
        self._main_window.color_grid.c_palegreen.clicked.connect(
            lambda: self.__slot_apply_color("palegreen")
        )
        self._main_window.color_grid.c_forest.clicked.connect(
            lambda: self.__slot_apply_color("forest")
        )

        self._main_window.color_grid.c_blue.clicked.connect(lambda: self.__slot_apply_color("blue"))
        self._main_window.color_grid.c_tv_blue.clicked.connect(
            lambda: self.__slot_apply_color("tv_blue")
        )
        self._main_window.color_grid.c_lightblue.clicked.connect(
            lambda: self.__slot_apply_color("lightblue")
        )
        self._main_window.color_grid.c_skyblue.clicked.connect(
            lambda: self.__slot_apply_color("skyblue")
        )

        self._main_window.color_grid.c_yellow.clicked.connect(
            lambda: self.__slot_apply_color("yellow")
        )
        self._main_window.color_grid.c_tv_yellow.clicked.connect(
            lambda: self.__slot_apply_color("tv_yellow")
        )
        self._main_window.color_grid.c_paleyellow.clicked.connect(
            lambda: self.__slot_apply_color("paleyellow")
        )
        self._main_window.color_grid.c_sand.clicked.connect(lambda: self.__slot_apply_color("sand"))

        self._main_window.color_grid.c_magenta.clicked.connect(
            lambda: self.__slot_apply_color("magenta")
        )
        self._main_window.color_grid.c_purple.clicked.connect(
            lambda: self.__slot_apply_color("purple")
        )
        self._main_window.color_grid.c_pink.clicked.connect(lambda: self.__slot_apply_color("pink"))
        self._main_window.color_grid.c_hotpink.clicked.connect(
            lambda: self.__slot_apply_color("hotpink")
        )

        self._main_window.color_grid.c_cyan.clicked.connect(lambda: self.__slot_apply_color("cyan"))
        self._main_window.color_grid.c_aquamarine.clicked.connect(
            lambda: self.__slot_apply_color("aquamarine")
        )
        self._main_window.color_grid.c_palecyan.clicked.connect(
            lambda: self.__slot_apply_color("palecyan")
        )
        self._main_window.color_grid.c_teal.clicked.connect(lambda: self.__slot_apply_color("teal"))

        self._main_window.color_grid.c_orange.clicked.connect(
            lambda: self.__slot_apply_color("orange")
        )
        self._main_window.color_grid.c_tv_orange.clicked.connect(
            lambda: self.__slot_apply_color("tv_orange")
        )
        self._main_window.color_grid.c_lightorange.clicked.connect(
            lambda: self.__slot_apply_color("lightorange")
        )
        self._main_window.color_grid.c_olive.clicked.connect(
            lambda: self.__slot_apply_color("olive")
        )

        self._main_window.color_grid.c_white.clicked.connect(
            lambda: self.__slot_apply_color("white")
        )
        self._main_window.color_grid.c_grey_70.clicked.connect(
            lambda: self.__slot_apply_color("grey70")
        )
        self._main_window.color_grid.c_grey_30.clicked.connect(
            lambda: self.__slot_apply_color("grey30")
        )
        self._main_window.color_grid.c_black.clicked.connect(
            lambda: self.__slot_apply_color("black")
        )
        # # </editor-fold>
        # # </editor-fold>
        # # <editor-fold desc="Selection slots">
        # self._main_window.show_sele_rb_panel_item.get_action().triggered.connect(
        #     self.__slot_show_sele
        # )
        # self._main_window.hide_sele_rb_panel_item.get_action().triggered.connect(
        #     self.__slot_hide_sele
        # )
        # self._main_window.clear_sele_rb_panel_item.get_action().triggered.connect(
        #     self.__slot_clear_sele
        # )
        # # </editor-fold>
        # # </editor-fold>

        # # </editor-fold>
        self._main_window.pyssa_objects_panel.tree_view.clicked.connect(
            self.__slot_select_pymol_object
        )
        # </editor-fold>
        # self.feedback_timer.setSingleShot(True)
        # self.feedback_timer.timeout.connect(self.update_protein_structure_tree_view) # TODO: Refactor the method update_protein_structure_tree_view() first!
        # </editor-fold>

    # </editor-fold>

    # <editor-fold desc="Public getters">
    def get_main_window(self) -> "main_window.MainWindow":
        """Returns the main window of the controller.

        Returns:
          The main window instance of the controller
        """
        return self._main_window

    # </editor-fold>

    # <editor-fold desc="Command runner bridge">
    def __slot_on_command_executed(self, command_name: str, line: str) -> None:
        """Dispatcher for special-processed commands coming from the Command Runner.

        Easily extend by adding entries into self._special_command_handlers mapping.
        """
        try:
            handler = getattr(self, "_special_command_handlers", {}).get(command_name)
            if handler is not None:
                handler(command_name, line)
        except Exception as e:
            logger.error(f"Error handling command '{command_name}': {e}")

    def __handle_refresh_molecule_objects(self, _command_name: str, _line: str) -> None:
        """Refresh Molecule Objects protein model from current PyMOL state.

        Uses get_model() and relies on ProteinModel.add_protein() to skip duplicates.
        """
        try:
            model = self._user_pymol.get_cmd_module().get_model()
            self._main_window.pyssa_objects_panel._model.add_protein(model)
        except Exception as e:
            logger.error(f"Failed to refresh molecule objects: {e}")
    # </editor-fold>

    # <editor-fold desc="Slot methods">

    # <editor-fold desc="Left side panel slots">
    def __slot_close_left_side_panel(self) -> None:
        """Closes the left side panel using the ToolWindowLayout."""
        try:
            self._main_window.tool_window_layout.set_left_panel_hidden(True)
        except Exception as e:
            logger.error(f"Failed to close left side panel: {e}")

    def __slot_show_protein_structure_side_panel(self) -> None:
        """Toggle the protein structure left side panel using ToolWindowLayout.

        - If the left panel is visible and already showing the protein structure page,
          hide the left panel.
        - Otherwise, show the left panel and switch to the protein structure page.
        """
        try:
            target_index = enums.LeftSidePanel.PROTEIN_STRUCTURE
            layout = self._main_window.tool_window_layout
            is_hidden = layout.is_left_panel_hidden
            current_index = layout.left_stack.currentIndex()
            if not is_hidden and current_index == target_index:
                layout.set_left_panel_hidden(True)
            else:
                layout.left_stack.setCurrentIndex(target_index)
                layout.set_left_panel_hidden(False)
        except Exception as e:
            logger.error(f"Failed to toggle left side panel: {e}")

    # </editor-fold>

    # <editor-fold desc="Project management">
    def __slot_create_project(self) -> None:
        """Initializes the CreateProjectViewController and shows the create project dialog."""
        try:
            logger.log(
                log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                "Menu entry 'Project/Create' clicked.",
            )
            self._external_controller = (
                create_project_view_controller.CreateProjectViewController(
                    self._interface_manager
                )
            )
            self._external_controller.user_input.connect(self._post_create_project)
            self._interface_manager.get_create_view().show()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message(
                "An unknown error occurred!"
            )

    def _post_create_project(self, user_input: tuple) -> None:
        """Creates a new project based on the user's input.

        Args:
            user_input: A tuple containing the user input for creating a new project. The tuple must have two elements:
                        - tmp_project_name: A string representing the name of the new project.
                        - tmp_protein_name: A string representing the name of the protein associated with the new project.

        Notes:
            This method is called after the user provides input to create a new project. It performs the following steps:
        """
        # <editor-fold desc="Checks">
        if user_input is None:
            logger.error("user_input is None.")
            self._interface_manager.status_bar_manager.show_error_message(
                "No data received!"
            )
            return

        # </editor-fold>

        try:
            self._interface_manager.block_gui(with_wait_cursor=True)
            tmp_project_name, tmp_protein_name = user_input
            tmp_project_database_filepath = str(
                pathlib.Path(
                    f"{self._interface_manager.get_application_settings().workspace_path}/{tmp_project_name}.db"
                )
            )
            with database_manager.DatabaseManager(
                    tmp_project_database_filepath
            ) as db_manager:
                db_manager.build_new_database()
            self._database_thread = database_thread.DatabaseThread(
                tmp_project_database_filepath
            )
            # self._database_thread.start()
            # self._active_task = tasks.LegacyTask(
            #     target=project_async.create_new_project,
            #     args=(
            #         tmp_project_name,
            #         self._interface_manager.get_application_settings().get_workspace_path(),
            #         self._interface_manager.watcher,
            #         self._interface_manager,
            #     ),
            #     post_func=self.__await_create_project,
            # )

            self._interface_manager.get_task_manager().append_task_result(
                task_result_factory.TaskResultFactory.run_task_result(
                    a_task_result=task_result.TaskResult.from_action(
                        an_action=action.Action(
                            a_target=project_async.create_new_project,
                            args=(
                                tmp_project_name,
                                self._interface_manager.get_application_settings().get_workspace_path(),
                                self._interface_manager.watcher,
                                self._interface_manager,
                            ),
                        ),
                        an_await_function=self.__await_create_project,
                    ),
                    a_task_scheduler=self._interface_manager.get_task_scheduler(),
                )
            )
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message(
                "An unknown error occurred!"
            )

    def __await_create_project(self, return_value: tuple[str, list[tuple[bool, tuple]]]) -> None:
        """Awaits the async method that creates the new project.

        Args:
            return_value (tuple[str, list[tuple[bool, tuple]]]): A tuple containing the return value of the async method.
        """
        # <editor-fold desc="Checks">
        if return_value is None:
            logger.error("return_value is None.")
            self._interface_manager.status_bar_manager.show_error_message(
                "No data received!"
            )
            return
        if return_value[0] == "":
            self._interface_manager.status_bar_manager.show_error_message(
                "Creating the project failed!"
            )
            self._interface_manager.refresh_main_view()
            self._interface_manager.stop_wait_cursor()
            return
        # </editor-fold>

        try:
            tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
            (
                _,
                tmp_project,
                self._interface_manager.watcher,
                self._interface_manager,
            ) = tmp_result
            self._interface_manager.set_new_project(tmp_project)
            self._interface_manager.refresh_workspace_model()
            self._interface_manager.pymol_session_manager.reinitialize_session()
            self._connect_sequence_selection_model()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message(
                "An unknown error occurred!"
            )
        finally:
            self._interface_manager.stop_wait_cursor()
            self._interface_manager.refresh_main_view()

    def __slot_open_project(self) -> None:
        """Opens the dialog to open a project."""
        try:
            logger.log(
                log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                "Menu entry 'Project/Open' clicked.",
            )
            self._external_controller = (
                open_project_view_controller.OpenProjectViewController(
                    self._interface_manager
                )
            )
            self._external_controller.return_value.connect(self._post_open_project)
            self._interface_manager.get_open_view().show()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message(
                "An unknown error occurred!"
            )

    def _post_open_project(self, return_value: tuple) -> None:
        """Post method that gets executed after the open project dialog closes.

        Args:
            return_value (tuple[str, list[tuple[bool, tuple]]]): The return value from the method.
        """
        # <editor-fold desc="Checks">
        if return_value is None:
            logger.error("return_value is None.")
            self._interface_manager.status_bar_manager.show_error_message(
                "No data received!"
            )
            return
        if return_value[1] is False:
            self._interface_manager.refresh_main_view()
            return

        # </editor-fold>

        try:
            self._interface_manager.status_bar_manager.show_temporary_message(
                enums.StatusMessages.OPENING_PROJECT.value,
                False,
            )
            self._interface_manager.block_gui(with_wait_cursor=True)
            tmp_project_name = return_value
            tmp_project_database_filepath = str(
                pathlib.Path(
                    f"{self._interface_manager.get_application_settings().workspace_path}/{tmp_project_name}.db",
                ),
            )
            self._database_thread = database_thread.DatabaseThread(
                tmp_project_database_filepath
            )
            # self._database_thread.start()
            self._database_manager.set_database_filepath(
                tmp_project_database_filepath
            )
            # self._active_task = tasks.LegacyTask(
            #     target=project_async.open_project,
            #     args=(
            #         tmp_project_name,
            #         tmp_project_database_filepath,
            #         self._interface_manager,
            #         self._interface_manager.pymol_session_manager,
            #         self.custom_progress_signal,
            #         self._interface_manager.watcher,
            #     ),
            #     post_func=self.__await_open_project,
            # )

            self._interface_manager.get_task_manager().append_task_result(
                task_result_factory.TaskResultFactory.run_task_result(
                    a_task_result=task_result.TaskResult.from_action(
                        an_action=action.Action(
                            a_target=project_async.open_project,
                            args=(
                                tmp_project_name,
                                tmp_project_database_filepath,
                                self._interface_manager,
                                self._interface_manager.pymol_session_manager,
                                self.custom_progress_signal,
                                self._interface_manager.watcher,
                            ),
                        ),
                        an_await_function=self.__await_open_project,
                    ),
                    a_task_scheduler=self._interface_manager.get_task_scheduler(),
                )
            )
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message(
                "An unknown error occurred!"
            )

    def __await_open_project(self, return_value: tuple[str, list[tuple[bool, tuple]]]) -> None:
        """Finishes the project opening process.

        Args:
            return_value (tuple[str, list[tuple[bool, tuple]]]): The return value from the async method.
        """
        self._interface_manager.status_bar_manager.hide_progress_bar()
        # <editor-fold desc="Checks">
        if return_value is None:
            logger.error("return_value is None.")
            self._interface_manager.status_bar_manager.show_error_message(
                "No data received!"
            )
            self._interface_manager.refresh_main_view()
            self._interface_manager.stop_wait_cursor()
            return
        if return_value[0] == "":
            self._interface_manager.status_bar_manager.show_error_message(
                enums.StatusMessages.OPENING_PROJECT_FAILED.value,
            )
            self._interface_manager.refresh_main_view()
            self._interface_manager.stop_wait_cursor()
            return

        # </editor-fold>

        try:
            tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
            _, tmp_project, tmp_interface_manager, tmp_watcher = tmp_result
            self._interface_manager = tmp_interface_manager
            self._interface_manager.watcher = tmp_watcher
            self._interface_manager.refresh_main_view()
            self._interface_manager.hide_progress_bar()
            self._interface_manager.status_bar_manager.show_temporary_message(
                enums.StatusMessages.OPENING_PROJECT_FINISHED.value,
            )
            self._connect_sequence_selection_model()
            # Expand all available proteins
            self.__slot_expand_all_proteins()
            # Expand all available protein pairs
            self.__slot_expand_all_protein_pairs()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message(
                "An unknown error occurred!"
            )
        finally:
            self._interface_manager.stop_wait_cursor()

    def __slot_use_project(self) -> None:
        """Opens the use project dialog."""
        try:
            logger.log(
                log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                "Menu entry 'Project/Use' clicked.",
            )
            self._external_controller = (
                use_project_view_controller.UseProjectViewController(
                    self._interface_manager
                )
            )
            self._external_controller.user_input.connect(self._post_use_project)
            self._interface_manager.get_use_project_view().show()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message(
                "An unknown error occurred!"
            )

    def _post_use_project(self, user_input: tuple) -> None:
        """Starts the use project process.

        Args:
            user_input (tuple): The user inputs from the use project dialog.
        """
        # <editor-fold desc="Checks">
        if user_input is None:
            logger.error("user_input is None.")
            self._interface_manager.status_bar_manager.show_error_message(
                "No data received!"
            )
            return

        # </editor-fold>

        try:
            tmp_project_database_filepath = str(
                pathlib.Path(
                    f"{self._interface_manager.get_application_settings().get_workspace_path()}/{user_input[0]}.db"
                )
            )
            with database_manager.DatabaseManager(
                    tmp_project_database_filepath
            ) as db_manager:
                db_manager.build_new_database()

            # self._active_task = tasks.LegacyTask(
            #   target=project_async.create_use_project,
            #   args=(
            #     user_input[0],
            #     self._interface_manager.get_application_settings().get_workspace_path(),
            #     user_input[1],
            #     self._interface_manager.watcher,
            #     self._interface_manager,
            #   ),
            #   post_func=self.__await_use_project,
            # )

            self._interface_manager.get_task_manager().append_task_result(
                task_result_factory.TaskResultFactory.run_task_result(
                    a_task_result=task_result.TaskResult.from_action(
                        an_action=action.Action(
                            a_target=project_async.create_use_project,
                            args=(
                                user_input[0],
                                self._interface_manager.get_application_settings().get_workspace_path(),
                                user_input[1],
                                self._interface_manager.watcher,
                                self._interface_manager,
                            ),
                        ),
                        an_await_function=self.__await_use_project,
                    ),
                    a_task_scheduler=self._interface_manager.get_task_scheduler(),
                )
            )
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message(
                "An unknown error occurred!"
            )
        else:
            self._interface_manager.block_gui(with_wait_cursor=True)
            # self._active_task.start()

    def __await_use_project(self, return_value: tuple[str, list[tuple[bool, tuple]]]) -> None:
        """Finishes the use project process.

        Args:
            return_value (tuple[str, list[tuple[bool, tuple]]]): The result data from the async method.
        """
        # <editor-fold desc="Checks">
        if return_value is None:
            logger.error("return_value is None.")
            self._interface_manager.status_bar_manager.show_error_message(
                "No data received!"
            )
            self._interface_manager.refresh_main_view()
            self._interface_manager.stop_wait_cursor()
            return
        if return_value[0] == "":
            self._interface_manager.status_bar_manager.show_error_message(
                "Using the project failed!"
            )
            self._interface_manager.refresh_main_view()
            self._interface_manager.stop_wait_cursor()
            return

        # </editor-fold>

        try:
            tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
            (
                _,
                tmp_project,
                self._interface_manager.watcher,
                self._interface_manager,
            ) = tmp_result
            self._interface_manager.set_new_project(tmp_project)
            self._interface_manager.add_project_to_workspace_model(
                tmp_project.get_project_name()
            )
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message(
                "An unknown error occurred!"
            )
        else:
            self._connect_sequence_selection_model()
            self._interface_manager.status_bar_manager.show_temporary_message(
                "Use process finished."
            )
        finally:
            self._interface_manager.pymol_session_manager.reinitialize_session()
            self._interface_manager.refresh_main_view()
            self._interface_manager.stop_wait_cursor()

    def __slot_delete_project(self) -> None:
        """Opens the delete project dialog."""
        try:
            logger.log(
                log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                "Menu entry 'Project/Delete' clicked.",
            )
            self._external_controller = (
                delete_project_view_controller.DeleteProjectViewController(
                    self._interface_manager
                )
            )
            self._interface_manager.get_delete_view().show()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message(
                "An unknown error occurred!"
            )

    def _post_delete_project(self) -> None:
        """Refreshes the main view after the delete project dialog closed."""
        self._interface_manager.refresh_main_view()

    def __slot_import_project(self) -> None:
        """Imports a project into the current workspace."""
        try:
            logger.log(
                log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                "Menu entry 'Project/Import' clicked.",
            )
            file_dialog = QtWidgets.QFileDialog()
            desktop_path = QtCore.QStandardPaths.standardLocations(
                QtCore.QStandardPaths.DesktopLocation
            )[0]
            file_dialog.setDirectory(desktop_path)
            file_path, _ = file_dialog.getOpenFileName(
                self._view,
                "Select a project file to import",
                "",
                "Project Database File (*.db)",
            )
            if not file_path:
                return
            file = QtCore.QFile(file_path)
            if not file.open(QtCore.QFile.ReadOnly | QtCore.QFile.Text):
                print("Error: Cannot open file for reading")
                return
            tmp_import_filepath = pathlib.Path(file_path)
            tmp_project_name_input_dialog = QtWidgets.QInputDialog()
            tmp_new_project_name = tmp_project_name_input_dialog.getText(
                self._view,
                "Project Name",
                "Enter A Project Name:",
                text=tmp_import_filepath.name.replace(".db", ""),
            )[0]
            if tmp_new_project_name == "":
                return
            # self._active_task = tasks.LegacyTask(
            #   target=project_async.import_project,
            #   args=(
            #     tmp_new_project_name,
            #     tmp_import_filepath,
            #     self._interface_manager,
            #   ),
            #   post_func=self._await__slot_import_project,
            # )

            self._interface_manager.get_task_manager().append_task_result(
                task_result_factory.TaskResultFactory.run_task_result(
                    a_task_result=task_result.TaskResult.from_action(
                        an_action=action.Action(
                            a_target=project_async.import_project,
                            args=(
                                tmp_new_project_name,
                                tmp_import_filepath,
                                self._interface_manager,
                            ),
                        ),
                        an_await_function=self._await__slot_import_project,
                    ),
                    a_task_scheduler=self._interface_manager.get_task_scheduler(),
                )
            )
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message(
                "An unknown error occurred!"
            )
        else:
            self._interface_manager.block_gui(with_wait_cursor=True)
            self._interface_manager.status_bar_manager.show_temporary_message(
                "Importing project ...", a_with_timeout_flag=False
            )
            # self._active_task.start()

    def _await__slot_import_project(self, return_value: tuple[str, list[tuple[bool, tuple]]]) -> None:
        """Finishes the import project process.

        Args:
            return_value (tuple[str, list[tuple[bool, tuple]]]): The result data from the async method.
        """
        # <editor-fold desc="Checks">
        if return_value is None:
            logger.error("return_value is None.")
            self._interface_manager.status_bar_manager.show_error_message(
                "No data received!"
            )
            self._interface_manager.stop_wait_cursor()
            self._interface_manager.refresh_main_view()
            return

        # </editor-fold>

        try:
            tmp_success_flag, tmp_result = task_result.TaskResult.get_single_action_result(return_value)
            self._database_thread = database_thread.DatabaseThread(tmp_result[1])
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message(
                "An unknown error occurred!"
            )
        finally:
            self._interface_manager.stop_wait_cursor()
            self._interface_manager.refresh_main_view()
            self._interface_manager.status_bar_manager.show_temporary_message(
                "Importing project finished."
            )

    def __slot_export_current_project(self) -> None:
        """Exports the current project to an importable format."""
        try:
            logger.log(
                log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                "Menu entry 'Project/Export' clicked.",
            )
            file_dialog = QtWidgets.QFileDialog()
            desktop_path = QtCore.QStandardPaths.standardLocations(
                QtCore.QStandardPaths.DesktopLocation
            )[0]
            file_dialog.setDirectory(desktop_path)
            file_path, _ = file_dialog.getSaveFileName(
                self._view,
                "Export current project",
                "",
                "Project Database File (*.db)",
            )
            if file_path:
                shutil.copyfile(
                    self._interface_manager.get_current_project().get_database_filepath(),
                    file_path,
                )
                # tmp_dialog = custom_message_box.CustomMessageBoxOk(
                #     "The project was successfully exported.", "Export Project",
                #     custom_message_box.CustomMessageBoxIcons.INFORMATION.value
                # )
                # tmp_dialog.exec_()
                self._interface_manager.status_bar_manager.show_temporary_message(
                    "The project was successfully exported."
                )
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message(
                "An unknown error occurred!"
            )

    def __slot_close_project(self) -> None:
        """Closes the current project."""
        try:
            logger.log(
                log_levels.SLOT_FUNC_LOG_LEVEL_VALUE,
                "Menu entry 'Project/Close' clicked.",
            )
            # self._active_task = tasks.LegacyTask(
            #     target=project_async.close_project,
            #     args=(
            #         self._database_thread,
            #         self._interface_manager.pymol_session_manager,
            #     ),
            #     post_func=self.__await_close_project,
            # )

            self._interface_manager.get_task_manager().append_task_result(
                task_result_factory.TaskResultFactory.run_task_result(
                    a_task_result=task_result.TaskResult.from_action(
                        an_action=action.Action(
                            a_target=project_async.close_project,
                            args=(
                                self._database_thread,
                                self._interface_manager.pymol_session_manager,
                            ),
                        ),
                        an_await_function=self.__await_close_project,
                    ),
                    a_task_scheduler=self._interface_manager.get_task_scheduler(),
                )
            )

            self._interface_manager.restore_default_main_view()
            self._interface_manager.close_job_notification_panel()
            self._interface_manager.close_job_overview_panel()
            self._disconnect_sequence_selection_model()
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message(
                "An unknown error occurred!"
            )
        else:
            self._interface_manager.block_gui(with_wait_cursor=True)
            self.update_status("Saving current project ...")
            # self._active_task.start()

    def __await_close_project(self, return_value: tuple[str, list[tuple[bool, tuple]]]) -> None:
        """Await the async closing process.

        Args:
            return_value (tuple[str, list[tuple[bool, tuple]]]): A tuple of the return value of the async closing process.
        """
        # <editor-fold desc="Checks">
        if return_value[0] == "":
            self._interface_manager.refresh_main_view()
            self._interface_manager.stop_wait_cursor()
            self._interface_manager.status_bar_manager.show_error_message(
                "Closing the project failed!"
            )
            return
        # </editor-fold>

        try:
            self._interface_manager.set_new_project(project.Project())
            self._view.ui.project_tab_widget.setCurrentIndex(0)
            self.update_status("Closing project finished.")
        except Exception as e:
            logger.error(f"An error occurred: {e}")
            self._interface_manager.status_bar_manager.show_error_message(
                "An unknown error occurred!"
            )
        finally:
            self._interface_manager.refresh_main_view()
            self._interface_manager.stop_wait_cursor()
    
    def __slot_exit_application(self) -> None:
        """Slot method for the Exit Application menu item"""
        self.shutdown_application_processes()
    
    # def __slot_close_application(self, an_event_signal: tuple[str, QtCore.QEvent]) -> None:
    #     """Closes all threads and process as well as the application itself.
    # 
    #     Args:
    #       an_event_signal: Signal that the dialog closes that consists of a string (default emtpy) and the QEvent
    #     """
    #     # <editor-fold desc="Checks">
    #     # psa_comm_api.macros.REQUIRE_NOT_NONE(an_event_signal)
    #     # </editor-fold>
    #     self.shutdown_application_processes()
    #     an_event_signal[1].accept()  # Closing QApplication

    # </editor-fold>

    # <editor-fold desc="Protein structure panel">
    def __slot_select_pymol_object(self) -> None:
        selected_indexes = self._main_window.pyssa_objects_panel.tree_view.selectionModel().selectedIndexes()
        selection_strings = []
        for idx in selected_indexes:
            selection_strings.append(
                self._main_window.pyssa_objects_panel.tree_view.model().construct_selection_string(
                    idx
                )
            )
        combined_selection_string = " or ".join(selection_strings)
        self._user_pymol.get_cmd_module().select("sele", combined_selection_string, enable=1)

    # </editor-fold>

    # <editor-fold desc="Session ribbon slots">
    # <editor-fold desc="Session slots">
    def __slot_open_session(self) -> None:
        """TODO: Change this implementation to correct session opening."""
        self._status_bar_manager.show_temporary_message("Session opened")
        self._user_pymol.get_cmd_module().save("test.pse")
        self._pymol_worker_connection.send(
            worker_command.WorkerCommand(
                "test.pse", "get_model", ("", ), True
            )
        )
        print(self._pymol_worker_connection.recv())

    # </editor-fold>

    # # <editor-fold desc="Scene slots">
    def _save_a_scene(self, scene_name: str):
        pass

    def __slot_save_scene(self) -> None:
      """Saves a pymol scene."""
      text, ok_pressed = QtWidgets.QInputDialog.getText(
        self._main_window,
        "Save Scene",
        "Enter scene name:",
        # QtWidgets.QLineEdit,
        # ""
      )
      if ok_pressed and text != '':
        self._user_pymol.get_cmd_module().scene(key=text, action="store")
        tmp_image_filepath = self._save_a_scene(text)
        self._main_window.side_panel_pymol_scenes.scenes_list.add_scene(text, pathlib.Path(tmp_image_filepath))

    def __slot_recall_scene(self, an_item) -> None:
      """Recalls an already created PyMOL scene."""
      tmp_widget = self._main_window.side_panel_pymol_scenes.scenes_list.list_widget.itemWidget(
        an_item)
      self._user_pymol.get_cmd_module().scene(tmp_widget.label.text(), "recall")

    def __slot_update_scene(self):
        """Update the currently selected PyMOL scene and refresh its thumbnail."""
        try:
            scene_name = self._main_window.side_panel_pymol_scenes.scenes_list.get_scene_name_from_item(
                self._main_window.side_panel_pymol_scenes.scenes_list.get_selected_item()
            )
            self._user_pymol.get_cmd_module().scene(key=scene_name, action="update")
            image_filepath = self._save_a_scene(scene_name)
            self._main_window.side_panel_pymol_scenes.scenes_list.update_scene_thumbnail(
                scene_name, pathlib.Path(image_filepath)
            )
        except Exception as error:
            logger.error(f"Failed to update scene: {error}")
            QtWidgets.QMessageBox.critical(
                self._main_window,
                "Error",
                f"Failed to update scene: {error}"
            )

    def __slot_delete_scene(self):
        """Deletes the currently selected PyMOL scene and removes it from the list."""
        try:
            scene_name = self._main_window.side_panel_pymol_scenes.scenes_list.get_scene_name_from_item(
                self._main_window.side_panel_pymol_scenes.scenes_list.get_selected_item()
            )
            self._user_pymol.get_cmd_module().scene(key=scene_name, action="clear")
            self._main_window.side_panel_pymol_scenes.scenes_list.remove_scene(scene_name)
        except Exception as error:
            logger.error(f"Failed to update scene: {error}")
            QtWidgets.QMessageBox.critical(
                self._main_window,
                "Error",
                f"Failed to update scene: {error}"
            )

    # # </editor-fold>

    # <editor-fold desc="Representation slots">
    # <editor-fold desc="Cartoon representation">
    def __slot_display_as_cartoon(self) -> None:
        """Opens a QMenu to show/hide the cartoon representation."""
        try:
            self._main_window.cartoon_show_hide_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("cartoon"))
            )
        except Exception as e:
            logger.error(e.__str__())

    def __slot_show_as_cartoon(self) -> None:
        """Shows the `pyssa_sele` selection in cartoon representation."""
        self._user_pymol.get_cmd_module().show("cartoon", "sele")

    def __slot_hide_cartoon(self) -> None:
        """Hides the cartoon representation of the `sele` selection."""
        self._user_pymol.get_cmd_module().hide("cartoon", "sele")

    # </editor-fold>

    # <editor-fold desc="Sticks representation">
    def __slot_display_as_sticks(self) -> None:
        """Opens a QMenu to show/hide the sticks representation."""
        try:
            self._main_window.sticks_show_hide_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("sticks"))
            )
        except Exception as e:
            logger.error(e.__str__())

    def __slot_show_as_sticks(self) -> None:
        """Shows the `pyssa_sele` selection in sticks representation."""
        self._user_pymol.get_cmd_module().show("sticks", "sele")

    def __slot_hide_sticks(self) -> None:
        """Hides the sticks representation of the `pyssa_sele` selection."""
        self._user_pymol.get_cmd_module().hide("sticks", "sele")

    # </editor-fold>

    # <editor-fold desc="Ribbon representation">
    def __slot_display_as_ribbon(self) -> None:
        """Opens a QMenu to show/hide the ribbon representation."""
        try:
            self._main_window.ribbon_show_hide_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("ribbon"))
            )
        except Exception as e:
            logger.error(e.__str__())

    def __slot_show_as_ribbon(self) -> None:
        """Shows the `sele` selection in ribbon representation."""
        self._user_pymol.get_cmd_module().show("ribbon", "sele")

    def __slot_hide_ribbon(self) -> None:
        """Hides the ribbon representation of the `sele` selection."""
        self._user_pymol.get_cmd_module().hide("ribbon", "sele")
    # </editor-fold>

    # <editor-fold desc="Lines representation">
    def __slot_display_as_lines(self) -> None:
        """Opens a QMenu to show/hide the lines representation."""
        try:
            self._main_window.lines_show_hide_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("lines"))
            )
        except Exception as e:
            logger.error(e.__str__())

    def __slot_show_as_lines(self) -> None:
        """Shows the `sele` selection in lines representation."""
        self._user_pymol.get_cmd_module().show("lines", "sele")

    def __slot_hide_lines(self) -> None:
        """Hides the lines representation of the `sele` selection."""
        self._user_pymol.get_cmd_module().hide("lines", "sele")
    # </editor-fold>

    # <editor-fold desc="Spheres representation">
    def __slot_display_as_spheres(self) -> None:
        """Opens a QMenu to show/hide the spheres representation."""
        try:
            self._main_window.spheres_show_hide_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("spheres"))
            )
        except Exception as e:
            logger.error(e.__str__())

    def __slot_show_as_spheres(self) -> None:
        """Shows the `pyssa_sele` selection in spheres representation."""
        self._user_pymol.get_cmd_module().show("spheres", "sele")

    def __slot_hide_spheres(self) -> None:
        """Hides the spheres representation of the `sele` selection."""
        self._user_pymol.get_cmd_module().hide("spheres", "sele")

    # </editor-fold>

    # <editor-fold desc="Dots representation">
    def __slot_display_as_dots(self) -> None:
        """Opens a QMenu to show/hide the dots representation."""
        try:
            self._main_window.dots_show_hide_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("dots"))
            )
        except Exception as e:
            logger.error(e.__str__())

    def __slot_show_as_dots(self) -> None:
        """Shows the `sele` selection in dots representation."""
        self._user_pymol.get_cmd_module().show("dots", "sele")

    def __slot_hide_dots(self) -> None:
        """Hides the dots representation of the `sele` selection."""
        self._user_pymol.get_cmd_module().hide("dots", "sele")
    # </editor-fold>

    # <editor-fold desc="Mesh representation">
    def __slot_display_as_mesh(self) -> None:
        """Opens a QMenu to show/hide the mesh representation."""
        try:
            self._main_window.mesh_show_hide_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("mesh"))
            )
        except Exception as e:
            logger.error(e.__str__())

    def __slot_show_as_mesh(self) -> None:
        """Shows the `sele` selection in mesh representation."""
        self._user_pymol.get_cmd_module().show("mesh", "sele")

    def __slot_hide_mesh(self) -> None:
        """Hides the mesh representation of the `sele` selection."""
        self._user_pymol.get_cmd_module().hide("mesh", "sele")
    # </editor-fold>

    # <editor-fold desc="Surface representation">
    def __slot_display_as_surface(self) -> None:
        """Opens a QMenu to show/hide the surface representation."""
        try:
            self._main_window.surface_show_hide_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("surface"))
            )
        except Exception as e:
            logger.error(e.__str__())

    def __slot_show_as_surface(self) -> None:
        """Shows the `pyssa_sele` selection in surface representation."""
        self._user_pymol.get_cmd_module().show("surface", "sele")

    def __slot_hide_surface(self) -> None:
        """Hides the surface representation of the `sele` selection."""
        self._user_pymol.get_cmd_module().hide("surface", "sele")

    # </editor-fold>
    # </editor-fold>

    # <editor-fold desc="Color slots">
    def __slot_display_color_grid(self) -> None:
        """Opens the color grid at the position of the 'color' quick access button."""
        try:
            # Prefer the viewer toolbar quick access button position
            self._main_window.color_grid_menu.exec(
                self._get_viewer_tool_bar_action_pos(self._main_window.viewer_toolbar_actions.get("color"))
            )
            # Fallback: open without explicit position
            self._main_window.color_grid_menu.exec()
        except Exception as e:
            logger.error(e.__str__())

    def __slot_apply_color(self, a_color_name) -> None:
        """Colors the default sele selection in the given color."""
        self._user_pymol.get_cmd_module().color(a_color_name, "sele")

    # </editor-fold>

    # <editor-fold desc="Selection slots">
    def __slot_show_sele(self) -> None:
        """Highlights the selection in PyMOL"""
        self._user_pymol.get_cmd_module().select("sele", enable=1)

    def __slot_hide_sele(self) -> None:
        """Hides the selection highlighting in PyMOL"""
        self._user_pymol.get_cmd_module().select("sele", enable=0)

    def __slot_clear_sele(self) -> None:
        self._user_pymol.get_cmd_module().select("sele", "none", enable=0)
        self.feedback_timer.start(100)

    # </editor-fold>
    # </editor-fold>

    # <editor-fold desc="PyMOL related">
    def __slot_pymol_single_left(self) -> None:
        """PyMOL single left click event."""
        try:
            self.feedback_timer.start(100)
        except Exception as e:
            print(e.__str__())

    # TODO: Refactor below
    def update_protein_structure_tree_view(self) -> None:
        """Updates the protein structure tree view in the side panel."""
        try:
            selection_strings = []
            self._user_pymol.get_cmd_module().select(
                "sele", enable=1
            )  # Highlights the selection even if clicked on the PyMOL "canvas"
            atoms = self._user_pymol.get_cmd_module().get_model("sele")
            for at in atoms.atom:
                selection_strings.append(
                    self.parse_selection_string(
                        f"/1nb1//{at.chain}/{at.resi}+{str(at.resn)}/{at.name}"  # TODO: The hard-coded 3bmp can be fixed by storing the active protein name in some sort of pymol manager
                    )
                )  # TODO: Needs more work!
            self._main_window.pyssa_objects_panel.tree_view.selectionModel().clearSelection()
            self.select_item(
                self._main_window.pyssa_objects_panel.tree_view,
                selection_strings,
            )
        except Exception as e:
            print(e.__str__())

    @staticmethod
    def parse_selection_string(selection_string):
        """Parse a PyMOL selection string into hierarchical components."""
        components = selection_string.strip("/").split("/")
        object_name = components[0]  # Object
        chain_id = components[2]  # Chain
        residue_info = components[3].replace("+", " - ")  # Residue (e.g., "22+ASP")
        atom_name = components[4]  # Atom
        return object_name, chain_id, residue_info, atom_name

    @staticmethod
    def select_item(tree_view, pymol_selection_strings):
        """Select an item in the tree view by its text."""
        model: "protein_model.ProteinModel" = tree_view.model()

        # --- Begin insert (should only run once!)
        tmp_numbering_map = (
            model.create_counts_map()
        )  # TODO: This needs to be moved in the loading routine of a protein!
        # --- end

        selection_model = tree_view.selectionModel()
        residue_atom_map = {}
        for tmp_pymol_selection_string in pymol_selection_strings:
            resi_value, atom_value, tmp_atom_index = model.search_single_selection(
                tmp_pymol_selection_string
            )
            if resi_value not in residue_atom_map.keys():
                residue_atom_map[resi_value] = []
            residue_atom_map[resi_value].append((atom_value, tmp_atom_index))
        for tmp_key in residue_atom_map.keys():
            # Checks if all atoms were selected
            if len(residue_atom_map[tmp_key]) == tmp_numbering_map["atom_counts"][tmp_key]:
                selection_model.select(
                    residue_atom_map[tmp_key][0][1].parent(),
                    selection_model.Select | selection_model.Rows,
                )
            else:
                for _, tmp_atom in residue_atom_map[tmp_key]:
                    selection_model.select(tmp_atom, selection_model.Select | selection_model.Rows)

    # </editor-fold>

    # </editor-fold>

    def shutdown_application_processes(self) -> None:
        """Closes all threads and process as well as the application itself."""
        # if not pyssa_constants.FRONTEND_ONLY:
        #   # TODO: Add correct pyssa_core logic here
        #   raise NotImplementedError()
        # self.aux_pymol_client.shutdown_service()
        self._pymol_worker_connection.send(worker_command.WorkerCommand("", "shutdown", ()))
        self._pymol_worker_process.join()
        self._main_window.close()

    def _get_viewer_tool_bar_action_pos(self, an_action) -> QtCore.QPoint:
        """Return a global point beneath the toolbar button for the given action.

        This prefers the visible viewer toolbar hosted by ToolWindowLayout. If unavailable,
        it falls back to the legacy MainWindow.viewer_toolbar. If no button can be located,
        it falls back to the current cursor position to ensure a valid QPoint is always returned.
        """
        try:
            # Prefer the active viewer toolbar inside the ToolWindowLayout
            viewer_toolbar = getattr(self._main_window.tool_window_layout, "viewer_toolbar", None)
        except Exception:
            viewer_toolbar = None
        if viewer_toolbar is None:
            # Fallback to legacy attribute if present
            viewer_toolbar = getattr(self._main_window, "viewer_toolbar", None)

        if an_action is not None and viewer_toolbar is not None:
            try:
                tmp_button = viewer_toolbar.get_tool_button_for_action(an_action)
                if tmp_button is not None and tmp_button.isVisible():
                    # Position just below the button
                    return tmp_button.mapToGlobal(tmp_button.rect().bottomLeft())
            except Exception:
                pass
        # Final fallback: cursor position
        return QtGui.QCursor.pos()
