import os

from pyssa.controller import interface_manager
from pyssa.gui.ui.styles import styles
from pyssa.gui.ui.views import delete_project_view
from pyssa.util import input_validator, gui_utils, constants, enums


class DeleteProjectViewController():
    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager", the_delete_project_view: "delete_project_view.DeleteProjectView") -> None:
        self._interface_manager = the_interface_manager
        self._view = the_delete_project_view
        self._fill_projects_list_view()
        self._connect_all_ui_elements_to_slot_functions()

    def _fill_projects_list_view(self) -> None:
        self._view.ui.list_delete_projects_view.setModel(self._interface_manager.get_workspace_model())

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.btn_delete_delete_project.clicked.connect(self.delete_project)
        self._view.ui.txt_delete_search.textChanged.connect(self.validate_delete_search)
        self._view.ui.list_delete_projects_view.clicked.connect(self.select_project_from_delete_list)
        self._view.ui.txt_delete_selected_projects.textChanged.connect(self.activate_delete_button)

    def delete_project(self) -> None:
        """Deletes an existing project."""
        # popup message which warns the user that the selected project gets deleted
        response: bool = gui_utils.warning_message_project_gets_deleted()
        tmp_index = self._view.ui.list_delete_projects_view.currentIndex()
        if response is True:
            os.remove(self._view.ui.list_delete_projects_view.model().data(tmp_index, enums.ModelEnum.FILEPATH_ROLE))
            self._view.ui.list_delete_projects_view.model().removeRow(tmp_index.row())  # removes item from model
        else:
            constants.PYSSA_LOGGER.info("No project has been deleted. No changes were made.")

    def select_project_from_delete_list(self) -> None:
        """Selects a project from the project list on the delete page."""
        try:
            self._view.ui.txt_delete_selected_projects.setText(self._view.ui.list_delete_projects_view.currentIndex().text())
        except AttributeError:
            self._view.ui.txt_delete_selected_projects.setText("")

    def activate_delete_button(self) -> None:
        """Activates the delete button."""
        if self._view.ui.txt_delete_selected_projects.text() == "":
            self._view.ui.btn_delete_delete_project.setEnabled(False)
            styles.color_button_not_ready(self._view.ui.btn_delete_delete_project)
        else:
            self._view.ui.btn_delete_delete_project.setEnabled(True)
            styles.color_button_ready(self._view.ui.btn_delete_delete_project)

    def validate_delete_search(self) -> None:
        """Validates the input of the project name in real-time."""
        if self._view.ui.list_delete_projects_view.currentIndex().isValid():
            self._view.ui.list_delete_projects_view.model().itemFromIndex(self._view.ui.list_delete_projects_view.currentIndex()).setSelected(False)
        # set color for lineEdit
        input_validator.InputValidator.validate_search_input(
            self._view.ui.list_delete_projects_view.model(),
            self._view.ui.txt_delete_search,
            self._view.ui.lbl_delete_status_search,
            self._view.ui.txt_delete_selected_projects,
        )


