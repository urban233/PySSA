from PyQt5 import QtGui
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt

from pyssa.gui.ui.styles import styles
from pyssa.gui.ui.views import main_view
from pyssa.internal.data_structures import project, settings, chain
from pyssa.util import enums, gui_utils


class InterfaceManager:
    """A manager for all views."""
    _current_project: "project.Project"
    _application_settings: "settings.Settings"
    _main_view: "main_view.MainView"

    def __init__(self, a_project, the_settings: "settings.Settings", the_main_view: "main_view.MainView"):
        self._current_project = a_project
        self._application_settings = the_settings
        self._main_view = the_main_view
    
    def set_new_project(self, the_current_project: "project.Project") -> None:
        """Sets the new current project into the interface manager."""
        self._current_project = the_current_project

    def get_current_project(self) -> "project.Project":
        """Returns the current project."""
        return self._current_project

    def refresh_main_view(self, a_protein_model: QtGui.QStandardItemModel):
        """Modifies the UI of the main view based on an app model."""
        self._main_view.ui.lbl_logo.hide()
        if self._current_project.get_project_name() != "":
            self._main_view.ui.lbl_project_name.show()
            self._main_view.ui.lbl_project_name.setText(f"Name: {self._current_project.get_project_name()}")
            self._main_view.ui.project_tab_widget.show()
            self._main_view.ui.action_new_project.setEnabled(False)
            self._main_view.ui.action_open_project.setEnabled(False)
            self._main_view.ui.action_delete_project.setEnabled(False)
            self._main_view.ui.actionImport.setEnabled(False)
            self._main_view.ui.actionExport.setEnabled(False)
            self._main_view.ui.action_close_project.setEnabled(True)
            self._main_view.ui.actionPreview.setEnabled(True)
            self._main_view.ui.actionCreate_ray_ta.setEnabled(True)
            self._main_view.ui.actionCreate_drawn.setEnabled(True)
            self._main_view.ui.actionPlaceholder.setEnabled(True)
            self._main_view.ui.project_tab_widget.setCurrentIndex(0)
        else:
            self._main_view.ui.lbl_project_name.hide()
            self._main_view.ui.project_tab_widget.hide()
            self._main_view.ui.lbl_logo.show()
            self._main_view.ui.action_new_project.setEnabled(True)
            self._main_view.ui.action_open_project.setEnabled(True)
            self._main_view.ui.action_delete_project.setEnabled(True)
            self._main_view.ui.actionImport.setEnabled(True)
            self._main_view.ui.actionExport.setEnabled(True)
            self._main_view.ui.action_close_project.setEnabled(False)
            self._main_view.ui.actionPreview.setEnabled(False)
            self._main_view.ui.actionCreate_ray_ta.setEnabled(False)
            self._main_view.ui.actionCreate_drawn.setEnabled(False)
            self._main_view.ui.actionPlaceholder.setEnabled(False)

        if len(self._current_project.proteins) > 0:
            tmp_root_item = a_protein_model.invisibleRootItem()
            for tmp_protein in self._current_project.proteins:
                tmp_protein_item = QtGui.QStandardItem(tmp_protein.get_molecule_object())
                tmp_protein_item.setData(tmp_protein, enums.ModelEnum.OBJECT_ROLE)
                tmp_protein_item.setData("protein", enums.ModelEnum.TYPE_ROLE)
                tmp_root_item.appendRow(tmp_protein_item)
                for tmp_chain in tmp_protein.chains:
                    tmp_chain_item = QtGui.QStandardItem(tmp_chain.chain_letter)
                    tmp_chain_item.setData(tmp_chain, enums.ModelEnum.OBJECT_ROLE)
                    tmp_chain_item.setData("chain", enums.ModelEnum.TYPE_ROLE)
                    tmp_protein_item.appendRow(tmp_chain_item)
            self._main_view.ui.proteins_tree_view.setModel(a_protein_model)
            self._main_view.ui.proteins_tree_view.setHeaderHidden(True)

    def show_chain_pymol_parameters(self, a_chain_item: QtGui.QStandardItem):
        tmp_chain: "chain.Chain" = a_chain_item.data(enums.ModelEnum.OBJECT_ROLE)
        self._main_view.setup_proteins_table(len(tmp_chain.pymol_parameters))
        i = 0
        for tmp_key in tmp_chain.pymol_parameters.keys():
            tmp_value_item = QtWidgets.QTableWidgetItem(tmp_chain.pymol_parameters[tmp_key])
            tmp_key_item = QtWidgets.QTableWidgetItem(str(tmp_key).replace("_", " "))
            tmp_key_item.setFlags(tmp_key_item.flags() & ~Qt.ItemIsEditable)
            self._main_view.ui.proteins_table_widget.setItem(i, 0, tmp_key_item)
            self._main_view.ui.proteins_table_widget.setItem(i, 1, tmp_value_item)
            i += 1
        self._main_view.cb_chain_color.setCurrentIndex(
            self._main_view.cb_chain_color.findText(tmp_chain.pymol_parameters["chain_color"])
        )
        self._main_view.cb_chain_representation.setCurrentIndex(
            self._main_view.cb_chain_representation.findText(tmp_chain.pymol_parameters["chain_representation"])
        )
        self._main_view.ui.proteins_table_widget.resizeColumnsToContents()

    def update_status_bar(self, message: str) -> None:
        """Sets a custom message into the status bar."""
        self._main_view.status_bar.showMessage(message)

    def start_wait_spinner(self) -> None:
        """Starts the spinner."""
        self._main_view.wait_spinner.start()

    def stop_wait_spinner(self) -> None:
        """Stops the spinner."""
        self._main_view.wait_spinner.stop()
