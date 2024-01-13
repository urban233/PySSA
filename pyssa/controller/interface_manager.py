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

    _workspace_model: QtGui.QStandardItemModel
    _sequence_model: QtGui.QStandardItemModel
    _protein_model: QtGui.QStandardItemModel
    _protein_pair_model: QtGui.QStandardItemModel

    def __init__(self, a_project, the_settings: "settings.Settings"):
        self._current_project = a_project
        self._application_settings = the_settings
        self._workspace_model = QtGui.QStandardItemModel()
        self._sequence_model = QtGui.QStandardItemModel()
        self._protein_model = QtGui.QStandardItemModel()
        self._protein_pair_model = QtGui.QStandardItemModel()
    
    def set_new_project(self, the_current_project: "project.Project") -> None:
        """Sets the new current project into the interface manager."""
        self._current_project = the_current_project
        self._build_proteins_model()

    def get_current_project(self) -> "project.Project":
        """Returns the current project."""
        return self._current_project

    def _build_proteins_model(self) -> None:
        if len(self._current_project.proteins) > 0:
            tmp_root_item = self._protein_model.invisibleRootItem()
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

    def refresh_main_view(self, the_main_view: "main_view.MainView"):
        """Modifies the UI of the main view based on an app model."""
        the_main_view.ui.lbl_logo.hide()
        if self._current_project.get_project_name() != "":
            the_main_view.ui.lbl_project_name.show()
            the_main_view.ui.lbl_project_name.setText(f"Name: {self._current_project.get_project_name()}")
            the_main_view.ui.project_tab_widget.show()
            the_main_view.ui.action_new_project.setEnabled(False)
            the_main_view.ui.action_open_project.setEnabled(False)
            the_main_view.ui.action_delete_project.setEnabled(False)
            the_main_view.ui.actionImport.setEnabled(False)
            the_main_view.ui.actionExport.setEnabled(False)
            the_main_view.ui.action_close_project.setEnabled(True)
            the_main_view.ui.actionPreview.setEnabled(True)
            the_main_view.ui.actionCreate_ray_ta.setEnabled(True)
            the_main_view.ui.actionCreate_drawn.setEnabled(True)
            the_main_view.ui.actionPlaceholder.setEnabled(True)
            the_main_view.ui.project_tab_widget.setCurrentIndex(0)
        else:
            the_main_view.ui.lbl_project_name.hide()
            the_main_view.ui.project_tab_widget.hide()
            the_main_view.ui.lbl_logo.show()
            the_main_view.ui.action_new_project.setEnabled(True)
            the_main_view.ui.action_open_project.setEnabled(True)
            the_main_view.ui.action_delete_project.setEnabled(True)
            the_main_view.ui.actionImport.setEnabled(True)
            the_main_view.ui.actionExport.setEnabled(True)
            the_main_view.ui.action_close_project.setEnabled(False)
            the_main_view.ui.actionPreview.setEnabled(False)
            the_main_view.ui.actionCreate_ray_ta.setEnabled(False)
            the_main_view.ui.actionCreate_drawn.setEnabled(False)
            the_main_view.ui.actionPlaceholder.setEnabled(False)

        if len(self._current_project.proteins) > 0:
            the_main_view.ui.proteins_tree_view.setModel(self._protein_model)
            the_main_view.ui.proteins_tree_view.setHeaderHidden(True)

    def show_chain_pymol_parameters(self, a_chain_item: QtGui.QStandardItem, the_main_view):
        tmp_chain: "chain.Chain" = a_chain_item.data(enums.ModelEnum.OBJECT_ROLE)
        the_main_view.setup_proteins_table(len(tmp_chain.pymol_parameters))
        i = 0
        for tmp_key in tmp_chain.pymol_parameters.keys():
            tmp_value_item = QtWidgets.QTableWidgetItem(tmp_chain.pymol_parameters[tmp_key])
            tmp_key_item = QtWidgets.QTableWidgetItem(str(tmp_key).replace("_", " "))
            tmp_key_item.setFlags(tmp_key_item.flags() & ~Qt.ItemIsEditable)
            the_main_view.ui.proteins_table_widget.setItem(i, 0, tmp_key_item)
            the_main_view.ui.proteins_table_widget.setItem(i, 1, tmp_value_item)
            i += 1
        the_main_view.cb_chain_color.setCurrentIndex(
            the_main_view.cb_chain_color.findText(tmp_chain.pymol_parameters["chain_color"])
        )
        the_main_view.cb_chain_representation.setCurrentIndex(
            the_main_view.cb_chain_representation.findText(tmp_chain.pymol_parameters["chain_representation"])
        )
        the_main_view.ui.proteins_table_widget.resizeColumnsToContents()

    def update_status_bar(self, message: str, the_main_view) -> None:
        """Sets a custom message into the status bar."""
        the_main_view.status_bar.showMessage(message)

    def start_wait_spinner(self, the_main_view) -> None:
        """Starts the spinner."""
        the_main_view.wait_spinner.start()

    def stop_wait_spinner(self, the_main_view) -> None:
        """Stops the spinner."""
        the_main_view.wait_spinner.stop()
