import glob
import os
import pathlib
import sys

from PyQt5 import QtGui, QtCore
from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt

from pyssa.controller import database_manager
from pyssa.gui.ui.custom_widgets import custom_line_edit
from pyssa.gui.ui.dialogs import dialog_startup
from pyssa.gui.ui.views import main_view, predict_monomer_view, distance_analysis_view, delete_project_view, \
    create_project_view, open_project_view, import_sequence_view, rename_protein_view, use_project_view
from pyssa.gui.ui.styles import styles
from pyssa.gui.ui.views import create_project_view, open_project_view, delete_project_view, import_sequence_view
from pyssa.gui.ui.views import main_view, predict_monomer_view, distance_analysis_view, results_view, add_protein_view
from pyssa.gui.ui.views import hotspots_protein_regions_view
from pyssa.gui.ui.styles import styles
from pyssa.internal.data_structures import project, settings, chain
from pyssa.internal.data_structures.data_classes import current_session
from pyssa.util import enums, constants, exception, main_window_util


class InterfaceManager:
    """A manager for all views."""
    string_model = QtCore.QStringListModel()

    _main_view: "main_view.MainView"
    _predict_monomer_view: "predict_monomer_view.PredictMonomerView"
    _distance_analysis_view: "distance_analysis_view.DistanceAnalysisView"
    _create_project_view: "create_project_view.CreateProjectView"
    _open_project_view: "open_project_view.OpenProjectView"
    _delete_project_view: "delete_project_view.DeleteProjectView"
    _results_view: "results_view.ResultsView"
    _add_protein_view: "add_protein_view.AddProteinView"
    _import_sequence_view: "import_sequence_view.ImportSequenceView"
    _rename_protein_view: "rename_protein_view.RenameProteinView"
    _use_project_view: "use_project_view.UseProjectView"
    _hotspots_protein_regions_view: "hotspots_protein_regions_view.HotspotsProteinRegionsView"

    _current_workspace: pathlib.Path
    _current_project: "project.Project"
    _current_pymol_session: "current_session.CurrentPymolSession"
    _application_settings: "settings.Settings"
    current_tab_index: int = 0

    _workspace_model: QtGui.QStandardItemModel
    _sequence_model: QtGui.QStandardItemModel
    _protein_model: QtGui.QStandardItemModel
    _protein_pair_model: QtGui.QStandardItemModel

    def __init__(self) -> None:
        # View definitions
        self._main_view = main_view.MainView()
        self._predict_monomer_view = predict_monomer_view.PredictMonomerView()
        self._distance_analysis_view = distance_analysis_view.DistanceAnalysisView()
        self._create_project_view = create_project_view.CreateProjectView()
        self._open_project_view = open_project_view.OpenProjectView()
        self._delete_project_view = delete_project_view.DeleteProjectView()
        self._hotspots_protein_regions_view = hotspots_protein_regions_view.HotspotsProteinRegionsView()
        self._results_view = results_view.ResultsView()
        self._add_protein_view = add_protein_view.AddProteinView()
        self._import_sequence_view: "import_sequence_view.ImportSequenceView" = import_sequence_view.ImportSequenceView()
        self._rename_protein_view = rename_protein_view.RenameProteinView()
        self._use_project_view = use_project_view.UseProjectView()

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

        # General attributes definitions
        self._current_workspace = self._application_settings.workspace_path
        self._current_project = project.Project()
        self._current_pymol_session = current_session.CurrentPymolSession("", "")
        # Model definitions
        self._workspace_model = QtGui.QStandardItemModel()
        self._sequence_model = QtGui.QStandardItemModel()
        self._protein_model = QtGui.QStandardItemModel()
        self._protein_pair_model = QtGui.QStandardItemModel()
        self._build_workspace_model()

    def get_main_view(self) -> "main_view.MainView":
        return self._main_view

    def get_create_view(self) -> "create_project_view.CreateProjectView":
        return self._create_project_view

    def get_open_view(self) -> "open_project_view.OpenProjectView":
        return self._open_project_view

    def get_delete_view(self) -> "delete_project_view.DeleteProjectView":
        return self._delete_project_view

    def get_import_sequence_view(self):
        return self._import_sequence_view

    def get_predict_monomer_view(self) -> "predict_monomer_view.PredictMonomerView":
        return self._predict_monomer_view

    def get_distance_analysis_view(self) -> "distance_analysis_view.DistanceAnalysisView":
        return self._distance_analysis_view

    def get_results_view(self):
        return self._results_view

    def get_add_protein_view(self):
        return self._add_protein_view

    def get_hotspots_protein_regions_view(self) -> "hotspots_protein_regions_view.HotspotsProteinRegionsView":
        return self._hotspots_protein_regions_view

    def get_rename_protein_view(self):
        return self._rename_protein_view

    def get_use_project_view(self):
        return self._use_project_view

    def get_application_settings(self) -> "settings.Settings":
        return self._application_settings

    def set_new_project(self, the_current_project: "project.Project") -> None:
        """Sets the new current project into the interface manager."""
        self._current_project = the_current_project
        self._sequence_model.clear()
        self._build_sequences_model()
        self._protein_model.clear()
        self._build_proteins_model()
        self._protein_pair_model.clear()
        self._build_protein_pairs_model()

    def get_current_project(self) -> "project.Project":
        """Returns the current project."""
        return self._current_project

    def get_workspace_path(self):
        return self._current_workspace

    def set_new_workspace(self, the_current_workspace) -> None:
        """Sets the new current workspace into the interface manager."""
        self._current_workspace = the_current_workspace
        self._workspace_model.clear()
        self._build_workspace_model()

    def get_workspace_model(self) -> QtGui.QStandardItemModel:
        """Returns the current workspace model"""
        return self._workspace_model

    def get_workspace_projects(self):
        """Returns the workspace projects."""
        db_pattern = os.path.join(self._application_settings.get_workspace_path(), '*.db')
        self.string_model.setStringList(
            # Filters the workspace for all project files based on the xml extension
            [os.path.basename(file).replace(".db", "") for file in glob.glob(db_pattern)]
        )
        return self.string_model

    def get_workspace_projects_as_list(self) -> list:
        db_pattern = os.path.join(self._application_settings.get_workspace_path(), '*.db')
        return [os.path.basename(file).replace(".db", "") for file in glob.glob(db_pattern)]

    def get_information_about_current_session(self):
        return self._current_pymol_session.session_name, self._current_pymol_session.object_type

    def set_new_session_information(self, a_session_name: str, an_object_name: str) -> None:
        self._current_pymol_session.session_name = a_session_name
        self._current_pymol_session.object_type = an_object_name

    def update_settings(self):
        """Deserializes the settings json file."""
        self._application_settings = self._application_settings.deserialize_settings()

    def refresh_protein_model(self):
        self._protein_model.clear()
        self._build_proteins_model()

    def refresh_protein_pair_model(self):
        self._protein_pair_model.clear()
        self._build_protein_pairs_model()

    def refresh_sequence_model(self):
        self._sequence_model.clear()
        self._build_sequences_model()

    def _build_workspace_model(self) -> None:
        tmp_workspace = self._current_workspace
        xml_pattern = os.path.join(tmp_workspace, '*.xml')
        tmp_root_item = self._workspace_model.invisibleRootItem()
        for tmp_filename in [os.path.basename(file).replace(".xml", "") for file in glob.glob(xml_pattern)]:
            tmp_project_item = QtGui.QStandardItem(tmp_filename)
            tmp_filepath = pathlib.Path(f"{tmp_workspace}/{tmp_filename}.xml")
            tmp_project_item.setData(tmp_filepath, enums.ModelEnum.FILEPATH_ROLE)
            tmp_root_item.appendRow(tmp_project_item)

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

    def _build_protein_pairs_model(self) -> None:
        if len(self._current_project.protein_pairs) > 0:
            tmp_root_item = self._protein_pair_model.invisibleRootItem()
            for tmp_protein_pair in self._current_project.protein_pairs:
                tmp_protein_pair_item = QtGui.QStandardItem(tmp_protein_pair.name)
                tmp_protein_pair_item.setData(tmp_protein_pair, enums.ModelEnum.OBJECT_ROLE)
                tmp_protein_pair_item.setData("protein_pair", enums.ModelEnum.TYPE_ROLE)
                # Create protein 1 item
                tmp_protein_item_1 = QtGui.QStandardItem(tmp_protein_pair.protein_1.get_molecule_object())
                tmp_protein_item_1.setData(tmp_protein_pair.protein_1, enums.ModelEnum.OBJECT_ROLE)
                tmp_protein_item_1.setData("protein", enums.ModelEnum.TYPE_ROLE)
                for tmp_chain in tmp_protein_pair.protein_1.chains:
                    tmp_chain_item = QtGui.QStandardItem(tmp_chain.chain_letter)
                    tmp_chain_item.setData(tmp_chain, enums.ModelEnum.OBJECT_ROLE)
                    tmp_chain_item.setData("chain", enums.ModelEnum.TYPE_ROLE)
                    tmp_protein_item_1.appendRow(tmp_chain_item)
                # Create protein 2 item
                tmp_protein_item_2 = QtGui.QStandardItem(tmp_protein_pair.protein_2.get_molecule_object())
                tmp_protein_item_2.setData(tmp_protein_pair.protein_2, enums.ModelEnum.OBJECT_ROLE)
                tmp_protein_item_2.setData("protein", enums.ModelEnum.TYPE_ROLE)
                for tmp_chain in tmp_protein_pair.protein_2.chains:
                    tmp_chain_item = QtGui.QStandardItem(tmp_chain.chain_letter)
                    tmp_chain_item.setData(tmp_chain, enums.ModelEnum.OBJECT_ROLE)
                    tmp_chain_item.setData("chain", enums.ModelEnum.TYPE_ROLE)
                    tmp_protein_item_2.appendRow(tmp_chain_item)
                tmp_protein_pair_item.appendRow(tmp_protein_item_1)
                tmp_protein_pair_item.appendRow(tmp_protein_item_2)
                tmp_root_item.appendRow(tmp_protein_pair_item)

    def _build_sequences_model(self):
        if len(self._current_project.sequences) > 0:
            tmp_root_item = self._sequence_model.invisibleRootItem()
            for tmp_sequence in self._current_project.sequences:
                tmp_sequence_item = QtGui.QStandardItem(tmp_sequence.name)
                tmp_sequence_item.setData(tmp_sequence, enums.ModelEnum.OBJECT_ROLE)
                tmp_root_item.appendRow(tmp_sequence_item)

    def refresh_main_view(self):
        """Modifies the UI of the main view based on an app model."""
        self._main_view.ui.lbl_logo.hide()
        if self._current_project.get_project_name() != "":
            # A project is open
            self._main_view.ui.lbl_project_name.show()
            self._main_view.ui.lbl_project_name.setText(f"Name: {self._current_project.get_project_name()}")
            self._main_view.ui.project_tab_widget.show()
            self._main_view.ui.action_new_project.setEnabled(False)
            self._main_view.ui.action_open_project.setEnabled(False)
            self._main_view.ui.action_use_project.setEnabled(True)
            self._main_view.ui.action_delete_project.setEnabled(False)
            self._main_view.ui.action_import_project.setEnabled(False)
            self._main_view.ui.action_export_project.setEnabled(True)
            self._main_view.ui.action_close_project.setEnabled(True)
            self._main_view.ui.action_predict_monomer.setEnabled(True)
            self._main_view.ui.action_distance_analysis.setEnabled(True)
            self._main_view.ui.action_preview_image.setEnabled(True)
            self._main_view.ui.action_ray_tracing_image.setEnabled(True)
            self._main_view.ui.action_simple_image.setEnabled(True)
            self._main_view.ui.action_protein_regions.setEnabled(False)
            self._main_view.ui.project_tab_widget.setCurrentIndex(self.current_tab_index)
        else:
            # No project is open
            self._main_view.ui.lbl_project_name.hide()
            self._main_view.ui.project_tab_widget.hide()
            self._main_view.ui.lbl_logo.show()
            self._main_view.ui.action_new_project.setEnabled(True)
            self._main_view.ui.action_open_project.setEnabled(True)
            self._main_view.ui.action_use_project.setEnabled(False)
            self._main_view.ui.action_delete_project.setEnabled(True)
            self._main_view.ui.action_import_project.setEnabled(True)
            self._main_view.ui.action_export_project.setEnabled(False)
            self._main_view.ui.action_close_project.setEnabled(False)
            self._main_view.ui.action_preview_image.setEnabled(False)
            self._main_view.ui.action_ray_tracing_image.setEnabled(False)
            self._main_view.ui.action_simple_image.setEnabled(False)
            self._main_view.ui.action_protein_regions.setEnabled(False)

        if len(self._current_project.proteins) > 0:
            self._main_view.ui.proteins_tree_view.setModel(self._protein_model)
            self._main_view.ui.proteins_tree_view.setHeaderHidden(True)

        if len(self._current_project.protein_pairs) > 0:
            self._main_view.ui.protein_pairs_tree_view.setModel(self._protein_pair_model)
            self._main_view.ui.protein_pairs_tree_view.setHeaderHidden(True)

        if len(self._current_project.sequences) > 0:
            self._main_view.ui.seqs_list_view.setModel(self._sequence_model)

    def restore_default_main_view(self):
        # Restore sequences table
        self._main_view.ui.seqs_table_widget.clear()
        self._main_view.ui.seqs_table_widget.setRowCount(0)
        # Restore proteins table
        self._main_view.ui.proteins_table_widget.clear()
        self._main_view.ui.proteins_table_widget.setRowCount(0)
        # Restore protein pairs table
        self._main_view.ui.protein_pairs_table_widget.clear()
        self._main_view.ui.protein_pairs_table_widget.setRowCount(0)
        # Initialize UI
        self._main_view.initialize_ui()

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

    def show_chain_pymol_parameter_for_protein_pairs(self,
                                                     a_chain_item,
                                                     a_protein_pair_id: int,
                                                     a_protein_id: int):
        tmp_chain: "chain.Chain" = a_chain_item.data(enums.ModelEnum.OBJECT_ROLE)
        self._main_view.setup_protein_pairs_table(len(tmp_chain.pymol_parameters))
        with database_manager.DatabaseManager(str(self._current_project.get_database_filepath())) as db_manager:
            db_manager.open_project_database()
            tmp_color = db_manager.get_pymol_parameter_for_certain_protein_chain_in_protein_pair(
                a_protein_pair_id, a_protein_id, tmp_chain.chain_letter, enums.PymolParameterEnum.COLOR.value
            )
            tmp_representation = db_manager.get_pymol_parameter_for_certain_protein_chain_in_protein_pair(
                a_protein_pair_id, a_protein_id, tmp_chain.chain_letter, enums.PymolParameterEnum.REPRESENTATION.value
            )
            db_manager.close_project_database()
        tmp_pymol_parameters = {
            enums.PymolParameterEnum.COLOR.value: tmp_color[0],
            enums.PymolParameterEnum.REPRESENTATION.value: tmp_representation[0],
        }
        i = 0
        for tmp_key in tmp_pymol_parameters.keys():
            tmp_value_item = QtWidgets.QTableWidgetItem(tmp_pymol_parameters[tmp_key])
            tmp_key_item = QtWidgets.QTableWidgetItem(str(tmp_key).replace("_", " "))
            tmp_key_item.setFlags(tmp_key_item.flags() & ~Qt.ItemIsEditable)
            self._main_view.ui.protein_pairs_table_widget.setItem(i, 0, tmp_key_item)
            self._main_view.ui.protein_pairs_table_widget.setItem(i, 1, tmp_value_item)
            i += 1
        self._main_view.cb_chain_color_protein_pair.setCurrentIndex(
            self._main_view.cb_chain_color_protein_pair.findText(tmp_pymol_parameters[enums.PymolParameterEnum.COLOR.value])
        )
        self._main_view.cb_chain_representation_protein_pair.setCurrentIndex(
            self._main_view.cb_chain_representation_protein_pair.findText(tmp_pymol_parameters[enums.PymolParameterEnum.REPRESENTATION.value])
        )
        self._main_view.ui.protein_pairs_table_widget.resizeColumnsToContents()

    def update_status_bar(self, message: str) -> None:
        """Sets a custom message into the status bar."""
        self._main_view.status_bar.showMessage(message)

    def show_sequence_parameters(self, a_sequence_item: QtGui.QStandardItem):
        self._main_view.setup_sequences_table(2)
        tmp_sequence = a_sequence_item.data(enums.ModelEnum.OBJECT_ROLE)
        # Table label items
        tmp_name_label_item = QtWidgets.QTableWidgetItem("Name")
        self._main_view.ui.seqs_table_widget.setItem(0, 0, tmp_name_label_item)
        tmp_sequence_label_item = QtWidgets.QTableWidgetItem("Sequence")
        self._main_view.ui.seqs_table_widget.setItem(1, 0, tmp_sequence_label_item)
        # Table value items
        tmp_seq_name_item = QtWidgets.QTableWidgetItem(tmp_sequence.name)
        tmp_seq_name_item.setToolTip("Click to edit name")
        self._main_view.ui.seqs_table_widget.setItem(0, 1, tmp_seq_name_item)
        tmp_sequence_item = QtWidgets.QTableWidgetItem(f"{tmp_sequence.seq[:15]} ...")
        tmp_sequence_item.setToolTip("Click to view complete sequence")
        tmp_sequence_item.setData(enums.ModelEnum.OBJECT_ROLE, tmp_sequence)
        self._main_view.ui.seqs_table_widget.setItem(1, 1, tmp_sequence_item)
        # Table item flags
        tmp_sequence_item.setFlags(tmp_sequence_item.flags() & ~Qt.ItemIsEditable)
        tmp_name_label_item.setFlags(tmp_name_label_item.flags() & ~Qt.ItemIsEditable)
        tmp_sequence_label_item.setFlags(tmp_sequence_label_item.flags() & ~Qt.ItemIsEditable)
        # Extra cell widgets

        self._main_view.line_edit_seq_name.setText(tmp_seq_name_item.data(Qt.DisplayRole))
        self._main_view.ui.seqs_table_widget.resizeColumnsToContents()

    def manage_buttons_for_proteins_tab(self,
                                        an_object_type: str,
                                        is_protein_in_pair: bool):
        if an_object_type == "protein":
            self._main_view.ui.btn_save_protein.setEnabled(True)
            self._main_view.ui.btn_open_protein_session.setEnabled(True)
            self._main_view.ui.btn_create_protein_scene.setEnabled(True)
            self._main_view.ui.btn_update_protein_scene.setEnabled(True)
        elif an_object_type == "chain":
            self._main_view.ui.btn_save_protein.setEnabled(False)
            self._main_view.ui.btn_open_protein_session.setEnabled(False)
            self._main_view.ui.btn_create_protein_scene.setEnabled(True)
            self._main_view.ui.btn_update_protein_scene.setEnabled(True)
        else:
            constants.PYSSA_LOGGER.warning("Unknown object type on proteins tab selected.")

        if is_protein_in_pair:
            self._main_view.ui.btn_delete_protein.setEnabled(False)
        else:
            self._main_view.ui.btn_delete_protein.setEnabled(True)

    def manage_buttons_for_protein_pairs_tab(self, an_object_type: str):
        if an_object_type == "protein_pair":
            self._main_view.ui.btn_delete_protein_pair.setEnabled(True)
            self._main_view.ui.btn_open_protein_pair_session.setEnabled(True)
            self._main_view.ui.btn_create_protein_pair_scene.setEnabled(True)
            self._main_view.ui.btn_update_protein_pair_scene.setEnabled(True)
        else:
            self._main_view.ui.btn_delete_protein_pair.setEnabled(False)
            self._main_view.ui.btn_open_protein_pair_session.setEnabled(False)
            self._main_view.ui.btn_create_protein_pair_scene.setEnabled(True)
            self._main_view.ui.btn_update_protein_pair_scene.setEnabled(True)

    # <editor-fold desc="Getter methods for protein tab in main view">
    def get_current_protein_tree_index(self):
        return self._main_view.ui.proteins_tree_view.currentIndex()

    def get_child_index_of_get_current_protein_tree_index(self):
        return self._main_view.ui.proteins_tree_view.currentIndex().child(0, 0)

    def get_current_protein_tree_index_type(self):
        return self._main_view.ui.proteins_tree_view.model().data(
            self.get_current_protein_tree_index(), enums.ModelEnum.TYPE_ROLE
        )

    def get_current_protein_tree_index_object(self):
        return self.get_current_protein_tree_index().data(enums.ModelEnum.OBJECT_ROLE)

    def get_parent_index_object_of_current_protein_tree_index(self):
        return self.get_current_protein_tree_index().parent().data(enums.ModelEnum.OBJECT_ROLE)

    # </editor-fold>

    # <editor-fold desc="Getter methods for protein pairs tab in main view">
    def get_current_protein_pair_tree_index(self):
        return self._main_view.ui.protein_pairs_tree_view.currentIndex()

    def get_child_index_of_get_current_protein_pair_tree_index(self):
        return self._main_view.ui.protein_pairs_tree_view.currentIndex().child(0, 0)

    def get_current_protein_pair_tree_index_type(self):
        return self._main_view.ui.protein_pairs_tree_view.model().data(
            self.get_current_protein_pair_tree_index(), enums.ModelEnum.TYPE_ROLE
        )

    def get_current_protein_pair_tree_index_object(self):
        return self.get_current_protein_pair_tree_index().data(enums.ModelEnum.OBJECT_ROLE)

    def get_parent_index_object_of_current_protein_pair_tree_index(self):
        return self.get_current_protein_pair_tree_index().parent().data(enums.ModelEnum.OBJECT_ROLE)

    def get_grand_parent_index_object_of_current_protein_pair_tree_index(self):
        return self.get_current_protein_pair_tree_index().parent().parent().data(enums.ModelEnum.OBJECT_ROLE)
    # </editor-fold>

    def start_wait_spinner(self) -> None:
        """Starts the spinner."""
        self._main_view.wait_spinner.start()

    def stop_wait_spinner(self) -> None:
        """Stops the spinner."""
        self._main_view.wait_spinner.stop()
