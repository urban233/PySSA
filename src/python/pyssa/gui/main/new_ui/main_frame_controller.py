import logging
import pathlib
import time
from typing import Optional

from Bio.PDB import PDBParser
from Bio import SeqIO

from PyQt6 import QtCore, QtWidgets
from tea.concurrent import task_manager, action
from tea.concurrent import task_scheduler
from tea.concurrent import task_result

from pyssa.gui.control import create_project_controller
from pyssa.gui.dialog import dialog_create_project
from pyssa.model.util.gui_style import icons
from pyssa.model.database import database_handler, databases_manager
from pyssa.gui.main import main_frame
from pyssa.gui.util import status_bar_manager
from pyssa.model.qmodel import project_model
from pyssa.model.qmodel import sequence_model
from pyssa.model.qmodel import protein_model
from pyssa.model.qmodel import ligand_model
from pyssa.model.central_objects import protein
from pyssa.model.central_objects import project
from pyssa.model.util import component_manager
from pyssa.model.preference import model_definitions
from pyssa.gui.util import import_sequence_component
from pyssa.gui.util import import_protein_component
from pyssa.gui.util import import_ligand_component
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging
from pyssa.extension.mv_wrappers.abstract_wrapper import imv_client

logger = default_logging.setup_logger(__file__)

__docformat__ = "google"


class MainFrameController:
  """Controller class for MainFrame"""

  def __init__(self, a_main_frame: "main_frame.MainFrame") -> None:
    """Constructor.

    Args:
      a_main_frame: The main frame.

    Raises:
      exception.NoneValueError: If `a_main_frame` is None.
    """
    # <editor-fold desc="Checks">
    if a_main_frame is None:
      default_logging.append_to_log_file(logger, "a_main_frame is None.", logging.ERROR)
      raise exception.NoneValueError("a_main_frame is None.")
    # </editor-fold>
    self.main_frame: "main_frame.MainFrame" = a_main_frame
    self.basic_controllers: dict = {
      "New": create_project_controller.CreateProjectController(dialog_create_project.DialogCreateProject())
    }
    self.status_bar_manager = status_bar_manager.StatusBarManager(self.main_frame)
    self.mv_client: Optional["imv_client.IMvClient"] = None
    # <editor-fold desc="Tasks">
    self.task_manager = task_manager.TaskManager()
    self.task_scheduler = task_scheduler.TaskScheduler()
    # </editor-fold>

    self.databases_manager = databases_manager.DatabasesManager()

    # <editor-fold desc="Construct component manager instance">
    self.component_manager = component_manager.ComponentManager(
      {
        model_definitions.ComponentsEnum.IMPORT_SEQUENCE: import_sequence_component.ImportSequenceComponent(
          self.schedule_tool_task_result_object),
        model_definitions.ComponentsEnum.IMPORT_PROTEIN: import_protein_component.ImportProteinComponent(
          self.schedule_tool_task_result_object),
        model_definitions.ComponentsEnum.IMPORT_LIGAND: import_ligand_component.ImportLigandComponent(
          self.schedule_tool_task_result_object),
      }
    )
    # </editor-fold>

    # <editor-fold desc="Models">
    self.project_model: "project_model.ProjectModel" = project_model.ProjectModel.from_workspace_path(
      model_definitions.ModelDefinitions.DEFAULT_WORKSPACE_PATH
    )
    self.sequence_model: "sequence_model.SequenceModel" = sequence_model.SequenceModel.from_a_list_of_sequences([])
    self.protein_model: "protein_model.ProteinModel" = protein_model.ProteinModel.from_a_list_of_proteins([])
    self.ligand_model: "ligand_model.LigandModel" = ligand_model.LigandModel.from_a_list_of_ligands([])
    # </editor-fold>

    self.project: Optional["project.Project"] = None

    # self.protein_protein_complex_model
    # self.protein_ligand_complex_model
    self.set_all_models_to_their_views()
    self.main_frame.fill_open_project_scroll_area(
      self.project_model, self.__slot_open_project
    )
    self.main_frame.fill_home_project_scroll_area(
      self.project_model, self.__slot_open_project
    )
    self.connect_all_slots()
    # TODO: below is just for testing purposes until project management is working
    #self.temp_project_setup()
    self.update_main_frame_gui()
    self.main_frame.ribbon_bar.setCurrentIndex(0)
    self.main_frame.ui.btn_project_page_back.hide()

  def connect_all_slots(self) -> None:
    """Connects all signals to their slot methods."""
    self.main_frame.dialogClosed.connect(self.close_application)
    self.main_frame.ribbon_bar.currentChanged.connect(self.switch_application_page)
    self.main_frame.ui.project_tab_widget.currentChanged.connect(self.__slot_switch_object_tab)
    # <editor-fold desc="Menu items">
    self.main_frame.ui.action_new_project.triggered.connect(self.slot_create_project_dialog)
    self.main_frame.ui.action_open_project.triggered.connect(self.slot_open_project_dialog)
    self.main_frame.ui.action_close_project.triggered.connect(self.close_project)
    self.main_frame.ui.action_exit_application.triggered.connect(self.exit_application)
    # </editor-fold>
    # <editor-fold desc="Sequences tab">
    #self.main_frame.ui.btn_import_seq.clicked.connect(self.slot_import_sequence_dialog)
    self.main_frame.ui.btn_import_seq.clicked.connect(self.main_frame.show_import_popup)
    # </editor-fold>
    # <editor-fold desc="Proteins tab">
    self.main_frame.ui.btn_import_protein.clicked.connect(self.slot_import_protein_dialog)
    self.main_frame.ui.lineEdit.textChanged.connect(self.select_based_on_protein_tree_selection_string)
    self.main_frame.ui.proteins_tree_view.clicked.connect(self.display_selection_string_for_the_current_protein)
    # </editor-fold>
    # <editor-fold desc="Ligands tab">
    self.main_frame.ui.btn_import_ligand.clicked.connect(self.slot_import_ligand_dialog)
    # </editor-fold>

    # <editor-fold desc="Project page">
    self.main_frame.ui.btn_project_page_back.clicked.connect(self.display_main_page)
    self.main_frame.ui.tab_widget_project_page.currentChanged.connect(self.display_project_tab)

    self.main_frame.ui.btn_new_blank_project.clicked.connect(self.open_create_project_dialog)
    self.main_frame.btn_blank_project.clicked.connect(self.open_create_project_dialog)
    self.basic_controllers["New"].component_task.connect(self.schedule_tool_task_result_object)
    # <editor-fold desc="Open project tab">
    # </editor-fold>
    # </editor-fold>

  def temp_project_setup(self) -> None:
    """Sets up a temporary project."""
    self.project: "project.Project" = project.Project("dummy")
    self.build_all_project_related_models()
    #self.main_frame.update_project_name(self.project.name)
    # Sequences
    records = SeqIO.read(r"C:\Users\student\Downloads\rcsb_pdb_6OMN.fasta", "fasta")
    print(records)
    # Proteins
    tmp_pdb_parser = PDBParser()
    tmp_filepath: pathlib.Path = pathlib.Path(r"C:\Users\student\github_repos\PyDD_Prototype\assets\test_files\6omn.pdb")
    tmp_structure = tmp_pdb_parser.get_structure(tmp_filepath.name.replace(".pdb", ""), str(tmp_filepath))
    tmp_protein = protein.Protein(tmp_filepath.name.replace(".pdb", ""), a_structure=tmp_structure)
    self.databases_manager.get_handler("dummy").insert_new_protein(tmp_protein, direct_write=True)
    tmp_protein_from_db = self.databases_manager.get_handler("dummy").awaits_get_protein_by_name(tmp_protein.name)
    self.protein_model.add_protein(tmp_protein_from_db)

  # <editor-fold desc="Util methods">
  # <editor-fold desc="GUI related">
  def update_main_frame_gui(self) -> None:
    """Updates the entire gui of the main frame."""
    self.main_frame.ui.tabWidget.show()
    self.main_frame.ui.menubar.hide()
    # <editor-fold desc="Always displayed and enabled">
    self.main_frame.ui.menuSettings.setEnabled(True)
    self.main_frame.ui.menuAbout.setEnabled(True)
    self.main_frame.ui.action_exit_application.setEnabled(True)
    self.main_frame.ui.action_edit_settings.setEnabled(True)
    self.main_frame.ui.action_restore_settings.setEnabled(True)
    self.main_frame.ui.action_documentation.setEnabled(True)
    self.main_frame.ui.action_arrange_windows.setEnabled(True)
    self.main_frame.ui.action_tutorials.setVisible(False)  # Hidden until tutorials are available
    self.main_frame.ui.action_show_log_in_explorer.setEnabled(True)
    self.main_frame.ui.action_clear_logs.setEnabled(True)
    self.main_frame.ui.action_restart_pymol.setEnabled(True)
    self.main_frame.ui.action_about.setEnabled(True)
    self.main_frame.ui.lbl_placeholder.hide()
    # </editor-fold>
    if self.project is None:
      # No project is open
      self.main_frame.ribbon_bar.setCurrentIndex(0)
      self.main_frame.ui.btn_project_page_back.hide()
      # <editor-fold desc="Menu items">
      self.main_frame.ui.action_new_project.setEnabled(True)
      self.main_frame.ui.action_import_project.setEnabled(True)
      self.main_frame.ui.action_export_project.setEnabled(True)
      self.main_frame.ui.action_close_project.setEnabled(False)
      if self.project_model.is_empty():
        # Empty workspace
        self.main_frame.ui.action_open_project.setEnabled(False)
        self.main_frame.ui.action_delete_project.setEnabled(False)
        self.main_frame.ui.action_use_project.setEnabled(False)
      else:
        # There are projects in the workspace
        self.main_frame.ui.action_open_project.setEnabled(True)
        self.main_frame.ui.action_delete_project.setEnabled(True)
        self.main_frame.ui.action_use_project.setEnabled(False)
      self.main_frame.ui.menuPrediction.setEnabled(False)
      self.main_frame.ui.menuAnalysis.setEnabled(False)
      self.main_frame.ui.menuDocking.setEnabled(False)
      self.main_frame.ui.menuResults.setEnabled(False)
      self.main_frame.ui.menuImage.setEnabled(False)
      self.main_frame.ui.menuHotspots.setEnabled(False)
      self.main_frame.ui.action_get_demo_projects.setEnabled(True)
      # </editor-fold>
      # <editor-fold desc="Macro UI elements">
      self.main_frame.ui.lbl_project_name.hide()
      self.main_frame.ui.lbl_logo.show()
      self.main_frame.ui.project_tab_widget.hide()
      self.main_frame.ui.tabWidget.hide()
      self.main_frame.ui.tabWidget_2.hide()
      # </editor-fold>
    else:
      # A project is open
      self.main_frame.setWindowTitle(f"PyDD - {self.project.name}")
      self.main_frame.ui.btn_project_page_back.show()
      # <editor-fold desc="Check project name display">
      if self.main_frame.ui.lbl_project_name.text() != "":
        self.main_frame.ui.lbl_project_name.setText(f"Project Name: {self.project.name}")
      # </editor-fold>
      # <editor-fold desc="Menu items">
      self.main_frame.ui.action_new_project.setEnabled(False)
      self.main_frame.ui.action_open_project.setEnabled(False)
      self.main_frame.ui.action_delete_project.setEnabled(False)
      self.main_frame.ui.action_use_project.setEnabled(True)
      self.main_frame.ui.action_import_project.setEnabled(False)
      self.main_frame.ui.action_export_project.setEnabled(False)
      self.main_frame.ui.action_close_project.setEnabled(True)
      self.main_frame.ui.action_get_demo_projects.setEnabled(False)
      # </editor-fold>
      # <editor-fold desc="Macro UI elements">
      self.main_frame.ui.lbl_project_name.show()
      self.main_frame.ui.lbl_logo.hide()
      self.main_frame.ui.project_tab_widget.show()
      self.main_frame.ui.tabWidget.show()
      self.main_frame.ui.tabWidget_2.show()
      # </editor-fold>
      # <editor-fold desc="Sequences">
      # <editor-fold desc="Contains sequences">
      if len(self.project.sequences) > 0:
        # TODO: Check which type of sequences are there!
        self.main_frame.ui.menuPrediction.setEnabled(True)
        self.main_frame.ui.action_predict_monomer.setEnabled(True)
        self.main_frame.ui.action_predict_multimer.setEnabled(True)
        self.main_frame.ui.btn_import_seq.setEnabled(True)
        self.main_frame.ui.btn_add_sequence.setEnabled(True)
        if self.main_frame.ui.project_tab_widget.currentIndex() == 0:
          # Tab widget is on the sequences tab
          self.main_frame.ui.tabWidget.setTabEnabled(1, False)
          self.main_frame.ui.tabWidget.setCurrentIndex(0)
          if len(self.main_frame.ui.seqs_list_view.selectedIndexes()) == 0:
            # NO sequence is selected
            self.main_frame.ui.btn_save_sequence.setEnabled(False)
            self.main_frame.ui.btn_delete_sequence.setEnabled(False)
          else:
            # A sequence is selected
            self.main_frame.ui.btn_save_sequence.setEnabled(True)
            self.main_frame.ui.btn_delete_sequence.setEnabled(True)
      # </editor-fold>
      # <editor-fold desc="NO sequences">
      else:
        self.main_frame.ui.menuPrediction.setEnabled(False)
        if self.main_frame.ui.project_tab_widget.currentIndex() == 0:
          # Tab widget is on the sequences tab
          self.main_frame.ui.tabWidget.setTabEnabled(1, False)
          self.main_frame.ui.tabWidget.hide()
          self.main_frame.ui.btn_import_seq.setEnabled(True)
          self.main_frame.ui.btn_add_sequence.setEnabled(True)
          self.main_frame.ui.btn_save_sequence.setEnabled(False)
          self.main_frame.ui.btn_delete_sequence.setEnabled(False)
      # </editor-fold>
      # </editor-fold>
      # <editor-fold desc="Proteins">
      # <editor-fold desc="Contains proteins">
      if len(self.project.proteins) > 0:
        self.main_frame.ui.menuAnalysis.setEnabled(True)
        self.main_frame.ui.action_distance_analysis.setEnabled(True)
        if self.main_frame.ui.project_tab_widget.currentIndex() == 1:
          self.main_frame.ui.tabWidget.setTabEnabled(1, True)
          # Tab widget is on proteins tab
          self.main_frame.ui.btn_protein_tree_view_expand.setEnabled(True)
          self.main_frame.ui.btn_protein_tree_view_collapse.setEnabled(True)
          self.main_frame.ui.btn_import_protein.setEnabled(True)
          if len(self.main_frame.ui.proteins_tree_view.selectedIndexes()) > 0:
            # A protein is selected
            self.main_frame.ui.btn_save_protein.setEnabled(True)
            self.main_frame.ui.btn_delete_protein.setEnabled(True)  # TODO: Here is a check needed if the protein is anywhere cross referenced!
          else:
            # NO protein is selected
            self.main_frame.ui.btn_save_protein.setEnabled(False)
            self.main_frame.ui.btn_delete_protein.setEnabled(False)
      # </editor-fold>
      # <editor-fold desc="NO proteins">
      else:
        self.main_frame.ui.menuAnalysis.setEnabled(False)
        if self.main_frame.ui.project_tab_widget.currentIndex() == 1:
          # Tab widget is on proteins tab
          self.main_frame.ui.tabWidget.setTabEnabled(1, False)
          self.main_frame.ui.btn_protein_tree_view_expand.setEnabled(False)
          self.main_frame.ui.btn_protein_tree_view_collapse.setEnabled(False)
          self.main_frame.ui.btn_import_protein.setEnabled(True)
          self.main_frame.ui.btn_save_protein.setEnabled(False)
          self.main_frame.ui.btn_delete_protein.setEnabled(False)
      # </editor-fold>
      # </editor-fold>
      # <editor-fold desc="Ligands">
      # <editor-fold desc="Contains ligands">
      if len(self.project.ligands) > 0:
        if self.main_frame.ui.project_tab_widget.currentIndex() == 2:
          self.main_frame.ui.tabWidget.setTabEnabled(1, False)
          # Tab widget is on ligands tab
          self.main_frame.ui.btn_ligand_tree_view_expand.setEnabled(True)
          self.main_frame.ui.btn_ligand_tree_view_collapse.setEnabled(True)
          self.main_frame.ui.btn_import_ligand.setEnabled(True)
          if len(self.main_frame.ui.ligand_tree_view.selectedIndexes()) > 0:
            # A protein is selected
            self.main_frame.ui.btn_save_ligand.setEnabled(True)
            self.main_frame.ui.btn_delete_ligand.setEnabled(
              True)  # TODO: Here is a check needed if the protein is anywhere cross referenced!
          else:
            # NO protein is selected
            self.main_frame.ui.btn_save_ligand.setEnabled(False)
            self.main_frame.ui.btn_delete_ligand.setEnabled(False)
      # </editor-fold>
      # <editor-fold desc="NO ligands">
      else:
        if self.main_frame.ui.project_tab_widget.currentIndex() == 2:
          # Tab widget is on ligands tab
          self.main_frame.ui.tabWidget.setTabEnabled(1, False)
          self.main_frame.ui.btn_ligand_tree_view_expand.setEnabled(False)
          self.main_frame.ui.btn_ligand_tree_view_collapse.setEnabled(False)
          self.main_frame.ui.btn_import_ligand.setEnabled(True)
          self.main_frame.ui.btn_save_ligand.setEnabled(False)
          self.main_frame.ui.btn_delete_ligand.setEnabled(False)
      # </editor-fold>
      # </editor-fold>
      # <editor-fold desc="Complexes">
      tmp_number_of_prot_prot_complexes = len(self.project.protein_protein_complexes)
      tmp_number_of_prot_lig_complexes = len(self.project.protein_ligand_complexes)
      if tmp_number_of_prot_prot_complexes == 0 and tmp_number_of_prot_lig_complexes == 0:
        # NO complexes
        self.main_frame.ui.project_tab_widget.setTabEnabled(3, False)
      else:
        # There are complexes
        # <editor-fold desc="Protein-protein complexes">
        if tmp_number_of_prot_prot_complexes > 0:
          # Contains protein-protein complexes
          self.main_frame.ui.tabWidget.setTabEnabled(1, True)
          if len(self.main_frame.ui.prot_prot_complex_tree_view.selectedIndexes()) > 0:
            # A protein-protein complex is selected
            self.main_frame.ui.btn_delete_prot_prot_complex.setEnabled(True)
          else:
            # NO protein-protein complex is selected
            self.main_frame.ui.btn_delete_prot_prot_complex.setEnabled(False)
        else:
          # NO protein-protein complexes
          self.main_frame.ui.tabWidget.setTabEnabled(1, False)
          self.main_frame.ui.tabWidget_2.setTabEnabled(0, False)
        # </editor-fold>
        # <editor-fold desc="Protein-ligand complexes">
        if tmp_number_of_prot_lig_complexes > 0:
          # Contains protein-ligand complexes
          if len(self.main_frame.ui.prot_prot_complex_tree_view.selectedIndexes()) > 0:
            # A protein-ligand complex is selected
            self.main_frame.ui.btn_delete_prot_lig_complex.setEnabled(True)
          else:
            # NO protein-ligand complex is selected
            self.main_frame.ui.btn_delete_prot_lig_complex.setEnabled(False)
        else:
          # NO protein-ligand complexes
          self.main_frame.ui.tabWidget_2.setTabEnabled(1, False)
        # </editor-fold>
      # </editor-fold>

  def stop_wait_cursor(self) -> None:
    """Stops the cursor."""
    QtWidgets.QApplication.restoreOverrideCursor()
    self.enable_job_panels()

  def block_gui(self, with_wait_cursor: bool = False) -> None:
    """Starts the wait cursor.

    Args:
        with_wait_cursor (bool): A boolean indicating whether to display a wait cursor.

    Raises:
        exception.IllegalArgumentError: If `with_wait_cursor` is None.
    """
    # <editor-fold desc="Checks">
    if with_wait_cursor is None:
      logger.error("with_wait_cursor is None.")
      raise exception.IllegalArgumentError("with_wait_cursor is None.")
    # </editor-fold>
    if with_wait_cursor is True:
      QtWidgets.QApplication.setOverrideCursor(QtCore.Qt.CursorShape.WaitCursor)
    self.disable_menu_bar_without_exit_application()
    self.main_frame.ui.project_tab_widget.setEnabled(False)
    self.main_frame.ui.tabWidget.setEnabled(False)
    self.main_frame.ui.frame_job_overview.setEnabled(False)
    self.main_frame.ui.frame_job_notification.setEnabled(False)
    self.main_frame.btn_open_job_overview.setEnabled(False)
    self.main_frame.btn_open_job_notification.setEnabled(False)

  def unblock_gui(self) -> None:
    """Unblocks the tabs of the main frame."""
    self.main_frame.ui.project_tab_widget.setEnabled(True)
    self.main_frame.ui.tabWidget.setEnabled(True)
    self.update_main_frame_gui()

  def disable_menu_bar_without_exit_application(self) -> None:
    """Disables the menu entries but not 'Exit Application'."""
    self.main_frame.ui.menuProject.setEnabled(True)
    self.main_frame.ui.action_new_project.setEnabled(False)
    self.main_frame.ui.action_open_project.setEnabled(False)
    self.main_frame.ui.action_close_project.setEnabled(False)
    self.main_frame.ui.action_use_project.setEnabled(False)
    self.main_frame.ui.action_delete_project.setEnabled(False)
    self.main_frame.ui.action_export_project.setEnabled(False)
    self.main_frame.ui.action_import_project.setEnabled(False)
    self.main_frame.ui.menuPrediction.setEnabled(False)
    self.main_frame.ui.menuAnalysis.setEnabled(False)
    self.main_frame.ui.menuResults.setEnabled(False)
    self.main_frame.ui.menuImage.setEnabled(False)
    self.main_frame.ui.menuHotspots.setEnabled(False)
    self.main_frame.ui.menuSettings.setEnabled(False)
    self.main_frame.ui.menuAbout.setEnabled(False)

  def enable_job_panels(self) -> None:
    """Enables the job overview panel, job notification panel,and the corresponding buttons."""
    self.main_frame.ui.frame_job_overview.setEnabled(True)
    self.main_frame.ui.frame_job_notification.setEnabled(True)
    self.main_frame.btn_open_job_overview.setEnabled(True)
    self.main_frame.btn_open_job_notification.setEnabled(True)

  def switch_application_page(self) -> None:
    """Displays the project page."""
    if self.main_frame.ribbon_bar.currentIndex() == 0:
      # Switched to Project page
      self.main_frame.ui.stackedWidget.setCurrentIndex(1)
      self.main_frame.ribbon_bar.hide()
      self.main_frame.restore_new_project_tab()
      self.main_frame.restore_open_project_tab()
      self.main_frame.restore_delete_project_tab()
    else:
      # Switched to Main page
      self.main_frame.ui.stackedWidget.setCurrentIndex(0)
      self.main_frame.ribbon_bar.show()

  def display_main_page(self) -> None:
    """Displays the main page."""
    self.main_frame.ui.stackedWidget.setCurrentIndex(0)
    self.main_frame.ribbon_bar.setCurrentIndex(1)
    self.main_frame.ribbon_bar.show()
  # </editor-fold>

  # <editor-fold desc="Model related">
  def build_all_project_related_models(self) -> None:
    """Builds all project related models."""
    self.sequence_model = sequence_model.SequenceModel().from_a_list_of_sequences(self.project.sequences)
    self.main_frame.ui.seqs_list_view.setModel(self.sequence_model)
    self.protein_model = protein_model.ProteinModel().from_a_list_of_proteins(self.project.proteins)
    self.main_frame.ui.proteins_tree_view.setModel(self.protein_model)
    self.ligand_model = ligand_model.LigandModel().from_a_list_of_ligands(self.project.ligands)
    self.main_frame.ui.ligand_tree_view.setModel(self.ligand_model)

  def set_all_models_to_their_views(self) -> None:
    """Sets all models to their respective views."""
    self.main_frame.ui.seqs_list_view.setModel(self.sequence_model)
    self.main_frame.ui.proteins_tree_view.setModel(self.protein_model)
    self.main_frame.ui.ligand_tree_view.setModel(self.ligand_model)
  # </editor-fold>

  def close_project(self) -> None:
    """Closes the currently opened project."""
    self.databases_manager.get_handler(self.project.name).is_current_project = False
    self.project = None
    self.update_main_frame_gui()

  def shutdown_application_processes(self) -> None:
    """Closes all threads and process as well as the application itself."""
    self.main_frame.close()
    # Waiting for the database manager to finish the workload
    self.databases_manager.stop_all_handlers()  # TODO: It could be beneficial to use a separate thread for this
    # TODO: It might be good to show a message box that the app is still saving?!
    self.mv_client.stop_viewer()
    #time.sleep(5)  # This is for debug purposes
    self.mv_client.process.terminate()

  def exit_application(self) -> None:
    """Slot method for the Exit Application menu item"""
    self.shutdown_application_processes()

  def close_application(self, value: tuple[str, QtCore.QEvent]):
    """Closes all threads and process as well as the application itself."""
    self.shutdown_application_processes()
    # Closing QApplication
    value[1].accept()

  def schedule_tool_task_result_object(
          self,
          a_task_result: tuple[bool, task_result.TaskResult]
  ) -> None:
    """Receives the task result object of the tool and schedules it.

    Args:
      a_task_result: The task result object of the 'create project' tool

    Raises:
      exception.NoneValueError: If `a_task_result` is None.
    """
    # <editor-fold desc="Checks">
    if a_task_result is None:
      default_logging.append_to_log_file(logger, "a_task_result is None.", logging.ERROR)
      raise exception.NoneValueError("a_task_result is None.")
    # </editor-fold>
    if a_task_result[0]:
      self.task_manager.append_task_result(a_task_result[1])
      self.task_scheduler.schedule(a_task_result[1])

  def __slot_switch_object_tab(self) -> None:
    """Slot method for changing the tab widget with the project objects."""
    self.update_main_frame_gui()
  # </editor-fold>

  # <editor-fold desc="Proteins tab">
  def display_selection_string_for_the_current_protein(self) -> None:
    """Displays the current selection as selection string on the protein tab."""
    if len(self.main_frame.ui.proteins_tree_view.selectionModel().selectedIndexes()) == 1:
      self.main_frame.ui.lineEdit.setText(
        self.protein_model.get_selection_string(
          self.main_frame.ui.proteins_tree_view.currentIndex()
        )
      )
    else:
      self.main_frame.ui.lineEdit.setText("")
    self.main_frame.ui.information_table.verticalHeader().hide()
    self.main_frame.ui.information_table.horizontalHeader().hide()
    self.main_frame.ui.information_table.setModel(self.protein_model.get_info_as_table_model(self.main_frame.ui.proteins_tree_view.currentIndex()))
    self.main_frame.ui.information_table.resizeColumnsToContents()
    self.update_main_frame_gui()

  def select_based_on_protein_tree_selection_string(self) -> None:
    """Make a selection based on the entered selection string."""
    tmp_indexes = self.protein_model.get_model_index_based_on_name(self.main_frame.ui.lineEdit.text())
    if len(tmp_indexes) > 0:
      self.main_frame.ui.proteins_tree_view.selectionModel().clearSelection()
      for tmp_index in tmp_indexes:
        self.main_frame.ui.proteins_tree_view.selectionModel().select(
          tmp_index,
          QtCore.QItemSelectionModel.SelectionFlag.Select
        )
  # </editor-fold>

  # <editor-fold desc="Handling of the 'create project' component.">
  def slot_create_project_dialog(self) -> None:
    """Opens the 'create project' tool."""
    try:
      default_logging.append_to_log_file(logger, "Run tool 'Create project'.", logging.INFO)
      self.component_manager.components[model_definitions.ComponentsEnum.CREATE_PROJECT].connect_project_model_to_component(
        self.project_model
      )
      self.component_manager.components[model_definitions.ComponentsEnum.CREATE_PROJECT].run_tool(
        self.__await_component_create_project_task
      )
    except Exception as e:
      default_logging.append_to_log_file(
        logger,
        f"An error occurred: {e}",
        logging.ERROR
      )
      self.status_bar_manager.show_error_message("An error occurred.")

  def __await_component_create_project_task(self, value) -> None:
    """Adds the opened project to the main frame controller.

    Args:
      value: Result tuple of the finished task
    """
    # <editor-fold desc="Checks">
    if value is None:
      self.main_frame.ui.statusbar.showMessage(
        "The result value is None!")  # TODO: Needs to be replaced by a statusbar manager like in PySSA
      default_logging.append_to_log_file(logger, "The result value is None!", logging.ERROR)
      return
    if value[1][0] is False:
      self.main_frame.ui.statusbar.showMessage("Creating project failed.")  # TODO: Needs to be replaced by a statusbar manager like in PySSA
      default_logging.append_to_log_file(logger, "Creating project failed!", logging.ERROR)
      return
    # </editor-fold>
    try:
      tmp_results = task_result.TaskResult.get_single_action_result(value)
      if tmp_results[0]:
        if self.project is not None:
          self.close_project()
        self.project: "project.Project" = tmp_results[1][0]
        self.databases_manager.add_handler(
          self.project.name,
          database_handler.DatabaseHandler(
            pathlib.Path(
              f"{model_definitions.ModelDefinitions.DEFAULT_WORKSPACE_PATH}/{self.project.name}.db"
            )
          )
        )
        self.databases_manager.get_handler(self.project.name).is_current_project = True
        tmp_database_work_task = self.databases_manager.get_handler(self.project.name).build_new_database_and_change_to_it(
          pathlib.Path(
            f"{model_definitions.ModelDefinitions.DEFAULT_WORKSPACE_PATH}/{self.project.name}.db"
          )
        )
        self.task_manager.append_task(tmp_database_work_task)
        self.task_scheduler.schedule_task(tmp_database_work_task)
        self.build_all_project_related_models()
    except Exception as e:
      default_logging.append_to_log_file(
        logger,
        f"An error occurred: {e}",
        logging.ERROR
      )
      self.status_bar_manager.show_error_message("An error occurred.")
    finally:
      self.update_main_frame_gui()
      self.main_frame.ribbon_bar.setCurrentIndex(1)
  # </editor-fold>

  # <editor-fold desc="Handling of the 'open project' component">
  def slot_open_project_dialog(self) -> None:
    """Slot method for the open project menu item."""
    default_logging.append_to_log_file(logger, "Run tool 'open project'.", logging.INFO)
    try:
      # This is a special case method for the open project component which extends the IComponent interface!
      self.component_manager.components[model_definitions.ComponentsEnum.OPEN_PROJECT].connect_project_model_to_component(
        self.project_model
      )
      self.component_manager.components[model_definitions.ComponentsEnum.OPEN_PROJECT].run_tool(
        self.__await_component_open_project_task
      )
    except Exception as e:
      default_logging.append_to_log_file(
        logger,
        f"An error occurred: {e}",
        logging.ERROR
      )
      self.status_bar_manager.show_error_message("An error occurred.")

  def __await_component_open_project_task(self, value) -> None:
    """Adds the opened project to the main frame controller.

    Args:
      value: Result tuple of the finished task
    """
    # <editor-fold desc="Checks">
    if value is None:
      self.main_frame.ui.statusbar.showMessage(
        "The result value is None!")  # TODO: Needs to be replaced by a statusbar manager like in PySSA
      default_logging.append_to_log_file(logger, "The result value is None!", logging.ERROR)
      return
    if value[1][0] is False:
      self.main_frame.ui.statusbar.showMessage("Opening project failed.")  # TODO: Needs to be replaced by a statusbar manager like in PySSA
      default_logging.append_to_log_file(logger, "Opening project failed!", logging.ERROR)
      return
    # </editor-fold>
    try:
      tmp_results = task_result.TaskResult.get_single_action_result(value)
      if tmp_results[0]:
        if self.project is not None:
          self.databases_manager.get_handler(self.project.name).is_current_project = False
        self.project: "project.Project" = tmp_results[1][0]
        if self.project.name not in self.databases_manager.handlers.keys():
          """
          Checks if the databases_manager contains already the handler 
          for the given project name.
          The handler could already exist because it may had some inserts 
          to do while not being the active project and therefore still 
          in the databases manager as handler. 
          """
          self.databases_manager.add_handler(
            self.project.name,
            database_handler.DatabaseHandler(
              pathlib.Path(f"{model_definitions.ModelDefinitions.DEFAULT_WORKSPACE_PATH}/{self.project.name}.db")
            )
          )
        self.databases_manager.get_handler(self.project.name).is_current_project = True
        if self.databases_manager.get_handler(self.project.name).is_working is False:
          """
          Checks if the handler is working. 
          If its an old handler it might be not in the do_work loop and 
          needs therefore be restarted.
          Or the handler is already in the do_work loop and without checking
          this, multiple threads for the same handler could coexist 
          which may interfere each other.
          """
          tmp_database_work_task = self.databases_manager.get_handler(self.project.name).create_do_work_task()
          self.task_manager.append_task(tmp_database_work_task)
          self.task_scheduler.schedule_task(tmp_database_work_task)
        tmp_task = task_result.TaskResult.from_action(
          an_action=action.Action(
            a_target=self.__async_fetch_project_from_database,
            args=(
              self.project.name, self.databases_manager.get_handler(self.project.name)
            )
          ),
          an_await_function=self.__await_fetch_project_from_database
        )
        self.task_manager.append_task_result(tmp_task)
        self.task_scheduler.schedule_task_result(tmp_task)
        self.block_gui(with_wait_cursor=True)
      else:
        self.update_main_frame_gui()
    except Exception as e:
      default_logging.append_to_log_file(
        logger,
        f"An error occurred: {e}",
        logging.ERROR
      )
      self.status_bar_manager.show_error_message("An error occurred.")

  def __async_fetch_project_from_database(self, a_project_name: str, a_database_handler: "database_handler.DatabaseHandler") -> "project.Project":
    """ASYNC method that fetches all data for the project from the database.

    Args:
      a_project_name: The name of the project
      a_database_handler: The database handler of the project to fetch

    Raises:
      exception.NoneValueError: If any of the arguments are None.
      exception.IllegalArgumentError: If `a_project_name` is an empty string.
    """
    # <editor-fold desc="Checks">
    if a_project_name is None:
      default_logging.append_to_log_file(logger, "a_project_name is None.", logging.ERROR)
      raise exception.NoneValueError("a_project_name is None.")
    if a_project_name == "":
      default_logging.append_to_log_file(logger, "a_project_name is an empty string.", logging.ERROR)
      raise exception.IllegalArgumentError("a_project_name is an empty string.")
    if a_database_handler is None:
      default_logging.append_to_log_file(logger, "a_database_handler is None.", logging.ERROR)
      raise exception.NoneValueError("a_database_handler is None.")
    # </editor-fold>
    return a_database_handler.get_complete_project(a_project_name)

  def __await_fetch_project_from_database(self, a_result: tuple[str, list]) -> None:
    """Finishes up the open project process."""
    try:
      if a_result[0]:
        self.project = task_result.TaskResult.get_single_action_result(a_result)[1]
        self.build_all_project_related_models()
        self.databases_manager.remove_unused_handlers()  # This must run after a GUI-reorganization method!
    except Exception as e:
      default_logging.append_to_log_file(
        logger,
        f"An error occurred: {e}",
        logging.ERROR
      )
      self.status_bar_manager.show_error_message("An error occurred.")
    finally:
      self.main_frame.ribbon_bar.setCurrentIndex(1)
      self.update_main_frame_gui()
      self.stop_wait_cursor()
      self.unblock_gui()

  # </editor-fold>

  # <editor-fold desc="Main page">
  # <editor-fold desc="Handling of the 'import sequence' component.">
  def slot_import_sequence_dialog(self):
    """Slot method for the import sequence button."""
    try:
      self.component_manager.components[model_definitions.ComponentsEnum.IMPORT_SEQUENCE].run_tool(
        self.__await_component_import_sequence_task
      )
    except Exception as e:
      default_logging.append_to_log_file(
        logger,
        f"An error occurred: {e}",
        logging.ERROR
      )
      self.status_bar_manager.show_error_message("An error occurred.")

  def __await_component_import_sequence_task(self, value) -> None:
    """Adds the imported sequence to the project object and model.

    Args:
      value: Result tuple of the finished task
    """
    # <editor-fold desc="Checks">
    if value is None:
      self.main_frame.ui.statusbar.showMessage(
        "The result value is None!")  # TODO: Needs to be replaced by a statusbar manager like in PySSA
      default_logging.append_to_log_file(logger, "The result value is None!", logging.ERROR)
      return
    if value[1][0] is False:
      self.main_frame.ui.statusbar.showMessage(
        "Opening project failed.")  # TODO: Needs to be replaced by a statusbar manager like in PySSA
      default_logging.append_to_log_file(logger, "Opening project failed!", logging.ERROR)
      return
    # </editor-fold>
    try:
      tmp_results = task_result.TaskResult.get_single_action_result(value)
      if tmp_results[0]:
        self.sequence_model.add_sequence(tmp_results[1])  # Can be ignored, because it's a project object not a tuple type
    except Exception as e:
      default_logging.append_to_log_file(
        logger,
        f"An error occurred: {e}",
        logging.ERROR
      )
      self.status_bar_manager.show_error_message("An error occurred.")
    finally:
      self.update_main_frame_gui()
  # </editor-fold>

  # <editor-fold desc="Handling of the 'import protein' component.">
  def slot_import_protein_dialog(self) -> None:
    """Slot method for the import protein button."""
    try:
      self.component_manager.components[model_definitions.ComponentsEnum.IMPORT_PROTEIN].run_tool(
        self.__await_component_import_protein_task
      )
    except Exception as e:
      default_logging.append_to_log_file(
        logger,
        f"An error occurred: {e}",
        logging.ERROR
      )
      self.status_bar_manager.show_error_message("An error occurred.")

  def __await_component_import_protein_task(self, value) -> None:
    """Adds the imported protein to the project object and model.

    Args:
      value: Result tuple of the finished task
    """
    # <editor-fold desc="Checks">
    if value is None:
      self.main_frame.ui.statusbar.showMessage(
        "The result value is None!")  # TODO: Needs to be replaced by a statusbar manager like in PySSA
      default_logging.append_to_log_file(logger, "The result value is None!", logging.ERROR)
      return
    if value[1][0] is False:
      self.main_frame.ui.statusbar.showMessage(
        "Opening project failed.")  # TODO: Needs to be replaced by a statusbar manager like in PySSA
      default_logging.append_to_log_file(logger, "Opening project failed!", logging.ERROR)
      return
    # </editor-fold>
    try:
      tmp_results = task_result.TaskResult.get_single_action_result(value)
      if tmp_results[0]:
        tmp_protein = tmp_results[1][0]
        self.project.add_protein(tmp_protein)
        self.protein_model.add_protein(tmp_protein)
        self.databases_manager.get_handler(self.project.name).insert_new_protein(
          tmp_protein
        )
    except Exception as e:
      default_logging.append_to_log_file(
        logger,
        f"An error occurred: {e}",
        logging.ERROR
      )
      self.status_bar_manager.show_error_message("An error occurred.")
    finally:
      self.update_main_frame_gui()
  # </editor-fold>

  # <editor-fold desc="Handling of the 'import ligand' component.">
  def slot_import_ligand_dialog(self) -> None:
    """Slot method for the import ligand button."""
    try:
      self.component_manager.components[model_definitions.ComponentsEnum.IMPORT_LIGAND].run_tool(
        self.__await_component_import_ligand_task
      )
    except Exception as e:
      default_logging.append_to_log_file(
        logger,
        f"An error occurred: {e}",
        logging.ERROR
      )
      self.status_bar_manager.show_error_message("An error occurred.")

  def __await_component_import_ligand_task(self, value) -> None:
    """Adds the imported ligand to the project object and model.

    Args:
      value: Result tuple of the finished task
    """
    # <editor-fold desc="Checks">
    if value is None:
      self.main_frame.ui.statusbar.showMessage(
        "The result value is None!")  # TODO: Needs to be replaced by a statusbar manager like in PySSA
      default_logging.append_to_log_file(logger, "The result value is None!", logging.ERROR)
      return
    if value[1][0] is False:
      self.main_frame.ui.statusbar.showMessage(
        "Opening project failed.")  # TODO: Needs to be replaced by a statusbar manager like in PySSA
      default_logging.append_to_log_file(logger, "Opening project failed!", logging.ERROR)
      return
    # </editor-fold>
    try:
      tmp_results = task_result.TaskResult.get_single_action_result(value)
      if tmp_results[0]:
        self.ligand_model.add_ligand(tmp_results[1])  # TODO: Cannot be ignored!
    except Exception as e:
      default_logging.append_to_log_file(
        logger,
        f"An error occurred: {e}",
        logging.ERROR
      )
      self.status_bar_manager.show_error_message("An error occurred.")
    finally:
      self.update_main_frame_gui()
  # </editor-fold>
  # </editor-fold>

  # <editor-fold desc="Project Page">
  def display_project_tab(self):
    """Displays a tab on the project page."""
    if self.main_frame.ui.tab_widget_project_page.currentIndex() == 1:
      self.main_frame.restore_new_project_tab()
    elif self.main_frame.ui.tab_widget_project_page.currentIndex() == 2:
      self.main_frame.restore_open_project_tab()
    elif self.main_frame.ui.tab_widget_project_page.currentIndex() == 3:
      self.main_frame.restore_delete_project_tab()
    elif self.main_frame.ui.tab_widget_project_page.currentIndex() == 4:
      self.close_project()
    elif self.main_frame.ui.tab_widget_project_page.currentIndex() == 8:
      self.exit_application()

  # <editor-fold desc="New project tab">
  def open_create_project_dialog(self):
    self.basic_controllers["New"].restore_ui()
    self.basic_controllers["New"].get_dialog().exec()
    tmp_task_result = task_result.TaskResult.from_action(
      action.Action(
        a_target=self.async_create_project
      ),
      self.__await_component_create_project_task
    )
    if self.basic_controllers["New"].was_canceled is False:
      self.basic_controllers["New"].component_task.emit((True, tmp_task_result))
    else:
      self.basic_controllers["New"].component_task.emit((False, tmp_task_result))

  def async_create_project(self) -> tuple["project.Project"]:
    """Creates the project based on the entered project name."""
    tmp_project = project.Project(self.basic_controllers["New"].entered_project_name)
    return tmp_project,
  # </editor-fold>

  # <editor-fold desc="Open project tab">
  def __slot_open_project(self, a_value: tuple) -> None:
    """Slot method for the open project button."""
    tmp_project_name, tmp_date_modified = a_value
    if self.project is not None:
      self.close_project()
    self.selected_project_name = tmp_project_name
    tmp_task_result = task_result.TaskResult.from_action(
      action.Action(a_target=self.async_open_project),
      self.__await_component_open_project_task
    )
    self.task_manager.append_task_result(tmp_task_result)
    self.task_scheduler.schedule(tmp_task_result)

  def async_open_project(self) -> tuple["project.Project"]:
    """Opens the selected project"""
    try:
      return project.Project(self.selected_project_name),
    except Exception as e:
      default_logging.append_to_log_file(logger, e.__str__(), logging.ERROR)
      raise exception.AsyncOperationFailedError(e.__str__())
    finally:
      default_logging.append_to_log_file(logger, "'open_project' method finished.", logging.DEBUG)
  # </editor-fold>

  # </editor-fold>
