import logging
from typing import Optional

from PyQt6 import QtCore, QtGui, QtWidgets

from pyssa.gui.custom_widgets import ribbon_bar, inline_menu, side_tabs, \
  custom_button, line_edit_search
from pyssa.model.data_classes import workspace_project
from pyssa.model.preference import model_definitions
from pyssa.model.util.gui_style import icons
from pyssa.model.util.gui_style import styles_utils
from pyssa.gui.main.forms.auto import auto_main_frame
from pyssa.gui.custom_widgets import color_grid
from pyssa.gui.custom_widgets import popup_dialog
from pyssa.model.util import exception
from pyssa.model.pyssa_logging import default_logging

logger = default_logging.setup_logger(__file__)


class MainFrame(QtWidgets.QMainWindow):
  """Main frame class."""

  dialogClosed = QtCore.pyqtSignal(tuple)
  """A signal indicating that the dialog is closed."""

  def __init__(self) -> None:
    """Constructor."""
    super().__init__()
    # build ui object
    self.ui = auto_main_frame.Ui_MainWindow()
    self.ui.setupUi(self)
    # Custom widgets
    self.color_gird = color_grid.ColorGrid()
    self.ribbon_bar = ribbon_bar.RibbonBar()
    self.btn_blank_project = custom_button.BigCardButton("Blank Project")
    self.ui.btn_scene_name.clicked.connect(self.show_menu)
    # Init gui
    self.init_ui()
    # Add experiment widgets
    self.line_edit_search = line_edit_search.LineEditSearch(self)
    self.ui.verticalLayout_55.insertWidget(0, self.line_edit_search)
    self.installEventFilter(self)

  def eventFilter(self, obj, event):
    if event.type() == QtCore.QEvent.Type.MouseButtonPress:
      # Close the dialog if the click is outside the dialog and QLineEdit
      if self.line_edit_search.search_dialog.isVisible():
        self.line_edit_search.search_dialog.hide()
    return super().eventFilter(obj, event)

  # <editor-fold desc="Basic">
  def closeEvent(self, event) -> None:  # noqa: ANN001
    """Overrides the closeEvent of the QMainWindow class.

    Args:
        event: The event object representing the close event.
    """
    # <editor-fold desc="Checks">
    # if event is None:
    #   logger.error("event is None.")
    #   raise exception.IllegalArgumentError("event is None.")

    # </editor-fold>
    # Emit the custom signal when the window is closed
    self.dialogClosed.emit(("", event))

  def init_ui(self) -> None:
    """Initialize the UI elements."""
    # <editor-fold desc="Hide job related content">
    self.ui.frame_job_overview.hide()
    self.ui.frame_job_notification.hide()
    # </editor-fold>

    self.ui.proteins_tree_view.setHeaderHidden(True)
    self.ui.proteins_tree_view.setEditTriggers(
      QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
    )
    self.ui.proteins_tree_view.setSelectionMode(
      QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection
    )
    self.add_custom_job_panels()
    self.add_custom_color_grid()
    self.add_custom_pop_up_dialogs()
    self.add_custom_ribbon()
    self.add_blank_project_button()
    self.set_icons()
    self.set_tab_texts_for_project_page(
      [
        ("Home", model_definitions.IconsEnum.HOME),
        ("New", model_definitions.IconsEnum.NEW),
        ("Open", model_definitions.IconsEnum.OPEN),
        ("Delete", model_definitions.IconsEnum.DELETE),
        ("Close", model_definitions.IconsEnum.CLOSE),
        ("Import", model_definitions.IconsEnum.IMPORT),
        ("Export", model_definitions.IconsEnum.EXPORT),
        ("Settings", None),
        ("Exit", None)
      ]
    )
    self.set_custom_stylesheet_for_project_page_tab_widget()
    styles_utils.set_stylesheet(self)
  # </editor-fold>

  # <editor-fold desc="Setup">
  def set_tab_texts_for_project_page(self, tab_texts: list[tuple[str, Optional["model_definitions.IconsEnum"]]]) -> None:
    """Sets the tab text for all tabs on the project page."""
    for i, tmp_tab_text in enumerate(tab_texts):
      self.ui.tab_widget_project_page.tabBar().setTabText(i, "")
      self.ui.tab_widget_project_page.tabBar().setTabButton(
        i, QtWidgets.QTabBar.ButtonPosition.LeftSide,
        side_tabs.SideTabs(tmp_tab_text[0], tmp_tab_text[1])
      )

  def set_custom_stylesheet_for_project_page_tab_widget(self):
    """Sets custom stylesheet for the project page tab widget."""
    self.ui.tab_widget_project_page.setStyleSheet(
      """
      QTabWidget::pane { /* The tab widget frame */
          border: 1px solid #DCDBE3;
          background: #f5f5f5;
          margin: 0px;
          padding: 0px;
      }

      QTabWidget::tab-bar {
          left: 0px; /* move to the right by 5px */
      }

      QTabBar::tab {
          border: none;
          border-radius: 0px;
          background: #f0f0f0;
          min-width: 70px;
          max-width: 70px;
          font-size: 13px;
          padding: 0px;
          padding-top: 5px;
          padding-bottom: 5px;
          padding-right: 30px;
      }

      QTabBar::tab:selected, QTabBar::tab:hover {
          background: #e6e6e6;
          border: 2px solid #B8B8B8;
          border-top: none;
          border-bottom: none;
          border-right: none;
      }

      QTabBar::tab:selected {
          border-color: #367AF6;
          background: #e6e6e6;
          font: bold;
      }

      QTabBar::tab:!selected {
          margin-top: 0px
      }

      QToolButton {
          background-color: white;
          /*min-width: 1.5em;*/
          min-width: 24px;
          min-height: 24px;
          max-height: 64px;
          padding: 0px;
          border: none;
      }

      QToolButton::hover {
          background: #f5f5f5;
          color: black;
      }
      """
    )

  def set_icons(self) -> None:
    """Sets all icons for the main frame."""
    # <editor-fold desc="Help">
    icons.set_icon(self.ui.btn_help, model_definitions.IconsEnum.HELP)
    icons.set_icon(self.ui.btn_help_2, model_definitions.IconsEnum.HELP)
    icons.set_icon(self.ui.btn_help_3, model_definitions.IconsEnum.HELP)
    icons.set_icon(self.ui.btn_help_4, model_definitions.IconsEnum.HELP)
    icons.set_icon(self.ui.btn_help_5, model_definitions.IconsEnum.HELP)
    icons.set_icon(self.ui.btn_help_6, model_definitions.IconsEnum.HELP)
    icons.set_icon(self.ui.btn_help_7, model_definitions.IconsEnum.HELP)
    icons.set_icon(self.ui.btn_help_8, model_definitions.IconsEnum.HELP)
    icons.set_icon(self.ui.btn_help_info, model_definitions.IconsEnum.HELP)
    icons.set_icon(self.ui.btn_help_session, model_definitions.IconsEnum.HELP)
    # </editor-fold>
    # <editor-fold desc="Sequence tab">
    # import
    icons.set_icon(
      self.ui.btn_import_seq,
      model_definitions.IconsEnum.IMPORT_SEQUENCE,
      model_definitions.IconsEnum.IMPORT_SEQUENCE_DISABLED,
    )
    # add
    icons.set_icon(
      self.ui.btn_add_sequence,
      model_definitions.IconsEnum.ADD_SEQUENCE,
      model_definitions.IconsEnum.ADD_SEQUENCE_DISABLED,
    )
    # save sequence
    icons.set_icon(
      self.ui.btn_save_sequence,
      model_definitions.IconsEnum.SAVE_SEQUENCE,
      model_definitions.IconsEnum.SAVE_SEQUENCE_DISABLED,
    )

    # delete sequence
    icons.set_icon(
      self.ui.btn_delete_sequence,
      model_definitions.IconsEnum.DELETE_SEQUENCE,
      model_definitions.IconsEnum.DELETE_SEQUENCE_DISABLED,
    )
    # </editor-fold>
    # <editor-fold desc="Proteins tab">
    # expand all
    icons.set_icon(
      self.ui.btn_protein_tree_view_expand,
      model_definitions.IconsEnum.EXPAND_ALL,
      size=QtCore.QSize(20, 20)
    )
    # collapse all
    icons.set_icon(
      self.ui.btn_protein_tree_view_collapse,
      model_definitions.IconsEnum.COLLAPSE_ALL,
      size=QtCore.QSize(20, 20)
    )
    # import
    icons.set_icon(
      self.ui.btn_import_protein,
      model_definitions.IconsEnum.IMPORT_PROTEIN,
      model_definitions.IconsEnum.IMPORT_PROTEIN_DISABLED,
    )
    # save
    icons.set_icon(
      self.ui.btn_save_protein,
      model_definitions.IconsEnum.SAVE_PROTEIN,
      model_definitions.IconsEnum.SAVE_PROTEIN_DISABLED,
    )
    # delete
    icons.set_icon(
      self.ui.btn_delete_protein,
      model_definitions.IconsEnum.DELETE_PROTEIN,
      model_definitions.IconsEnum.DELETE_PROTEIN_DISABLED,
    )
    # </editor-fold>
    # <editor-fold desc="Ligand tab">
    # expand all
    icons.set_icon(
      self.ui.btn_ligand_tree_view_expand,
      model_definitions.IconsEnum.EXPAND_ALL,
      size=QtCore.QSize(20, 20)
    )
    # collapse all
    icons.set_icon(
      self.ui.btn_ligand_tree_view_collapse,
      model_definitions.IconsEnum.COLLAPSE_ALL,
      size=QtCore.QSize(20, 20)
    )
    # import
    icons.set_icon(
      self.ui.btn_import_ligand,
      model_definitions.IconsEnum.IMPORT_PROTEIN,
      model_definitions.IconsEnum.IMPORT_PROTEIN_DISABLED,
    )
    # save
    icons.set_icon(
      self.ui.btn_save_ligand,
      model_definitions.IconsEnum.SAVE_PROTEIN,
      model_definitions.IconsEnum.SAVE_PROTEIN_DISABLED,
    )
    # delete
    icons.set_icon(
      self.ui.btn_delete_ligand,
      model_definitions.IconsEnum.DELETE_PROTEIN,
      model_definitions.IconsEnum.DELETE_PROTEIN_DISABLED,
    )
    # </editor-fold>
    # <editor-fold desc="Protein-Protein complex tab">
    # expand all
    icons.set_icon(
      self.ui.btn_prot_prot_complex_expand,
      model_definitions.IconsEnum.EXPAND_ALL,
      size=QtCore.QSize(20, 20)
    )
    # collapse all
    icons.set_icon(
      self.ui.btn_prot_prot_complex_collapse,
      model_definitions.IconsEnum.COLLAPSE_ALL,
      size=QtCore.QSize(20, 20)
    )
    # delete
    icons.set_icon(
      self.ui.btn_delete_prot_prot_complex,
      model_definitions.IconsEnum.DELETE_PROTEIN_PAIR,
      model_definitions.IconsEnum.DELETE_PROTEIN_PAIR_DISABLED,
    )
    # </editor-fold>
    # <editor-fold desc="Protein-Ligand complex tab">
    # expand all
    icons.set_icon(
      self.ui.btn_prot_lig_complex_expand,
      model_definitions.IconsEnum.EXPAND_ALL,
      size=QtCore.QSize(20, 20)
    )
    # collapse all
    icons.set_icon(
      self.ui.btn_prot_lig_complex_collapse,
      model_definitions.IconsEnum.COLLAPSE_ALL,
      size=QtCore.QSize(20, 20)
    )
    # delete
    icons.set_icon(
      self.ui.btn_delete_prot_lig_complex,
      model_definitions.IconsEnum.DELETE_PROTEIN_PAIR,
      model_definitions.IconsEnum.DELETE_PROTEIN_PAIR_DISABLED,
    )
    # </editor-fold>
    # <editor-fold desc="Session tab">
    # open
    icons.set_icon(
      self.ui.btn_open_session,
      model_definitions.IconsEnum.OPEN_SESSION,
      model_definitions.IconsEnum.OPEN_SESSION_DISABLED,
    )
    # create
    icons.set_icon(
      self.ui.btn_create_scene,
      model_definitions.IconsEnum.CREATE_SESSION_SCENE,
      model_definitions.IconsEnum.CREATE_SESSION_SCENE_DISABLED,
    )
    # update
    icons.set_icon(
      self.ui.btn_update_scene,
      model_definitions.IconsEnum.UPDATE_SESSION_SCENE,
      model_definitions.IconsEnum.UPDATE_SESSION_SCENE_DISABLED,
    )
    # delete scene
    icons.set_icon(
      self.ui.btn_delete_scene,
      model_definitions.IconsEnum.DELETE_SESSION_SCENE,
      model_definitions.IconsEnum.DELETE_SESSION_SCENE_DISABLED,
    )
    # </editor-fold>

    icons.set_icon(self.ui.btn_project_page_back, model_definitions.IconsEnum.BACK)
    self.ui.btn_project_page_back.setStyleSheet(
      """
      QPushButton#btn_project_page_back {
          background-color: white;
          border: 1px solid #5f6368;
          border-radius: 13px;
          border-color: #5f6368;
          padding: 2px;
          margin-top: 10px;
          margin-left: 10px;
          margin-right: 20px;
          min-width: 20px;
          max-width: 20px;
          min-height: 20px;
          max-height: 20px
      }
      """
    )

  # </editor-fold>

  # <editor-fold desc="Add custom widgets">
  def add_custom_job_panels(self) -> None:
    """Add custom job panels.

    Creates and configures the GUI elements for the job panels, including buttons, labels, and icons.
    """
    # job panel related gui elements
    self.btn_open_job_overview = QtWidgets.QPushButton()
    self.lbl_job_overview = QtWidgets.QLabel("No running jobs.")
    self.btn_open_job_notification = QtWidgets.QPushButton()
    self.lbl_job_notification = QtWidgets.QLabel("No new notifications.")
    icons.set_icon(self.btn_open_job_overview, model_definitions.IconsEnum.JOBS)
    icons.set_icon(self.btn_open_job_notification, model_definitions.IconsEnum.NOTIFY)
    self.btn_open_job_overview.setStyleSheet(
      """
          QPushButton {
              background-color: rgba(220, 219, 227, 0.01);
              border: none;
              border-radius: 4px;
              min-width: 20px;
              max-width: 20px;
              min-height: 20px;
              max-height: 20px;
          }
          QPushButton::hover {
              background-color: rgba(220, 219, 227, 0.5);
              border: none;
              min-width: 24px;
              max-width: 24px;
              min-height: 24px;
              max-height: 24px;
          }
      """
    )
    self.lbl_job_overview.setStyleSheet("""padding-top: 30px;""")
    self.btn_open_job_notification.setStyleSheet(
      """
          QPushButton {
              background-color: rgba(220, 219, 227, 0.01);
              border: none;
              border-radius: 4px;
              min-width: 20px;
              max-width: 20px;
              min-height: 20px;
              max-height: 20px;
          }
          QPushButton::hover {
              background-color: rgba(220, 219, 227, 0.5);
              border: none;
              min-width: 24px;
              max-width: 24px;
              min-height: 24px;
              max-height: 24px;
          }
      """
    )
    """QPushButton {
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
            }"""
    self.lbl_job_notification.setStyleSheet("""padding-top: 30px;""")
    self.ui.job_overview_layout.insertWidget(
      0, self.lbl_job_overview
    )  # After inserting the widget count is 2
    self.ui.job_overview_layout.setAlignment(
      self.lbl_job_overview, QtCore.Qt.AlignmentFlag.AlignHCenter
    )
    self.ui.job_notification_layout.insertWidget(
      0, self.lbl_job_notification
    )  # After inserting the widget count is 2
    self.ui.job_notification_layout.setAlignment(
      self.lbl_job_notification, QtCore.Qt.AlignmentFlag.AlignHCenter
    )

  def add_custom_color_grid(self) -> None:
    """Adds the custom color grid widget to the main frame."""
    self.ui.verticalLayout_46.insertWidget(1, self.color_gird)

  def add_custom_pop_up_dialogs(self):
    """Adds custom pop-up dialogs to the main frame."""
    # Add custom scene dropdown
    self.drop_down = popup_dialog.ScenesPopUpDialog(self)
    self.import_popup = popup_dialog.ImportPopUpDialog(self)

  def add_custom_ribbon(self):
    """Adds a custom Microsoft style ribbon widget to the main frame."""
    self.ui.verticalLayout.insertWidget(0, self.ribbon_bar)

  def add_blank_project_button(self):
    self.ui.scroll_area_new_project_layout.insertWidget(0, self.btn_blank_project)

  def fill_open_project_scroll_area(self, the_project_model: QtGui.QStandardItemModel, a_callable):
    for i in range(the_project_model.rowCount()):
      tmp_workspace_project: "workspace_project.WorkspaceProject" = the_project_model.data(the_project_model.index(i, 0), model_definitions.RolesEnum.OBJECT_ROLE)
      self.ui.scroll_area_open_project_layout.insertWidget(
        i + 1,
        tmp_workspace_project.get_as_button(a_callable)
      )
    self.ui.lbl_scroll_area_open_project_name.setStyleSheet("margin-left: 10px;")
    self.ui.lbl_scroll_area_open_project_date_modified.setStyleSheet("margin-right: 25px;")

  def fill_home_project_scroll_area(self, the_project_model: QtGui.QStandardItemModel, a_callable):
    for i in range(the_project_model.rowCount()):
      if i < 11:
        tmp_workspace_project: "workspace_project.WorkspaceProject" = the_project_model.data(the_project_model.index(i, 0), model_definitions.RolesEnum.OBJECT_ROLE)
        self.ui.scroll_area_home_project_layout.insertWidget(
          i + 1,
          tmp_workspace_project.get_as_button(a_callable)
        )
    self.ui.lbl_scroll_area_home_project_name.setStyleSheet("margin-left: 10px;")
    self.ui.lbl_scroll_area_home_project_date_modified.setStyleSheet("margin-right: 25px;")
  # </editor-fold>

  # <editor-fold desc="Pop-up related">
  def show_menu(self) -> None:
    """Shows the dropdown menu for the scene selection."""
    # Determine the position to show the dialog (beneath the button)
    button_pos = self.ui.btn_scene_name.mapToGlobal(QtCore.QPoint(0, self.ui.btn_scene_name.height()))
    self.drop_down.move(button_pos)
    self.drop_down.show()

  def show_import_popup(self) -> None:
    """Shows the dropdown menu for the scene selection."""
    # Show the dialog momentarily to ensure its size is calculated
    self.import_popup.adjustSize()
    # Get the button's position in global coordinates
    button_pos = self.ui.btn_import_seq.mapToGlobal(QtCore.QPoint(0, 0))
    # Subtract the dialog's height to position it above the button
    dialog_height = self.import_popup.height()
    adjusted_pos = button_pos - QtCore.QPoint(0, dialog_height)
    # Move the dialog to the adjusted position
    self.import_popup.move(adjusted_pos)
    self.import_popup.show()
  # </editor-fold>

  def restore_new_project_tab(self):
    """Restores the 'New' project tab on the project page."""
    icons.set_icon(self.ui.btn_help_5, model_definitions.IconsEnum.HELP)

  def restore_open_project_tab(self):
    """Restores the 'Open' project tab on the project page."""
    icons.set_icon(self.ui.btn_help_6, model_definitions.IconsEnum.HELP)

  def restore_delete_project_tab(self):
    """Restores the 'Delete' project tab on the project page."""
    self.ui.lbl_delete_status_search.setText("")
    self.ui.btn_delete_delete_project.setEnabled(False)
    self.ui.list_delete_projects_view.setEditTriggers(
      QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers
    )
    self.ui.list_delete_projects_view.setSelectionMode(
      QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection
    )
    self.ui.btn_help_6.setIcon(QtGui.QIcon(":/icons/help_w200.png"))
    self.ui.btn_help_6.setIconSize(
      self.ui.btn_help_6.icon().actualSize(QtCore.QSize(30, 30))
    )
    self.ui.btn_help_6.setText("")
