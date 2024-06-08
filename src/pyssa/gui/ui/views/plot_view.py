#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/zielesny/PySSA>
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
"""Module for the plot view."""
import copy
import logging
from typing import Optional

import numpy as np
from PyQt5.QtCore import Qt
from PyQt5 import QtGui
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5.QtWidgets import QVBoxLayout, QWidget, QScrollArea
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backend_bases import MouseButton
from matplotlib.figure import Figure
from matplotlib import ticker

from src.pyssa.controller import help_manager
from src.pyssa.gui.ui.styles import styles
from src.pyssa.gui.ui.views import histogram_properties_view
from src.pyssa.internal.data_structures import protein_pair
from src.pyssa.internal.thread import tasks
from src.pyssa.internal.thread.async_pyssa import util_async
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import pyssa_keys, enums
from src.pyssa.util import constants
from src.pyssa.util import exception

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)
__docformat__ = "google"


class PlotWidget(QWidget):
  """Class for a custom QWidget for plotting."""

  def __init__(self, parent=None) -> None:  # noqa: ANN001
    """Constructor.

    Args:
        parent: the parent
    """
    super(PlotWidget, self).__init__(parent)
    self.figure = Figure(figsize=(9, 5.25))
    self.canvas = FigureCanvas(self.figure)
    layout = QVBoxLayout()
    layout.addWidget(self.canvas)
    self.setLayout(layout)

  def set_figure_size(self, width: float, height: float) -> None:
    """Set the size of the figure.

    Args:
        width (float): The width of the figure.
        height (float): The height of the figure.

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None or have a value less than zero.
    """
    # <editor-fold desc="Checks">
    if width is None or width < 0:
      logger.error("width is either None or a value less than 0.")
      raise exception.IllegalArgumentError(
          "width is either None or a value less than 0."
      )
    if height is None or height < 0:
      logger.error("height is either None or a value less than 0.")
      raise exception.IllegalArgumentError(
          "height is either None or a value less than 0."
      )

    # </editor-fold>

    self.figure.set_size_inches(width, height)


class PlotView(QtWidgets.QDialog):
  """Dialog window for the distance data visualizer."""

  def __init__(
      self,
      protein_pair_from_project: "protein_pair.ProteinPair",
      a_project: "project.Project",
      the_protein_pair: "protein_pair.ProteinPair",
      the_help_manager: "help_manager.HelpManager",
      parent=None,
  ) -> None:  # noqa: ANN001
    """Constructor.

    Args:
        protein_pair_from_project (protein_pair.ProteinPair): The protein pair object obtained from the project.
        a_project (project.Project): The project object.
        the_protein_pair (protein_pair.ProteinPair): Another protein pair object.
        parent: The parent widget (optional).

    Raises:
        exception.IllegalArgumentError: If any of the arguments are None.
    """
    # <editor-fold desc="Checks">
    if protein_pair_from_project is None:
      logger.error("protein_pair_from_project is None.")
      raise exception.IllegalArgumentError("protein_pair_from_project is None.")
    if a_project is None:
      logger.error("a_project is None.")
      raise exception.IllegalArgumentError("a_project is None.")
    if the_protein_pair is None:
      logger.error("the_protein_pair is None.")
      raise exception.IllegalArgumentError("the_protein_pair is None.")

    # </editor-fold>

    QtWidgets.QDialog.__init__(self, parent)
    self._protein_pair = the_protein_pair
    self._current_project = a_project
    self._help_manager = the_help_manager
    self.clicked_point_scatter = None
    self.highlighted_bin_index = None
    self.active_row_information = None
    self._sync_with_pymol_flag = False
    self.bars = None
    self._histogram_properties = {
        enums.HistogramPropertiesEnum.X_AXIS_UNITS: 10,
        enums.HistogramPropertiesEnum.DISTANCE_INTERVAL: 1.0,
    }

    # self.resizeEvent = self.handle_resize
    # Create a timer for delayed updates
    self.resize_timer = QtCore.QTimer(self)
    self.resize_timer.setInterval(250)  # Set the interval in milliseconds
    self.resize_timer.setSingleShot(True)

    self.protein_pair_for_analysis: protein_pair.ProteinPair = (
        protein_pair_from_project
    )
    # custom_pyssa_styles.set_stylesheet(self)

    self._initialize_ui()

    self.setModal(False)
    # --BEGIN

    # <editor-fold desc="Distance table logic">
    self.csv_model = QtGui.QStandardItemModel()
    self.csv_model.setColumnCount(7)
    labels = [
        "Residue pair no.",
        "Protein 1 Chain",
        "Protein 1 Position",
        "Protein 1 Residue",
        "Protein 2 Chain",
        "Protein 2 Position",
        "Protein 2 Residue",
        "Distance in Å",
    ]
    self.csv_model.setHorizontalHeaderLabels(labels)
    self.table_view.setModel(self.csv_model)

    # csv_filepath = pathlib.Path(f"{constants.CACHE_CSV_DIR}/{self._protein_pair.name}.csv")
    # if not os.path.exists(constants.CACHE_CSV_DIR):
    #     os.mkdir(constants.CACHE_CSV_DIR)
    #
    # distance_data = self._protein_pair.distance_analysis.analysis_results.distance_data
    # distance_data_array = np.array(
    #     [
    #         distance_data[pyssa_keys.ARRAY_DISTANCE_INDEX],
    #         distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN],
    #         distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION],
    #         distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI],
    #         distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN],
    #         distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION],
    #         distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI],
    #         distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES],
    #     ],
    # )
    # distance_data_array_transpose = distance_data_array.transpose()
    # with open(csv_filepath, mode="w", newline="") as file:
    #     writer = csv.writer(file, delimiter=",")
    #     writer.writerows(distance_data_array_transpose)
    #
    # file.close()
    import decimal

    # with open(csv_filepath, "r", encoding="utf-8") as csv_file:
    #     i = 0
    #     for line in csv_file:
    #         tmp_list = line.split(",")
    #         # tmp_list.pop(0)
    #         standard_item_list = []
    #         pair_no_item = QtGui.QStandardItem()
    #         pair_no_item.setData(int(tmp_list[0]), role=QtCore.Qt.DisplayRole)
    #         ref_chain_item = QtGui.QStandardItem()
    #         ref_chain_item.setData(str(tmp_list[1]), role=QtCore.Qt.DisplayRole)
    #         ref_pos_item = QtGui.QStandardItem()
    #         ref_pos_item.setData(int(tmp_list[2]), role=QtCore.Qt.DisplayRole)
    #         ref_resi_item = QtGui.QStandardItem()
    #         ref_resi_item.setData(str(tmp_list[3]), role=QtCore.Qt.DisplayRole)
    #         model_chain_item = QtGui.QStandardItem()
    #         model_chain_item.setData(str(tmp_list[4]), role=QtCore.Qt.DisplayRole)
    #         model_pos_item = QtGui.QStandardItem()
    #         model_pos_item.setData(int(tmp_list[5]), role=QtCore.Qt.DisplayRole)
    #         model_resi_item = QtGui.QStandardItem()
    #         model_resi_item.setData(str(tmp_list[6]), role=QtCore.Qt.DisplayRole)
    #         distance_item = QtGui.QStandardItem()
    #         tmp_distance_value = decimal.Decimal(float(tmp_list[7]))
    #         distance_item.setData(str(tmp_distance_value.quantize(decimal.Decimal('0.00'))), role=QtCore.Qt.DisplayRole)
    #         standard_item_list.append(pair_no_item)
    #         standard_item_list.append(ref_chain_item)
    #         standard_item_list.append(ref_pos_item)
    #         standard_item_list.append(ref_resi_item)
    #         standard_item_list.append(model_chain_item)
    #         standard_item_list.append(model_pos_item)
    #         standard_item_list.append(model_resi_item)
    #         standard_item_list.append(distance_item)
    #         self.csv_model.insertRow(i, standard_item_list)
    #     i += 1
    # csv_file.close()

    print(self._protein_pair.distance_analysis.analysis_results.distance_data)

    tmp_distance_data: dict[str, np.ndarray] = (
        self._protein_pair.distance_analysis.analysis_results.distance_data
    )
    for i in range(len(tmp_distance_data[pyssa_keys.ARRAY_DISTANCE_INDEX])):
      standard_item_list = []
      pair_no_item = QtGui.QStandardItem()
      pair_no_item.setData(
          int(tmp_distance_data[pyssa_keys.ARRAY_DISTANCE_INDEX][i]),
          role=QtCore.Qt.DisplayRole,
      )
      ref_chain_item = QtGui.QStandardItem()
      ref_chain_item.setData(
          str(tmp_distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_CHAIN][i]),
          role=QtCore.Qt.DisplayRole,
      )
      ref_pos_item = QtGui.QStandardItem()
      ref_pos_item.setData(
          int(tmp_distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_POSITION][i]),
          role=QtCore.Qt.DisplayRole,
      )
      ref_resi_item = QtGui.QStandardItem()
      ref_resi_item.setData(
          str(tmp_distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_1_RESI][i]),
          role=QtCore.Qt.DisplayRole,
      )
      model_chain_item = QtGui.QStandardItem()
      model_chain_item.setData(
          str(tmp_distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_CHAIN][i]),
          role=QtCore.Qt.DisplayRole,
      )
      model_pos_item = QtGui.QStandardItem()
      model_pos_item.setData(
          int(tmp_distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_POSITION][i]),
          role=QtCore.Qt.DisplayRole,
      )
      model_resi_item = QtGui.QStandardItem()
      model_resi_item.setData(
          str(tmp_distance_data[pyssa_keys.ARRAY_DISTANCE_PROT_2_RESI][i]),
          role=QtCore.Qt.DisplayRole,
      )
      distance_item = QtGui.QStandardItem()
      tmp_distance_value = decimal.Decimal(
          float(tmp_distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES][i])
      )
      distance_item.setData(
          str(tmp_distance_value.quantize(decimal.Decimal("0.00"))),
          role=QtCore.Qt.DisplayRole,
      )
      standard_item_list.append(pair_no_item)
      standard_item_list.append(ref_chain_item)
      standard_item_list.append(ref_pos_item)
      standard_item_list.append(ref_resi_item)
      standard_item_list.append(model_chain_item)
      standard_item_list.append(model_pos_item)
      standard_item_list.append(model_resi_item)
      standard_item_list.append(distance_item)
      self.csv_model.insertRow(i, standard_item_list)

    self.csv_model.removeRow(0)
    self.table_view.setAlternatingRowColors(True)
    self.table_view.resizeColumnsToContents()
    self.table_view.verticalHeader().setVisible(False)
    self.table_view.setSortingEnabled(True)
    self.table_view.sortByColumn(0, QtCore.Qt.AscendingOrder)
    self.table_view.setEditTriggers(
        QtWidgets.QAbstractItemView.NoEditTriggers
    )  # disables editing of cells
    self.table_view.setLocale(QtCore.QLocale("English"))
    # </editor-fold>

    # --END

    self.create_all_graphics()
    self.delayed_resize()
    self.setup_context_menu()
    self._connect_all_signals()

  def _initialize_ui(self) -> None:
    """Initializes the user interface by defining and setting up all the required UI elements such as scroll areas, plot widgets, menu bar, labels, and actions."""
    # <editor-fold desc="Define basic ui elements">
    self.scroll_area = QScrollArea()
    self.plot_widget_dplot = PlotWidget()
    self.plot_widget_dhistogram = PlotWidget()
    self.menubar = QtWidgets.QMenuBar(self)
    self.lbl_status = QtWidgets.QLabel()
    self.lbl_status.setText("Status: ")

    # </editor-fold>

    # <editor-fold desc="Create a menu and add it to the menu bar">
    # plots_menu = QtWidgets.QMenu('Views', self)  # TODO: raises an memory access violation error
    table_menu = QtWidgets.QMenu("Table", self)
    hide_column_menu = QtWidgets.QMenu("Hide Column", self)
    show_column_menu = QtWidgets.QMenu("Show Column", self)
    help_menu = QtWidgets.QMenu("Help", self)
    # self.menubar.addMenu(plots_menu)
    self.menubar.addMenu(table_menu)
    table_menu.addMenu(hide_column_menu)
    table_menu.addMenu(show_column_menu)
    self.menubar.addMenu(help_menu)

    # </editor-fold>

    # <editor-fold desc="Create actions and add them to the menu">
    self.action_plot = QtWidgets.QAction("Plot", self)
    self.action_plot.setCheckable(True)
    self.action_plot.setChecked(True)
    self.action_histogram = QtWidgets.QAction("Histogram", self)
    self.action_histogram.setCheckable(True)
    self.action_histogram.setChecked(True)

    self.action_table = QtWidgets.QAction("Table", self)
    self.action_table.setCheckable(True)
    self.action_table.setChecked(True)

    self.action_sync_with_pymol = QtWidgets.QAction("Sync With PyMOL", self)
    self.action_sync_with_pymol.setCheckable(True)
    self.action_sync_with_pymol.setChecked(False)

    # plots_menu.addAction(self.action_plot)
    # plots_menu.addAction(self.action_histogram)
    # plots_menu.addAction(self.action_table)
    # plots_menu.addAction(self.action_sync_with_pymol)

    self.action_hide_residue_pair_no = QtWidgets.QAction("Residue Pair No.")
    self.action_hide_protein_1_chain = QtWidgets.QAction("Protein 1 Chain")
    self.action_hide_protein_1_position = QtWidgets.QAction(
        "Protein 1 Position"
    )
    self.action_hide_protein_1_residue = QtWidgets.QAction("Protein 1 Residue")
    self.action_hide_protein_2_chain = QtWidgets.QAction("Protein 2 Chain")
    self.action_hide_protein_2_position = QtWidgets.QAction(
        "Protein 2 Position"
    )
    self.action_hide_protein_2_residue = QtWidgets.QAction("Protein 2 Residue")
    self.action_hide_distance = QtWidgets.QAction("Distance")

    hide_column_menu.addAction(self.action_hide_residue_pair_no)
    hide_column_menu.addAction(self.action_hide_protein_1_chain)
    hide_column_menu.addAction(self.action_hide_protein_1_position)
    hide_column_menu.addAction(self.action_hide_protein_1_residue)
    hide_column_menu.addAction(self.action_hide_protein_2_chain)
    hide_column_menu.addAction(self.action_hide_protein_2_position)
    hide_column_menu.addAction(self.action_hide_protein_2_residue)
    hide_column_menu.addAction(self.action_hide_distance)

    self.action_show_all_columns = QtWidgets.QAction("All Columns")
    self.action_show_residue_pair_no = QtWidgets.QAction("Residue Pair No.")
    self.action_show_protein_1_chain = QtWidgets.QAction("Protein 1 Chain")
    self.action_show_protein_1_position = QtWidgets.QAction(
        "Protein 1 Position"
    )
    self.action_show_protein_1_residue = QtWidgets.QAction("Protein 1 Residue")
    self.action_show_protein_2_chain = QtWidgets.QAction("Protein 2 Chain")
    self.action_show_protein_2_position = QtWidgets.QAction(
        "Protein 2 Position"
    )
    self.action_show_protein_2_residue = QtWidgets.QAction("Protein 2 Residue")
    self.action_show_distance = QtWidgets.QAction("Distance")

    show_column_menu.addAction(self.action_show_all_columns)
    show_column_menu.addSeparator()
    show_column_menu.addAction(self.action_show_residue_pair_no)
    show_column_menu.addAction(self.action_show_protein_1_chain)
    show_column_menu.addAction(self.action_show_protein_1_position)
    show_column_menu.addAction(self.action_show_protein_1_residue)
    show_column_menu.addAction(self.action_show_protein_2_chain)
    show_column_menu.addAction(self.action_show_protein_2_position)
    show_column_menu.addAction(self.action_show_protein_2_residue)
    show_column_menu.addAction(self.action_show_distance)

    self.action_docs = QtWidgets.QAction("PySSA Documentation", self)
    self.action_docs.triggered.connect(self._open_help_center)
    self.action_help = QtWidgets.QAction("Get Help", self)
    self.action_help.triggered.connect(self._open_distance_data_visualizer_help)

    help_menu.addAction(self.action_docs)
    help_menu.addAction(self.action_help)

    # </editor-fold>

    # self.toolbar = NavigationToolbar2QT(self.plot_widget.canvas, self)
    # Find the action you want to remove
    # items_to_remove = ["Home", "Back", "Forward", "Pan", "Zoom", "Subplots"]
    # for tmp_item in items_to_remove:
    #     for action in self.toolbar.actions():
    #         print(action.text())
    #         if action.text() == tmp_item:
    #             self.toolbar.removeAction(action)

    # <editor-fold desc="Set layouts">
    self.scroll_area.setWidget(self.plot_widget_dhistogram)

    # Create labels
    self.lbl_status1 = QtWidgets.QLabel(
        f"Protein 1: {self.protein_pair_for_analysis.protein_1.get_molecule_object()}"
    )
    self.lbl_status2 = QtWidgets.QLabel(
        f"Protein 2: {self.protein_pair_for_analysis.protein_2.get_molecule_object()}"
    )

    # Create a QTableView
    self.table_view = QtWidgets.QTableView()
    self.table_view.setMinimumWidth(450)

    # Create a QWidget to hold labels and QTableView
    self.container_widget = QtWidgets.QWidget()
    container_layout = QtWidgets.QVBoxLayout(self.container_widget)
    container_layout.addWidget(self.lbl_status1)
    container_layout.addWidget(self.lbl_status2)
    container_layout.addWidget(self.table_view)

    stylesheet = """
                    QLabel {
                        background-color: white;
                        font-size: 12px;
                        padding: 5px;
                        border-style: solid;
                        border-width: 2px;
                        border-radius: 6px;
                        border-color: #DCDBE3;
                    }
                    """
    self.container_widget.setStyleSheet(stylesheet)
    # # Create a QHBoxLayout for the scroll area
    # self.scroll_area_layout = QtWidgets.QHBoxLayout()
    #
    # # Assuming self.plot_widget and self.scroll_area are already defined
    # self.scroll_area_layout.addWidget(self.plot_widget)
    # self.scroll_area_layout.addWidget(self.container_widget)
    #
    # # Set up scroll area
    # self.scroll_area = QtWidgets.QScrollArea()
    # self.scroll_area.setLayout(self.scroll_area_layout)
    #
    # # Create main layout
    # self.main_Layout = QtWidgets.QVBoxLayout()
    # self.main_Layout.addWidget(self.scroll_area)
    # self.main_Layout.addWidget(self.lbl_status)
    # self.main_Layout.setMenuBar(self.menubar)
    # self.setLayout(self.main_Layout)

    # Create the first splitter to divide the dialog into two sections
    self.vertical_splitter = QtWidgets.QSplitter()
    self.main_layout = QtWidgets.QVBoxLayout(self)
    self.main_layout.addWidget(self.vertical_splitter)

    # Left side (plot area)
    plot_area = QtWidgets.QFrame()
    self.vertical_splitter.addWidget(plot_area)
    # Right side (table)
    self.vertical_splitter.addWidget(self.container_widget)
    # Second splitter within the plot area to split it horizontally
    self.horizontal_splitter = QtWidgets.QSplitter()
    plot_area_layout = QVBoxLayout()
    plot_area.setLayout(plot_area_layout)
    plot_area_layout.addWidget(self.horizontal_splitter)
    # Left part of the plot area
    self.horizontal_splitter.addWidget(self.plot_widget_dplot)
    # Right part of the plot area
    self.horizontal_splitter.addWidget(self.scroll_area)
    self.horizontal_splitter.setOrientation(0)  # Set orientation to horizontal

    self.vertical_splitter.setCollapsible(0, False)
    self.main_layout.setMenuBar(self.menubar)
    self.main_layout.addWidget(self.lbl_status)
    self.setLayout(self.main_layout)

    # </editor-fold>

    # <editor-fold desc="Setup subplots">
    self._ax_plot = self.plot_widget_dplot.figure.add_subplot()
    self._ax_hist = self.plot_widget_dhistogram.figure.add_subplot()

    # </editor-fold>

    # <editor-fold desc="Set window styles">
    styles.set_stylesheet(self)
    stylesheet = """
                QDialog {background-color: #F6F4F8;}
                QTableWidget {background-color: white;}
                """
    self.setStyleSheet(stylesheet)
    self.resize(1400, 800)
    self.setWindowFlag(QtCore.Qt.WindowMaximizeButtonHint, True)
    self.setWindowFlag(QtCore.Qt.WindowCloseButtonHint, True)
    # self.setWindowFlags(self.windowFlags() ^ QtCore.Qt.WindowContextHelpButtonHint)
    self.setWindowIcon(QtGui.QIcon(constants.PLUGIN_LOGO_FILEPATH))
    self.setWindowTitle("Distance Data Visualizer")

    # </editor-fold>

    self.change_selection_color()

  # <editor-fold desc="Help related methods">
  def _open_help_center(self) -> None:
    """Opens the help dialog."""
    self._help_manager.open_general_help_page()

  def _open_distance_data_visualizer_help(self) -> None:
    """Opens the help dialog on the distance_data_visualizer page."""
    self._help_manager.open_distance_data_visualizer_page()

  # </editor-fold>

  def _connect_all_signals(self) -> None:
    """Connects all UI elements to their corresponding slot functions in the class."""
    self.action_table.triggered.connect(self.hide_distance_table)

    self.action_hide_residue_pair_no.triggered.connect(
        self.__action_hide_residue_pair_no
    )
    self.action_hide_protein_1_chain.triggered.connect(
        self.__action_hide_protein_1_chain
    )
    self.action_hide_protein_1_position.triggered.connect(
        self.__action_hide_protein_1_position
    )
    self.action_hide_protein_1_residue.triggered.connect(
        self.__action_hide_protein_1_residue
    )
    self.action_hide_protein_2_chain.triggered.connect(
        self.__action_hide_protein_2_chain
    )
    self.action_hide_protein_2_position.triggered.connect(
        self.__action_hide_protein_2_position
    )
    self.action_hide_protein_2_residue.triggered.connect(
        self.__action_hide_protein_2_residue
    )
    self.action_hide_distance.triggered.connect(self.__action_hide_distances)

    self.action_show_residue_pair_no.triggered.connect(
        self.__action_show_residue_pair_no
    )
    self.action_show_protein_1_chain.triggered.connect(
        self.__action_show_protein_1_chain
    )
    self.action_show_protein_1_position.triggered.connect(
        self.__action_show_protein_1_position
    )
    self.action_show_protein_1_residue.triggered.connect(
        self.__action_show_protein_1_residue
    )
    self.action_show_protein_2_chain.triggered.connect(
        self.__action_show_protein_2_chain
    )
    self.action_show_protein_2_position.triggered.connect(
        self.__action_show_protein_2_position
    )
    self.action_show_protein_2_residue.triggered.connect(
        self.__action_show_protein_2_residue
    )
    self.action_show_distance.triggered.connect(self.__action_show_distances)
    self.action_show_all_columns.triggered.connect(
        self.__action_show_all_columns
    )

    self.table_view.clicked.connect(self.highlight_table_selection_in_plot)
    self.action_plot.triggered.connect(self.move_horizontal_splitter)
    self.action_histogram.triggered.connect(self.move_horizontal_splitter)
    self.action_sync_with_pymol.triggered.connect(self.__slot_sync_with_pymol)

    self.resize_timer.timeout.connect(self.actual_resize)
    self.vertical_splitter.splitterMoved.connect(self.delayed_resize)
    self.horizontal_splitter.splitterMoved.connect(self.delayed_resize)

    self.plot_widget_dplot.canvas.mpl_connect(
        "button_press_event", self.on_canvas_click
    )
    # self.plot_widget_dhistogram.canvas.mpl_connect('button_press_event', self.__slot_open_context_menu_for_histogram)
    # self.plot_widget_dplot.canvas.mpl_connect('motion_notify_event', self._on_move)

  # <editor-fold desc="QSplitter related methods">
  def hide_distance_table(self) -> None:
    """Hide or show the distance table widget based on the state of the action_table checkbox.

    If the action_table checkbox is checked, the distance table widget is hidden by moving the splitter
    of the vertical splitter to a lower position. The position is calculated by subtracting 450 from the
    highest position of the vertical splitter.

    If the action_table checkbox is unchecked, the distance table widget is shown by moving the splitter
    of the vertical splitter to the highest position.
    """
    if self.action_table.isChecked():
      tmp_index = self.vertical_splitter.indexOf(self.container_widget)
      tmp_lowest, tmp_highest = self.vertical_splitter.getRange(tmp_index)
      self.vertical_splitter.moveSplitter(tmp_highest - 450, tmp_index)
    else:
      tmp_index = self.vertical_splitter.indexOf(self.container_widget)
      tmp_lowest, tmp_highest = self.vertical_splitter.getRange(tmp_index)
      self.vertical_splitter.moveSplitter(tmp_highest, tmp_index)

  def move_horizontal_splitter(self) -> None:
    """Moves the horizontal splitter based on the state of the plot and histogram actions.

    This method checks the state of the plot and histogram actions and moves the horizontal splitter accordingly.
    It determines the index and range of the horizontal splitter, and then moves the splitter based on the checked actions.
    """
    if self.action_plot.isChecked() and self.action_histogram.isChecked():
      # Both should be displayed
      tmp_index = self.horizontal_splitter.indexOf(self.plot_widget_dhistogram)
      tmp_lowest, tmp_highest = self.horizontal_splitter.getRange(tmp_index)
      self.horizontal_splitter.moveSplitter(tmp_lowest + 400, tmp_index)
    elif self.action_plot.isChecked() and not self.action_histogram.isChecked():
      # Only the plot should be displayed
      tmp_index = self.horizontal_splitter.indexOf(self.plot_widget_dhistogram)
      tmp_lowest, tmp_highest = self.horizontal_splitter.getRange(tmp_index)
      self.horizontal_splitter.moveSplitter(tmp_highest, tmp_index)
    elif not self.action_plot.isChecked() and self.action_histogram.isChecked():
      # Only the histogram should be displayed
      tmp_index = self.horizontal_splitter.indexOf(self.plot_widget_dhistogram)
      tmp_lowest, tmp_highest = self.horizontal_splitter.getRange(tmp_index)
      self.horizontal_splitter.moveSplitter(tmp_lowest, tmp_index)
    elif (
        not self.action_plot.isChecked()
        and not self.action_histogram.isChecked()
        and self.action_table.isChecked()
    ):
      # No plots should be displayed
      tmp_index = self.vertical_splitter.indexOf(self.container_widget)
      tmp_lowest, tmp_highest = self.vertical_splitter.getRange(tmp_index)
      self.vertical_splitter.moveSplitter(tmp_lowest, tmp_index)
    elif (
        not self.action_plot.isChecked()
        and not self.action_histogram.isChecked()
        and not self.action_table.isChecked()
    ):
      tmp_index = self.vertical_splitter.indexOf(self.container_widget)
      tmp_lowest, tmp_highest = self.vertical_splitter.getRange(tmp_index)
      self.vertical_splitter.moveSplitter(tmp_highest, tmp_index)

  def delayed_resize(self) -> None:
    """Starts a delay timer and clears the figures.

    This method starts or restarts a timer when the splitter is moved.
    When the timer is triggered, the plot widgets are cleared and the resize timer is started.
    """
    # Start or restart the timer when the splitter is moved
    self.plot_widget_dplot.figure.clear()
    self.plot_widget_dhistogram.figure.clear()
    self.resize_timer.start()

  def actual_resize(self) -> None:
    """Performs the actual resizing operation.

    This method will be called after the timer interval has elapsed.
    """
    logger.debug("Resizing starts now ...")
    if self.container_widget.size().width() == 0:
      self.action_table.setChecked(False)
    else:
      self.action_table.setChecked(True)

    if self.plot_widget_dplot.size().height() == 0:
      self.action_plot.setChecked(False)
    else:
      self.action_plot.setChecked(True)

    if self.plot_widget_dhistogram.size().height() == 0:
      self.action_histogram.setChecked(False)
    elif self.plot_widget_dhistogram.size().width() == 0:
      self.action_plot.setChecked(False)
      self.action_histogram.setChecked(False)
    else:
      self.action_histogram.setChecked(True)

    self.toggle_graphics_visibility()

  def toggle_graphics_visibility(self) -> None:
    """Toggle the visibility of the graphics components.

    This method is used to toggle the visibility of the graphics components in the application.
    The visibility is determined by the states of the 'Plot' and 'Histogram' actions.
    If both actions are checked, the distance plot and histogram will be displayed.
    If only the 'Plot' action is checked, only the distance plot will be displayed.
    If only the 'Histogram' action is checked, only the distance histogram will be displayed.
    If neither action is checked, both components will be hidden.
    """
    self.table_view.selectionModel().clearSelection()
    self.plot_widget_dplot.figure.clear()
    self.plot_widget_dhistogram.figure.clear()
    if self.action_plot.isChecked() and self.action_histogram.isChecked():
      logger.debug(self.scroll_area.size())
      logger.debug(self.scroll_area.width() / 100)

      self.plot_widget_dplot.show()
      self.plot_widget_dhistogram.show()

      self._ax_plot = self.plot_widget_dplot.figure.subplots()
      self._ax_hist = self.plot_widget_dhistogram.figure.subplots()
      self.create_distance_plot()
      self.setup_plot_defaults()
      self.create_distance_histogram()
      self.setup_histogram_defaults()

      self.plot_widget_dplot.figure.tight_layout()
      self.plot_widget_dplot.canvas.draw()
      try:  # TODO: this is not an ideal way, but I couldn't find anything better
        self.plot_widget_dhistogram.figure.tight_layout()
        self.plot_widget_dhistogram.canvas.draw()
      except np.linalg.LinAlgError:
        logger.debug("A layout cannot be applied to the histogram.")

    elif self.action_plot.isChecked() and not self.action_histogram.isChecked():
      self.plot_widget_dplot.show()

      self._ax_plot = self.plot_widget_dplot.figure.subplots()
      self.create_distance_plot()
      self.setup_plot_defaults()

      self.plot_widget_dplot.figure.tight_layout()
      self.plot_widget_dplot.canvas.draw()

    elif not self.action_plot.isChecked() and self.action_histogram.isChecked():
      self.plot_widget_dhistogram.show()
      self._ax_hist = self.plot_widget_dhistogram.figure.subplots()
      self.create_distance_histogram()
      self.setup_histogram_defaults()

      try:  # TODO: this is not an ideal way, but I couldn't find anything better
        self.plot_widget_dhistogram.figure.tight_layout()
        self.plot_widget_dhistogram.canvas.draw()
      except np.linalg.LinAlgError:
        logger.debug("A layout cannot be applied to the histogram.")
    else:
      tmp_index = self.vertical_splitter.indexOf(self.container_widget)
      tmp_lowest, tmp_highest = self.vertical_splitter.getRange(tmp_index)
      self.vertical_splitter.moveSplitter(tmp_lowest, tmp_index)
      self.plot_widget_dplot.hide()
      self.plot_widget_dhistogram.hide()

  # </editor-fold>

  def __slot_sync_with_pymol(self) -> None:
    """This method was implemented in an earlier version of PySSA but was removed.

    Use the code below as starting point if the feature should come back.
    """
    # if self.action_sync_with_pymol.isChecked():
    #     self._sync_with_pymol_flag = True
    # else:
    #     self._sync_with_pymol_flag = False
    raise NotImplementedError()

  # <editor-fold desc="Plotting related methods">
  def create_all_graphics(self) -> None:
    """Create all graphics for the plot widgets.

    This method shows the plot widgets and creates the distance plot and histogram.
    It also sets up the default settings for the plots.
    """
    self.plot_widget_dplot.show()
    self.plot_widget_dhistogram.show()
    # num_subplots = 2  # Adjust this number based on your requirements
    # self.plot_widget_dplot.figure.clear()  # Clear existing figure
    # self.plot_widget_dhistogram.figure.clear()

    # # Adjust figure layout based on the number of subplots
    # self.plot_widget.figure.subplots(num_subplots, 1)
    #
    # # Create and set up the subplots
    # self._ax_plot = self.plot_widget.figure.axes[0]  # First subplot
    # self.create_distance_plot()
    # self.setup_plot_defaults()
    #
    # self._ax_hist = self.plot_widget.figure.axes[1]  # Second subplot
    # self.create_distance_histogram()
    # self.setup_histogram_defaults()
    #
    # # Adjust layout
    # self.plot_widget.figure.tight_layout()
    # self.plot_widget.canvas.draw()

    self.create_distance_plot()
    self.setup_plot_defaults()
    self.create_distance_histogram()
    self.setup_histogram_defaults()
    self.plot_widget_dplot.figure.tight_layout()
    self.plot_widget_dhistogram.figure.tight_layout()
    self.plot_widget_dplot.canvas.draw()
    self.plot_widget_dhistogram.canvas.draw()

  def setup_plot_defaults(self) -> None:
    """Set up the default plot settings.

    This method sets up the default settings for the plot, including the axis labels, title, x-axis limits, and tick settings.
    """
    # Set axis labels
    self._ax_plot.set_xlabel("Residue Pair No.")
    self._ax_plot.set_ylabel("Distance In Å")
    self._ax_plot.set_title("Distance Plot")
    # Set the x-axis limits with minimum value set to 0
    self._ax_plot.set_xlim(0)
    # Set the number of ticks for the x and y axes
    self._ax_plot.xaxis.set_major_locator(ticker.MultipleLocator(10))

  def setup_histogram_defaults(self) -> None:
    """Sets up default settings for the histogram plot.

    This method is used to customize the appearance of the histogram plot.
    It moves the entire x-axis to the top, sets the title, axis labels, and y-axis inversion.
    It also removes spines where no axis are present and sets the major locator for the x-axis.
    """
    # Move the entire x-axis to the top
    self._ax_hist.xaxis.tick_top()
    self._ax_hist.xaxis.set_label_position("top")
    self._ax_hist.set_title("Distance Histogram")
    # Set axis labels
    self._ax_hist.set_xlabel("Count")
    self._ax_hist.set_ylabel("Bins")
    # Invert the y-axis
    self._ax_hist.invert_yaxis()
    # Remove the spines where no ax is are present
    self._ax_hist.spines["right"].set_visible(False)
    self._ax_hist.spines["bottom"].set_visible(False)
    # Remove the spines where no axis are present
    self._ax_hist.spines["right"].set_visible(False)
    self._ax_hist.spines["bottom"].set_visible(False)
    # Set axis labels
    self._ax_hist.set_xlabel("Frequency Of C-α Distances")
    self._ax_hist.set_ylabel("Distance Interval In Å")
    self._ax_hist.xaxis.set_major_locator(
        ticker.MultipleLocator(
            self._histogram_properties[
                enums.HistogramPropertiesEnum.X_AXIS_UNITS
            ]
        )
    )

  def create_distance_plot(self) -> None:
    """Creates a distance plot by plotting the distance values against the index of the values."""
    # data for actual distance plot line
    distance_data = (
        self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data
    )
    distance_list = distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES]
    # data for cutoff line
    # cutoff_line = []
    # for i in range(len(distance_list)):
    #     cutoff_line.append(self.protein_pair_for_analysis.distance_analysis.cutoff)
    self._ax_plot.plot(distance_list, color="#367AF6")
    self._ax_plot.grid(True, linestyle="--", color="lightgray", linewidth=0.7)

  def create_distance_histogram(self) -> None:
    """Create a distance histogram based on the distance data."""
    distance_data: dict[
        str,
        np.ndarray,
    ] = (
        self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data
    )
    distance_list = copy.deepcopy(
        distance_data[pyssa_keys.ARRAY_DISTANCE_DISTANCES]
    )
    distance_list.sort()
    length = len(distance_list)
    max_distance = distance_list[length - 1]
    self.bins = np.arange(
        0,
        max_distance + 1,
        self._histogram_properties[
            enums.HistogramPropertiesEnum.DISTANCE_INTERVAL
        ],
    )
    hist, bin_edges = np.histogram(distance_list, bins=self.bins)

    frequencies = hist.tolist()
    self.freqs_without_zeros = []
    self.bins_without_zeros_label = []
    self.bins_without_zeros = []
    for i, freq in enumerate(frequencies):
      if freq != 0:
        self.freqs_without_zeros.append(freq)
        tmp_str_bin = str(float(bin_edges[i]))
        tmp_str_bin_2 = str(
            float(bin_edges[i])
            + self._histogram_properties[
                enums.HistogramPropertiesEnum.DISTANCE_INTERVAL
            ]
        )
        self.bins_without_zeros.append(bin_edges[i])
        self.bins_without_zeros_label.append(f"[{tmp_str_bin},{tmp_str_bin_2}]")

    self.bars = self._ax_hist.barh(
        self.bins_without_zeros_label,
        self.freqs_without_zeros,
        color="#367AF6",
        height=0.6,
    )  # 0.6 default
    self._ax_hist.bar_label(self.bars, padding=4)
    # Extra canvas and figure settings
    tmp_histogram_width = self.scroll_area.width() - 20
    tmp_histogram_height = int(((5 / 6) * len(self.bars)) * 100)
    self.plot_widget_dhistogram.resize(
        tmp_histogram_width, tmp_histogram_height
    )
    self.plot_widget_dhistogram.set_figure_size(
        tmp_histogram_width / 100, tmp_histogram_height / 100
    )

  # </editor-fold>

  def on_canvas_click(self, event) -> None:  # noqa: ANN001
    """Modifies the plots based on the on_canvas_click event.

    Args:
        event: The event object that contains information about the click event.
    """
    if event.button is MouseButton.LEFT:
      x_clicked, y_clicked = event.xdata, event.ydata
      # Clear the previous clicked point
      if self.clicked_point_scatter is not None:
        try:
          self.clicked_point_scatter.remove()
        except Exception as e:
          constants.PYSSA_LOGGER.warning(f"Something went wrong: {e}")

      # Find the nearest point on the line
      distance_data = self.protein_pair_for_analysis.distance_analysis.analysis_results.distance_data[
          pyssa_keys.ARRAY_DISTANCE_DISTANCES
      ]
      x_line = range(len(distance_data))
      y_line = distance_data
      distances = np.sqrt((x_line - x_clicked) ** 2 + (y_line - y_clicked) ** 2)
      index_nearest = np.argmin(distances)
      x_nearest, y_nearest = x_line[index_nearest], y_line[index_nearest]
      self.clicked_point_scatter = self._ax_plot.scatter(
          x_nearest, y_nearest, color="red", marker="o", label="Clicked Point"
      )
      self.plot_widget_dplot.canvas.draw()
      self.active_row_information = self.search_point_by_id(x_nearest)
      if self.active_row_information is not None:
        tmp_prot_1_name = (
            self.protein_pair_for_analysis.protein_1.get_molecule_object()
        )
        tmp_prot_2_name = (
            self.protein_pair_for_analysis.protein_2.get_molecule_object()
        )
        msg: str = (
            f"Status: Residue pair no. = {x_nearest}, Distance (Å) = {y_nearest} Between {tmp_prot_1_name} "
            f"Chain {self.active_row_information[1]} Position {self.active_row_information[2]} Residue {self.active_row_information[3]} and "
            f"{tmp_prot_2_name} Chain {self.active_row_information[4]} Position {self.active_row_information[5]} Residue {self.active_row_information[6]}"
        )
        self.lbl_status.setText(msg)
        # if self._sync_with_pymol_flag:
        #     self.highlight_residue_in_pymol(x_nearest, tmp_prot_1_name)
      else:
        self.lbl_status.setText(
            f"Status: Residue pair no. = {x_nearest}, Distance (Å) = {y_nearest}"
        )
      # Change the selection color
      self.highlight_histogram_bar(y_nearest)
    elif event.button is MouseButton.RIGHT:
      print("Open Context Menu ...")

  def change_selection_color(self) -> None:
    """Changes the color palette of the table view."""
    palette = self.table_view.palette()
    palette.setColor(QtGui.QPalette.Highlight, QtGui.QColor("red"))
    palette.setColor(QtGui.QPalette.HighlightedText, QtGui.QColor("white"))
    self.table_view.setPalette(palette)

  def search_point_by_id(self, target_id: int) -> Optional[tuple]:
    """Searches a point by its ID.

    Args:
        target_id (int): The id of the point to search for.

    Returns:
        A tuple containing the information of the point if found, otherwise None.

    Raises:
        exception.IllegalArgumentError: If `target_id` is None.
    """
    # <editor-fold desc="Checks">
    if target_id is None:
      logger.error("target_id is None.")
      raise exception.IllegalArgumentError("target_id is None.")

    # </editor-fold>

    for row in range(self.csv_model.rowCount()):
      id_item = self.csv_model.item(
          row, 0
      )  # Assuming running id is in the first column
      if id_item and id_item.text() == str(target_id):
        # Select the entire row in the QTableView
        self.table_view.selectRow(row)
        self.active_row_information = (
            self.table_view.model().index(row, 0).data(Qt.DisplayRole),
            self.table_view.model().index(row, 1).data(Qt.DisplayRole),
            self.table_view.model().index(row, 2).data(Qt.DisplayRole),
            self.table_view.model().index(row, 3).data(Qt.DisplayRole),
            self.table_view.model().index(row, 4).data(Qt.DisplayRole),
            self.table_view.model().index(row, 5).data(Qt.DisplayRole),
            self.table_view.model().index(row, 6).data(Qt.DisplayRole),
        )
        return self.active_row_information  # Return True if the id is found
    return None

  def highlight_table_selection_in_plot(self) -> None:
    """Highlights the selected row in the table by plotting a scatter point in the corresponding position on the plot.

    Also updates the status label with information about the selected row.
    """
    tmp_x_value = (
        self.table_view.model()
        .index(
            self.table_view.currentIndex().row(),
            0,
        )
        .data(Qt.DisplayRole)
    )
    tmp_y_value = float(
        self.table_view.model()
        .index(
            self.table_view.currentIndex().row(),
            7,
        )
        .data(Qt.DisplayRole)
    )

    # Clear the previous clicked point
    if self.clicked_point_scatter is not None:
      try:
        self.clicked_point_scatter.remove()
      except Exception as e:
        constants.PYSSA_LOGGER.warning(f"Something went wrong: {e}")
    self.clicked_point_scatter = self._ax_plot.scatter(
        tmp_x_value, tmp_y_value, color="red", marker="o", label="Clicked Point"
    )
    self.plot_widget_dplot.canvas.draw()

    tmp_row = self.table_view.currentIndex().row()
    self.table_view.selectRow(tmp_row)
    self.active_row_information = (
        self.table_view.model().index(tmp_row, 0).data(Qt.DisplayRole),
        self.table_view.model().index(tmp_row, 1).data(Qt.DisplayRole),
        self.table_view.model().index(tmp_row, 2).data(Qt.DisplayRole),
        self.table_view.model().index(tmp_row, 3).data(Qt.DisplayRole),
        self.table_view.model().index(tmp_row, 4).data(Qt.DisplayRole),
        self.table_view.model().index(tmp_row, 5).data(Qt.DisplayRole),
        self.table_view.model().index(tmp_row, 6).data(Qt.DisplayRole),
    )
    if self.active_row_information is not None:
      tmp_prot_1_name = (
          self.protein_pair_for_analysis.protein_1.get_molecule_object()
      )
      tmp_prot_2_name = (
          self.protein_pair_for_analysis.protein_2.get_molecule_object()
      )
      msg: str = (
          f"Status: Residue pair no. = {tmp_x_value}, Distance (Å) = {tmp_y_value} Between {tmp_prot_1_name} "
          f"Chain {self.active_row_information[1]} Position {self.active_row_information[2]} Residue {self.active_row_information[3]} and "
          f"{tmp_prot_2_name} Chain {self.active_row_information[4]} Position {self.active_row_information[5]} Residue {self.active_row_information[6]}"
      )
      self.lbl_status.setText(msg)
    else:
      self.lbl_status.setText(
          f"Status: Residue pair no. = {tmp_x_value}, Distance (Å) = {tmp_y_value}"
      )
    self.highlight_histogram_bar(float(tmp_y_value))

  def highlight_histogram_bar(self, point_to_highlight: float) -> None:
    """Highlights the histogram bar that lies in the interval of the point in the plot.

    Args:
        point_to_highlight (float): The data point to be highlighted in the histogram.

    Raises:
        exception.IllegalArgumentError: If `point_to_highlight` is None.
    """
    # <editor-fold desc="Checks">
    if point_to_highlight is None:
      logger.error("point_to_highlight is None.")
      raise exception.IllegalArgumentError("point_to_highlight is None.")

    # </editor-fold>

    # Highlight a specific data point (example: point_to_highlight)
    if self.highlighted_bin_index is not None:
      self.bars[self.highlighted_bin_index].set_facecolor("#367AF6")

    for tmp_bin in self.bins:
      tmp_lower_bin = tmp_bin - 1
      tmp_upper_bin = tmp_bin
      if (
          tmp_lower_bin <= point_to_highlight <= tmp_upper_bin
      ):
        self.highlighted_bin_index = int(tmp_lower_bin)
        break
    self.highlighted_bin_index = self.bins_without_zeros.index(
        self.highlighted_bin_index
    )  # TODO: Could be cleaner
    if 0 <= self.highlighted_bin_index < len(self.bars):
      bar = self.bars[self.highlighted_bin_index]
      bar.set_facecolor("red")

    # Update the plot
    self.plot_widget_dhistogram.canvas.draw()

  def highlight_residue_in_pymol(
      self, the_current_id: int, tmp_prot_1_name: str
  ) -> None:
    """This method was implemented in an earlier version of PySSA but was removed.

    Use the code below as starting point if the feature should come back.
    """
    # self.table_view.currentIndex()
    # tmp_id, tmp_chain_1, tmp_pos_1, tmp_residue_1, _, _, _ = self.search_point_by_id(the_current_id)
    # tmp_pymol_selection = selection.Selection(tmp_prot_1_name)
    # tmp_pymol_selection.set_single_selection("", tmp_chain_1, tmp_residue_1, "")
    # # TODO: this needs to be changed to new architecture
    # # cmd.zoom(tmp_pymol_selection.selection_string)
    # # cmd.show("sticks", tmp_pymol_selection.selection_string)
    raise NotImplementedError()

  # <editor-fold desc="Distance table related methods for show and hide colums">
  def __action_hide_residue_pair_no(self) -> None:
    """Hides the residue pair no. column in the table view."""
    self.table_view.hideColumn(0)

  def __action_hide_protein_1_chain(self) -> None:
    """Hides the protein 1 chain column in the table view."""
    self.table_view.hideColumn(1)

  def __action_hide_protein_1_position(self) -> None:
    """Hides the protein 1 position column in the table view."""
    self.table_view.hideColumn(2)

  def __action_hide_protein_1_residue(self) -> None:
    """Hides the protein 1 residue column in the table view."""
    self.table_view.hideColumn(3)

  def __action_hide_protein_2_chain(self) -> None:
    """Hides the protein 2 chain column in the table view."""
    self.table_view.hideColumn(4)

  def __action_hide_protein_2_position(self) -> None:
    """Hides the protein 2 position column in the table view."""
    self.table_view.hideColumn(5)

  def __action_hide_protein_2_residue(self) -> None:
    """Hides the protein 2 residue column in the table view."""
    self.table_view.hideColumn(6)

  def __action_hide_distances(self) -> None:
    """Hides the distances column in the table view."""
    self.table_view.hideColumn(7)

  def __action_show_residue_pair_no(self) -> None:
    """Shows the residue pair no. column in the table view."""
    self.table_view.showColumn(0)

  def __action_show_protein_1_chain(self) -> None:
    """Shows the protein 1 chain column in the table view."""
    self.table_view.showColumn(1)

  def __action_show_protein_1_position(self) -> None:
    """Shows the protein 1 position column in the table view."""
    self.table_view.showColumn(2)

  def __action_show_protein_1_residue(self) -> None:
    """Shows the protein 1 residue column in the table view."""
    self.table_view.showColumn(3)

  def __action_show_protein_2_chain(self) -> None:
    """Shows the protein 2 chain column in the table view."""
    self.table_view.showColumn(4)

  def __action_show_protein_2_position(self) -> None:
    """Shows the protein 2 position column in the table view."""
    self.table_view.showColumn(5)

  def __action_show_protein_2_residue(self) -> None:
    """Shows the protein 2 residue column in the table view."""
    self.table_view.showColumn(6)

  def __action_show_distances(self) -> None:
    """Shows the distances column in the table view."""
    self.table_view.showColumn(7)

  def __action_show_all_columns(self) -> None:
    """Hides all columns in the table view."""
    self.table_view.showColumn(0)
    self.table_view.showColumn(1)
    self.table_view.showColumn(2)
    self.table_view.showColumn(3)
    self.table_view.showColumn(4)
    self.table_view.showColumn(5)
    self.table_view.showColumn(6)
    self.table_view.showColumn(7)

  # </editor-fold>

  def __slot_open_context_menu_for_histogram(self) -> None:
    raise NotImplementedError()

  def setup_context_menu(self) -> None:
    """Sets up the context menu for the distance data visualizer."""
    self.context_menu = QtWidgets.QMenu()
    self.hide_selected_column = QtWidgets.QAction(self.tr("Hide Column"))
    self.context_menu.addAction(self.hide_selected_column)
    self.hide_selected_column.triggered.connect(self._hide_selected_column)

    # Set the context menu for the buttons
    self.table_view.setContextMenuPolicy(3)
    self.table_view.customContextMenuRequested.connect(self._show_context_menu)

    # <editor-fold desc="Context menu setup for histogram">
    # for matplotlib histogram
    self.context_menu_hist = QtWidgets.QMenu()
    self.properties_hist = QtWidgets.QAction(self.tr("Properties"))
    self.context_menu_hist.addAction(self.properties_hist)
    self.properties_hist.triggered.connect(self._open_properties_view)

    self.properties_hist_restore = QtWidgets.QAction(
        self.tr("Restore Defaults")
    )
    self.context_menu_hist.addAction(self.properties_hist_restore)
    self.properties_hist_restore.triggered.connect(
        self._restore_default_histogram_properties
    )

    # Set the context menu for the buttons
    # fixme: This context menu makes problems, needs to be better implemented in the future
    # self.plot_widget_dhistogram.setContextMenuPolicy(3)
    # self.plot_widget_dhistogram.customContextMenuRequested.connect(
    #     self._show_context_menu_for_histogram
    # )
    # </editor-fold>

  def _show_context_menu(self, a_point: QtCore.QPoint) -> None:
    """Shows the context menu.

    Args:
        a_point (QtCore.QPoint): A QPoint object representing the position where the context menu should be shown.
    """
    self.hide_selected_column.triggered.disconnect()
    self.hide_selected_column.triggered.connect(self._hide_selected_column)
    self.context_menu.exec_(self.table_view.mapToGlobal(a_point))

  def _show_context_menu_for_histogram(self, a_point: QtCore.QPoint) -> None:
    """Shows the context menu for the histogram figure.

    Args:
        a_point (QtCore.QPoint): A QPoint object representing the position where the context menu should be shown.
    """
    self.properties_hist.triggered.disconnect()
    self.properties_hist.triggered.connect(self._open_properties_view)
    self.properties_hist_restore.triggered.disconnect()
    self.properties_hist_restore.triggered.connect(
        self._restore_default_histogram_properties
    )
    # add here more action connections
    self.context_menu_hist.exec_(
        self.plot_widget_dhistogram.mapToGlobal(a_point)
    )

  def _open_properties_view(self) -> None:
    """Opens the histogram properties dialog."""
    # self.properties_view =
    self._histogram_properties_view = (
        histogram_properties_view.HistogramPropertiesView(
            self._histogram_properties
        )
    )
    self._histogram_properties_view.new_properties.connect(
        self.post_open_properties_view
    )
    self._histogram_properties_view.show()

  def post_open_properties_view(self, new_properties: tuple) -> None:
    """Updates the histogram properties view after coming back from the properties' dialog.

    Args:
        new_properties (tuple): A tuple containing the new values for the x-axis units and distance interval.

    Raises:
        exception.IllegalArgumentError: If `new_properties` is None.
    """
    # <editor-fold desc="Checks">
    if new_properties is None:
      logger.error("new_properties is None.")
      raise exception.IllegalArgumentError("new_properties is None.")

    # </editor-fold>

    tmp_current_x_axis_units = self._histogram_properties[
        enums.HistogramPropertiesEnum.X_AXIS_UNITS
    ]
    tmp_current_distance_interval = self._histogram_properties[
        enums.HistogramPropertiesEnum.DISTANCE_INTERVAL
    ]
    tmp_x_axis_units = int(new_properties[0])
    tmp_distance_interval = float(new_properties[1])
    if (
        tmp_current_x_axis_units != tmp_x_axis_units
        or tmp_current_distance_interval != tmp_distance_interval
    ):
      # at least one property changed
      self._histogram_properties[enums.HistogramPropertiesEnum.X_AXIS_UNITS] = (
          tmp_x_axis_units
      )
      self._histogram_properties[
          enums.HistogramPropertiesEnum.DISTANCE_INTERVAL
      ] = tmp_distance_interval

      self.actual_resize()

  def _restore_default_histogram_properties(self) -> None:
    """Restores the histogram properties to the default values."""
    self._histogram_properties[enums.HistogramPropertiesEnum.X_AXIS_UNITS] = (
        constants.DEFAULT_HISTOGRAM_PROPERTIES[
            enums.HistogramPropertiesEnum.X_AXIS_UNITS
        ]
    )
    self._histogram_properties[
        enums.HistogramPropertiesEnum.DISTANCE_INTERVAL
    ] = constants.DEFAULT_HISTOGRAM_PROPERTIES[
        enums.HistogramPropertiesEnum.DISTANCE_INTERVAL
    ]
    self.actual_resize()

  def _hide_selected_column(self) -> None:
    """Hides a selected column from the table view."""
    self.table_view.hideColumn(self.table_view.currentIndex().column())
