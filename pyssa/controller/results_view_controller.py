import os
from PyQt5 import QtWidgets
from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal

from pyssa.controller import interface_manager
from pyssa.gui.ui.dialogs import dialog_advanced_prediction_configurations
from pyssa.gui.ui.messageboxes import basic_boxes
from pyssa.gui.ui.styles import styles
from pyssa.gui.ui.views import results_view, plot_view
from pyssa.internal.data_structures import protein_pair
from pyssa.internal.data_structures.data_classes import prediction_protein_info, prediction_configuration
from pyssa.internal.thread import tasks
from pyssa.io_pyssa import safeguard
from pyssa.presenter import main_presenter_async
from pyssa.util import gui_utils, tools, constants, exit_codes, prediction_util


class ResultsViewController(QtCore.QObject):

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager", the_protein_pair: protein_pair.ProteinPair):
        super().__init__()
        self._interface_manager = the_interface_manager
        self._the_protein_pair = the_protein_pair
        self._view: "results_view.ResultsView" = the_interface_manager.get_results_view()
        self._connect_all_ui_elements_to_slot_functions()

    def _connect_all_ui_elements_to_slot_functions(self):
        self._view.ui.btn_view_plots.clicked.connect(self._open_plot_view)
        #self._view.ui.btn_save_data.clicked.connect()

    def _open_plot_view(self) -> None:
        tmp_dialog = plot_view.PlotView(self._the_protein_pair,
                                        self._interface_manager.get_current_project(),
                                        self._the_protein_pair)
        tmp_dialog.exec_()
