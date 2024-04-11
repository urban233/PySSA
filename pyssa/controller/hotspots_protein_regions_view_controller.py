#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2022
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
"""Module for the Hotspots Dialog."""
import logging
import os
import subprocess
import typing

import pygetwindow
from PyQt5 import QtCore
from pymol import cmd

from pyssa.controller import interface_manager
from pyssa.gui.ui.custom_dialogs import custom_message_box
from pyssa.internal.data_structures import protein_pair
from pyssa.internal.thread import tasks
from pyssa.internal.thread.async_pyssa import util_async
from pyssa.util import session_util, constants
from pyssa.logging_pyssa import log_levels, log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class HotspotsProteinRegionsViewController(QtCore.QObject):
    """Class for the Hotspots Protein Regions View Controller."""
    return_value = QtCore.pyqtSignal(str)

    def __init__(self, the_interface_manager: "interface_manager.InterfaceManager") -> None:
        super().__init__()
        self._interface_manager = the_interface_manager
        self._view = the_interface_manager.get_hotspots_protein_regions_view()
        self._connect_all_ui_elements_to_slot_functions()
        self._protein_names: tuple[str, typing.Union[str, int]] = self._get_protein_names_from_tree_view()

    def open_help(self, a_page_name: str):
        """Opens the pyssa documentation window if it's not already open.

        Args:
            a_page_name (str): a name of a documentation page to display
        """
        self._interface_manager.status_bar_manager.show_temporary_message("Opening help center ...")
        if len(pygetwindow.getWindowsWithTitle(constants.WINDOW_TITLE_OF_HELP_CENTER)) != 1:
            self._interface_manager.documentation_window = None
        self._active_task = tasks.Task(
            target=util_async.open_documentation_on_certain_page,
            args=(a_page_name, self._interface_manager.documentation_window),
            post_func=self.__await_open_help,
        )
        self._active_task.start()

    def __await_open_help(self, return_value):
        self._interface_manager.documentation_window = return_value[2]
        if not os.path.exists(constants.HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH):
            tmp_dialog = custom_message_box.CustomMessageBoxOk(
                "The script for bringing the documentation window in front could not be found!", "Documentation",
                custom_message_box.CustomMessageBoxIcons.ERROR.value
            )
            tmp_dialog.exec_()
        else:
            self._interface_manager.documentation_window.restore()
            subprocess.run([constants.HELP_CENTER_BRING_TO_FRONT_EXE_FILEPATH])
            self._interface_manager.status_bar_manager.show_temporary_message("Opening help center finished.")

    def _open_help_for_dialog(self):
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Help' button was clicked.")
        self.open_help("help/hotspots/protein_regions/")

    def _get_protein_names_from_tree_view(self):
        """Checks if the object is a protein or a chain."""
        if self._interface_manager.current_tab_index == 1:
            # proteins tab
            if self._interface_manager.get_current_protein_tree_index_type() == "protein":
                return self._interface_manager.get_current_active_protein_object().get_molecule_object(), 0
            else:
                return self._interface_manager.get_current_active_protein_object().get_molecule_object(), 0
        elif self._interface_manager.current_tab_index == 2:
            # protein pairs tab
            tmp_type: str = self._interface_manager.get_current_protein_pair_tree_index_type()
            if tmp_type == "protein_pair" or tmp_type == "protein" or tmp_type == "chain" or tmp_type == "scene" or tmp_type == "header":
                tmp_protein_pair: "protein_pair.ProteinPair" = self._interface_manager.get_current_active_protein_pair_object()
                return tmp_protein_pair.protein_1.get_molecule_object(), tmp_protein_pair.protein_2.get_molecule_object()
            else:
                raise ValueError("Invalid tree view selection.")
        else:
            raise ValueError("Invalid tab index for this operation.")  # TODO: Add logger message

    def _connect_all_ui_elements_to_slot_functions(self) -> None:
        self._view.ui.btn_help.clicked.connect(self._open_help_for_dialog)
        self._view.ui.btn_sticks_show.clicked.connect(self.__slot_show_resi_sticks)
        self._view.ui.btn_sticks_hide.clicked.connect(self.__slot_hide_resi_sticks)
        self._view.ui.btn_disulfide_bonds_show.clicked.connect(self.__slot_show_disulfide_bonds)
        self._view.ui.btn_disulfide_bonds_hide.clicked.connect(self.__slot_hide_disulfide_bonds)
        self._view.ui.btn_position_zoom.clicked.connect(self.__slot_zoom_resi_position)

    def __slot_show_resi_sticks(self) -> None:
        """Shows the pymol selection as sticks."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Show' sticks button was clicked.")
        if session_util.check_if_sele_is_empty():
            return
        cmd.show(representation="sticks", selection="sele and not hydrogens")
        cmd.select(name="sele", selection="sele and not hydrogens")
        cmd.color(color="atomic", selection="sele and not elem C")
        cmd.set("valence", 0)  # this needs to be better implemented

    def __slot_hide_resi_sticks(self) -> None:
        """Hides the balls and sticks representation of the pymol selection."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Hide' sticks button was clicked.")
        if session_util.check_if_sele_is_empty():
            return
        cmd.hide(representation="sticks", selection="sele")

    def __slot_show_disulfide_bonds(self) -> None:
        """Shows all disulfid bonds within the pymol session."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Show' disulfide bonds button was clicked.")
        tmp_pymol_selection_option: str = "byres (resn CYS and name SG) within 2 of (resn CYS and name SG)"
        for tmp_protein_name in self._protein_names:
            if tmp_protein_name != 0:
                cmd.select(
                    name="disulfides",
                    selection=f"{tmp_protein_name} & {tmp_pymol_selection_option}",
                )
                cmd.color(color="atomic", selection="disulfides and not elem C")
                cmd.set("valence", 0)  # this needs to be better implemented
                cmd.show("sticks", "disulfides")
                cmd.hide("sticks", "elem H")

    def __slot_hide_disulfide_bonds(self) -> None:
        """Hides all disulfid bonds within the pymol session."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Hide' disulfide bonds button was clicked.")
        tmp_pymol_selection_option: str = "byres (resn CYS and name SG) within 2 of (resn CYS and name SG)"
        for tmp_protein_name in self._protein_names:
            if tmp_protein_name != 0:
                cmd.select(
                    name="disulfides",
                    selection=f"{tmp_protein_name} & {tmp_pymol_selection_option}",
                )
                cmd.hide("sticks", "disulfides")

    def __slot_zoom_resi_position(self) -> None:
        """Zooms to the pymol selection."""
        logger.log(log_levels.SLOT_FUNC_LOG_LEVEL_VALUE, "'Zoom' button was clicked.")
        session_util.check_if_sele_is_empty()
        cmd.zoom(selection="sele", buffer=8.0, state=0, complete=0)
