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
"""Module that functions as client for the User PyMOL interface."""
import logging

import zmq

from application_process import application_process_manager
from pyssa.internal.data_structures.data_classes import pymol_command
from pyssa.logging_pyssa import log_handlers
from pyssa.util import pyssa_exception
from pyssa_pymol import pymol_enums


logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


class PyMOLInterface:
    """Functions as interface to connect PySSA to PyMOL.

    Notes:
        This class is used within PySSA to communicate PyMOL commands.
    """
    def __init__(self, an_app_process_manager: "application_process_manager.ApplicationProcessManager"):
        self._app_process_manager: "application_process_manager.ApplicationProcessManager" = an_app_process_manager
        context = zmq.Context()
        self._main_socket = context.socket(zmq.REQ)
        self._main_socket.setsockopt(zmq.REQ_RELAXED, 1)
        self._main_socket.connect("tcp://127.0.0.1:9070")
        self._poller = zmq.Poller()
        self._poller.register(self._main_socket, zmq.POLLIN)

    def reinitialize_session(self):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.REINITIALIZE_SESSION,
            (0, 0)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def load_pymol_session(self, a_session_filepath) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.LOAD_PYMOL_SESSION,
            (0, str(a_session_filepath))
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def save_pymol_session(self, a_session_filepath) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.SAVE_PYMOL_SESSION,
            (0, str(a_session_filepath))
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def get_all_object_names(self) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.GET_ALL_OBJECT_NAMES,
            (0, 0)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def get_model(self, a_selection_string: str) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.GET_MODEL,
            (0, a_selection_string)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def select(self, a_name, a_selection_string):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.SELECT,
            (a_name, a_selection_string)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def scene(self, a_key, an_action):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.SCENE,
            (a_key, an_action)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def load_scene(self, a_scene_name) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.LOAD_SCENE,
            (0, str(a_scene_name))
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def get_scene_list(self) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.GET_SCENE_LIST,
            (0, 0)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def set_custom_setting(self, a_setting_name, a_value) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.SET_CUSTOM_SETTINGS,
            (a_setting_name, a_value)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def get_residue_colors(self, a_selection_string):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.GET_RESIDUE_COLORS,
            (0, a_selection_string)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def get_chain_color(self, a_selection_string, a_chain_letter):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.GET_CHAIN_COLOR,
            (a_selection_string, a_chain_letter)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def get_residue_color_config(self, a_selection_string, a_chain_letter):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.GET_RESIDUE_COLOR_CONFIG,
            (a_selection_string, a_chain_letter)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def get_chain_repr_state(self, a_selection_string, a_chain_letter):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.GET_CHAIN_REPR_STATE,
            (a_selection_string, a_chain_letter)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def show_custom_representation(self, a_representation, a_selection_string):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.SHOW_CUSTOM_REPRESENTATION,
            (a_representation, a_selection_string)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def hide_custom_representation(self, a_representation, a_selection_string):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.HIDE_CUSTOM_REPRESENTATION,
            (a_representation, a_selection_string)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def zoom_with_custom_parameters(self, a_selection_string, a_buffer_size=8.0, a_state=0, a_complete_flag=0):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.ZOOM_WITH_CUSTOM_PARAMETERS,
            (a_selection_string, a_buffer_size, a_state, a_complete_flag)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def color_selection(self, a_pymol_color, a_selection_string):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.COLOR_SELECTION,
            (a_pymol_color, a_selection_string)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def set_background_color(self, a_pymol_color):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.SET_BACKGROUND_COLOR,
            (0, a_pymol_color)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def set_default_graphic_settings(self):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.SET_DEFAULT_GRAPHIC_SETTINGS,
            (0, 0)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def ray(self, a_width: int, a_height: int, a_renderer: int) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.RAY,
            (a_width, a_height, a_renderer)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply


    def draw(self, a_width: int, a_height: int, an_antialias_value: int) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.DRAW,
            (a_width, a_height, an_antialias_value)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

    def png(self, an_image_filepath: str, a_dpi_value: int) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.PNG,
            (an_image_filepath, a_dpi_value)
        )
        try:
            tmp_reply = send_command_to_pymol(self._main_socket, tmp_pymol_command, self._poller, self._app_process_manager)
        except pyssa_exception.PyMOLNotRespondingError as e:
            logger.error(e)
            return {}
        else:
            return tmp_reply

def send_command_to_pymol(
        the_main_socket,
        a_pymol_command: "pymol_command.PyMOLCommand",
        the_poller,
        the_app_process_manager: "application_process_manager.ApplicationProcessManager",
        a_timeout=3000,
) -> dict:
    """Sends a job request to the auxiliary pymol process."""
    # First check if main socket can receive messages (fixme: There is a better way to test this)
    the_main_socket.send_string("Check availability ...")
    the_main_socket.recv_string()
    the_main_socket.send_json(a_pymol_command.get_command())
    tmp_reply_is_ready = False
    while True:
        tmp_events = the_poller.poll(a_timeout)
        if tmp_events:
            tmp_reply_is_ready = True
            break
        if the_app_process_manager.pymol_crashed():
            break
    if tmp_reply_is_ready:
        return the_main_socket.recv_json()
    else:
        raise pyssa_exception.PyMOLNotRespondingError()
