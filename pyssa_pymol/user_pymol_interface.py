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
"""Module that provides an interface to the User PyMOL."""
from typing import Optional

import zmq
from PyQt5 import QtCore
from pyssa.internal.thread import tasks
from pyssa_pymol import pymol_enums, commands, local_logging

logger = local_logging.setup_logger(__file__)
__docformat__ = "google"


class UserPyMOLInterface(QtCore.QObject):
    """Functions as interface to connect PySSA to User PyMOL.

    Notes:
        This class is used to receive command requests form PySSA and send it to User PyMOL.
    """

    def __init__(self) -> None:
        """Constructor."""
        super().__init__()
        self.commands = {
            pymol_enums.CommandEnum.REINITIALIZE_SESSION: commands.reinitialize_session,
            pymol_enums.CommandEnum.LOAD_PYMOL_SESSION: commands.load_pymol_session,
            pymol_enums.CommandEnum.SAVE_PYMOL_SESSION: commands.save_pymol_session,
            pymol_enums.CommandEnum.GET_ALL_OBJECT_NAMES: commands.get_all_object_names,
            pymol_enums.CommandEnum.GET_MODEL: commands.get_model,
            pymol_enums.CommandEnum.LOAD_SCENE: commands.load_scene,
            pymol_enums.CommandEnum.GET_SCENE_LIST: commands.get_scene_list,
            pymol_enums.CommandEnum.SET_CUSTOM_SETTINGS: commands.set_custom_setting,
            pymol_enums.CommandEnum.GET_RESIDUE_COLORS: commands.get_residue_colors,
            pymol_enums.CommandEnum.GET_CHAIN_COLOR: commands.get_chain_color,
            pymol_enums.CommandEnum.GET_RESIDUE_COLOR_CONFIG: commands.get_residue_color_config,
            pymol_enums.CommandEnum.GET_CHAIN_REPR_STATE: commands.get_chain_repr_state,
            pymol_enums.CommandEnum.SHOW_CUSTOM_REPRESENTATION: commands.show_custom_representation,
            pymol_enums.CommandEnum.HIDE_CUSTOM_REPRESENTATION: commands.hide_custom_representation,
            pymol_enums.CommandEnum.ZOOM_WITH_CUSTOM_PARAMETERS: commands.zoom_with_custom_parameters,
            pymol_enums.CommandEnum.COLOR_SELECTION: commands.color_selection,
            pymol_enums.CommandEnum.SET_BACKGROUND_COLOR: commands.set_background_color,
            pymol_enums.CommandEnum.SCENE: commands.scene,
            pymol_enums.CommandEnum.SELECT: commands.select,
            pymol_enums.CommandEnum.RAY: commands.ray,
            pymol_enums.CommandEnum.DRAW: commands.draw,
            pymol_enums.CommandEnum.PNG: commands.png,
        }
        self._should_restart_service = False
        self._task = tasks.LegacyTask(
            target=self.main_loop,
            args=(0, 0),
            #post_func=self.finished_main_loop
        )
        self._task.start()

    def main_loop(self, placeholder_1: int, placeholder_2: int) -> tuple[int, int]:
        """Main loop method for handling socket communication and executing commands.

        Args:
            placeholder_1: Placeholder argument 1 for use with LegacyTask class.
            placeholder_2: Placeholder argument 2 for use with LegacyTask class.
        """
        # <editor-fold desc="Setup socket connection">
        try:
            logger.info("Establishing connection...")
            context = zmq.Context()
            main_socket = context.socket(zmq.REP)
            main_socket.bind("tcp://127.0.0.1:9070")
        except Exception as e:
            logger.error(f"Exception raised during startup! {e}")
            exit(1)
        else:
            logger.info("Binding socket was successful.")
        
        # </editor-fold>

        tmp_user_pymol_should_be_closed = False
        try:
            while tmp_user_pymol_should_be_closed is False:
                logger.info("Waiting for data from the client ...")
                data = main_socket.recv_json()
                print("Received data: ", data)
                logger.info(f"Received arguments: {data['arguments']}")

                if data["command"] == pymol_enums.CommandEnum.CLOSE_USER_PYMOL.value:
                    logger.info("Set boolean to close user pymol ...")
                    tmp_user_pymol_should_be_closed = True
                elif data["command"] == pymol_enums.CommandEnum.REINITIALIZE_SESSION.value:
                    logger.info("Running command: reinitialize_session")
                    tmp_result: tuple[bool, str] = commands.reinitialize_session()
                    logger.info(f"Command result: {tmp_result}")
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])},
                    )
                elif data["command"] == pymol_enums.CommandEnum.LOAD_PYMOL_SESSION.value:
                    logger.info("Running command: load_pymol_session")
                    _, tmp_filepath = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.load_pymol_session(tmp_filepath)
                    logger.info(f"Command result: {tmp_result}")
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])},
                    )
                elif data["command"] == pymol_enums.CommandEnum.SAVE_PYMOL_SESSION.value:
                    logger.info("Running command: save_pymol_session")
                    _, tmp_filepath = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.save_pymol_session(tmp_filepath)
                    logger.info(f"Command result: {tmp_result}")
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])},
                    )
                elif data["command"] == pymol_enums.CommandEnum.GET_ALL_OBJECT_NAMES.value:
                    logger.info("Running command: get_all_object_names")
                    tmp_result: tuple[bool, str, Optional[list]] = commands.get_all_object_names()
                    logger.info(f"Command result: {tmp_result}")
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1]), "data": tmp_result[2]},
                    )
                elif data["command"] == pymol_enums.CommandEnum.GET_MODEL.value:
                    logger.info("Running command: get_model")
                    _, tmp_selection_string = data["arguments"]
                    tmp_result: tuple = commands.get_model(tmp_selection_string)
                    logger.info(f"Command result: {tmp_result}")
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1]), "data": tmp_result[2]},
                    )
                elif data["command"] == pymol_enums.CommandEnum.LOAD_SCENE.value:
                    logger.info("Running command: load_scene")
                    _, tmp_scene_name = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.load_scene(tmp_scene_name)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])},
                    )
                elif data["command"] == pymol_enums.CommandEnum.GET_SCENE_LIST.value:
                    logger.info("Running command: get_scene_list")
                    tmp_result: tuple[bool, str, Optional[list]] = commands.get_scene_list()
                    logger.info(f"Command result: {tmp_result}")
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1]), "data": tmp_result[2]},
                    )
                elif data["command"] == pymol_enums.CommandEnum.SET_CUSTOM_SETTINGS.value:
                    logger.info("Running command: set_custom_settings")
                    tmp_setting_name, tmp_value = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.set_custom_setting(tmp_setting_name, tmp_value)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])},
                    )
                elif data["command"] == pymol_enums.CommandEnum.GET_RESIDUE_COLORS.value:
                    logger.info("Running command: get_residue_colors")
                    _, tmp_selection_string = data["arguments"]
                    tmp_result: tuple[bool, str, Optional[dict]] = commands.get_residue_colors(tmp_selection_string)
                    logger.info(f"Command result: {tmp_result}")
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1]), "data": tmp_result[2]},
                    )
                elif data["command"] == pymol_enums.CommandEnum.GET_CHAIN_COLOR.value:
                    logger.info("Running command: get_chain_color")
                    tmp_selection_string, tmp_chain_letter = data["arguments"]
                    tmp_result: tuple[bool, str, Optional[tuple]] = commands.get_chain_color(tmp_selection_string, tmp_chain_letter)
                    logger.info(f"Command result: {tmp_result}")
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1]), "data": tmp_result[2]},
                    )
                elif data["command"] == pymol_enums.CommandEnum.GET_RESIDUE_COLOR_CONFIG.value:
                    logger.info("Running command: get_residue_color_config")
                    tmp_selection_string, tmp_chain_letter = data["arguments"]
                    tmp_result: tuple[bool, str, Optional[list]] = commands.get_residue_color_config(tmp_selection_string, tmp_chain_letter)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1]), "data": tmp_result[2]},
                    )
                elif data["command"] == pymol_enums.CommandEnum.GET_CHAIN_REPR_STATE.value:
                    logger.info("Running command: get_chain_repr_state")
                    tmp_selection_string, tmp_chain_letter = data["arguments"]
                    tmp_result: tuple[bool, str, Optional[dict]] = commands.get_chain_repr_state(tmp_selection_string, tmp_chain_letter)
                    logger.info(f"Command result: {tmp_result}")
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1]), "data": tmp_result[2]},
                    )
                elif data["command"] == pymol_enums.CommandEnum.SHOW_CUSTOM_REPRESENTATION.value:
                    logger.info("Running command: show_custom_representation")
                    tmp_representation, tmp_selection_string = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.show_custom_representation(tmp_representation, tmp_selection_string)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])},
                    )
                elif data["command"] == pymol_enums.CommandEnum.HIDE_CUSTOM_REPRESENTATION.value:
                    logger.info("Running command: hide_custom_representation")
                    tmp_representation, tmp_selection_string = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.hide_custom_representation(tmp_representation, tmp_selection_string)
                    logger.info(f"Command result: {tmp_result}")
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])},
                    )
                elif data["command"] == pymol_enums.CommandEnum.ZOOM_WITH_CUSTOM_PARAMETERS.value:
                    logger.info("Running command: zoom_with_custom_parameters")
                    tmp_selection_string, tmp_buffer_size, tmp_state, tmp_complete_flag = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.zoom_with_custom_parameters(
                        tmp_selection_string, tmp_buffer_size, tmp_state, tmp_complete_flag,
                    )
                    logger.info(f"Command result: {tmp_result}")
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])},
                    )
                elif data["command"] == pymol_enums.CommandEnum.COLOR_SELECTION.value:
                    logger.info("Running command: color_selection")
                    tmp_pymol_color, tmp_selection_string = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.color_selection(tmp_pymol_color, tmp_selection_string)
                    logger.info(f"Command result: {tmp_result}")
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])},
                    )
                elif data["command"] == pymol_enums.CommandEnum.SET_BACKGROUND_COLOR.value:
                    logger.info("Running command: set_background_color")
                    _, tmp_pymol_color = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.set_background_color(tmp_pymol_color)
                    logger.info(f"Command result: {tmp_result}")
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])},
                    )
                elif data["command"] == pymol_enums.CommandEnum.SET_DEFAULT_GRAPHIC_SETTINGS.value:
                    logger.info("Running command: set_default_graphic_settings")
                    tmp_result: tuple[bool, str] = commands.set_default_graphics_settings()
                    logger.info(f"Command result: {tmp_result}")
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])},
                    )
                elif data["command"] == pymol_enums.CommandEnum.SCENE.value:
                    logger.info("Running command: scene")
                    a_key, an_action = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.scene(a_key, an_action)
                    logger.info(f"Command result: {tmp_result}")
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])},
                    )
                elif data["command"] == pymol_enums.CommandEnum.SELECT.value:
                    logger.info("Running command: select")
                    a_name, a_selection_string = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.select(a_name, a_selection_string)
                    logger.info(f"Command result: {tmp_result}")
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])},
                    )
                elif data["command"] == pymol_enums.CommandEnum.RAY.value:
                    logger.info("Running command: ray")
                    tmp_width, tmp_height, tmp_renderer = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.ray(tmp_width, tmp_height, tmp_renderer)
                    logger.info(f"Command result: {tmp_result}")
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])},
                    )
                elif data["command"] == pymol_enums.CommandEnum.DRAW.value:
                    logger.info("Running command: draw")
                    tmp_width, tmp_height, tmp_antialias_value = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.draw(tmp_width, tmp_height, tmp_antialias_value)
                    logger.info(f"Command result: {tmp_result}")
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])},
                    )
                elif data["command"] == pymol_enums.CommandEnum.PNG.value:
                    logger.info("Running command: png")
                    tmp_image_filepath, tmp_dpi_value = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.png(tmp_image_filepath, tmp_dpi_value)
                    logger.info(f"Command result: {tmp_result}")
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])},
                    )
        except Exception as e:
            main_socket.send_json({"success": False, "message": str(e)})
        finally:
            return 0, 0

    def finished_main_loop(self, a_placeholder_1: int, a_placeholder_2: int) -> None:
        """Logs the message "Finished main loop" using the logger.

        Args:
            a_placeholder_1: Placeholder argument 1 for use with LegacyTask class.
            a_placeholder_2: Placeholder argument 2 for use with LegacyTask class.
        """
        logger.info("Finished main loop")
