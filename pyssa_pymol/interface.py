import zmq
from PyQt5 import QtCore
from pymol import cmd

from pyssa.internal.thread import tasks
from pyssa_pymol import pymol_enums, commands


class Interface(QtCore.QObject):
    """Functions as interface to connect PySSA to PyMOL.

    Notes:
        This class is used to receive command requests form PySSA and send it to PyMOL.
    """
    def __init__(self):
        super().__init__()
        self.command_map = {
            pymol_enums.CommandEnum.LOAD_PYMOL_SESSION: commands.load_pymol_session
        }
        self._task = tasks.Task(
            target=self.main_loop,
            args=(0, 0),
            post_func=self.finished_main_loop
        )
        self._task.start()

    def main_loop(self, placeholder_1, placeholder_2):
        # <editor-fold desc="Setup socket connection">
        try:
            print("Establishing connection...")
            context = zmq.Context()
            main_socket = context.socket(zmq.REP)
            main_socket.bind("tcp://127.0.0.1:9070")  # Binding to localhost on port 9070
        except Exception as e:
            print("Exception raised during startup!")
            print(e)
            exit(1)
        else:
            print("Binding socket was successful.")
            pass
        # </editor-fold>

        print("Entering while loop ...")
        while True:
            # Wait for a request from the client
            try:
                print("Waiting to receive message ...")
                message = main_socket.recv_string()
                print("Received message: ", message)
                main_socket.send_string("I am ready to receive data.")
                print("Waiting for data from the client ...")
                data = main_socket.recv_json()
                print(data)

                if data["command"] == pymol_enums.CommandEnum.REINITIALIZE_SESSION.value:
                    tmp_result: tuple[bool, str] = commands.reinitialize_session()
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])}
                    )
                elif data["command"] == pymol_enums.CommandEnum.LOAD_PYMOL_SESSION.value:
                    _, tmp_filepath = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.load_pymol_session(tmp_filepath)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])}
                    )
                elif data["command"] == pymol_enums.CommandEnum.SAVE_PYMOL_SESSION.value:
                    _, tmp_filepath = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.save_pymol_session(tmp_filepath)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])}
                    )
                elif data["command"] == pymol_enums.CommandEnum.GET_ALL_OBJECT_NAMES.value:
                    tmp_result: tuple[bool, str] = commands.get_all_object_names()
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1]), "data": tmp_result[2]}
                    )
                elif data["command"] == pymol_enums.CommandEnum.GET_MODEL.value:
                    _, tmp_selection_string = data["arguments"]
                    tmp_result: tuple = commands.get_model(tmp_selection_string)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1]), "data": tmp_result[2]}
                    )
                elif data["command"] == pymol_enums.CommandEnum.LOAD_SCENE.value:
                    _, tmp_scene_name = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.load_scene(tmp_scene_name)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])}
                    )
                elif data["command"] == pymol_enums.CommandEnum.GET_SCENE_LIST.value:
                    tmp_result: tuple[bool, str] = commands.get_scene_list()
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1]), "data": tmp_result[2]}
                    )
                elif data["command"] == pymol_enums.CommandEnum.SET_CUSTOM_SETTINGS.value:
                    tmp_setting_name, tmp_value = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.set_custom_setting(tmp_setting_name, tmp_value)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])}
                    )
                elif data["command"] == pymol_enums.CommandEnum.GET_RESIDUE_COLORS.value:
                    _, tmp_selection_string = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.get_residue_colors(tmp_selection_string)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1]), "data": tmp_result[2]}
                    )
                elif data["command"] == pymol_enums.CommandEnum.GET_CHAIN_COLOR.value:
                    tmp_selection_string, tmp_chain_letter = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.get_chain_color(tmp_selection_string, tmp_chain_letter)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1]), "data": tmp_result[2]}
                    )
                elif data["command"] == pymol_enums.CommandEnum.GET_RESIDUE_COLOR_CONFIG.value:
                    tmp_selection_string, tmp_chain_letter = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.get_residue_color_config(tmp_selection_string, tmp_chain_letter)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1]), "data": tmp_result[2]}
                    )
                elif data["command"] == pymol_enums.CommandEnum.GET_CHAIN_REPR_STATE.value:
                    tmp_selection_string, tmp_chain_letter = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.get_chain_repr_state(tmp_selection_string, tmp_chain_letter)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1]), "data": tmp_result[2]}
                    )
                elif data["command"] == pymol_enums.CommandEnum.SHOW_CUSTOM_REPRESENTATION.value:
                    tmp_representation, tmp_selection_string = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.show_custom_representation(tmp_representation, tmp_selection_string)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])}
                    )
                elif data["command"] == pymol_enums.CommandEnum.HIDE_CUSTOM_REPRESENTATION.value:
                    tmp_representation, tmp_selection_string = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.hide_custom_representation(tmp_representation, tmp_selection_string)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])}
                    )
                elif data["command"] == pymol_enums.CommandEnum.ZOOM_WITH_CUSTOM_PARAMETERS.value:
                    tmp_selection_string, tmp_buffer_size, tmp_state, tmp_complete_flag = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.zoom_with_custom_parameters(
                        tmp_selection_string, tmp_buffer_size, tmp_state, tmp_complete_flag
                    )
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])}
                    )
                elif data["command"] == pymol_enums.CommandEnum.COLOR_SELECTION.value:
                    tmp_pymol_color, tmp_selection_string = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.color_selection(tmp_pymol_color, tmp_selection_string)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])}
                    )
                elif data["command"] == pymol_enums.CommandEnum.SET_BACKGROUND_COLOR.value:
                    _, tmp_pymol_color = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.set_background_color(tmp_pymol_color)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])}
                    )
                elif data["command"] == pymol_enums.CommandEnum.SET_DEFAULT_GRAPHIC_SETTINGS.value:
                    tmp_result: tuple[bool, str] = commands.set_default_graphics_settings()
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])}
                    )
                elif data["command"] == pymol_enums.CommandEnum.SCENE.value:
                    a_key, an_action = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.scene(a_key, an_action)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])}
                    )
                elif data["command"] == pymol_enums.CommandEnum.SELECT.value:
                    a_name, a_selection_string = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.select(a_name, a_selection_string)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])}
                    )
                elif data["command"] == pymol_enums.CommandEnum.RAY.value:
                    tmp_width, tmp_height, tmp_renderer = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.ray(tmp_width, tmp_height, tmp_renderer)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])}
                    )
                elif data["command"] == pymol_enums.CommandEnum.DRAW.value:
                    tmp_width, tmp_height, tmp_antialias_value = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.draw(tmp_width, tmp_height, tmp_antialias_value)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])}
                    )
                elif data["command"] == pymol_enums.CommandEnum.PNG.value:
                    tmp_image_filepath, tmp_dpi_value = data["arguments"]
                    tmp_result: tuple[bool, str] = commands.png(tmp_image_filepath, tmp_dpi_value)
                    main_socket.send_json(
                        {"success": tmp_result[0], "message": str(tmp_result[1])}
                    )
            except Exception as e:
                main_socket.send_json({"success": False, "message": str(e)})

    def finished_main_loop(self):
        print("Finished main loop")
