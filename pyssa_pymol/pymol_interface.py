import zmq

from pyssa.internal.data_structures.data_classes import pymol_command
from pyssa_pymol import pymol_enums


class PyMOLInterface:
    """Functions as interface to connect PySSA to PyMOL.

    Notes:
        This class is used within PySSA to communicate PyMOL commands.
    """
    def __init__(self):
        context = zmq.Context()
        self._main_socket = context.socket(zmq.REQ)
        self._main_socket.connect("tcp://127.0.0.1:9070")

    def reinitialize_session(self):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.REINITIALIZE_SESSION,
            (0, 0)
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def load_pymol_session(self, a_session_filepath) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.LOAD_PYMOL_SESSION,
            (0, str(a_session_filepath))
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def save_pymol_session(self, a_session_filepath) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.SAVE_PYMOL_SESSION,
            (0, str(a_session_filepath))
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def get_all_object_names(self) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.GET_ALL_OBJECT_NAMES,
            (0, 0)
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def get_model(self, a_selection_string: str) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.GET_MODEL,
            (0, a_selection_string)
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def select(self, a_name, a_selection_string):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.SELECT,
            (a_name, a_selection_string)
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def scene(self, a_key, an_action):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.SCENE,
            (a_key, an_action)
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def load_scene(self, a_scene_name) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.LOAD_SCENE,
            (0, str(a_scene_name))
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def get_scene_list(self) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.GET_SCENE_LIST,
            (0, 0)
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def set_custom_setting(self, a_setting_name, a_value) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.SET_CUSTOM_SETTINGS,
            (a_setting_name, a_value)
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def get_residue_colors(self, a_selection_string):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.GET_RESIDUE_COLORS,
            (0, a_selection_string)
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def get_chain_color(self, a_selection_string, a_chain_letter):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.GET_CHAIN_COLOR,
            (a_selection_string, a_chain_letter)
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def get_chain_repr_state(self, a_selection_string, a_chain_letter):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.GET_CHAIN_REPR_STATE,
            (a_selection_string, a_chain_letter)
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def show_custom_representation(self, a_representation, a_selection_string):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.SHOW_CUSTOM_REPRESENTATION,
            (a_representation, a_selection_string)
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def hide_custom_representation(self, a_representation, a_selection_string):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.HIDE_CUSTOM_REPRESENTATION,
            (a_representation, a_selection_string)
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def zoom_with_custom_parameters(self, a_selection_string, a_buffer_size=8.0, a_state=0, a_complete_flag=0):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.ZOOM_WITH_CUSTOM_PARAMETERS,
            (a_selection_string, a_buffer_size, a_state, a_complete_flag)
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def color_selection(self, a_pymol_color, a_selection_string):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.COLOR_SELECTION,
            (a_pymol_color, a_selection_string)
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def set_background_color(self, a_pymol_color):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.SET_BACKGROUND_COLOR,
            (0, a_pymol_color)
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def set_default_graphic_settings(self):
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.SET_DEFAULT_GRAPHIC_SETTINGS,
            (0, 0)
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def ray(self, a_width: int, a_height: int, a_renderer: int) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.RAY,
            (a_width, a_height, a_renderer)
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def draw(self, a_width: int, a_height: int, an_antialias_value: int) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.DRAW,
            (a_width, a_height, an_antialias_value)
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

    def png(self, an_image_filepath: str, a_dpi_value: int) -> dict:
        tmp_pymol_command = pymol_command.PyMOLCommand(
            pymol_enums.CommandEnum.PNG,
            (an_image_filepath, a_dpi_value)
        )
        return send_command_to_pymol(self._main_socket, tmp_pymol_command)

def send_command_to_pymol(
        the_main_socket,
        a_pymol_command: "pymol_command.PyMOLCommand"
) -> dict:
    """Sends a job request to the auxiliary pymol process."""
    # First check if main socket can receive messages (fixme: There is a better way to test this)
    the_main_socket.send_string("Check availability ...")
    the_main_socket.recv_string()
    the_main_socket.send_json(a_pymol_command.get_command())
    return the_main_socket.recv_json()
