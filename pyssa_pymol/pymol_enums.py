import enum


class CommandEnum(enum.Enum):
    REINITIALIZE_SESSION = "reinitialize_session"
    GET_ALL_OBJECT_NAMES = "get_all_object_names"
    GET_MODEL = "get_model"

    LOAD_PYMOL_SESSION = "load_pymol_session"
    SAVE_PYMOL_SESSION = "save_pymol_session"

    SELECT = "select"
    SCENE = "scene"
    LOAD_SCENE = "load_scene"
    GET_SCENE_LIST = "get_scene_list"

    SHOW_CUSTOM_REPRESENTATION = "show_custom_representation"
    HIDE_CUSTOM_REPRESENTATION = "hide_custom_representation"
    ZOOM_WITH_CUSTOM_PARAMETERS = "zoom_with_custom_parameters"
    SET_CUSTOM_SETTINGS = "set_custom_settings"
    SET_BACKGROUND_COLOR = "set_background_color"
    SET_DEFAULT_GRAPHIC_SETTINGS = "set_default_graphic_settings"

    GET_RESIDUE_COLORS = "get_residue_colors"
    GET_CHAIN_COLOR = "get_chain_color"
    GET_CHAIN_REPR_STATE = "get_chain_repr_state"

    COLOR_SELECTION = "color_selection"

    # Image commands
    RAY = "ray"
    DRAW = "draw"
    PNG = "png"
