import pymol
from pymol import cmd
from pymol import CmdException

from pyssa.util import constants


def reinitialize_session():
    try:
        cmd.reinitialize()
    except CmdException as e:
        return False, e
    else:
        return True, ""


def load_pymol_session(a_filepath):
    """Loads a pymol session from an existing .pse file."""
    try:
        cmd.load(a_filepath)
    except CmdException as e:
        return False, e
    else:
        return True, ""


def save_pymol_session(a_filepath):
    """Loads a pymol session from an existing .pse file."""
    try:
        cmd.save(a_filepath)
    except CmdException as e:
        return False, e
    else:
        return True, ""


def get_all_object_names():
    try:
        tmp_names = cmd.get_names()
    except CmdException as e:
        return False, e, None
    else:
        return True, "", tmp_names


def get_model(a_selection_string):
    try:
        tmp_model = cmd.get_model(a_selection_string)
        tmp_resi = tmp_model.atom[0].resi
    except CmdException as e:
        return True, e, None
    except IndexError as e:
        return True, e, None
    else:
        return True, "", 0


def select(a_name, a_selection_string):
    try:
        print(cmd.select(name=a_name, selection=a_selection_string))
    except CmdException as e:
        return False, e
    else:
        return True, ""


def scene(a_key, an_action):
    try:
        cmd.scene(key=a_key, action=an_action)
    except CmdException as e:
        return False, e
    else:
        return True, ""


def load_scene(a_scene_name):
    try:
        cmd.scene(a_scene_name, "recall")
    except CmdException as e:
        return False, e
    else:
        return True, ""


def get_scene_list():
    try:
        tmp_scene_names = cmd.get_scene_list()
    except CmdException as e:
        return False, e, None
    else:
        return True, "", tmp_scene_names


def set_custom_setting(a_setting_name, a_value):
    try:
        cmd.set(a_setting_name, a_value)
    except CmdException as e:
        return False, e
    else:
        return True, ""


def get_residue_colors(a_selection_string):
    pymol.stored.colors = []
    cmd.iterate(a_selection_string, "stored.colors.append((chain, resi, name, color))")
    res_colors = {}
    for chain, resi, name, color_index in pymol.stored.colors:
        if name == 'CA':  # c-alpha atom
            res_colors[(chain, resi, name)] = constants.PYMOL_COLORS_WITH_INDICES[color_index]
    return True, "", res_colors


def get_chain_color(a_selection_string, chain_letter):
    pymol.stored.colors = []
    cmd.iterate(a_selection_string, "stored.colors.append((chain, resi, name, color))")
    tmp_ca_atom_colors = []
    tmp_residue_atom_colors = []
    tmp_is_colored_by_elements = False
    for chain, resi, name, color_index in pymol.stored.colors:
        if chain == chain_letter and name == 'CA':  # c-alpha atom
            try:
                tmp_ca_atom_colors.append(constants.PYMOL_COLORS_WITH_INDICES[color_index])
            except KeyError:
                tmp_ca_atom_colors.append("By RMSD")
        else:
            try:
                tmp_ca_atom_colors.append(constants.PYMOL_COLORS_WITH_INDICES[color_index])
            except KeyError:
                tmp_is_colored_by_elements = True
    tmp_chain_color = list(set(tmp_ca_atom_colors))[0]
    return True, "", (tmp_chain_color, tmp_ca_atom_colors, tmp_is_colored_by_elements)


def get_residue_color_config(a_selection_string, chain_letter):
    """Gets the colors of C-, N-, and O-atoms for the first residue of the given selection."""
    pymol.stored.colors = []
    cmd.iterate(a_selection_string, "stored.colors.append((chain, resi, name, color, elem))")
    tmp_residue_element_colors = ["", "", ""]
    for chain, resi, name, color_index, element in pymol.stored.colors:
        if chain == chain_letter and element == "C":  # C atom
            tmp_residue_element_colors[0] = constants.PYMOL_COLORS_WITH_INDICES[color_index]
        elif chain == chain_letter and element == "N":  # N atom
            tmp_residue_element_colors[1] = constants.PYMOL_COLORS_WITH_INDICES[color_index]
        elif chain == chain_letter and element == "O":  # O atom
            tmp_residue_element_colors[2] = constants.PYMOL_COLORS_WITH_INDICES[color_index]
    return True, "", tmp_residue_element_colors


def get_chain_repr_state(a_selection_string, chain_letter):
    """Returns the index for the current representation state.

    Note:
        If no match occurs, then -1 is returned
    """
    pymol.stored.reps = []
    cmd.iterate(a_selection_string, "stored.reps.append((chain, resi, name, reps))")
    for chain, resi, name, reps_index in pymol.stored.reps:
        if chain == chain_letter:  # c-alpha atom
            constants.PYSSA_LOGGER.debug("Representation index: " + str(reps_index))
            return True, "", constants.PYMOL_REPR_STATES_WITH_INDICES[reps_index]
    return True, "", {-1: ""}


def show_custom_representation(a_representation, a_selection_string):
    try:
        cmd.show(representation=a_representation, selection=a_selection_string)
    except CmdException as e:
        return False, e
    else:
        return True, ""


def hide_custom_representation(a_representation, a_selection_string):
    try:
        cmd.hide(representation=a_representation, selection=a_selection_string)
    except CmdException as e:
        return False, e
    else:
        return True, ""


def zoom_with_custom_parameters(a_selection_string, a_buffer_size=8.0, a_state=0, a_complete_flag=0):
    try:
        cmd.zoom(selection=a_selection_string, buffer=a_buffer_size, state=a_state, complete=a_complete_flag)
    except CmdException as e:
        return False, e
    else:
        return True, ""


def color_selection(a_pymol_color, a_selection_string):
    try:
        cmd.color(a_pymol_color, a_selection_string)
    except CmdException as e:
        return False, e
    else:
        return True, ""


def set_background_color(a_pymol_color):
    try:
        cmd.bg_color(a_pymol_color)
    except CmdException as e:
        return False, e
    else:
        return True, ""


def set_default_graphics_settings():
    try:
        cmd.set("valence", 0)
        cmd.set("scene_buttons", 0)
        cmd.set("ray_trace_mode", constants.PYMOL_DEFAULT_RAY_TRACE_MODE)
        cmd.set("antialias", constants.PYMOL_DEFAULT_ANTIALIAS)
        cmd.set("ambient", constants.PYMOL_DEFAULT_AMBIENT)
        cmd.set("cartoon_fancy_helices", constants.PYMOL_DEFAULT_FANCY_HELICES)
        cmd.set("cartoon_discrete_colors", constants.PYMOL_DEFAULT_CARTOON_DISCRETE_COLORS)
        cmd.set("cartoon_sampling", constants.PYMOL_DEFAULT_CARTOON_SAMPLING)
        cmd.set("spec_power", constants.PYMOL_DEFAULT_SPEC_POWER)
        cmd.set("spec_reflect", constants.PYMOL_DEFAULT_SPEC_REFLECT)
        cmd.set("ray_transparency_contrast", constants.PYMOL_DEFAULT_RAY_TRANSPARENCY_CONTRAST)
        cmd.set("ray_transparency_oblique", constants.PYMOL_DEFAULT_RAY_TRANSPARENCY_OBLIQUE)  # noqa: E501
        cmd.set("ray_transparency_oblique_power", constants.PYMOL_DEFAULT_RAY_OBLIQUE_POWER)
        cmd.set("ray_trace_color", constants.PYMOL_DEFAULT_RAY_COLOR)
        cmd.unset("depth_cue")
    except CmdException as e:
        return False, e
    else:
        return True, ""


def ray(a_width: int, a_height: int, a_renderer: int):
    try:
        cmd.ray(width=a_width, height=a_height, renderer=a_renderer)
    except CmdException as e:
        return False, e
    else:
        return True, ""


def draw(a_width: int, a_height: int, an_antialias_value: int):
    try:
        cmd.draw(width=a_width, height=a_height, antialias=an_antialias_value)
    except CmdException as e:
        return False, e
    else:
        return True, ""


def png(an_image_filepath: str, a_dpi_value: int):
    try:
        cmd.png(filename=an_image_filepath, dpi=a_dpi_value)
    except CmdException as e:
        return False, e
    else:
        return True, ""
