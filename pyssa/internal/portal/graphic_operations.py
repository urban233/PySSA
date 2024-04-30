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
"""Module for graphic operations in pymol."""
from pyssa.internal.data_structures import protein_pair
from pyssa.util import constants, protein_pair_util


# def get_residue_colors(a_selection_string: str):
#     pymol.stored.colors = []
#     cmd.iterate(a_selection_string, "stored.colors.append((chain, resi, name, color))")
#     res_colors = {}
#     for chain, resi, name, color_index in pymol.stored.colors:
#         if name == 'CA':  # c-alpha atom
#             res_colors[(chain, resi, name)] = constants.PYMOL_COLORS_WITH_INDICES[color_index]
#     return res_colors
#
#
# def get_chain_color(a_selection_string: str, chain_letter: str):
#     pymol.stored.colors = []
#     cmd.iterate(a_selection_string, "stored.colors.append((chain, resi, name, color))")
#     tmp_ca_atom_colors = []
#     tmp_residue_atom_colors = []
#     tmp_is_colored_by_elements = False
#     for chain, resi, name, color_index in pymol.stored.colors:
#         if chain == chain_letter and name == 'CA':  # c-alpha atom
#             #print(color_index)
#             try:
#                 tmp_ca_atom_colors.append(constants.PYMOL_COLORS_WITH_INDICES[color_index])
#             except KeyError:
#                 tmp_ca_atom_colors.append("By RMSD")
#         else:
#             try:
#                 tmp_ca_atom_colors.append(constants.PYMOL_COLORS_WITH_INDICES[color_index])
#             except KeyError:
#                 tmp_is_colored_by_elements = True
#     tmp_chain_color = list(set(tmp_ca_atom_colors))[0]
#     return tmp_chain_color, tmp_ca_atom_colors, tmp_is_colored_by_elements
#
#
# def get_chain_repr_state(a_selection_string: str, chain_letter: str) -> dict:
#     """Returns the index for the current representation state.
#
#     Note:
#         If no match occurs, then -1 is returned
#     """
#     pymol.stored.reps = []
#     cmd.iterate(a_selection_string, "stored.reps.append((chain, resi, name, reps))")
#     for chain, resi, name, reps_index in pymol.stored.reps:
#         if chain == chain_letter:  # c-alpha atom
#             constants.PYSSA_LOGGER.debug("Representation index: " + str(reps_index))
#             return constants.PYMOL_REPR_STATES_WITH_INDICES[reps_index]
#     return {-1: ""}
#
#
# def show_protein_selection_as_balls_and_sticks(selection: str) -> None:
#     """Shows the protein as balls and sticks in representation mode.
#
#     Args:
#         selection: a pymol selection string.
#     """
#     if not pymol_safeguard.PymolSafeguard.check_if_protein_in_session():
#         raise pymol.CmdException("No protein is in pymol session.")
#     try:
#         cmd.show(representation="sticks", selection=selection)
#     except pymol.CmdException:
#         print("No sticks and balls can be shown in protein.")
#
#
# def hide_protein_selection_as_balls_and_sticks(selection: str) -> None:
#     """Hides the protein from balls and sticks.
#
#     Args:
#         selection: a pymol selection string.
#
#     """
#     if not pymol_safeguard.PymolSafeguard.check_if_protein_in_session():
#         raise pymol.CmdException("No protein is in pymol session.")
#     try:
#         cmd.hide(representation="sticks", selection=selection)
#     except pymol.CmdException:
#         print("No sticks and balls can be hidden in protein.")
#
#
# def zoom_to_residue_in_protein_position(selection: str) -> None:
#     """Zooms to the residue in protein.
#
#     Args:
#         selection: a pymol selection string
#
#     """
#     if not pymol_safeguard.PymolSafeguard.check_if_protein_in_session():
#         raise pymol.CmdException("No protein is in pymol session.")
#     try:
#         cmd.zoom(selection=selection, buffer=8.0, state=0, complete=0)
#     except pymol.CmdException:
#         print("No residue can be shown in protein.")


# def color_protein(pymol_color: str, a_selection_string: str) -> None:
#     """Colors a specific protein selection with a given PyMOL color.
#
#     Args:
#         pymol_color: a color which is available in PyMOL
#         a_selection_string: a PyMOL conform selection string
#
#     """
#     if pymol_color == "":
#         return
#     if pymol_color not in constants.PYMOL_COLORS:
#         raise ValueError(f"An illegal color argument. {pymol_color}")
#     if not pymol_safeguard.PymolSafeguard.check_if_protein_in_session():
#         raise pymol.CmdException("No protein is in pymol session.")
#     try:
#         cmd.color(pymol_color, a_selection_string)
#     except pymol.CmdException:
#         print("Color process was unsuccessful.")


# def color_protein_pair_by_rmsd(a_protein_pair: "protein_pair.ProteinPair") -> None:
#     """Colors a specific protein pair based on their rmsd value.
#
#     Args:
#         a_protein_pair: the protein pair to color.
#     """
#     cutoff_1 = 0.5
#     cutoff_2 = 1.0
#     cutoff_3 = 2
#     cutoff_4 = 4
#     cutoff_5 = 6
#
#     color_1 = "br0"
#     color_2 = "br2"
#     color_3 = "br4"
#     color_4 = "br6"
#     color_5 = "br8"
#     color_6 = "red"
#
#     cmd.color("hydrogen", a_protein_pair.protein_2.get_molecule_object())
#
#     i: int = 0
#     for distance_value in a_protein_pair.distance_analysis.analysis_results.distance_data.get("distance"):
#         if distance_value <= cutoff_1:
#             atom_info = protein_pair_util.get_chain_and_position(
#                 a_protein_pair.distance_analysis.analysis_results.distance_data,
#                 i,
#             )
#             # create two atoms for the get_distance command
#             atom1: str = (
#                 f"/{a_protein_pair.protein_1.get_molecule_object()}//" f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
#             )
#             atom2: str = (
#                 f"/{a_protein_pair.protein_2.get_molecule_object()}//" f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
#             )
#             # coloring
#             cmd.color(color_1, atom1)
#             cmd.color(color_1, atom2)
#             i += 1
#
#         elif distance_value <= cutoff_2:
#             atom_info = protein_pair_util.get_chain_and_position(
#                 a_protein_pair.distance_analysis.analysis_results.distance_data,
#                 i,
#             )
#             # create two atoms for the get_distance command
#             atom1: str = (
#                 f"/{a_protein_pair.protein_1.get_molecule_object()}//" f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
#             )
#             atom2: str = (
#                 f"/{a_protein_pair.protein_2.get_molecule_object()}//" f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
#             )
#             # coloring
#             cmd.color(color_2, atom1)
#             cmd.color(color_2, atom2)
#             i += 1
#
#         elif distance_value <= cutoff_3:
#             atom_info = protein_pair_util.get_chain_and_position(
#                 a_protein_pair.distance_analysis.analysis_results.distance_data,
#                 i,
#             )
#             # create two atoms for the get_distance command
#             atom1: str = (
#                 f"/{a_protein_pair.protein_1.get_molecule_object()}//"
#                 f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}/CA"
#             )
#             atom2: str = (
#                 f"/{a_protein_pair.protein_2.get_molecule_object()}//"
#                 f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}/CA"
#             )
#             # coloring
#             cmd.color(color_3, atom1)
#             cmd.color(color_3, atom2)
#             i += 1
#
#         elif distance_value <= cutoff_4:
#             atom_info = protein_pair_util.get_chain_and_position(
#                 a_protein_pair.distance_analysis.analysis_results.distance_data,
#                 i,
#             )
#             # create two atoms for the get_distance command
#             atom1: str = (
#                 f"/{a_protein_pair.protein_1.get_molecule_object()}//" f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
#             )
#             atom2: str = (
#                 f"/{a_protein_pair.protein_2.get_molecule_object()}//" f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
#             )
#             # coloring
#             cmd.color(color_4, atom1)
#             cmd.color(color_4, atom2)
#             i += 1
#
#         elif distance_value <= cutoff_5:
#             atom_info = protein_pair_util.get_chain_and_position(
#                 a_protein_pair.distance_analysis.analysis_results.distance_data,
#                 i,
#             )
#             # create two atoms for the get_distance command
#             atom1: str = (
#                 f"/{a_protein_pair.protein_1.get_molecule_object()}//" f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
#             )
#             atom2: str = (
#                 f"/{a_protein_pair.protein_2.get_molecule_object()}//" f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
#             )
#             # coloring
#             cmd.color(color_5, atom1)
#             cmd.color(color_5, atom2)
#             i += 1
#
#         elif distance_value > cutoff_5:
#             atom_info = protein_pair_util.get_chain_and_position(
#                 a_protein_pair.distance_analysis.analysis_results.distance_data,
#                 i,
#             )
#             # create two atoms for the get_distance command
#             atom1: str = (
#                 f"/{a_protein_pair.protein_1.get_molecule_object()}//" f"{atom_info[0]}/{atom_info[2]}`{atom_info[1]}"
#             )
#             atom2: str = (
#                 f"/{a_protein_pair.protein_2.get_molecule_object()}//" f"{atom_info[3]}/{atom_info[5]}`{atom_info[4]}"
#             )
#             # coloring
#             cmd.color(color_6, f"({atom1})")
#             cmd.color(color_6, f"({atom2})")
#             i += 1


# def setup_default_session_graphic_settings() -> None:
#     """This functions modifies the pymol session to look fancy."""
#     cmd.bg_color(constants.PYMOL_DEFAULT_BACKGROUND_COLOR)
#     cmd.set("valence", 0)
#     cmd.set("scene_buttons", 0)
#     cmd.set("ray_trace_mode", constants.PYMOL_DEFAULT_RAY_TRACE_MODE)
#     cmd.set("antialias", constants.PYMOL_DEFAULT_ANTIALIAS)
#     cmd.set("ambient", constants.PYMOL_DEFAULT_AMBIENT)
#     cmd.set("cartoon_fancy_helices", constants.PYMOL_DEFAULT_FANCY_HELICES)
#     cmd.set("cartoon_discrete_colors", constants.PYMOL_DEFAULT_CARTOON_DISCRETE_COLORS)
#     cmd.set("cartoon_sampling", constants.PYMOL_DEFAULT_CARTOON_SAMPLING)
#     cmd.set("spec_power", constants.PYMOL_DEFAULT_SPEC_POWER)
#     cmd.set("spec_reflect", constants.PYMOL_DEFAULT_SPEC_REFLECT)
#     cmd.set("ray_transparency_contrast", constants.PYMOL_DEFAULT_RAY_TRANSPARENCY_CONTRAST)
#     cmd.set("ray_transparency_oblique", constants.PYMOL_DEFAULT_RAY_TRANSPARENCY_OBLIQUE)  # noqa: E501
#     cmd.set("ray_transparency_oblique_power", constants.PYMOL_DEFAULT_RAY_OBLIQUE_POWER)
#     cmd.set("ray_trace_color", constants.PYMOL_DEFAULT_RAY_COLOR)
#     cmd.unset("depth_cue")


# def setup_default_image_graphic_settings(ray_shadows: bool, opaque_background: int = 0) -> None:
#     """Sets up the default image graphic settings for PyMOL.
#
#     Args:
#         ray_shadows (bool): false if no shadows, true if shadows should be displayed.
#         opaque_background (int, optional): 0 for a transparent background and 1 for a white background.
#     """
#     if not ray_shadows:
#         opt_ray_shadows: str = "off"
#     else:
#         opt_ray_shadows: str = "on"
#     cmd.bg_color(constants.PYMOL_DEFAULT_BACKGROUND_COLOR)
#     cmd.set("ray_trace_mode", constants.PYMOL_DEFAULT_RAY_TRACE_MODE)
#     cmd.set("antialias", constants.PYMOL_DEFAULT_ANTIALIAS)
#     cmd.set("ray_shadows", opt_ray_shadows)
#     cmd.set("ray_opaque_background", opaque_background)
#
#
# def setup_default_graphic_settings_for_interesting_regions() -> None:
#     """Sets up the default graphic settings for interesting regions."""
#     cmd.set("label_size", 14)
#     cmd.set("label_font_id", 13)
#     cmd.set("label_color", "hotpink")
#     cmd.set("depth_cue", 0)
#     cmd.bg_color(constants.PYMOL_DEFAULT_BACKGROUND_COLOR)
#     # interacts directly with molecule objects in the session
#     cmd.hide("cartoon", "all")
#     cmd.show("ribbon", "all")
