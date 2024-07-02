#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
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
"""Module contains an enumeration for pymol commands."""
import enum

__docformat__ = "google"


class CommandEnum(enum.Enum):
  """Contains enumeration for pymol commands."""

  CLOSE_USER_PYMOL = "close_user_pymol"

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
  GET_RESIDUE_COLOR_CONFIG = "get_residue_color_config"

  COLOR_SELECTION = "color_selection"

  # Image commands
  RAY = "ray"
  DRAW = "draw"
  PNG = "png"

  # Session commands
  RESET = "reset"
