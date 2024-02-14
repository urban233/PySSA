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
"""Module for handling everything related to the settings.xml."""
import json
import os
import pathlib

from pyssa.io_pyssa import safeguard
from pyssa.util import constants


class Settings:
    """This class is for the handling of the settings.xml file.

    Var:
        workspace_path:
            path of the workspace
        model_path:
            path where models are stored
        prediction_path:
            path where the models from predictions are initially stored
        cycles:
            number of structure alignment cycles
        cutoff:
            cutoff of the structure alignment
    """

    def __init__(self, dir_settings: str, filename: str) -> None:
        """Constructor.

        Args:
            dir_settings:
                directory where the settings.xml is stored
            filename:
                name of the settings.xml
        """
        if not os.path.exists(constants.SETTINGS_DIR):
            os.mkdir(constants.SETTINGS_DIR)
        if not os.path.exists(constants.DEFAULT_WORKSPACE_PATH):
            os.mkdir(constants.DEFAULT_WORKSPACE_PATH)

        self.workspace_path = constants.DEFAULT_WORKSPACE_PATH
        self.cycles: int = 0
        self.cutoff: float = 1.0
        self.app_launch = 1
        self.dir_settings: str = dir_settings
        self.filename: str = filename
        self.local_colabfold: int = 1
        self.wsl_install: int = 1
        self.color_vision_mode: str = "normal"
        self.ask_save_pymol_session: int = 0
        self.image_background_color: str = "black"
        self.image_renderer: str = "-1"  # or "0"
        self.image_ray_trace_mode: int = 1  # ranges from 0 to 3
        self.image_ray_texture: int = 0  # ranges from 0 to 5
        self.start_help_at_startup: int = 1

    def serialize_settings(self) -> None:
        """This function serialize the protein object."""
        if not os.path.exists(constants.SETTINGS_DIR):
            os.mkdir(constants.SETTINGS_DIR)
        settings_dict = self.__dict__
        update = {
            "workspace_path": str(self.workspace_path),
            "dir_settings": str(self.dir_settings),
        }
        settings_dict.update(update)
        settings_file = open(f"{constants.SETTINGS_DIR}\\{constants.SETTINGS_FILENAME}", "w", encoding="utf-8")
        json.dump(settings_dict, settings_file, indent=4)

    @staticmethod
    def deserialize_settings() -> "Settings":
        """This function constructs the protein object from the json file.

        Returns:
            a complete project object deserialized from a json file
        """
        try:
            settings_obj_file = open(
                pathlib.Path(f"{constants.SETTINGS_DIR}\\{constants.SETTINGS_FILENAME}"),
                "r",
                encoding="utf-8",
            )
        except FileNotFoundError:
            # TODO: log message is needed
            return  # noqa: RET502 #TODO: needs to be fixed
        settings_dict = json.load(settings_obj_file)

        try:
            tmp_settings: Settings = Settings(settings_dict.get("dir_settings"), settings_dict.get("filename"))

            tmp_settings.workspace_path = Settings._check_integrity_of_workspace_path(
                settings_dict.get("workspace_path")
            )
            tmp_settings.cycles = Settings._check_integrity_of_cycles(int(settings_dict.get("cycles")))
            tmp_settings.cutoff = Settings._check_integrity_of_cutoff(float(settings_dict.get("cutoff")))
            tmp_settings.app_launch = Settings._check_integrity_of_app_start_value(settings_dict.get("app_launch"))
            tmp_settings.wsl_install = Settings._check_integrity_of_wsl_install_flag(
                int(settings_dict.get("wsl_install"))
            )
            tmp_settings.local_colabfold = Settings._check_integrity_of_colabfold_install_flag(
                int(settings_dict.get("local_colabfold"))
            )
            tmp_settings.color_vision_mode = Settings._check_integrity_of_color_blindness_value(
                settings_dict.get("color_vision_mode")
            )
            tmp_settings.ask_save_pymol_session = Settings._check_integrity_of_save_pymol_session_flag(
                int(settings_dict.get("ask_save_pymol_session"))
            )
            tmp_settings.image_background_color = Settings._check_integrity_of_bg_color_value(
                settings_dict.get("image_background_color")
            )
            tmp_settings.image_renderer = Settings._check_integrity_of_renderer_value(
                settings_dict.get("image_renderer")
            )
            tmp_settings.image_ray_trace_mode = Settings._check_integrity_of_ray_trace_mode_value(
                settings_dict.get("image_ray_trace_mode")
            )
            tmp_settings.image_ray_texture = Settings._check_integrity_of_ray_texture_value(
                settings_dict.get("image_ray_texture")
            )
            tmp_settings.start_help_at_startup = Settings._check_integrity_of_start_help_at_startup_flag(
                int(settings_dict.get("start_help_at_startup"))
            )

        except ValueError as e:
            raise AttributeError(f"An error occurred during the deserialization of the settings: {e}")
        return tmp_settings

    # <editor-fold desc="Integrity checks">
    @staticmethod
    def _check_integrity_of_workspace_path(a_workspace_path: str) -> str:
        if safeguard.Safeguard.check_filepath(pathlib.Path(a_workspace_path)):
            return a_workspace_path
        else:
            raise ValueError

    @staticmethod
    def _check_integrity_of_cycles(a_cycles_value: int) -> int:
        if safeguard.Safeguard.check_if_number_is_positive(a_cycles_value):
            return a_cycles_value
        else:
            raise ValueError

    @staticmethod
    def _check_integrity_of_cutoff(a_cutoff_value: float) -> float:
        if safeguard.Safeguard.check_if_number_is_positive(a_cutoff_value):
            return a_cutoff_value
        else:
            raise ValueError

    @staticmethod
    def _check_integrity_of_app_start_value(a_app_start_value: int) -> int:
        if a_app_start_value == 0 or a_app_start_value == 1:
            return a_app_start_value
        else:
            raise ValueError

    @staticmethod
    def _check_integrity_of_wsl_install_flag(a_wsl_install_flag: int) -> int:
        if a_wsl_install_flag == 0 or a_wsl_install_flag == 1:
            return a_wsl_install_flag
        else:
            raise ValueError

    @staticmethod
    def _check_integrity_of_colabfold_install_flag(a_colabfold_install_flag: int) -> int:
        if a_colabfold_install_flag == 0 or a_colabfold_install_flag == 1:
            return a_colabfold_install_flag
        else:
            raise ValueError

    @staticmethod
    def _check_integrity_of_color_blindness_value(a_color_blindness_value: str) -> str:
        if a_color_blindness_value != "":
            return a_color_blindness_value
        else:
            raise ValueError

    @staticmethod
    def _check_integrity_of_save_pymol_session_flag(save_pymol_session_flag: int) -> int:
        if save_pymol_session_flag == 0 or save_pymol_session_flag == 1:
            return save_pymol_session_flag
        else:
            raise ValueError

    @staticmethod
    def _check_integrity_of_bg_color_value(a_bg_color_value: str) -> str:
        if a_bg_color_value == "white" or a_bg_color_value == "black":
            return a_bg_color_value
        else:
            raise ValueError(f"Invalid background color! The value is not black or white: {a_bg_color_value}")

    @staticmethod
    def _check_integrity_of_renderer_value(a_renderer_value: str) -> str:
        if a_renderer_value == "0" or a_renderer_value == "-1":
            return a_renderer_value
        else:
            raise ValueError(f"Invalid renderer! The value is not 0 or -1: {a_renderer_value}")

    @staticmethod
    def _check_integrity_of_ray_trace_mode_value(a_ray_trace_mode_value: int) -> int:
        tmp_possible_values: list = list(range(3))
        if a_ray_trace_mode_value in tmp_possible_values:
            return a_ray_trace_mode_value
        else:
            raise ValueError(f"Invalid ray trace mode! The value is not between 0 and 3: {a_ray_trace_mode_value}")

    @staticmethod
    def _check_integrity_of_ray_texture_value(a_ray_texture_value: int) -> int:
        tmp_possible_values: list = list(range(6))
        if a_ray_texture_value in tmp_possible_values:
            return a_ray_texture_value
        else:
            raise ValueError(f"Invalid ray trace mode! The value is not between 0 and 3: {a_ray_texture_value}")

    @staticmethod
    def _check_integrity_of_start_help_at_startup_flag(start_help_at_startup_flag: int) -> int:
        if start_help_at_startup_flag == 0 or start_help_at_startup_flag == 1:
            return start_help_at_startup_flag
        else:
            raise ValueError

    # </editor-fold>

    def get_workspace_path(self) -> pathlib.Path:
        """This function gets the value of the workspace_path variable.

        Returns (str):
            workspace_path
        """
        return self.workspace_path

    def set_workspace_path(self, value: pathlib.Path) -> None:
        """This function gets the value of the workspace_path variable."""
        self.workspace_path = value

    def get_prediction_path(self) -> pathlib.Path:
        """This function gets the value of the prediction_path variable.

        Returns (str):
            prediction_path
        """
        return self.prediction_path

    def get_cycles(self) -> int:
        """This function gets the value of the cycles variable.

        Returns (int):
            cycles
        """
        return int(self.cycles)

    def set_cycles(self, value: int) -> None:
        """This function gets the value of the cycles variable."""
        self.cycles = value

    def get_cutoff(self) -> float:
        """This function gets the value of the cutoff variable.

        Returns (float):
            cutoff
        """
        return float(self.cutoff)

    def set_cutoff(self, value: float) -> None:
        """This function gets the value of the cutoff variable."""
        self.cutoff = value

    def get_app_launch(self) -> int:
        """This function gets the value of the app_launch variable.

        Returns (int):
           app_launch
        """
        return int(self.app_launch)

    def set_app_launch(self, value: int) -> None:
        """This function gets the value of the app_launch variable."""
        self.app_launch = value

    def restore_settings(self, dir_settings: str, filename: str) -> None:
        """Resets the settings to the default values."""
        self.workspace_path = constants.DEFAULT_WORKSPACE_PATH
        self.cycles: int = 0
        self.cutoff: float = 1.0
        self.app_launch = 1
        self.dir_settings: str = dir_settings
        self.filename: str = filename
        self.wsl_install: int = 0
        self.local_colabfold: int = 0
        self.ask_save_pymol_session: int = 0
        self.image_background_color: str = "black"
        self.image_renderer: str = "-1"  # or "0"
        self.image_ray_trace_mode: int = 1  # ranges from 0 to 3
        self.image_ray_texture: int = 0  # ranges from 0 to 5
        self.start_help_at_startup: bool = True

        self.serialize_settings()
