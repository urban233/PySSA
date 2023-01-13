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
"""Module for handling everything related to the settings.xml"""
import json
import os
import pathlib
from pyssa.gui.utilities import constants
from pyssa.gui.data_structures import safeguard


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
        """Constructor

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
        self.prediction_path = pathlib.Path(f"{os.path.expanduser('~')}/Downloads")
        self.cycles: int = 0
        self.cutoff: float = 1.0
        self.app_launch = 1
        self.dir_settings: str = dir_settings
        self.filename: str = filename
        self.wsl_install: int = 0
        self.local_colabfold: int = 0
        self.wsl_username: str = "no_user_name"

    def serialize_settings(self) -> None:
        """This function serialize the protein object

        """
        if not os.path.exists(constants.SETTINGS_DIR):
            os.mkdir(constants.SETTINGS_DIR)
        settings_dict = self.__dict__
        update = {
            'workspace_path': str(self.workspace_path),
            'prediction_path': str(self.prediction_path),
            'dir_settings': str(self.dir_settings),
        }
        settings_dict.update(update)
        settings_file = open(f"{constants.SETTINGS_DIR}\\{constants.SETTINGS_FILENAME}", "w", encoding="utf-8")
        json.dump(settings_dict, settings_file, indent=4)

    @staticmethod
    def deserialize_settings():
        """This function constructs the protein object from
        the json file

        Returns:
            a complete project object deserialized from a json file
        """
        try:
            settings_obj_file = open(pathlib.Path(f"{constants.SETTINGS_DIR}\\{constants.SETTINGS_FILENAME}"), "r", encoding="utf-8")
        except FileNotFoundError:
            print(f"There is no valid json file under: {constants.SETTINGS_DIR}\\{constants.SETTINGS_FILENAME}. Please restore the settings!")
            return
        settings_dict = json.load(settings_obj_file)
        tmp_settings: Settings = Settings(settings_dict.get("dir_settings"), settings_dict.get("filename"))
        if safeguard.Safeguard.check_filepath(settings_dict.get("workspace_path")):
            tmp_settings.workspace_path = settings_dict.get("workspace_path")
        else:
            raise ValueError
        if safeguard.Safeguard.check_filepath(settings_dict.get("prediction_path")):
            tmp_settings.prediction_path = settings_dict.get("prediction_path")
        else:
            raise ValueError
        if safeguard.Safeguard.check_if_number_is_positive(int(settings_dict.get("cycles"))):
            tmp_settings.cycles = settings_dict.get("cycles")
        else:
            raise ValueError
        if safeguard.Safeguard.check_if_number_is_positive(float(settings_dict.get("cutoff"))):
            tmp_settings.cutoff = settings_dict.get("cutoff")
        else:
            raise ValueError
        if int(settings_dict.get("app_launch")) == 0 or int(settings_dict.get("app_launch")) == 1:
            tmp_settings.app_launch = settings_dict.get("app_launch")
        else:
            raise ValueError
        if not safeguard.Safeguard.check_filepath(settings_dict.get("dir_settings")):
            raise ValueError
        if not settings_dict.get("filename") == "settings.json":
            raise ValueError
        if int(settings_dict.get("wsl_install")) == 0 or int(settings_dict.get("wsl_install")) == 1:
            tmp_settings.wsl_install = int(settings_dict.get("wsl_install"))
        else:
            raise ValueError
        if int(settings_dict.get("local_colabfold")) == 0 or int(settings_dict.get("local_colabfold")) == 1:
            tmp_settings.local_colabfold = int(settings_dict.get("local_colabfold"))
        else:
            raise ValueError
        tmp_settings.wsl_username = settings_dict.get("wsl_username")
        return tmp_settings

    def get_workspace_path(self) -> pathlib.Path:
        """This function gets the value of the workspace_path variable

        Returns (str):
            workspace_path
        """
        return self.workspace_path

    def set_workspace_path(self, value) -> None:
        """This function gets the value of the workspace_path variable

        """
        self.workspace_path = value

    def get_prediction_path(self) -> pathlib.Path:
        """This function gets the value of the prediction_path variable

        Returns (str):
            prediction_path
        """
        return self.prediction_path

    def set_prediction_path(self, value) -> None:
        """This function gets the value of the prediction_path variable

        """
        self.prediction_path = value

    def get_cycles(self) -> int:
        """This function gets the value of the cycles variable

        Returns (int):
            cycles
        """
        return int(self.cycles)

    def set_cycles(self, value) -> None:
        """This function gets the value of the cycles variable

        """
        self.cycles = value

    def get_cutoff(self) -> float:
        """This function gets the value of the cutoff variable

        Returns (float):
            cutoff
        """
        return float(self.cutoff)

    def set_cutoff(self, value) -> None:
        """This function gets the value of the cutoff variable

        """
        self.cutoff = value

    def get_app_launch(self) -> int:
        """This function gets the value of the app_launch variable

       Returns (int):
           app_launch
       """
        return int(self.app_launch)

    def set_app_launch(self, value) -> None:
        """This function gets the value of the app_launch variable

        """
        self.app_launch = value

    def restore_settings(self, dir_settings, filename):
        """This function resets the settings to the default values

        """
        self.workspace_path = constants.DEFAULT_WORKSPACE_PATH
        self.prediction_path = pathlib.Path(f"{os.path.expanduser('~')}/Downloads")
        self.cycles: int = 0
        self.cutoff: float = 1.0
        self.app_launch = 1
        self.dir_settings: str = dir_settings
        self.filename: str = filename
        self.wsl_install: int = 0
        self.local_colabfold: int = 0
        self.serialize_settings()
