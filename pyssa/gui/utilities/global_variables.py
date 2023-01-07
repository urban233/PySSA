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
"""Module for global variables which will be used across this project"""

import os
from pyssa.gui.data_structures import settings
from pyssa.gui.utilities import constants


global_var_root_dir = os.path.realpath(os.path.join(os.path.dirname(__file__), "..", ".."))
# global_var is used to access the settings project-wide
#global_var_settings_obj = settings.Settings(constants.SETTINGS_DIR, constants.SETTINGS_FILE)
#global_var_settings_obj.deserialize_settings()
global_var_tmp_project_info = []
global_var_workspace_proteins: dict = {}
