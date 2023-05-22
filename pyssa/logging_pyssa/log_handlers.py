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
"""This module contains the handlers used by the loggers."""
import logging
import datetime
from pyssa.util import constants

current_time = datetime.datetime.now()
filename = f"{current_time.year}-{current_time.month:02d}-{current_time.day:02d}_{current_time.hour:02d}-{current_time.minute:02d}.log"

log_formatter = logging.Formatter("%(asctime)s: %(name)s %(levelname)s - %(message)s")
log_file_handler = logging.FileHandler(constants.LOG_FILEPATH)
log_file_handler.setLevel(logging.DEBUG)
log_file_handler.setFormatter(log_formatter)
