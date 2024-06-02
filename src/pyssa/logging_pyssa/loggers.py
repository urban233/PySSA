#
# PySSA - Python-Plugin for Sequence-to-Structure Analysis
# Copyright (C) 2024
# Martin Urban (martin.urban@studmail.w-hs.de)
# Hannah Kullik (hannah.kullik@studmail.w-hs.de)
#
# Source code is available at <https://github.com/zielesny/PySSA>
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
"""Module contains all loggers used in the pyssa project."""
from src.pyssa.logging_pyssa import log_handlers
from src.pyssa.util import constants


# adding necessary handlers
constants.PYSSA_LOGGER.addHandler(log_handlers.log_file_handler)
constants.PREDICTION_WORKER_LOGGER.addHandler(log_handlers.log_file_handler)
constants.ANALYSIS_WORKER_LOGGER.addHandler(log_handlers.log_file_handler)
