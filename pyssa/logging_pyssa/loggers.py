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
"""This module contains all loggers used in the pyssa project."""
import logging
from pyssa.logging_pyssa import log_handlers

# all loggers used in the pyssa project
pyssa = logging.getLogger("PySSA-Logger")
prediction_worker = logging.getLogger("PredictionWorker")
analysis_worker = logging.getLogger("AnalysisWorker")

# adding necessary handlers
pyssa.addHandler(log_handlers.log_file_handler)
prediction_worker.addHandler(log_handlers.log_file_handler)
analysis_worker.addHandler(log_handlers.log_file_handler)


def log_multiple_messages(logger, level, messages: list[str]):
    for message in messages:
        logger.log(level=level, msg=message)
