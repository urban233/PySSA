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
"""Module contains job summary dataclasses."""
from dataclasses import dataclass
from pyssa.internal.data_structures.data_classes import prediction_protein_info
from pyssa.util import enums


@dataclass
class JobBaseInformation:
    """Class which contains all general information about a job."""
    job_type: "enums.JobType"
    project_name: str
    protein_names: list[str]
    protein_pair_names: list[str]
    job_progress: "enums.JobProgress"

    def add_image_filepath(self, a_filepath: str):
        self._filepath = a_filepath

    def get_image_filepath(self):
        if self._filepath is None:
            return ""
        return self._filepath


@dataclass
class PredictionJobSummary:
    """Class which contains all information about a prediction job."""

    prediction_protein_infos: list[prediction_protein_info.PredictionProteinInfo]

    def get_protein_names(self):
        tmp_protein_names = []
        for tmp_prediction_info in self.prediction_protein_infos:
            tmp_protein_names.append(tmp_prediction_info.name)
        return tmp_protein_names


@dataclass
class DistanceAnalysisJobSummary:
    """Class which contains all information about a distance analysis job."""

    analysis_names: list[str]
