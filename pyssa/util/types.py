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
from PyQt5 import QtWidgets
from pyssa.internal.data_structures import protein
from pyssa.internal.data_structures import protein_pair
from pyssa.internal.data_structures import project
from pyssa.internal.data_structures import settings
from internal.data_structures.data_classes import protein_analysis_info
from internal.data_structures.data_classes import prediction_configuration
from pyssa.internal.data_structures import sequence
from pyssa.internal.data_structures import chain
from pyssa.internal.data_structures import selection
from pyssa.internal.analysis_types import distance_analysis

PROTEIN = protein.Protein
PROTEIN_PAIR = protein_pair.ProteinPair
CHAIN = chain.Chain
SELECTION = selection.Selection

PROTEIN_ANALYSIS_INFO = protein_analysis_info.ProteinAnalysisInfo
PROTEIN_SEQUENCE = sequence.ProteinSequence

PREDICTION_CONFIG = prediction_configuration.PredictionConfiguration

PROJECT = project.Project
SETTINGS = settings.Settings

TABLE_WIDGET = QtWidgets.QTableWidget
LIST_WIDGET = QtWidgets.QListWidget
CHECKBOX_WIDGET = QtWidgets.QCheckBox
STATUS_BAR_WIDGET = QtWidgets.QStatusBar

DISTANCE_ANALYSIS = distance_analysis.DistanceAnalysis
