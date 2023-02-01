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
"""This module contains helper function for the analysis process."""
import pathlib
import logging
from pyssa.util import constants
from pyssa.logging_pyssa import log_handlers

logger = logging.getLogger(__file__)
logger.addHandler(log_handlers.log_file_handler)


def get_analysis_runs(self):
    """This function creates a data format which is used for the analysis

    """
    batch_analysis = []
    for row_no in range(self.list_analysis_overview.count()):
        tmp_batch_analysis = self.list_analysis_overview.item(row_no).text()
        separator_index = tmp_batch_analysis.find("_vs_")
        prot_1 = tmp_batch_analysis[:separator_index]

        if prot_1.find(";") != -1:
            prot_1_name = prot_1[:prot_1.find(";")]
            prot_1_chains = prot_1[prot_1.find(";") + 1:].split(",")
        else:
            prot_1_name = prot_1
            prot_1_chains = None
        prot_2 = tmp_batch_analysis[separator_index + 4:]
        if prot_2.find(";") != -1:
            prot_2_name = prot_2[:prot_2.find(";")]
            prot_2_chains = prot_2[prot_2.find(";") + 1:].split(",")
        else:
            prot_2_name = prot_2
            prot_2_chains = None
        tmp_prot_1 = protein_analysis_info.ProteinAnalysisInfo(prot_1_name, prot_1_chains, tmp_batch_analysis)
        tmp_prot_2 = protein_analysis_info.ProteinAnalysisInfo(prot_2_name, prot_2_chains, tmp_batch_analysis)
        batch_analysis.append((tmp_prot_1, tmp_prot_2))

    transformer = data_transformer.DataTransformer(self.app_project, batch_analysis)
    # contains analysis-ready data format: list(tuple(prot_1, prot_2, export_dir, name), ...)
    return transformer.transform_data_for_analysis()