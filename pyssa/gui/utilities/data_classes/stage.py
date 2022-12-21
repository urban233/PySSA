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
from dataclasses import dataclass


@dataclass
class Stage:
    """This class stores all gui elements of one stage in a specific form

    Notes:
        List of valid keys for primary_elements dict:
            - label_protein_name
            - text_field_protein_name
            - label_protein_name_status
            - label_protein_sequence
            - text_field_protein_sequence
            - label_protein_sequence_status
            - button_add_sequence
            - label_protein_sequence_overview
            - table_protein_sequence_overview
            - label_advanced_config
            - button_advanced_config
            - label_prediction_mode
            - label_protein_sequence_batch
            - button_add_sequence_batch
            - label_advanced_config_batch
            - button_advanced_config_batch

        List of valid keys for control_elements dict:
            - back_button
            - next_button
            - predict_button
            - single_button
            - batch_button
    """
    primary_elements: dict
    control_elements: dict

