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
"""Module for the gui page management class"""
from pyssa.util import tools, gui_utils
from pyssa.gui.ui.styles import styles
from pyssa.internal.data_structures.data_classes import stage


class GuiPageManagement:
    """This class is used to control gui elements on one page"""

    """Var: elements_to_show is used to show gui elements after they where used -> shown after pressing a next button"""
    # this dict can to be expanded according to further needs
    elements_to_show: dict = {
        "text_field_protein_name": "label_protein_name",
        "text_field_protein_sequence": "label_protein_sequence",
        "table_protein_sequence_overview": "label_protein_sequence_overview",
        "label_protein_structure_1": "label_protein_structure_1",
        "label_protein_structure_2": "label_protein_structure_2",
        "label_vs": "label_vs",
        "list_protein_structure_1_chains": "label_protein_structure_1_chains",
        "list_protein_structure_2_chains": "label_protein_structure_2_chains",
    }

    """Var: status_labels is used to set an empty string to these specific labels"""
    status_labels: dict = {
        "label_protein_name_status": None,
        "label_protein_sequence_status": None,
    }

    def __init__(self, stages: list[stage.Stage]) -> None:
        """Constructor

        Args:
            stages:
                is a list of Stage objects
        """
        self.stages: list[stage.Stage] = stages

    def show_stage_x(self, x: int):
        primary_elements_to_show = []
        control_elements_to_show = []
        primary_elements_to_hide = []
        control_elements_to_hide = []
        i = 0
        for tmp_stage in self.stages:
            if i < x:
                # iterate over txt_lbl_dict
                for key in self.elements_to_show:
                    if tmp_stage.primary_elements.__contains__(key):
                        if key.find("text") != -1:
                            gui_utils.disable_text_box(
                                tmp_stage.primary_elements[key],
                                tmp_stage.primary_elements[self.elements_to_show[key]]
                            )
                        primary_elements_to_show.append(tmp_stage.primary_elements[key])
                        primary_elements_to_show.append(tmp_stage.primary_elements[self.elements_to_show[key]])

                # iterate over stages dict primary elements
                for key in tmp_stage.primary_elements:
                    if tmp_stage.primary_elements[key] not in primary_elements_to_show:
                        primary_elements_to_hide.append(tmp_stage.primary_elements[key])
                # iterate over stages dict control elements
                for key in tmp_stage.control_elements:
                    control_elements_to_hide.append(tmp_stage.control_elements[key])
            elif i == x:
                for key in self.elements_to_show:
                    if tmp_stage.primary_elements.__contains__(key):
                        gui_utils.enable_text_box(
                            tmp_stage.primary_elements[key],
                            tmp_stage.primary_elements[self.elements_to_show[key]]
                        )
                for key in tmp_stage.primary_elements:
                    primary_elements_to_show.append(tmp_stage.primary_elements[key])
                for key in tmp_stage.control_elements:
                    control_elements_to_show.append(tmp_stage.control_elements[key])
            else:
                for key in tmp_stage.primary_elements:
                    if tmp_stage.primary_elements[key] not in primary_elements_to_show:
                        primary_elements_to_hide.append(tmp_stage.primary_elements[key])
                for key in tmp_stage.control_elements:
                    if tmp_stage.control_elements[key] not in control_elements_to_show:
                        control_elements_to_hide.append(tmp_stage.control_elements[key])
            i += 1
        gui_utils.show_gui_elements(primary_elements_to_show)
        gui_utils.show_gui_elements(control_elements_to_show)
        gui_utils.hide_gui_elements(primary_elements_to_hide)
        gui_utils.hide_gui_elements(control_elements_to_hide)

    def show_gui_elements_stage_x(self, stages_to_show: list, stages_to_hide: list,
                                  hide_specific_elements=[], show_specific_elements=[]):
        """Shows and hides precisely specific stages

        Args:
            stages_to_show:
                list which contains the stage numbers of the stages which should get displayed
            stages_to_hide:
                 list which contains the stage numbers of the stages which should get hidden
            hide_specific_elements (optional):
                list of pyqt widgets which should get hidden
            show_specific_elements (optional):
                list of pyqt widgets which should get shown
        """
        tmp_stages_to_show = []
        tmp_stages_to_hide = []
        for stage_number in stages_to_show:
            tmp_stages_to_show.append(self.stages[stage_number])
        for stage_number in stages_to_hide:
            tmp_stages_to_hide.append(self.stages[stage_number])
        gui_elements_to_show = []
        gui_elements_to_hide = []
        i = 1
        for tmp_stage in tmp_stages_to_show:
            if i < len(stages_to_show):
                for key in tmp_stage.primary_elements:
                    gui_elements_to_show.append(tmp_stage.primary_elements[key])
            elif i == len(stages_to_show):
                for key in tmp_stage.primary_elements:
                    gui_elements_to_show.append(tmp_stage.primary_elements[key])
                for key in tmp_stage.control_elements:
                    gui_elements_to_show.append(tmp_stage.control_elements[key])
            i += 1
        for tmp_stage in tmp_stages_to_hide:
            for key in tmp_stage.primary_elements:
                if tmp_stage.primary_elements[key] not in gui_elements_to_show:
                    gui_elements_to_hide.append(tmp_stage.primary_elements[key])
            for key in tmp_stage.control_elements:
                if tmp_stage.control_elements[key] not in gui_elements_to_show:
                    gui_elements_to_hide.append(tmp_stage.control_elements[key])
        if len(show_specific_elements) != 0:
            for element in show_specific_elements:
                gui_elements_to_show.append(element)
        if len(hide_specific_elements) != 0:
            for element in hide_specific_elements:
                gui_elements_to_hide.append(element)
        gui_utils.show_gui_elements(gui_elements_to_show)
        gui_utils.hide_gui_elements(gui_elements_to_hide)

    def enable_text_boxes_stage_x(self, stages_to_enable: list):
        tmp_stages_to_enable = []
        for stage_number in stages_to_enable:
            tmp_stages_to_enable.append(self.stages[stage_number])
        for tmp_stage in tmp_stages_to_enable:
            for key in self.elements_to_show:
                if tmp_stage.primary_elements.__contains__(key):
                    gui_utils.enable_text_box(
                        tmp_stage.primary_elements[key],
                        tmp_stage.primary_elements[self.elements_to_show[key]]
                    )

    def disable_text_boxes_stage_x(self, stages_to_disable: list):
        tmp_stages_to_disable = []
        for stage_number in stages_to_disable:
            tmp_stages_to_disable.append(self.stages[stage_number])
        for tmp_stage in tmp_stages_to_disable:
            for key in self.elements_to_show:
                if tmp_stage.primary_elements.__contains__(key):
                    gui_utils.disable_text_box(
                        tmp_stage.primary_elements[key],
                        tmp_stage.primary_elements[self.elements_to_show[key]]
                    )

    def disable_all_next_buttons(self):
        for tmp_stage in self.stages:
            if tmp_stage.control_elements.keys().__contains__('next_button') is True:
                tmp_stage.control_elements['next_button'].setEnabled(False)
            else:
                if tmp_stage.control_elements.keys().__contains__('predict_button') is True:
                    tmp_stage.control_elements['predict_button'].setEnabled(False)
                else:
                    print("No relevant button!")

    def create_validation(self):
        for tmp_stage in self.stages:
            if tmp_stage.primary_elements.keys().__contains__('text_field_protein_name') and tmp_stage.primary_elements['text_field_protein_name'].isEnabled():
                tools.validate_protein_name(
                    tmp_stage.primary_elements['text_field_protein_name'],
                    tmp_stage.primary_elements['label_protein_name_status'],
                    tmp_stage.control_elements['next_button'],
                )
            elif tmp_stage.primary_elements.keys().__contains__('text_field_protein_sequence') and tmp_stage.primary_elements['text_field_protein_sequence'].isEnabled():
                tools.validate_protein_sequence(
                    tmp_stage.primary_elements['text_field_protein_sequence'],
                    tmp_stage.primary_elements['label_protein_sequence_status'],
                    tmp_stage.control_elements['next_button'],
                )

    def clear_all_text_boxes(self):
        for tmp_stage in self.stages:
            for key in self.elements_to_show:
                if tmp_stage.primary_elements.__contains__(key):
                    if key.find("text") != -1:
                        tmp_stage.primary_elements[key].clear()

    def set_empty_string_in_label(self):
        for tmp_stage in self.stages:
            for key in self.status_labels:
                if tmp_stage.primary_elements.__contains__(key):
                    tmp_stage.primary_elements[key].clear()

    @staticmethod
    def activate_specific_button(button):
        styles.color_button_ready(button)
        button.setEnabled(True)

    @staticmethod
    def deactivate_specific_button(button):
        styles.color_button_not_ready(button)
        button.setEnabled(False)

    # def __show_stage_x(self, x):
    #     if len(self.stages) < x+1:
    #         raise ValueError(f"The number of stages ({x}) is to low for this kind of operation!")
    #     gui_elements_to_show = []
    #
    #     gui_utils.show_gui_elements(self.stages[x])
    #     for stage in self.stages:
    #         if stage is not self.stages[x]:
    #             gui_utils.hide_gui_elements(stage)


