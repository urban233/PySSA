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
from pymol import Qt

from pyssa.gui.ui.forms.auto_generated.auto_dialog_advanced_prediction_configurations import Ui_Dialog
from pyssa.util import gui_utils
from pyssa.internal.data_structures.data_classes import prediction_configuration


class DialogAdvancedPredictionConfigurations(Qt.QtWidgets.QDialog):

    def __init__(self, prediction_config, parent=None):
        """Constructor

        Args:
            args
            kwargs
        """
        Qt.QtWidgets.QDialog.__init__(self, parent)
        # build ui object
        self.ui = Ui_Dialog()
        self.ui.setupUi(self)

        self.prediction_config: prediction_configuration.PredictionConfiguration = prediction_config

        self.ui.cb_amber.setChecked(self.prediction_config.amber_force_field)
        item_list_templates = [
            "none",
            "pdb70",
            "custom",
        ]
        gui_utils.fill_combo_box(self.ui.combo_box_template, item_list_templates)
        self.ui.combo_box_template.setCurrentIndex(self.ui.combo_box_template.findText(self.prediction_config.templates))

        self.ui.btn_cancel.clicked.connect(self.close_dialog)
        self.ui.btn_save.clicked.connect(self.save_config)
        self.ui.cb_amber.stateChanged.connect(self.change_amber_force_field)
        self.ui.combo_box_template.currentIndexChanged.connect(self.change_template_mode)
        self.setWindowTitle("Advanced prediction configuration")

    # @SLOT
    def close_dialog(self):
        self.close()

    def save_config(self):
        self.close()
        return self.prediction_config

    def change_amber_force_field(self):
        if self.ui.cb_amber.checkState() == 0:
            self.prediction_config.amber_force_field = False
        else:
            self.prediction_config.amber_force_field = True

    def change_template_mode(self):
        self.prediction_config.templates = self.ui.combo_box_template.currentText()
